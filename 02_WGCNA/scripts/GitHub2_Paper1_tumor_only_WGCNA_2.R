
############################################################
# Paper 1: Tumor-only WGCNA analysis
#
# Workflow:
#   Workflow-04-3A to Workflow-04-2R
#
# Research objective:
#   To construct a weighted gene co-expression network using
#   primary TCGA-LUAD tumor samples, identify co-expression
#   modules, and determine the module assignment and
#   module membership of AOC3, MS4A7, and IL7R.
#
# Input data:
#   Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData
#   Workflow-04-2D_WGCNA_input_checkpoint_DONE.txt
#
# Input data characteristics:
#   Primary TCGA-LUAD tumor samples: 513
#   Expression filtering:
#     FPKM-UQ > 1 in at least 103 of 513 samples
#   Expression transformation:
#     log2(FPKM-UQ + 1)
#
# Main analysis:
#   1. Recover and validate the WGCNA input checkpoint.
#   2. Select the soft-thresholding power.
#   3. Construct unsigned co-expression modules.
#   4. Calculate module eigengenes.
#   5. Annotate gene-module assignments.
#   6. Evaluate module membership for AOC3, MS4A7, and IL7R.
#
# Software:
#   R 4.5.1
#   WGCNA
#   ggplot2
#   here
#   AnnotationDbi
#   org.Hs.eg.db
#
# Reproducibility:
#   Random seed: 20261106
#   WGCNA computation: single-threaded
#
# Important interpretation:
#   Module membership and module assignment describe
#   co-expression relationships and do not establish
#   causal biological relationships.
############################################################


# ==========================================================
# Step 1. Project initialization and software environment
# ==========================================================

# ----------------------------------------------------------
# 1.1 Project directories
# ----------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required. Please install it before running this script.")
}

project_root <- here::here()

input_dir <- here::here("data")

result_dir <- here::here(
  "results",
  "02_tumor_only_WGCNA"
)

figure_dir <- here::here(
  "figures",
  "02_tumor_only_WGCNA"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# The checkpoint files are expected in the project data directory.
# Change this one directory only if the actual files are stored elsewhere.

checkpoint_dir <- input_dir

candidate_checkpoint <- file.path(
  checkpoint_dir,
  "Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData"
)

done_file <- file.path(
  checkpoint_dir,
  "Workflow-04-2D_WGCNA_input_checkpoint_DONE.txt"
)

# ----------------------------------------------------------
# 1.2 Reproducibility settings
# ----------------------------------------------------------

set.seed(20261106)

options(
  stringsAsFactors = FALSE
)

# Windows-compatible conservative setting:
# WGCNA will run without parallel threading.

if (!requireNamespace("WGCNA", quietly = TRUE)) {
  stop("Package 'WGCNA' is required.")
}

WGCNA::disableWGCNAThreads()

# ----------------------------------------------------------
# 1.3 Required packages and package versions
# ----------------------------------------------------------

required_packages <- c(
  "here",
  "WGCNA",
  "ggplot2",
  "AnnotationDbi",
  "org.Hs.eg.db"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "The following required packages are missing: ",
    paste(missing_packages, collapse = ", ")
  )
}

package_versions <- data.frame(
  Package = required_packages,
  Version = vapply(
    required_packages,
    function(pkg) {
      as.character(utils::packageVersion(pkg))
    },
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("Project and software environment\n")
cat("========================================\n")

cat("Project root:", project_root, "\n")
cat("Input directory:", input_dir, "\n")
cat("Result directory:", result_dir, "\n")
cat("Figure directory:", figure_dir, "\n")
cat("R version:", R.version.string, "\n")
cat("Random seed: 20261106\n")
cat("WGCNA threading: disabled\n\n")

print(
  package_versions,
  row.names = FALSE
)

write.table(
  package_versions,
  file = file.path(
    result_dir,
    "Workflow-04-3A_software_versions.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==========================================================
# Step 2. Verify the Workflow-04-2D checkpoint files
# ==========================================================

cat("\n========================================\n")
cat("Step 2. Checkpoint file verification\n")
cat("========================================\n")

if (!file.exists(candidate_checkpoint)) {
  stop(
    "Candidate checkpoint file was not found:\n",
    candidate_checkpoint
  )
}

if (!file.exists(done_file)) {
  stop(
    "Checkpoint completion file was not found:\n",
    done_file
  )
}

cat("Candidate checkpoint: FOUND\n")
cat("DONE file: FOUND\n")

checkpoint_file_info <- data.frame(
  File = c(
    candidate_checkpoint,
    done_file
  ),
  Size_bytes = file.info(
    c(candidate_checkpoint, done_file)
  )$size,
  stringsAsFactors = FALSE
)

print(
  checkpoint_file_info,
  row.names = FALSE
)


# ==========================================================
# Step 3. Read the checkpoint completion record
# ==========================================================

cat("\n========================================\n")
cat("Step 3. Checkpoint completion record\n")
cat("========================================\n")

done_content <- readLines(
  done_file,
  warn = FALSE
)

cat(
  paste(
    done_content,
    collapse = "\n"
  ),
  "\n"
)

if (length(done_content) == 0) {
  stop("The checkpoint DONE file is empty.")
}

writeLines(
  done_content,
  con = file.path(
    result_dir,
    "Workflow-04-3A_input_checkpoint_provenance.txt"
  )
)


# ==========================================================
# Step 4. Recover the WGCNA input checkpoint
# ==========================================================

cat("\n========================================\n")
cat("Step 4. Recover checkpoint objects\n")
cat("========================================\n")

# Load the checkpoint into an isolated environment first.
# This avoids unintentionally overwriting unrelated objects.

checkpoint_env <- new.env(parent = globalenv())

loaded_objects <- load(
  candidate_checkpoint,
  envir = checkpoint_env
)

cat(
  "Objects found in checkpoint:\n",
  paste(loaded_objects, collapse = ", "),
  "\n"
)

required_checkpoint_objects <- c(
  "datExpr_MAD",
  "clinical_aligned",
  "gene_MAD"
)

missing_checkpoint_objects <- setdiff(
  required_checkpoint_objects,
  loaded_objects
)

if (length(missing_checkpoint_objects) > 0) {
  stop(
    "Required objects are missing from the checkpoint: ",
    paste(
      missing_checkpoint_objects,
      collapse = ", "
    )
  )
}

# Recover only the objects required by the current workflow.

datExpr_MAD <- get(
  "datExpr_MAD",
  envir = checkpoint_env
)

clinical_aligned <- get(
  "clinical_aligned",
  envir = checkpoint_env
)

gene_MAD <- get(
  "gene_MAD",
  envir = checkpoint_env
)

rm(checkpoint_env)


# ==========================================================
# Step 5. Validate the recovered WGCNA input data
# ==========================================================

cat("\n========================================\n")
cat("Step 5. Input checkpoint validation\n")
cat("========================================\n")

# ----------------------------------------------------------
# 5.1 Expression matrix structure
# ----------------------------------------------------------

if (!is.matrix(datExpr_MAD)) {
  datExpr_MAD <- as.matrix(datExpr_MAD)
}

if (!is.numeric(datExpr_MAD)) {
  stop("datExpr_MAD must be a numeric expression matrix.")
}

if (nrow(datExpr_MAD) != 513) {
  stop(
    "Expected 513 primary tumor samples, but found ",
    nrow(datExpr_MAD),
    ". Please verify the input checkpoint."
  )
}

if (is.null(rownames(datExpr_MAD)) ||
    is.null(colnames(datExpr_MAD))) {
  stop("Expression matrix must have sample and gene identifiers.")
}

if (anyDuplicated(rownames(datExpr_MAD)) > 0) {
  stop("Duplicated sample identifiers detected in datExpr_MAD.")
}

if (anyDuplicated(colnames(datExpr_MAD)) > 0) {
  stop("Duplicated gene identifiers detected in datExpr_MAD.")
}

if (anyNA(datExpr_MAD) ||
    any(!is.finite(datExpr_MAD))) {
  stop(
    "Missing or non-finite values detected in datExpr_MAD."
  )
}

# ----------------------------------------------------------
# 5.2 Clinical sample structure
# ----------------------------------------------------------

if (!is.data.frame(clinical_aligned)) {
  stop("clinical_aligned must be a data frame.")
}

if (!"Sample" %in% colnames(clinical_aligned)) {
  stop("Column 'Sample' is missing from clinical_aligned.")
}

if (anyNA(clinical_aligned$Sample) ||
    anyDuplicated(clinical_aligned$Sample) > 0) {
  stop(
    "Missing or duplicated sample identifiers detected ",
    "in clinical_aligned."
  )
}

if (!identical(
  rownames(datExpr_MAD),
  as.character(clinical_aligned$Sample)
)) {
  stop(
    "Sample order mismatch between datExpr_MAD ",
    "and clinical_aligned."
  )
}

# ----------------------------------------------------------
# 5.3 MAD vector validation
# ----------------------------------------------------------

if (!is.numeric(gene_MAD)) {
  stop("gene_MAD must be a numeric vector.")
}

if (length(gene_MAD) != ncol(datExpr_MAD)) {
  stop(
    "gene_MAD length does not match the number of genes ",
    "in datExpr_MAD."
  )
}

if (is.null(names(gene_MAD))) {
  stop(
    "gene_MAD has no gene identifiers. ",
    "Gene-level correspondence cannot be verified."
  )
}

if (anyDuplicated(names(gene_MAD)) > 0) {
  stop("Duplicated gene identifiers detected in gene_MAD.")
}

if (!setequal(
  names(gene_MAD),
  colnames(datExpr_MAD)
)) {
  stop(
    "Gene identifiers in gene_MAD do not match ",
    "the expression matrix."
  )
}

# Align gene_MAD to the expression matrix explicitly.

gene_MAD <- gene_MAD[
  match(
    colnames(datExpr_MAD),
    names(gene_MAD)
  )
]

names(gene_MAD) <- colnames(datExpr_MAD)

if (anyNA(gene_MAD) ||
    any(!is.finite(gene_MAD))) {
  stop("Missing or non-finite values detected in gene_MAD.")
}

# ----------------------------------------------------------
# 5.4 Summary of recovered data
# ----------------------------------------------------------

cat(
  "Expression matrix:",
  nrow(datExpr_MAD),
  "samples ×",
  ncol(datExpr_MAD),
  "genes\n"
)

cat(
  "Clinical data:",
  nrow(clinical_aligned),
  "samples ×",
  ncol(clinical_aligned),
  "columns\n"
)

cat(
  "Missing expression values:",
  sum(is.na(datExpr_MAD)),
  "\n"
)

cat(
  "Non-finite expression values:",
  sum(!is.finite(datExpr_MAD)),
  "\n"
)

cat(
  "Duplicated samples:",
  anyDuplicated(rownames(datExpr_MAD)),
  "\n"
)

cat(
  "Duplicated genes:",
  anyDuplicated(colnames(datExpr_MAD)),
  "\n"
)

cat(
  "Sample order check:",
  identical(
    rownames(datExpr_MAD),
    as.character(clinical_aligned$Sample)
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

if (sum(gene_MAD > 0) != ncol(datExpr_MAD)) {
  warning(
    "The number of MAD-positive genes differs from ",
    "the number of columns in datExpr_MAD. ",
    "Please verify the upstream filtering record."
  )
}

# Save an input validation report.

input_validation <- data.frame(
  Metric = c(
    "Number of samples",
    "Number of genes",
    "Missing expression values",
    "Non-finite expression values",
    "Duplicated samples",
    "Duplicated genes",
    "MAD-positive genes",
    "Sample order matched"
  ),
  Value = c(
    nrow(datExpr_MAD),
    ncol(datExpr_MAD),
    sum(is.na(datExpr_MAD)),
    sum(!is.finite(datExpr_MAD)),
    anyDuplicated(rownames(datExpr_MAD)),
    anyDuplicated(colnames(datExpr_MAD)),
    sum(gene_MAD > 0),
    identical(
      rownames(datExpr_MAD),
      as.character(clinical_aligned$Sample)
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  input_validation,
  file = file.path(
    result_dir,
    "Workflow-04-3A_input_validation.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==========================================================
# Step 6. WGCNA soft-threshold selection
# ==========================================================

cat("\n========================================\n")
cat("Step 6. WGCNA soft-threshold selection\n")
cat("========================================\n")

# ----------------------------------------------------------
# 6.1 Analysis settings
# ----------------------------------------------------------

networkType <- "unsigned"

correlation_method <- "Pearson"

RsquaredCut <- 0.85

powers <- c(
  1:10,
  seq(12, 30, by = 2)
)

cat(
  "Network type:",
  networkType,
  "\n"
)

cat(
  "Correlation:",
  correlation_method,
  "\n"
)

cat(
  "Scale-free topology fit threshold:",
  RsquaredCut,
  "\n"
)

cat(
  "Candidate powers:",
  paste(powers, collapse = ", "),
  "\n"
)

# ----------------------------------------------------------
# 6.2 Run pickSoftThreshold
# ----------------------------------------------------------

sft <- WGCNA::pickSoftThreshold(
  datExpr_MAD,
  dataIsExpr = TRUE,
  corFnc = "cor",
  corOptions = list(
    use = "p",
    method = correlation_method
  ),
  networkType = networkType,
  powerVector = powers,
  RsquaredCut = RsquaredCut,
  moreNetworkConcepts = TRUE,
  verbose = 5
)

fit <- as.data.frame(
  sft$fitIndices
)

if (nrow(fit) == 0) {
  stop("No soft-threshold fit indices were returned.")
}

required_fit_columns <- c(
  "Power",
  "SFT.R.sq",
  "slope",
  "mean.k.",
  "median.k.",
  "max.k."
)

missing_fit_columns <- setdiff(
  required_fit_columns,
  colnames(fit)
)

if (length(missing_fit_columns) > 0) {
  stop(
    "Missing expected columns in sft$fitIndices: ",
    paste(
      missing_fit_columns,
      collapse = ", "
    )
  )
}

# ----------------------------------------------------------
# 6.3 Record the prespecified soft-threshold power
# ----------------------------------------------------------

# The previous analysis selected beta = 7.
# This value is retained only after verifying that it is
# present in the pickSoftThreshold results.

selected_power <- 7

selected_row <- which(
  fit$Power == selected_power
)

if (length(selected_row) != 1) {
  stop(
    "Selected power 7 is not uniquely represented ",
    "in the soft-threshold results."
  )
}

selected_fit <- fit[
  selected_row,
  ,
  drop = FALSE
]

cat("\nSelected soft-threshold power:\n")

cat(
  "Power (beta):",
  selected_power,
  "\n"
)

cat(
  "Scale-free topology fit (R-squared):",
  selected_fit$SFT.R.sq,
  "\n"
)

cat(
  "Mean connectivity:",
  selected_fit$mean.k.,
  "\n"
)

cat(
  "Median connectivity:",
  selected_fit$median.k.,
  "\n"
)

cat(
  "Maximum connectivity:",
  selected_fit$max.k.,
  "\n"
)

# Do not automatically change the previously selected power
# based only on an R-squared cutoff. The selected power and
# its connectivity characteristics should be interpreted
# together.

# ----------------------------------------------------------
# 6.4 Save soft-threshold results
# ----------------------------------------------------------

WGCNA_Step6 <- list(
  sft = sft,
  fitIndices = fit,
  candidate_powers = powers,
  selected_power = selected_power,
  networkType = networkType,
  correlation = correlation_method,
  RsquaredCut = RsquaredCut,
  n_samples = nrow(datExpr_MAD),
  n_genes = ncol(datExpr_MAD),
  random_seed = 20261106,
  package_version = as.character(
    utils::packageVersion("WGCNA")
  ),
  R_version = R.version.string
)

saveRDS(
  WGCNA_Step6,
  file = file.path(
    result_dir,
    "Workflow-04-2E_WGCNA_soft_threshold_checkpoint.rds"
  )
)

write.table(
  fit,
  file = file.path(
    result_dir,
    "Workflow-04-2E_WGCNA_soft_threshold_fitIndices.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 6.5 Generate soft-threshold plots
# ----------------------------------------------------------

signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq

grDevices::pdf(
  file = file.path(
    figure_dir,
    "Workflow-04-2E_soft_threshold_selection.pdf"
  ),
  width = 8,
  height = 6
)

plot(
  fit$Power,
  signed_R2,
  type = "n",
  xlab = "Soft-threshold power (beta)",
  ylab = "Scale-free topology fit (signed R-squared)",
  main = "Soft-threshold selection"
)

text(
  fit$Power,
  signed_R2,
  labels = fit$Power,
  cex = 0.8
)

abline(
  h = RsquaredCut,
  lty = 2
)

abline(
  v = selected_power,
  lty = 2
)

plot(
  fit$Power,
  fit$mean.k.,
  type = "n",
  xlab = "Soft-threshold power (beta)",
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

grDevices::dev.off()

# ----------------------------------------------------------
# 6.6 Verify the saved checkpoint
# ----------------------------------------------------------

soft_threshold_file <- file.path(
  result_dir,
  "Workflow-04-2E_WGCNA_soft_threshold_checkpoint.rds"
)

if (!file.exists(soft_threshold_file)) {
  stop("Soft-threshold checkpoint was not created.")
}

sft_check <- readRDS(
  soft_threshold_file
)

stopifnot(
  identical(
    sft_check$selected_power,
    selected_power
  )
)

cat(
  "\nStep 6 completed. Checkpoint verification passed.\n"
)


# ==========================================================
# Step 7. Tumor-only WGCNA module construction
# ==========================================================

cat("\n========================================\n")
cat("Step 7. Tumor-only WGCNA module construction\n")
cat("========================================\n")

# ----------------------------------------------------------
# 7.1 Confirm the expression matrix
# ----------------------------------------------------------

cat(
  "Expression matrix:",
  nrow(datExpr_MAD),
  "samples ×",
  ncol(datExpr_MAD),
  "genes\n"
)

stopifnot(
  identical(
    rownames(datExpr_MAD),
    as.character(clinical_aligned$Sample)
  )
)

stopifnot(
  all(is.finite(datExpr_MAD))
)

# ----------------------------------------------------------
# 7.2 WGCNA parameters
# ----------------------------------------------------------

softPower <- selected_power

networkType <- "unsigned"

TOMType <- "unsigned"

minModuleSize <- 30

deepSplit <- 3

mergeCutHeight <- 0.25

numericLabels <- TRUE

cat("\nWGCNA parameters:\n")

cat("Soft-threshold power:", softPower, "\n")
cat("Network type:", networkType, "\n")
cat("TOM type:", TOMType, "\n")
cat("Minimum module size:", minModuleSize, "\n")
cat("Deep split:", deepSplit, "\n")
cat("Merge cut height:", mergeCutHeight, "\n")
cat("Numeric labels:", numericLabels, "\n")

# ----------------------------------------------------------
# 7.3 Keep WGCNA single-threaded
# ----------------------------------------------------------

WGCNA::disableWGCNAThreads()

# ----------------------------------------------------------
# 7.4 Construct co-expression modules
# ----------------------------------------------------------

set.seed(20261106)

net <- WGCNA::blockwiseModules(
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
  numericLabels = numericLabels,
  verbose = 3
)

# ----------------------------------------------------------
# 7.5 Extract module labels and colors
# ----------------------------------------------------------

if (length(net$colors) != ncol(datExpr_MAD)) {
  stop(
    "The number of module labels does not match ",
    "the number of input genes."
  )
}

moduleColors <- WGCNA::labels2colors(
  net$colors
)

names(moduleColors) <- colnames(datExpr_MAD)

if (anyNA(moduleColors)) {
  stop("Missing module colors detected.")
}

moduleSizes <- sort(
  table(moduleColors),
  decreasing = TRUE
)

cat("\nModule summary:\n")

cat(
  "Number of detected module colors:",
  length(unique(moduleColors)),
  "\n"
)

cat(
  "Genes assigned to non-grey modules:",
  sum(moduleColors != "grey"),
  "\n"
)

cat(
  "Grey/unassigned genes:",
  sum(moduleColors == "grey"),
  "\n"
)

print(moduleSizes)

# ----------------------------------------------------------
# 7.6 Plot module size distribution
# ----------------------------------------------------------

moduleSize_df <- data.frame(
  Module = names(moduleSizes),
  GeneCount = as.numeric(moduleSizes),
  stringsAsFactors = FALSE
)

moduleSize_df$Module <- factor(
  moduleSize_df$Module,
  levels = names(moduleSizes)
)

p_module_size <- ggplot2::ggplot(
  moduleSize_df,
  ggplot2::aes(
    x = Module,
    y = GeneCount
  )
) +
  ggplot2::geom_col() +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "WGCNA module",
    y = "Number of genes",
    title = "Tumor-only WGCNA module sizes"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

ggplot2::ggsave(
  filename = file.path(
    figure_dir,
    "Workflow-04-2F_WGCNA_module_size.pdf"
  ),
  plot = p_module_size,
  width = 9,
  height = 6
)

# ----------------------------------------------------------
# 7.7 Extract and validate module eigengenes
# ----------------------------------------------------------

MEs0 <- net$MEs

MEs <- WGCNA::orderMEs(MEs0)

if (is.null(rownames(MEs)) ||
    !identical(
      rownames(MEs),
      rownames(datExpr_MAD)
    )) {
  stop(
    "Module eigengene sample order does not match ",
    "the expression matrix."
  )
}

if (anyNA(MEs) ||
    any(!is.finite(MEs))) {
  stop(
    "Missing or non-finite values detected ",
    "in module eigengenes."
  )
}

cat(
  "\nModule eigengene matrix:",
  nrow(MEs),
  "samples ×",
  ncol(MEs),
  "modules\n"
)

# ----------------------------------------------------------
# 7.8 Create gene-module assignment table
# ----------------------------------------------------------

gene_module_assignment <- data.frame(
  Gene = colnames(datExpr_MAD),
  ModuleLabel = as.numeric(net$colors),
  ModuleColor = as.character(moduleColors),
  stringsAsFactors = FALSE
)

if (nrow(gene_module_assignment) != ncol(datExpr_MAD)) {
  stop(
    "Gene-module assignment table has an unexpected number of rows."
  )
}

write.table(
  gene_module_assignment,
  file = file.path(
    result_dir,
    "Workflow-04-2F_gene_module_assignment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7.9 Create numeric-label to module-color mapping
# ----------------------------------------------------------

module_label_color <- data.frame(
  ModuleLabel = sort(unique(as.numeric(net$colors))),
  stringsAsFactors = FALSE
)

module_label_color$ModuleColor <- WGCNA::labels2colors(
  module_label_color$ModuleLabel
)

module_label_color$ME <- paste0(
  "ME",
  module_label_color$ModuleLabel
)

module_label_color$ME_present <- (
  module_label_color$ME %in% colnames(MEs)
)

write.table(
  module_label_color,
  file = file.path(
    result_dir,
    "Workflow-04-2F_module_label_color_mapping.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7.10 Save the Step 7 checkpoint
# ----------------------------------------------------------

WGCNA_Step7 <- list(
  net = net,
  moduleColors = moduleColors,
  moduleSizes = moduleSizes,
  module_label_color = module_label_color,
  gene_module_assignment = gene_module_assignment,
  MEs = MEs,
  softPower = softPower,
  networkType = networkType,
  TOMType = TOMType,
  minModuleSize = minModuleSize,
  deepSplit = deepSplit,
  mergeCutHeight = mergeCutHeight,
  numericLabels = numericLabels,
  nSamples = nrow(datExpr_MAD),
  nGenes = ncol(datExpr_MAD),
  random_seed = 20261106,
  package_version = as.character(
    utils::packageVersion("WGCNA")
  ),
  R_version = R.version.string
)

module_checkpoint_file <- file.path(
  result_dir,
  "Workflow-04-2F_WGCNA_tumor_only_modules.rds"
)

saveRDS(
  WGCNA_Step7,
  file = module_checkpoint_file
)

if (!file.exists(module_checkpoint_file)) {
  stop("Step 7 module checkpoint was not created.")
}

cat(
  "\nStep 7 completed. Module checkpoint saved.\n"
)


# ==========================================================
# Step 7A. Annotation and module membership of target genes
# ==========================================================

cat("\n========================================\n")
cat("Step 7A. Target gene annotation and module membership\n")
cat("========================================\n")

# ----------------------------------------------------------
# 7A.1 Define the three prioritized study genes
# ----------------------------------------------------------

target_annotation <- data.frame(
  Gene = c(
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


# ----------------------------------------------------------
# 7A.2 Check target genes using version-independent IDs
# ----------------------------------------------------------

target_annotation$Ensembl_base <- sub(
  "\\.[0-9]+$",
  "",
  target_annotation$Ensembl
)

all_ensembl <- colnames(datExpr_MAD)

all_ensembl_base <- sub(
  "\\.[0-9]+$",
  "",
  all_ensembl
)

if (anyDuplicated(all_ensembl_base) > 0) {
  stop(
    "Duplicated Ensembl identifiers were detected ",
    "after removing version suffixes. ",
    "Target-gene mapping is ambiguous and requires ",
    "manual resolution before proceeding."
  )
}

target_annotation$WGCNA_Ensembl <- all_ensembl[
  match(
    target_annotation$Ensembl_base,
    all_ensembl_base
  )
]

target_annotation$In_WGCNA <- !is.na(
  target_annotation$WGCNA_Ensembl
)

target_annotation$ModuleLabel <- NA_real_

target_annotation$ModuleColor <- NA_character_

matched_target_idx <- match(
  target_annotation$WGCNA_Ensembl,
  gene_module_assignment$Gene
)

target_annotation$ModuleLabel <- (
  gene_module_assignment$ModuleLabel[
    matched_target_idx
  ]
)

target_annotation$ModuleColor <- (
  gene_module_assignment$ModuleColor[
    matched_target_idx
  ]
)

print(
  target_annotation,
  row.names = FALSE
)

write.table(
  target_annotation,
  file = file.path(
    result_dir,
    "Workflow-04-2F_target_gene_module_assignment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7A.3 Verify Ensembl identifiers and annotation package
# ----------------------------------------------------------

# Remove Ensembl version suffixes for annotation lookup,
# while retaining the original identifiers in the results.

all_ensembl <- colnames(datExpr_MAD)

all_ensembl_base <- sub(
  "\\.[0-9]+$",
  "",
  all_ensembl
)

if (anyDuplicated(all_ensembl_base) > 0) {
  warning(
    "Duplicate Ensembl identifiers were generated ",
    "after removing version suffixes. ",
    "Annotation mappings will be retained without ",
    "assuming one-to-one correspondence."
  )
}

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(all_ensembl_base),
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL"
  )
)

gene_annotation <- unique(
  gene_annotation
)

# Retain the original WGCNA Ensembl identifier.

gene_annotation$Original_Ensembl <- all_ensembl[
  match(
    gene_annotation$ENSEMBL,
    all_ensembl_base
  )
]

# Check mapping characteristics.

cat(
  "\nNumber of WGCNA genes:",
  length(all_ensembl),
  "\n"
)

cat(
  "Number of annotation records:",
  nrow(gene_annotation),
  "\n"
)

cat(
  "Number of unique annotated Ensembl IDs:",
  length(unique(gene_annotation$ENSEMBL)),
  "\n"
)

cat(
  "Number of records without gene symbols:",
  sum(
    is.na(gene_annotation$SYMBOL) |
      gene_annotation$SYMBOL == ""
  ),
  "\n"
)

write.table(
  gene_annotation,
  file = file.path(
    result_dir,
    "Workflow-04-2F_gene_Ensembl_symbol_annotation.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7A.4 Calculate gene-module membership correlations
# ----------------------------------------------------------

# Module membership (MM) is defined here as the Pearson
# correlation between gene expression and the corresponding
# module eigengene.

MM_all <- stats::cor(
  datExpr_MAD,
  MEs,
  use = "pairwise.complete.obs",
  method = "pearson"
)

if (is.null(rownames(MM_all)) ||
    is.null(colnames(MM_all))) {
  stop("MM_all does not have valid gene and module identifiers.")
}

# ----------------------------------------------------------
# 7A.5 Map target genes to the correct module eigengene
# ----------------------------------------------------------

target_MM <- target_annotation

target_MM$ME <- NA_character_

target_MM$MM <- NA_real_

for (i in seq_len(nrow(target_MM))) {
  
  if (!isTRUE(target_MM$In_WGCNA[i])) {
    next
  }
  
  current_label <- target_MM$ModuleLabel[i]
  
  current_mapping <- module_label_color[
    module_label_color$ModuleLabel == current_label,
    ,
    drop = FALSE
  ]
  
  if (nrow(current_mapping) != 1) {
    warning(
      "Could not uniquely map module label for ",
      target_MM$Gene[i]
    )
    next
  }
  
  current_ME <- current_mapping$ME[1]
  
  if (!current_ME %in% colnames(MM_all)) {
    warning(
      "Module eigengene ",
      current_ME,
      " is unavailable for ",
      target_MM$Gene[i],
      ". This may occur for genes assigned to the grey module."
    )
    next
  }
  
  target_MM$ME[i] <- current_ME
  
  target_MM$MM[i] <- MM_all[
    target_MM$WGCNA_Ensembl[i],
    current_ME
  ]
}

print(
  target_MM,
  row.names = FALSE
)

write.table(
  target_MM,
  file = file.path(
    result_dir,
    "Workflow-04-2F_target_gene_module_membership.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 7A.6 Final target-gene audit
# ----------------------------------------------------------

target_gene_audit <- data.frame(
  Gene = target_MM$Gene,
  Ensembl = target_MM$Ensembl,
  In_WGCNA = target_MM$In_WGCNA,
  ModuleLabel = target_MM$ModuleLabel,
  ModuleColor = target_MM$ModuleColor,
  ME = target_MM$ME,
  MM = target_MM$MM,
  stringsAsFactors = FALSE
)

cat("\nFinal target-gene audit:\n")

print(
  target_gene_audit,
  row.names = FALSE
)

# Save a consolidated Step 7A checkpoint.

WGCNA_Step7A <- list(
  target_annotation = target_annotation,
  gene_annotation = gene_annotation,
  module_label_color = module_label_color,
  target_MM = target_MM,
  MM_all = MM_all,
  target_gene_audit = target_gene_audit
)

saveRDS(
  WGCNA_Step7A,
  file = file.path(
    result_dir,
    "Workflow-04-2F_target_gene_annotation_MM_checkpoint.rds"
  )
)

cat(
  "\nStep 7A completed.\n"
)

cat(
  "All Step 1–7A outputs were saved under:\n",
  result_dir,
  "\n"
)


############################################################
# Paper 1: Tumor-only WGCNA analysis
#
# Part 2: Clinical trait preparation and module-trait
#          association analysis
#
# Workflow:
#   Step 8-1 to Step 8B-3B
#
# Research objectives:
#   1. Validate the clinical data aligned to the
#      513 primary TCGA-LUAD tumor samples.
#   2. Construct standardized clinical trait variables.
#   3. Extract biological module eigengenes, excluding
#      the grey/unassigned module.
#   4. Evaluate associations between module eigengenes
#      and age, gender, pathological stage, and T stage.
#
# Input objects:
#   clinical_aligned
#   datExpr_MAD
#   MEs
#   net
#   moduleColors
#   module_label_color
#   MM_all
#   target_MM
#
# Statistical methods:
#   Age: Pearson correlation
#   Gender: two-sided Wilcoxon rank-sum test
#   Pathological stage: Kruskal-Wallis test
#   Pathological T stage: Kruskal-Wallis test
#
# Multiple-testing correction:
#   Benjamini-Hochberg (BH), separately within
#   each prespecified clinical-trait analysis.
#
# Survival-time handling:
#   Documented OS.time = 0 observations are retained
#   in the clinical trait data and are not arbitrarily
#   modified or automatically excluded.
#
# Interpretation:
#   These analyses evaluate associations and do not
#   establish causal relationships.
############################################################


# ==========================================================
# Step 8-1. Inspect and validate clinical_aligned
# ==========================================================

cat("\n========================================\n")
cat("Step 8-1. Clinical data validation\n")
cat("========================================\n")

# ----------------------------------------------------------
# 8-1.1 Confirm required objects
# ----------------------------------------------------------

required_objects <- c(
  "clinical_aligned",
  "datExpr_MAD",
  "MEs",
  "net",
  "moduleColors",
  "module_label_color",
  "MM_all",
  "target_MM"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    envir = .GlobalEnv,
    inherits = FALSE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The following objects from Part 1 are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ----------------------------------------------------------
# 8-1.2 Check clinical data structure
# ----------------------------------------------------------

if (!is.data.frame(clinical_aligned)) {
  stop("clinical_aligned must be a data frame.")
}

required_clinical_columns <- c(
  "Sample",
  "Age",
  "Gender",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m",
  "OS.status",
  "OS.time"
)

missing_clinical_columns <- setdiff(
  required_clinical_columns,
  colnames(clinical_aligned)
)

if (length(missing_clinical_columns) > 0) {
  stop(
    "Required clinical columns are missing: ",
    paste(
      missing_clinical_columns,
      collapse = ", "
    )
  )
}

if (anyNA(clinical_aligned$Sample) ||
    anyDuplicated(clinical_aligned$Sample) > 0) {
  stop(
    "Missing or duplicated Sample identifiers ",
    "were detected in clinical_aligned."
  )
}

if (nrow(clinical_aligned) != 513) {
  stop(
    "Expected 513 clinical samples, but found ",
    nrow(clinical_aligned),
    ". Please verify the upstream cohort."
  )
}

if (!identical(
  as.character(rownames(datExpr_MAD)),
  as.character(clinical_aligned$Sample)
)) {
  stop(
    "Sample order mismatch between datExpr_MAD ",
    "and clinical_aligned."
  )
}

cat(
  "Clinical data dimensions:",
  nrow(clinical_aligned),
  "samples ×",
  ncol(clinical_aligned),
  "variables\n"
)

cat("Sample identifiers: unique\n")
cat("Sample order: matched to datExpr_MAD\n")

# ----------------------------------------------------------
# 8-1.3 Summarize clinical variables
# ----------------------------------------------------------

clinical_structure <- data.frame(
  Variable = colnames(clinical_aligned),
  Class = vapply(
    clinical_aligned,
    function(x) {
      paste(class(x), collapse = "/")
    },
    FUN.VALUE = character(1)
  ),
  Missing_N = vapply(
    clinical_aligned,
    function(x) {
      sum(is.na(x))
    },
    FUN.VALUE = integer(1)
  ),
  stringsAsFactors = FALSE
)

print(
  clinical_structure,
  row.names = FALSE
)

write.table(
  clinical_structure,
  file = file.path(
    result_dir,
    "Workflow-04-2G_clinical_data_structure.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 8-1.4 Inspect the original clinical category distributions
# ----------------------------------------------------------

clinical_categorical_variables <- c(
  "Gender",
  "pathologic_stage_2",
  "pathologic_t",
  "pathologic_n",
  "pathologic_m",
  "OS.status"
)

clinical_category_summary <- do.call(
  rbind,
  lapply(
    clinical_categorical_variables,
    function(variable) {
      
      x <- as.character(
        clinical_aligned[[variable]]
      )
      
      tab <- table(
        x,
        useNA = "ifany"
      )
      
      data.frame(
        Variable = variable,
        Category = names(tab),
        N = as.integer(tab),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  clinical_category_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2G_original_clinical_category_distribution.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(
  clinical_category_summary,
  row.names = FALSE
)

# ----------------------------------------------------------
# 8-1.5 Inspect continuous clinical variables
# ----------------------------------------------------------

clinical_continuous_summary <- data.frame(
  Variable = c(
    "Age",
    "OS.time"
  ),
  N_nonmissing = c(
    sum(!is.na(clinical_aligned$Age)),
    sum(!is.na(clinical_aligned$OS.time))
  ),
  Minimum = c(
    suppressWarnings(
      min(
        as.numeric(as.character(clinical_aligned$Age)),
        na.rm = TRUE
      )
    ),
    suppressWarnings(
      min(
        as.numeric(as.character(clinical_aligned$OS.time)),
        na.rm = TRUE
      )
    )
  ),
  Median = c(
    median(
      as.numeric(as.character(clinical_aligned$Age)),
      na.rm = TRUE
    ),
    median(
      as.numeric(as.character(clinical_aligned$OS.time)),
      na.rm = TRUE
    )
  ),
  Maximum = c(
    suppressWarnings(
      max(
        as.numeric(as.character(clinical_aligned$Age)),
        na.rm = TRUE
      )
    ),
    suppressWarnings(
      max(
        as.numeric(as.character(clinical_aligned$OS.time)),
        na.rm = TRUE
      )
    )
  ),
  stringsAsFactors = FALSE
)

print(
  clinical_continuous_summary,
  row.names = FALSE
)

write.table(
  clinical_continuous_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2G_clinical_continuous_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==========================================================
# Step 8-pre. Save the WGCNA and target-gene checkpoint
# ==========================================================

cat("\n========================================\n")
cat("Step 8-pre. Save pre-clinical-trait checkpoint\n")
cat("========================================\n")

# The first part defines MEs as the ordered module
# eigengene matrix. Do not refer to MEs_ordered here.

WGCNA_Step8_pre <- list(
  net = net,
  moduleColors = moduleColors,
  module_label_color = module_label_color,
  MEs = MEs,
  MM_all = MM_all,
  target_MM = target_MM,
  clinical_aligned = clinical_aligned,
  datExpr_MAD = datExpr_MAD,
  softPower = softPower,
  networkType = networkType,
  TOMType = TOMType,
  minModuleSize = minModuleSize,
  deepSplit = deepSplit,
  mergeCutHeight = mergeCutHeight,
  random_seed = 20261106,
  R_version = R.version.string,
  WGCNA_version = as.character(
    utils::packageVersion("WGCNA")
  )
)

saveRDS(
  WGCNA_Step8_pre,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_tumor_only_MM_checkpoint.rds"
  )
)

cat("Pre-clinical-trait checkpoint saved.\n")


# ==========================================================
# Step 8A. Construct the clinical trait matrix
# ==========================================================

cat("\n========================================\n")
cat("Step 8A. Clinical trait matrix construction\n")
cat("========================================\n")

clinical_WGCNA_traits <- clinical_aligned

# ----------------------------------------------------------
# 8A.1 Age
# ----------------------------------------------------------

clinical_WGCNA_traits$Age_trait <- suppressWarnings(
  as.numeric(
    as.character(
      clinical_WGCNA_traits$Age
    )
  )
)

if (any(
  !is.na(clinical_aligned$Age) &
  is.na(clinical_WGCNA_traits$Age_trait)
)) {
  warning(
    "Some non-missing Age values could not be converted ",
    "to numeric values."
  )
}

# ----------------------------------------------------------
# 8A.2 Gender
#
# female = 0
# male   = 1
# ----------------------------------------------------------

gender_raw <- tolower(
  trimws(
    as.character(
      clinical_WGCNA_traits$Gender
    )
  )
)

clinical_WGCNA_traits$Gender_trait <- NA_real_

clinical_WGCNA_traits$Gender_trait[
  gender_raw == "female"
] <- 0

clinical_WGCNA_traits$Gender_trait[
  gender_raw == "male"
] <- 1

# ----------------------------------------------------------
# 8A.3 Pathological stage
#
# Stage I:
#   Stage I, Stage IA, Stage IB
#
# Stage II:
#   Stage II, Stage IIA, Stage IIB
#
# Stage III:
#   Stage IIIA, Stage IIIB
#
# Stage IV:
#   Stage IV
# ----------------------------------------------------------

stage_raw <- trimws(
  as.character(
    clinical_WGCNA_traits$pathologic_stage_2
  )
)

clinical_WGCNA_traits$Stage_trait <- NA_character_

clinical_WGCNA_traits$Stage_trait[
  stage_raw %in% c(
    "Stage I",
    "Stage IA",
    "Stage IB"
  )
] <- "Stage I"

clinical_WGCNA_traits$Stage_trait[
  stage_raw %in% c(
    "Stage II",
    "Stage IIA",
    "Stage IIB"
  )
] <- "Stage II"

clinical_WGCNA_traits$Stage_trait[
  stage_raw %in% c(
    "Stage IIIA",
    "Stage IIIB"
  )
] <- "Stage III"

clinical_WGCNA_traits$Stage_trait[
  stage_raw == "Stage IV"
] <- "Stage IV"

clinical_WGCNA_traits$Stage_trait <- factor(
  clinical_WGCNA_traits$Stage_trait,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  ),
  ordered = TRUE
)

# ----------------------------------------------------------
# 8A.4 Pathological T stage
#
# T1: T1, T1a, T1b
# T2: T2, T2a, T2b
# T3: T3
# T4: T4
#
# TX and unrecognized categories are excluded.
# ----------------------------------------------------------

T_raw <- toupper(
  trimws(
    as.character(
      clinical_WGCNA_traits$pathologic_t
    )
  )
)

clinical_WGCNA_traits$T_trait <- NA_character_

clinical_WGCNA_traits$T_trait[
  T_raw %in% c(
    "T1",
    "T1A",
    "T1B"
  )
] <- "T1"

clinical_WGCNA_traits$T_trait[
  T_raw %in% c(
    "T2",
    "T2A",
    "T2B"
  )
] <- "T2"

clinical_WGCNA_traits$T_trait[
  T_raw == "T3"
] <- "T3"

clinical_WGCNA_traits$T_trait[
  T_raw == "T4"
] <- "T4"

clinical_WGCNA_traits$T_trait <- factor(
  clinical_WGCNA_traits$T_trait,
  levels = c(
    "T1",
    "T2",
    "T3",
    "T4"
  ),
  ordered = TRUE
)

# ----------------------------------------------------------
# 8A.5 Pathological N stage
#
# Retain N0, N1, N2, N3.
# NX and unrecognized categories are excluded.
# ----------------------------------------------------------

N_raw <- toupper(
  trimws(
    as.character(
      clinical_WGCNA_traits$pathologic_n
    )
  )
)

clinical_WGCNA_traits$N_trait <- NA_character_

clinical_WGCNA_traits$N_trait[
  N_raw %in% c(
    "N0",
    "N1",
    "N2",
    "N3"
  )
] <- N_raw[
  N_raw %in% c(
    "N0",
    "N1",
    "N2",
    "N3"
  )
]

clinical_WGCNA_traits$N_trait <- factor(
  clinical_WGCNA_traits$N_trait,
  levels = c(
    "N0",
    "N1",
    "N2",
    "N3"
  ),
  ordered = TRUE
)

# ----------------------------------------------------------
# 8A.6 Pathological M stage
#
# M0: M0
# M1: M1, M1a, M1b
#
# MX and unrecognized categories are excluded.
# ----------------------------------------------------------

M_raw <- toupper(
  trimws(
    as.character(
      clinical_WGCNA_traits$pathologic_m
    )
  )
)

clinical_WGCNA_traits$M_trait <- NA_character_

clinical_WGCNA_traits$M_trait[
  M_raw == "M0"
] <- "M0"

clinical_WGCNA_traits$M_trait[
  M_raw %in% c(
    "M1",
    "M1A",
    "M1B"
  )
] <- "M1"

clinical_WGCNA_traits$M_trait <- factor(
  clinical_WGCNA_traits$M_trait,
  levels = c(
    "M0",
    "M1"
  )
)

# ----------------------------------------------------------
# 8A.7 Overall survival variables
#
# OS_event:
#   Dead  = 1
#   Alive = 0
#
# OS_time:
#   Original OS.time, with zero values retained.
# ----------------------------------------------------------

OS_status_raw <- tolower(
  trimws(
    as.character(
      clinical_WGCNA_traits$OS.status
    )
  )
)

clinical_WGCNA_traits$OS_time <- suppressWarnings(
  as.numeric(
    as.character(
      clinical_WGCNA_traits$OS.time
    )
  )
)

clinical_WGCNA_traits$OS_event <- NA_integer_

clinical_WGCNA_traits$OS_event[
  OS_status_raw == "dead"
] <- 1L

clinical_WGCNA_traits$OS_event[
  OS_status_raw == "alive"
] <- 0L

# ----------------------------------------------------------
# 8A.8 Audit unrecognized clinical categories
# ----------------------------------------------------------

unrecognized_category_audit <- data.frame(
  Variable = character(0),
  Original_Category = character(0),
  N = integer(0),
  stringsAsFactors = FALSE
)

category_mapping_inputs <- list(
  Gender = gender_raw,
  Stage = stage_raw,
  T_stage = T_raw,
  N_stage = N_raw,
  M_stage = M_raw,
  OS_status = OS_status_raw
)

recognized_categories <- list(
  Gender = c("female", "male"),
  Stage = c(
    "Stage I",
    "Stage IA",
    "Stage IB",
    "Stage II",
    "Stage IIA",
    "Stage IIB",
    "Stage IIIA",
    "Stage IIIB",
    "Stage IV"
  ),
  T_stage = c(
    "T1", "T1A", "T1B",
    "T2", "T2A", "T2B",
    "T3", "T4"
  ),
  N_stage = c(
    "N0", "N1", "N2", "N3"
  ),
  M_stage = c(
    "M0", "M1", "M1A", "M1B"
  ),
  OS_status = c(
    "dead", "alive"
  )
)

for (variable in names(category_mapping_inputs)) {
  
  current_values <- category_mapping_inputs[[variable]]
  
  current_recognized <- recognized_categories[[variable]]
  
  unknown_values <- current_values[
    !is.na(current_values) &
      current_values != "" &
      !current_values %in% current_recognized
  ]
  
  if (length(unknown_values) > 0) {
    
    unknown_table <- table(
      unknown_values
    )
    
    unrecognized_category_audit <- rbind(
      unrecognized_category_audit,
      data.frame(
        Variable = variable,
        Original_Category = names(unknown_table),
        N = as.integer(unknown_table),
        stringsAsFactors = FALSE
      )
    )
  }
}

if (nrow(unrecognized_category_audit) > 0) {
  
  warning(
    "Unrecognized clinical categories were detected. ",
    "They were retained as missing values in the ",
    "derived trait variables and have been documented."
  )
  
  print(
    unrecognized_category_audit,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No unrecognized non-missing clinical categories detected.\n"
  )
}

write.table(
  unrecognized_category_audit,
  file = file.path(
    result_dir,
    "Workflow-04-2G_unrecognized_clinical_categories.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 8A.9 Check the derived clinical traits
# ----------------------------------------------------------

clinical_trait_summary <- data.frame(
  Variable = c(
    "Age_trait",
    "Gender_trait",
    "Stage_trait",
    "T_trait",
    "N_trait",
    "M_trait",
    "OS_time",
    "OS_event"
  ),
  N_nonmissing = c(
    sum(!is.na(clinical_WGCNA_traits$Age_trait)),
    sum(!is.na(clinical_WGCNA_traits$Gender_trait)),
    sum(!is.na(clinical_WGCNA_traits$Stage_trait)),
    sum(!is.na(clinical_WGCNA_traits$T_trait)),
    sum(!is.na(clinical_WGCNA_traits$N_trait)),
    sum(!is.na(clinical_WGCNA_traits$M_trait)),
    sum(!is.na(clinical_WGCNA_traits$OS_time)),
    sum(!is.na(clinical_WGCNA_traits$OS_event))
  ),
  N_missing = c(
    sum(is.na(clinical_WGCNA_traits$Age_trait)),
    sum(is.na(clinical_WGCNA_traits$Gender_trait)),
    sum(is.na(clinical_WGCNA_traits$Stage_trait)),
    sum(is.na(clinical_WGCNA_traits$T_trait)),
    sum(is.na(clinical_WGCNA_traits$N_trait)),
    sum(is.na(clinical_WGCNA_traits$M_trait)),
    sum(is.na(clinical_WGCNA_traits$OS_time)),
    sum(is.na(clinical_WGCNA_traits$OS_event))
  ),
  stringsAsFactors = FALSE
)

print(
  clinical_trait_summary,
  row.names = FALSE
)

write.table(
  clinical_trait_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2G_derived_clinical_trait_summary.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 8A.10 Audit zero-valued survival times
# ----------------------------------------------------------

zero_time_idx <- which(
  !is.na(clinical_WGCNA_traits$OS_time) &
    is.finite(clinical_WGCNA_traits$OS_time) &
    clinical_WGCNA_traits$OS_time == 0
)

zero_time_audit <- clinical_WGCNA_traits[
  zero_time_idx,
  c(
    "Sample",
    "OS.status",
    "OS.time",
    "OS_time",
    "OS_event"
  ),
  drop = FALSE
]

cat(
  "\nNumber of observations with OS_time = 0:",
  nrow(zero_time_audit),
  "\n"
)

if (nrow(zero_time_audit) > 0) {
  print(
    zero_time_audit,
    row.names = FALSE
  )
}

write.table(
  zero_time_audit,
  file = file.path(
    result_dir,
    "Workflow-04-2G_zero_OS_time_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Do not remove, modify, or automatically exclude
# zero-time observations at this stage.

# ----------------------------------------------------------
# 8A.11 Verify sample alignment
# ----------------------------------------------------------

if (!identical(
  as.character(rownames(MEs)),
  as.character(clinical_WGCNA_traits$Sample)
)) {
  stop(
    "Sample order mismatch between MEs ",
    "and clinical_WGCNA_traits."
  )
}

cat(
  "\nClinical trait matrix construction completed.\n"
)

cat(
  "Number of samples:",
  nrow(clinical_WGCNA_traits),
  "\n"
)


# ==========================================================
# Step 8B-1. Prepare biological module eigengenes
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-1. Prepare biological module eigengenes\n")
cat("========================================\n")

# ----------------------------------------------------------
# 8B-1.1 Confirm ME column naming
# ----------------------------------------------------------

if (is.null(colnames(MEs))) {
  stop("MEs has no module identifiers.")
}

ME_labels <- suppressWarnings(
  as.integer(
    sub(
      "^ME",
      "",
      colnames(MEs)
    )
  )
)

if (anyNA(ME_labels)) {
  stop(
    "ME column names do not follow the expected ",
    "numeric-label format (e.g. ME1, ME2)."
  )
}

# ----------------------------------------------------------
# 8B-1.2 Match ME labels to module colors
# ----------------------------------------------------------

if (!all(
  c(
    "ModuleLabel",
    "ModuleColor",
    "ME"
  ) %in% colnames(module_label_color)
)) {
  stop(
    "module_label_color does not contain the expected ",
    "ModuleLabel, ModuleColor, and ME columns."
  )
}

ME_colors <- module_label_color$ModuleColor[
  match(
    ME_labels,
    module_label_color$ModuleLabel
  )
]

if (anyNA(ME_colors)) {
  stop(
    "Some module eigengenes could not be matched ",
    "to their module colors."
  )
}

# ----------------------------------------------------------
# 8B-1.3 Exclude the grey/unassigned module
# ----------------------------------------------------------

keep_ME <- ME_colors != "grey"

if (!any(keep_ME)) {
  stop("No biological module eigengenes remain after grey exclusion.")
}

MEs_biological <- MEs[
  ,
  keep_ME,
  drop = FALSE
]

# Rename ME columns using module colors.

colnames(MEs_biological) <- paste0(
  "ME_",
  ME_colors[keep_ME]
)

if (anyDuplicated(colnames(MEs_biological)) > 0) {
  stop(
    "Duplicated biological module names detected."
  )
}

# ----------------------------------------------------------
# 8B-1.4 Verify sample alignment and data integrity
# ----------------------------------------------------------

if (!identical(
  rownames(MEs_biological),
  as.character(clinical_WGCNA_traits$Sample)
)) {
  stop(
    "Sample order mismatch between biological MEs ",
    "and clinical_WGCNA_traits."
  )
}

if (anyNA(MEs_biological) ||
    any(!is.finite(MEs_biological))) {
  stop(
    "Missing or non-finite values detected ",
    "in MEs_biological."
  )
}

cat(
  "Biological module eigengene matrix:",
  nrow(MEs_biological),
  "samples ×",
  ncol(MEs_biological),
  "modules\n"
)

cat(
  "Grey module excluded:",
  sum(!keep_ME),
  "module(s)\n"
)

cat(
  "Sample order check: PASSED\n"
)

# Save biological module eigengenes.

saveRDS(
  MEs_biological,
  file = file.path(
    result_dir,
    "Workflow-04-2G_biological_module_eigengenes.rds"
  )
)

write.table(
  MEs_biological,
  file = file.path(
    result_dir,
    "Workflow-04-2G_biological_module_eigengenes.txt"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ==========================================================
# Step 8B-2A. Module eigengene association with Age
#
# Statistical test:
#   Pearson correlation
#
# Multiple-testing correction:
#   BH across all biological modules
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-2A. Module association with age\n")
cat("========================================\n")

age_results <- data.frame(
  Module = colnames(MEs_biological),
  Correlation = NA_real_,
  P_value = NA_real_,
  FDR = NA_real_,
  N = NA_integer_,
  Test_status = "Not tested",
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  
  age_values <- clinical_WGCNA_traits$Age_trait
  
  complete_idx <- complete.cases(
    module_values,
    age_values
  )
  
  module_complete <- module_values[complete_idx]
  
  age_complete <- age_values[complete_idx]
  
  age_results$N[i] <- length(age_complete)
  
  if (length(age_complete) < 3) {
    age_results$Test_status[i] <- "Insufficient observations"
    next
  }
  
  if (length(unique(age_complete)) < 2 ||
      length(unique(module_complete)) < 2) {
    age_results$Test_status[i] <- "No variation"
    next
  }
  
  test <- tryCatch(
    stats::cor.test(
      module_complete,
      age_complete,
      method = "pearson"
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(test)) {
    age_results$Test_status[i] <- "Test failed"
    next
  }
  
  age_results$Correlation[i] <- unname(
    test$estimate
  )
  
  age_results$P_value[i] <- test$p.value
  
  age_results$Test_status[i] <- "OK"
}

# BH correction across the complete prespecified
# biological-module testing family.

age_results$FDR <- p.adjust(
  age_results$P_value,
  method = "BH"
)

age_results <- age_results[
  order(
    age_results$P_value,
    na.last = TRUE
  ),
]

rownames(age_results) <- NULL

print(
  age_results,
  row.names = FALSE
)

write.table(
  age_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_age_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  age_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_age_association.rds"
  )
)


# ==========================================================
# Step 8B-2B. Module eigengene association with Gender
#
# Statistical test:
#   Two-sided Wilcoxon rank-sum test
#
# Gender coding:
#   Female = 0
#   Male   = 1
#
# Multiple-testing correction:
#   BH across all biological modules
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-2B. Module association with gender\n")
cat("========================================\n")

gender_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  FDR = NA_real_,
  Female_median = NA_real_,
  Male_median = NA_real_,
  N_female = NA_integer_,
  N_male = NA_integer_,
  N = NA_integer_,
  Test_status = "Not tested",
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
  
  gender_results$N_female[i] <- length(female_values)
  
  gender_results$N_male[i] <- length(male_values)
  
  gender_results$N[i] <- (
    length(female_values) +
      length(male_values)
  )
  
  if (length(female_values) > 0) {
    gender_results$Female_median[i] <- median(
      female_values
    )
  }
  
  if (length(male_values) > 0) {
    gender_results$Male_median[i] <- median(
      male_values
    )
  }
  
  if (length(female_values) < 2 ||
      length(male_values) < 2) {
    gender_results$Test_status[i] <- "Insufficient group size"
    next
  }
  
  test <- tryCatch(
    stats::wilcox.test(
      female_values,
      male_values,
      alternative = "two.sided",
      exact = FALSE
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(test)) {
    gender_results$Test_status[i] <- "Test failed"
    next
  }
  
  gender_results$P_value[i] <- test$p.value
  
  gender_results$Test_status[i] <- "OK"
}

gender_results$FDR <- p.adjust(
  gender_results$P_value,
  method = "BH"
)

gender_results <- gender_results[
  order(
    gender_results$P_value,
    na.last = TRUE
  ),
]

rownames(gender_results) <- NULL

print(
  gender_results,
  row.names = FALSE
)

write.table(
  gender_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_gender_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  gender_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_gender_association.rds"
  )
)


# ==========================================================
# Step 8B-3A. Module eigengene association with
#             pathological stage
#
# Statistical test:
#   Kruskal-Wallis test
#
# Stage groups:
#   Stage I, Stage II, Stage III, Stage IV
#
# Multiple-testing correction:
#   BH across all biological modules
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-3A. Module association with pathological stage\n")
cat("========================================\n")

stage_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  FDR = NA_real_,
  N = NA_integer_,
  N_Stage_I = NA_integer_,
  N_Stage_II = NA_integer_,
  N_Stage_III = NA_integer_,
  N_Stage_IV = NA_integer_,
  Test_status = "Not tested",
  stringsAsFactors = FALSE
)

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
  
  stage_results$N[i] <- length(module_complete)
  
  stage_results$N_Stage_I[i] <- sum(
    stage_complete == "Stage I"
  )
  
  stage_results$N_Stage_II[i] <- sum(
    stage_complete == "Stage II"
  )
  
  stage_results$N_Stage_III[i] <- sum(
    stage_complete == "Stage III"
  )
  
  stage_results$N_Stage_IV[i] <- sum(
    stage_complete == "Stage IV"
  )
  
  if (nlevels(stage_complete) < 2) {
    stage_results$Test_status[i] <- "Fewer than two observed groups"
    next
  }
  
  test <- tryCatch(
    stats::kruskal.test(
      module_complete,
      stage_complete
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(test)) {
    stage_results$Test_status[i] <- "Test failed"
    next
  }
  
  stage_results$P_value[i] <- test$p.value
  
  stage_results$Test_status[i] <- "OK"
}

stage_results$FDR <- p.adjust(
  stage_results$P_value,
  method = "BH"
)

stage_results <- stage_results[
  order(
    stage_results$P_value,
    na.last = TRUE
  ),
]

rownames(stage_results) <- NULL

print(
  stage_results,
  row.names = FALSE
)

# Significant modules are defined by BH-adjusted FDR < 0.05.

stage_sig <- stage_results[
  !is.na(stage_results$FDR) &
    stage_results$FDR < 0.05,
  ,
  drop = FALSE
]

write.table(
  stage_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_stage_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  stage_results,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_stage_association.rds"
  )
)

write.table(
  stage_sig,
  file = file.path(
    result_dir,
    "Workflow-04-2G_WGCNA_stage_significant_modules.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nNumber of Stage-associated modules (FDR < 0.05):",
  nrow(stage_sig),
  "\n"
)


# ==========================================================
# Step 8B-3B. Module eigengene association with
#             pathological T stage
#
# Statistical test:
#   Kruskal-Wallis test
#
# T-stage groups:
#   T1, T2, T3, T4
#
# TX and missing values are excluded.
#
# Multiple-testing correction:
#   BH across all biological modules
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-3B. Module association with pathological T stage\n")
cat("========================================\n")

T_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  FDR = NA_real_,
  N = NA_integer_,
  N_T1 = NA_integer_,
  N_T2 = NA_integer_,
  N_T3 = NA_integer_,
  N_T4 = NA_integer_,
  Test_status = "Not tested",
  stringsAsFactors = FALSE
)

for (i in seq_len(ncol(MEs_biological))) {
  
  module_values <- MEs_biological[, i]
  
  T_values <- clinical_WGCNA_traits$T_trait
  
  complete_idx <- complete.cases(
    module_values,
    T_values
  )
  
  module_complete <- module_values[complete_idx]
  
  T_complete <- droplevels(
    T_values[complete_idx]
  )
  
  T_results$N[i] <- length(module_complete)
  
  T_results$N_T1[i] <- sum(
    T_complete == "T1"
  )
  
  T_results$N_T2[i] <- sum(
    T_complete == "T2"
  )
  
  T_results$N_T3[i] <- sum(
    T_complete == "T3"
  )
  
  T_results$N_T4[i] <- sum(
    T_complete == "T4"
  )
  
  if (nlevels(T_complete) < 2) {
    T_results$Test_status[i] <- "Fewer than two observed groups"
    next
  }
  
  test <- tryCatch(
    stats::kruskal.test(
      module_complete,
      T_complete
    ),
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(test)) {
    T_results$Test_status[i] <- "Test failed"
    next
  }
  
  T_results$P_value[i] <- test$p.value
  
  T_results$Test_status[i] <- "OK"
}

T_results$FDR <- p.adjust(
  T_results$P_value,
  method = "BH"
)

T_results <- T_results[
  order(
    T_results$P_value,
    na.last = TRUE
  ),
]

rownames(T_results) <- NULL

print(
  T_results,
  row.names = FALSE
)

# Significant modules after BH correction.

T_sig <- T_results[
  !is.na(T_results$FDR) &
    T_results$FDR < 0.05,
  ,
  drop = FALSE
]

write.table(
  T_results,
  file = file.path(
    result_dir,
    "Workflow-04-2H_WGCNA_T_stage_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  T_results,
  file = file.path(
    result_dir,
    "Workflow-04-2H_WGCNA_T_stage_association.rds"
  )
)

write.table(
  T_sig,
  file = file.path(
    result_dir,
    "Workflow-04-2H_WGCNA_T_stage_significant_modules.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nNumber of T-stage-associated modules (FDR < 0.05):",
  nrow(T_sig),
  "\n"
)


# ==========================================================
# Step 8B-3B. Final QC and checkpoint
# ==========================================================

cat("\n========================================\n")
cat("Step 8B-3B. Final QC\n")
cat("========================================\n")

cat(
  "Number of biological modules:",
  ncol(MEs_biological),
  "\n"
)

cat(
  "Number of modules tested for age:",
  sum(!is.na(age_results$P_value)),
  "\n"
)

cat(
  "Number of modules tested for gender:",
  sum(!is.na(gender_results$P_value)),
  "\n"
)

cat(
  "Number of modules tested for pathological stage:",
  sum(!is.na(stage_results$P_value)),
  "\n"
)

cat(
  "Number of modules tested for T stage:",
  sum(!is.na(T_results$P_value)),
  "\n"
)

# Consolidated checkpoint for downstream analyses.

WGCNA_Step8B <- list(
  clinical_WGCNA_traits = clinical_WGCNA_traits,
  MEs_biological = MEs_biological,
  age_results = age_results,
  gender_results = gender_results,
  stage_results = stage_results,
  stage_sig = stage_sig,
  T_results = T_results,
  T_sig = T_sig,
  zero_time_audit = zero_time_audit,
  unrecognized_category_audit = unrecognized_category_audit,
  n_samples = nrow(MEs_biological),
  n_biological_modules = ncol(MEs_biological),
  random_seed = 20261106,
  R_version = R.version.string
)

saveRDS(
  WGCNA_Step8B,
  file = file.path(
    result_dir,
    "Workflow-04-2H_WGCNA_clinical_trait_association_checkpoint.rds"
  )
)

cat(
  "\nStep 8B-3B completed.\n"
)

cat(
  "Clinical trait and module association results saved to:\n",
  result_dir,
  "\n"
)

############################################################
# Step 8B-3C
# Module eigengene association with pathological N stage
#
# N0 / N1 / N2 / N3
# NX and NA excluded
# Kruskal-Wallis test
# Benjamini-Hochberg correction
############################################################

# ----------------------------------------------------------
# 1. Input validation
# ----------------------------------------------------------

if (!exists("clinical_WGCNA_traits") ||
    !exists("MEs_biological") ||
    !exists("result_dir")) {
  stop("Required objects are missing.")
}

if (!"N_trait" %in% colnames(clinical_WGCNA_traits)) {
  stop("N_trait is missing from clinical_WGCNA_traits.")
}

if (!identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)) {
  stop("Sample order mismatch between MEs_biological and clinical_WGCNA_traits.")
}

if (!is.factor(clinical_WGCNA_traits$N_trait)) {
  stop("N_trait must be a factor generated in Step 8A.")
}

expected_N_levels <- c("N0", "N1", "N2", "N3")

if (!all(levels(clinical_WGCNA_traits$N_trait) %in%
         expected_N_levels)) {
  stop("Unexpected levels detected in N_trait.")
}

if (!dir.exists(result_dir)) {
  stop("result_dir does not exist.")
}

# ----------------------------------------------------------
# 2. Check N-stage distribution
# ----------------------------------------------------------

N_distribution <- table(
  clinical_WGCNA_traits$N_trait,
  useNA = "ifany"
)

print(N_distribution)

cat(
  "\nN-stage levels:",
  paste(levels(clinical_WGCNA_traits$N_trait),
        collapse = ", "),
  "\n"
)

# ----------------------------------------------------------
# 3. Initialize results table
# ----------------------------------------------------------

N_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  FDR = NA_real_,
  N = NA_integer_,
  N_N0 = NA_integer_,
  N_N1 = NA_integer_,
  N_N2 = NA_integer_,
  N_N3 = NA_integer_,
  Test_status = "Not tested",
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 4. Kruskal-Wallis test for each biological module
# ----------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_name <- colnames(MEs_biological)[i]
  
  module_values <- MEs_biological[, i]
  N_values <- clinical_WGCNA_traits$N_trait
  
  complete_idx <- complete.cases(
    module_values,
    N_values
  ) & is.finite(module_values)
  
  module_complete <- module_values[complete_idx]
  
  N_complete <- droplevels(
    N_values[complete_idx]
  )
  
  N_counts <- table(
    factor(
      N_complete,
      levels = expected_N_levels
    )
  )
  
  N_results$N[i] <- length(module_complete)
  
  N_results$N_N0[i] <- as.integer(N_counts["N0"])
  N_results$N_N1[i] <- as.integer(N_counts["N1"])
  N_results$N_N2[i] <- as.integer(N_counts["N2"])
  N_results$N_N3[i] <- as.integer(N_counts["N3"])
  
  # Require at least two represented groups,
  # with at least two observations in each represented group.
  represented_counts <- N_counts[N_counts > 0]
  
  if (length(represented_counts) < 2) {
    N_results$Test_status[i] <-
      "Not tested: fewer than two observed N-stage groups"
    next
  }
  
  if (any(represented_counts < 2)) {
    N_results$Test_status[i] <-
      "Not tested: at least one observed group has fewer than two samples"
    next
  }
  
  if (length(unique(module_complete)) < 2) {
    N_results$Test_status[i] <-
      "Not tested: module eigengene has no variation"
    next
  }
  
  test_result <- tryCatch(
    stats::kruskal.test(
      module_complete ~ N_complete
    ),
    error = function(e) e
  )
  
  if (inherits(test_result, "error")) {
    N_results$Test_status[i] <- paste(
      "Model error:",
      conditionMessage(test_result)
    )
    next
  }
  
  N_results$P_value[i] <- test_result$p.value
  N_results$Test_status[i] <- "OK"
}

# ----------------------------------------------------------
# 5. Benjamini-Hochberg correction
#    One correction across all biological modules
# ----------------------------------------------------------

N_results$FDR <- p.adjust(
  N_results$P_value,
  method = "BH"
)

# ----------------------------------------------------------
# 6. Sort results
# ----------------------------------------------------------

N_results <- N_results[
  order(N_results$P_value, na.last = TRUE),
  ,
  drop = FALSE
]

rownames(N_results) <- NULL

# ----------------------------------------------------------
# 7. Significant modules
# ----------------------------------------------------------

N_sig <- N_results[
  !is.na(N_results$FDR) &
    N_results$FDR < 0.05,
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 8. QC summary
# ----------------------------------------------------------

cat("\n========== N-stage association QC ==========\n")

cat("Number of biological modules:",
    nrow(N_results), "\n")

cat("Original N-stage distribution:\n")
print(N_distribution)

cat("Number of successfully tested modules:",
    sum(N_results$Test_status == "OK"), "\n")

cat("Number of FDR-significant modules:",
    nrow(N_sig), "\n")

cat("\nModule-specific sample sizes:\n")

print(
  N_results[
    ,
    c(
      "Module",
      "N",
      "N_N0",
      "N_N1",
      "N_N2",
      "N_N3",
      "P_value",
      "FDR",
      "Test_status"
    )
  ],
  row.names = FALSE
)

# ----------------------------------------------------------
# 9. Save N-stage association results
# ----------------------------------------------------------

write.table(
  N_results,
  file = file.path(
    result_dir,
    "Workflow-04-2I_WGCNA_N_stage_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  N_results,
  file = file.path(
    result_dir,
    "Workflow-04-2I_WGCNA_N_stage_association.rds"
  )
)

write.table(
  N_sig,
  file = file.path(
    result_dir,
    "Workflow-04-2I_WGCNA_N_stage_significant_modules.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nN-stage association analysis completed.\n")


############################################################
# Step 8B-3D
# Module eigengene association with pathological M stage
#
# M0 vs M1
# MX and NA excluded
# Two-sided Wilcoxon rank-sum test
# Benjamini-Hochberg correction
############################################################

# ----------------------------------------------------------
# 1. Input validation
# ----------------------------------------------------------

if (!"M_trait" %in% colnames(clinical_WGCNA_traits)) {
  stop("M_trait is missing from clinical_WGCNA_traits.")
}

if (!is.factor(clinical_WGCNA_traits$M_trait)) {
  stop("M_trait must be a factor generated in Step 8A.")
}

expected_M_levels <- c("M0", "M1")

if (!all(levels(clinical_WGCNA_traits$M_trait) %in%
         expected_M_levels)) {
  stop("Unexpected levels detected in M_trait.")
}

if (!identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)) {
  stop("Sample order mismatch between MEs_biological and clinical_WGCNA_traits.")
}

# ----------------------------------------------------------
# 2. Check M-stage distribution
# ----------------------------------------------------------

M_distribution <- table(
  clinical_WGCNA_traits$M_trait,
  useNA = "ifany"
)

print(M_distribution)

cat(
  "\nM-stage levels:",
  paste(levels(clinical_WGCNA_traits$M_trait),
        collapse = ", "),
  "\n"
)

# ----------------------------------------------------------
# 3. Initialize results table
# ----------------------------------------------------------

M_results <- data.frame(
  Module = colnames(MEs_biological),
  P_value = NA_real_,
  FDR = NA_real_,
  N = NA_integer_,
  N_M0 = NA_integer_,
  N_M1 = NA_integer_,
  Median_M0 = NA_real_,
  Median_M1 = NA_real_,
  Test_status = "Not tested",
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 4. Two-sided Wilcoxon rank-sum test
# ----------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_name <- colnames(MEs_biological)[i]
  
  module_values <- MEs_biological[, i]
  M_values <- clinical_WGCNA_traits$M_trait
  
  complete_idx <- complete.cases(
    module_values,
    M_values
  ) & is.finite(module_values)
  
  module_complete <- module_values[complete_idx]
  
  M_complete <- droplevels(
    M_values[complete_idx]
  )
  
  M_counts <- table(
    factor(
      M_complete,
      levels = expected_M_levels
    )
  )
  
  M_results$N[i] <- length(module_complete)
  
  M_results$N_M0[i] <- as.integer(M_counts["M0"])
  M_results$N_M1[i] <- as.integer(M_counts["M1"])
  
  if (M_counts["M0"] > 0) {
    M_results$Median_M0[i] <- median(
      module_complete[M_complete == "M0"]
    )
  }
  
  if (M_counts["M1"] > 0) {
    M_results$Median_M1[i] <- median(
      module_complete[M_complete == "M1"]
    )
  }
  
  # Require at least two observations in each group.
  if (any(M_counts < 2)) {
    M_results$Test_status[i] <-
      "Not tested: fewer than two samples in at least one M-stage group"
    next
  }
  
  if (length(unique(module_complete)) < 2) {
    M_results$Test_status[i] <-
      "Not tested: module eigengene has no variation"
    next
  }
  
  test_result <- tryCatch(
    stats::wilcox.test(
      module_complete ~ M_complete,
      alternative = "two.sided",
      exact = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(test_result, "error")) {
    M_results$Test_status[i] <- paste(
      "Model error:",
      conditionMessage(test_result)
    )
    next
  }
  
  M_results$P_value[i] <- test_result$p.value
  M_results$Test_status[i] <- "OK"
}

# ----------------------------------------------------------
# 5. Benjamini-Hochberg correction
#    One correction across all biological modules
# ----------------------------------------------------------

M_results$FDR <- p.adjust(
  M_results$P_value,
  method = "BH"
)

# ----------------------------------------------------------
# 6. Sort results
# ----------------------------------------------------------

M_results <- M_results[
  order(M_results$P_value, na.last = TRUE),
  ,
  drop = FALSE
]

rownames(M_results) <- NULL

# ----------------------------------------------------------
# 7. Significant modules
# ----------------------------------------------------------

M_sig <- M_results[
  !is.na(M_results$FDR) &
    M_results$FDR < 0.05,
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 8. QC summary
# ----------------------------------------------------------

cat("\n========== M-stage association QC ==========\n")

cat("Number of biological modules:",
    nrow(M_results), "\n")

cat("Original M-stage distribution:\n")
print(M_distribution)

cat("Number of successfully tested modules:",
    sum(M_results$Test_status == "OK"), "\n")

cat("Number of FDR-significant modules:",
    nrow(M_sig), "\n")

cat("\nModule-specific sample sizes and results:\n")

print(
  M_results[
    ,
    c(
      "Module",
      "N",
      "N_M0",
      "N_M1",
      "Median_M0",
      "Median_M1",
      "P_value",
      "FDR",
      "Test_status"
    )
  ],
  row.names = FALSE
)

# ----------------------------------------------------------
# 9. Save M-stage association results
# ----------------------------------------------------------

write.table(
  M_results,
  file = file.path(
    result_dir,
    "Workflow-04-2J_WGCNA_M_stage_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  M_results,
  file = file.path(
    result_dir,
    "Workflow-04-2J_WGCNA_M_stage_association.rds"
  )
)

write.table(
  M_sig,
  file = file.path(
    result_dir,
    "Workflow-04-2J_WGCNA_M_stage_significant_modules.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nM-stage association analysis completed.\n")


############################################################
# Step 8C
# Module eigengene association with overall survival
#
# Univariable Cox proportional hazards regression
# Module eigengene treated as a continuous variable
# Benjamini-Hochberg correction
#
# OS_time = 0 is retained.
# Negative and non-finite OS times are excluded.
############################################################

# ----------------------------------------------------------
# 1. Input validation
# ----------------------------------------------------------

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required but is not installed.")
}

required_survival_vars <- c(
  "Sample",
  "OS_time",
  "OS_event"
)

missing_survival_vars <- setdiff(
  required_survival_vars,
  colnames(clinical_WGCNA_traits)
)

if (length(missing_survival_vars) > 0) {
  stop(
    "Missing survival variables: ",
    paste(missing_survival_vars, collapse = ", ")
  )
}

if (!identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)) {
  stop("Sample order mismatch between MEs_biological and clinical_WGCNA_traits.")
}

if (!is.numeric(clinical_WGCNA_traits$OS_time)) {
  stop("OS_time must be numeric.")
}

if (!is.numeric(clinical_WGCNA_traits$OS_event)) {
  stop("OS_event must be numeric and coded as 0/1.")
}

# ----------------------------------------------------------
# 2. Basic survival data QC
# ----------------------------------------------------------

OS_time <- clinical_WGCNA_traits$OS_time
OS_event <- clinical_WGCNA_traits$OS_event

cat("\n========== Survival data QC ==========\n")

cat("Total samples:", length(OS_time), "\n")
cat("Missing OS_time:", sum(is.na(OS_time)), "\n")
cat("Missing OS_event:", sum(is.na(OS_event)), "\n")

cat("Non-finite OS_time:",
    sum(!is.na(OS_time) & !is.finite(OS_time)), "\n")

cat("Negative OS_time:",
    sum(!is.na(OS_time) & OS_time < 0), "\n")

cat("OS_time = 0:",
    sum(!is.na(OS_time) & OS_time == 0), "\n")

cat("\nOS_event distribution:\n")
print(table(OS_event, useNA = "ifany"))

# Strict event coding check
invalid_event_values <- unique(
  OS_event[
    !is.na(OS_event) &
      !OS_event %in% c(0, 1)
  ]
)

if (length(invalid_event_values) > 0) {
  stop(
    "Invalid OS_event values detected: ",
    paste(invalid_event_values, collapse = ", ")
  )
}

# ----------------------------------------------------------
# 3. Audit OS_time = 0 observations
# ----------------------------------------------------------

zero_time_audit <- clinical_WGCNA_traits[
  !is.na(OS_time) & OS_time == 0,
  c("Sample", "OS_time", "OS_event"),
  drop = FALSE
]

cat("\nNumber of OS_time = 0 observations:",
    nrow(zero_time_audit), "\n")

print(zero_time_audit, row.names = FALSE)

write.table(
  zero_time_audit,
  file = file.path(
    result_dir,
    "Workflow-04-2K_OS_time_zero_audit.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 4. Define eligible survival observations
#
# Retain OS_time = 0.
# Exclude NA, negative, and non-finite times.
# Exclude missing event status.
# ----------------------------------------------------------

survival_eligible <- (
  !is.na(OS_time) &
    is.finite(OS_time) &
    OS_time >= 0 &
    !is.na(OS_event) &
    OS_event %in% c(0, 1)
)

cat("\nEligible survival observations:",
    sum(survival_eligible), "\n")

cat("Eligible events:",
    sum(OS_event[survival_eligible] == 1), "\n")

cat("Eligible censored observations:",
    sum(OS_event[survival_eligible] == 0), "\n")

cat("Eligible OS_time = 0 observations:",
    sum(OS_time[survival_eligible] == 0), "\n")

if (sum(survival_eligible) == 0) {
  stop("No eligible survival observations remain.")
}

if (sum(OS_event[survival_eligible] == 1) == 0) {
  stop("No death events are present in the eligible survival cohort.")
}

# ----------------------------------------------------------
# 5. Initialize Cox results table
# ----------------------------------------------------------

Cox_results <- data.frame(
  Module = colnames(MEs_biological),
  N = NA_integer_,
  Events = NA_integer_,
  Censored = NA_integer_,
  HR = NA_real_,
  CI_lower = NA_real_,
  CI_upper = NA_real_,
  P_value = NA_real_,
  FDR = NA_real_,
  Model_status = "Not tested",
  Model_warning = NA_character_,
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------
# 6. Univariable Cox regression for each module
# ----------------------------------------------------------

for (i in seq_len(ncol(MEs_biological))) {
  
  module_name <- colnames(MEs_biological)[i]
  
  module_values <- MEs_biological[, i]
  
  complete_idx <- (
    survival_eligible &
      !is.na(module_values) &
      is.finite(module_values)
  )
  
  module_complete <- module_values[complete_idx]
  
  OS_time_complete <- OS_time[complete_idx]
  OS_event_complete <- OS_event[complete_idx]
  
  n_complete <- length(module_complete)
  
  n_events <- sum(OS_event_complete == 1)
  n_censored <- sum(OS_event_complete == 0)
  
  Cox_results$N[i] <- n_complete
  Cox_results$Events[i] <- n_events
  Cox_results$Censored[i] <- n_censored
  
  # Require at least two observations and one event.
  if (n_complete < 2 || n_events < 1) {
    Cox_results$Model_status[i] <-
      "Not fitted: insufficient observations or no events"
    next
  }
  
  # Require variation in the module eigengene.
  if (length(unique(module_complete)) < 2) {
    Cox_results$Model_status[i] <-
      "Not fitted: module eigengene has no variation"
    next
  }
  
  cox_df <- data.frame(
    OS_time = OS_time_complete,
    OS_event = OS_event_complete,
    ME_value = module_complete
  )
  
  warning_messages <- character(0)
  
  cox_model <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        survival::Surv(OS_time, OS_event) ~ ME_value,
        data = cox_df,
        ties = "efron"
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (length(warning_messages) > 0) {
    Cox_results$Model_warning[i] <- paste(
      unique(warning_messages),
      collapse = " | "
    )
  }
  
  if (inherits(cox_model, "error")) {
    Cox_results$Model_status[i] <- paste(
      "Model error:",
      conditionMessage(cox_model)
    )
    next
  }
  
  model_coef <- stats::coef(cox_model)
  
  if (!"ME_value" %in% names(model_coef) ||
      !is.finite(model_coef["ME_value"])) {
    Cox_results$Model_status[i] <-
      "Not estimable: coefficient unavailable or non-finite"
    next
  }
  
  model_summary <- summary(cox_model)
  
  coef_table <- model_summary$coefficients
  
  if (!"ME_value" %in% rownames(coef_table)) {
    Cox_results$Model_status[i] <-
      "Not estimable: coefficient row unavailable"
    next
  }
  
  ci_result <- tryCatch(
    stats::confint(cox_model),
    error = function(e) e
  )
  
  if (inherits(ci_result, "error") ||
      !"ME_value" %in% rownames(ci_result)) {
    Cox_results$Model_status[i] <-
      "Not estimable: confidence interval unavailable"
    next
  }
  
  beta <- unname(model_coef["ME_value"])
  
  ci_lower_beta <- ci_result["ME_value", 1]
  ci_upper_beta <- ci_result["ME_value", 2]
  
  p_value <- coef_table["ME_value", "Pr(>|z|)"]
  
  if (!all(is.finite(c(
    beta,
    ci_lower_beta,
    ci_upper_beta,
    p_value
  )))) {
    Cox_results$Model_status[i] <-
      "Not estimable: non-finite model estimate"
    next
  }
  
  Cox_results$HR[i] <- exp(beta)
  Cox_results$CI_lower[i] <- exp(ci_lower_beta)
  Cox_results$CI_upper[i] <- exp(ci_upper_beta)
  Cox_results$P_value[i] <- p_value
  
  Cox_results$Model_status[i] <- "OK"
}

# ----------------------------------------------------------
# 7. BH correction
#    One correction across the module-level Cox tests
# ----------------------------------------------------------

Cox_results$FDR <- p.adjust(
  Cox_results$P_value,
  method = "BH"
)

# ----------------------------------------------------------
# 8. Sort results
# ----------------------------------------------------------

Cox_results <- Cox_results[
  order(Cox_results$P_value, na.last = TRUE),
  ,
  drop = FALSE
]

rownames(Cox_results) <- NULL

# ----------------------------------------------------------
# 9. Significant modules
# ----------------------------------------------------------

Cox_sig <- Cox_results[
  !is.na(Cox_results$FDR) &
    Cox_results$FDR < 0.05,
  ,
  drop = FALSE
]

# ----------------------------------------------------------
# 10. QC summary
# ----------------------------------------------------------

cat("\n========== Step 8C Cox QC ==========\n")

cat("Number of biological modules:",
    nrow(Cox_results), "\n")

cat("Number of successfully fitted models:",
    sum(Cox_results$Model_status == "OK"), "\n")

cat("Number of FDR-significant modules:",
    nrow(Cox_sig), "\n")

cat("\nModule-specific sample sizes and event counts:\n")

print(
  Cox_results[
    ,
    c(
      "Module",
      "N",
      "Events",
      "Censored",
      "HR",
      "CI_lower",
      "CI_upper",
      "P_value",
      "FDR",
      "Model_status"
    )
  ],
  row.names = FALSE
)

# ----------------------------------------------------------
# 11. Save Cox results
# ----------------------------------------------------------

write.table(
  Cox_results,
  file = file.path(
    result_dir,
    "Workflow-04-2K_WGCNA_OS_Cox_association.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  Cox_results,
  file = file.path(
    result_dir,
    "Workflow-04-2K_WGCNA_OS_Cox_association.rds"
  )
)

write.table(
  Cox_sig,
  file = file.path(
    result_dir,
    "Workflow-04-2K_WGCNA_OS_Cox_significant_modules.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C Cox analysis completed.\n")


############################################################
# Step 8C-1
# Proportional hazards assumption assessment
#
# Modules selected from Step 8C:
# FDR < 0.05
#
# Each module is assessed using its own complete-case cohort,
# consistent with its univariable Cox model in Step 8C.
############################################################

# ----------------------------------------------------------
# Step 8C-1A. Identify FDR-significant modules
# ----------------------------------------------------------

if (!exists("Cox_results")) {
  stop("Cox_results is missing. Run Step 8C first.")
}

ph_modules <- Cox_results$Module[
  !is.na(Cox_results$FDR) &
    Cox_results$FDR < 0.05 &
    Cox_results$Model_status == "OK"
]

ph_modules <- unique(ph_modules)

cat("\nModules selected for PH assessment:\n")
print(ph_modules)

cat("Number of modules selected:",
    length(ph_modules), "\n")

# ----------------------------------------------------------
# Step 8C-1B. Check module availability
# ----------------------------------------------------------

missing_ph_modules <- setdiff(
  ph_modules,
  colnames(MEs_biological)
)

if (length(missing_ph_modules) > 0) {
  stop(
    "The following PH modules are absent from MEs_biological: ",
    paste(missing_ph_modules, collapse = ", ")
  )
}

if (length(ph_modules) == 0) {
  warning(
    "No FDR-significant modules were identified in Step 8C. ",
    "PH assessment will not be performed."
  )
}

# ----------------------------------------------------------
# Step 8C-1C. Prepare common survival eligibility rule
# ----------------------------------------------------------

OS_time <- clinical_WGCNA_traits$OS_time
OS_event <- clinical_WGCNA_traits$OS_event

survival_eligible <- (
  !is.na(OS_time) &
    is.finite(OS_time) &
    OS_time >= 0 &
    !is.na(OS_event) &
    OS_event %in% c(0, 1)
)

# ----------------------------------------------------------
# Step 8C-1D. Run cox.zph() separately for each module
# ----------------------------------------------------------

ph_results <- list()
ph_summary_list <- list()

for (module_name in ph_modules) {
  
  module_values <- MEs_biological[, module_name]
  
  # Same module-specific complete-case rule as Step 8C
  complete_idx <- (
    survival_eligible &
      !is.na(module_values) &
      is.finite(module_values)
  )
  
  ph_df <- data.frame(
    OS_time = OS_time[complete_idx],
    OS_event = OS_event[complete_idx],
    ME_value = module_values[complete_idx]
  )
  
  n_ph <- nrow(ph_df)
  events_ph <- sum(ph_df$OS_event == 1)
  censored_ph <- sum(ph_df$OS_event == 0)
  
  cat(
    "\nPH assessment:",
    module_name,
    "\n"
  )
  
  cat("Complete cases:", n_ph, "\n")
  cat("Events:", events_ph, "\n")
  cat("Censored:", censored_ph, "\n")
  cat("OS_time = 0:", sum(ph_df$OS_time == 0), "\n")
  
  # Initialize a result row
  ph_row <- data.frame(
    Module = module_name,
    N = n_ph,
    Events = events_ph,
    Censored = censored_ph,
    OS_time_zero = sum(ph_df$OS_time == 0),
    PH_chisq = NA_real_,
    PH_df = NA_real_,
    PH_P = NA_real_,
    Global_chisq = NA_real_,
    Global_df = NA_real_,
    Global_P = NA_real_,
    PH_status = "Not tested",
    stringsAsFactors = FALSE
  )
  
  if (n_ph < 2 || events_ph < 1 ||
      length(unique(ph_df$ME_value)) < 2) {
    
    ph_row$PH_status <-
      "Not tested: insufficient observations, events, or eigengene variation"
    
    ph_results[[module_name]] <- list(
      module = module_name,
      data = ph_df,
      cox_model = NULL,
      zph = NULL,
      zph_table = NULL,
      status = ph_row$PH_status
    )
    
    ph_summary_list[[module_name]] <- ph_row
    next
  }
  
  # Refit the same univariable Cox model
  cox_model <- tryCatch(
    survival::coxph(
      survival::Surv(OS_time, OS_event) ~ ME_value,
      data = ph_df,
      ties = "efron"
    ),
    error = function(e) e
  )
  
  if (inherits(cox_model, "error")) {
    
    ph_row$PH_status <- paste(
      "Cox model error:",
      conditionMessage(cox_model)
    )
    
    ph_results[[module_name]] <- list(
      module = module_name,
      data = ph_df,
      cox_model = NULL,
      zph = NULL,
      zph_table = NULL,
      status = ph_row$PH_status
    )
    
    ph_summary_list[[module_name]] <- ph_row
    next
  }
  
  # Test proportional hazards assumption
  zph_test <- tryCatch(
    survival::cox.zph(cox_model),
    error = function(e) e
  )
  
  if (inherits(zph_test, "error")) {
    
    ph_row$PH_status <- paste(
      "cox.zph error:",
      conditionMessage(zph_test)
    )
    
    ph_results[[module_name]] <- list(
      module = module_name,
      data = ph_df,
      cox_model = cox_model,
      zph = NULL,
      zph_table = NULL,
      status = ph_row$PH_status
    )
    
    ph_summary_list[[module_name]] <- ph_row
    next
  }
  
  zph_table <- as.data.frame(zph_test$table)
  
  zph_table$Term <- rownames(zph_table)
  rownames(zph_table) <- NULL
  
  # Extract the module-specific row
  module_row <- zph_table[
    zph_table$Term == "ME_value",
    ,
    drop = FALSE
  ]
  
  # Extract the GLOBAL row
  global_row <- zph_table[
    zph_table$Term == "GLOBAL",
    ,
    drop = FALSE
  ]
  
  if (nrow(module_row) == 1) {
    ph_row$PH_chisq <- module_row$chisq
    ph_row$PH_df <- module_row$df
    ph_row$PH_P <- module_row$p
  }
  
  if (nrow(global_row) == 1) {
    ph_row$Global_chisq <- global_row$chisq
    ph_row$Global_df <- global_row$df
    ph_row$Global_P <- global_row$p
  }
  
  if (is.finite(ph_row$PH_P)) {
    ph_row$PH_status <- "OK"
  } else {
    ph_row$PH_status <-
      "PH test completed but module-specific result unavailable"
  }
  
  ph_results[[module_name]] <- list(
    module = module_name,
    data = ph_df,
    cox_model = cox_model,
    zph = zph_test,
    zph_table = zph_table,
    status = ph_row$PH_status
  )
  
  ph_summary_list[[module_name]] <- ph_row
}

# ----------------------------------------------------------
# Step 8C-1E. Combine PH results
# ----------------------------------------------------------

if (length(ph_summary_list) > 0) {
  
  ph_summary <- do.call(
    rbind,
    ph_summary_list
  )
  
  rownames(ph_summary) <- NULL
  
} else {
  
  ph_summary <- data.frame(
    Module = character(),
    N = integer(),
    Events = integer(),
    Censored = integer(),
    OS_time_zero = integer(),
    PH_chisq = numeric(),
    PH_df = numeric(),
    PH_P = numeric(),
    Global_chisq = numeric(),
    Global_df = numeric(),
    Global_P = numeric(),
    PH_status = character(),
    stringsAsFactors = FALSE
  )
}

# ----------------------------------------------------------
# Step 8C-1F. Interpret module-specific PH results
# ----------------------------------------------------------

ph_summary$PH_interpretation <- ifelse(
  is.na(ph_summary$PH_P),
  "Not assessed / unavailable",
  ifelse(
    ph_summary$PH_P < 0.05,
    "Potential PH violation",
    "No evidence of PH violation"
  )
)

# Keep GLOBAL test as supplementary information.
# The module-specific PH_P is the primary result
# for a univariable model containing one module.

print(ph_summary)

# ----------------------------------------------------------
# Step 8C-1G. Save PH assessment results
# ----------------------------------------------------------

write.table(
  ph_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2L_WGCNA_OS_PH_assessment.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(
  ph_results,
  file = file.path(
    result_dir,
    "Workflow-04-2L_WGCNA_OS_PH_models.rds"
  )
)

saveRDS(
  ph_summary,
  file = file.path(
    result_dir,
    "Workflow-04-2L_WGCNA_OS_PH_summary.rds"
  )
)

cat("\nPH assessment completed.\n")

############################################################
# Step 8C-1H
# Visual inspection of PH assumption
############################################################

ph_plot_file <- file.path(
  result_dir,
  "Workflow-04-2L_WGCNA_OS_PH_diagnostic_plots.pdf"
)

pdf(
  file = ph_plot_file,
  width = 7,
  height = 7,
  onefile = TRUE
)

for (module_name in names(ph_results)) {
  
  zph_object <- ph_results[[module_name]]$zph
  
  if (is.null(zph_object)) {
    next
  }
  
  plot(
    zph_object,
    var = 1,
    main = paste(
      "Proportional hazards assessment:",
      module_name
    )
  )
  
  abline(
    h = 0,
    lty = 2
  )
}

dev.off()

cat(
  "PH diagnostic plots saved to:\n",
  ph_plot_file,
  "\n"
)

############################################################
# Paper 1 LUAD tumor-only WGCNA
#
# Script:
# GitHub2_Paper1_tumor-only_WGCNA_1.R
#
# Analysis section:
# Step 8C-2
# Clinical covariate design and multivariable Cox analysis
#
# Research objective:
# To evaluate whether selected WGCNA module eigengenes
# are associated with overall survival after adjustment
# for prespecified clinical covariates.
#
# Input data:
# 1. clinical_aligned
#    Reconstructed and QC-checked TCGA-LUAD clinical data.
#
# 2. clinical_WGCNA_traits
#    Clinical variables generated in Step 8A.
#
# 3. MEs_biological
#    Biological module eigengenes generated in Step 8B.
#
# 4. Cox_results
#    Univariable module-level Cox results from Step 8C.
#
# 5. ph_summary
#    Module-specific PH assessment results from Step 8C-1.
#
# Input sample set:
# 513 primary LUAD tumor samples.
#
# Survival eligibility:
# - OS_time must be finite and >= 0.
# - OS_event must be coded as 0 or 1.
# - OS_time = 0 observations are retained.
# - Missing clinical covariates are handled by
#   model-specific complete-case analysis.
#
# Primary clinical adjustment:
# Age + Gender + pathological Stage.
#
# Sensitivity clinical adjustment:
# Age + Gender + T stage + N stage + M stage.
#
# Module selection:
# Based on Step 8C univariable Cox results and
# Step 8C-1 proportional hazards assessment.
#
# Statistical methods:
# - Multivariable Cox proportional hazards regression.
# - Schoenfeld residual-based PH assessment.
#
# Multiple testing:
# No additional BH correction is applied to the
# individual covariate coefficients in these models.
# Module-level selection was based on Step 8C.
#
# Software:
# R 4.5.1
# survival package (version recorded at runtime).
#
# Random seed:
# 20261106
# (No stochastic procedure is used in the Cox models.)
#
# Output:
# Project-relative result_dir.
#
# Reproducibility:
# - No setwd().
# - No absolute paths.
# - No save.image().
# - Explicit object checkpoints.
############################################################

# ----------------------------------------------------------
# 0. Software environment and setup
# ----------------------------------------------------------

set.seed(20261106)

if (!requireNamespace("survival", quietly = TRUE)) {
  stop("Package 'survival' is required.")
}

if (!requireNamespace("here", quietly = TRUE)) {
  stop("Package 'here' is required.")
}

if (!exists("result_dir") || !dir.exists(result_dir)) {
  stop("result_dir is missing or does not exist.")
}

cat("\n========== Software environment ==========\n")

cat("R version:",
    as.character(getRversion()), "\n")

cat("survival version:",
    as.character(utils::packageVersion("survival")), "\n")

cat("here version:",
    as.character(utils::packageVersion("here")), "\n")

cat("Random seed: 20261106\n")

cat("Output directory:", result_dir, "\n")

cat("==========================================\n")


############################################################
# Step 8C-2A
# Clinical covariate design and sample alignment
#
# No multivariable Cox model is fitted at this step.
############################################################

# ----------------------------------------------------------
# 1. Check required input objects
# ----------------------------------------------------------

required_objects <- c(
  "clinical_WGCNA_traits",
  "MEs_biological",
  "Cox_results",
  "ph_summary",
  "result_dir"
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
    "Missing objects: ",
    paste(missing_objects, collapse = ", ")
  )
}

required_clinical_vars <- c(
  "Sample",
  "Age_trait",
  "Gender_trait",
  "Stage_trait",
  "T_trait",
  "N_trait",
  "M_trait",
  "OS_time",
  "OS_event"
)

missing_clinical_vars <- setdiff(
  required_clinical_vars,
  colnames(clinical_WGCNA_traits)
)

if (length(missing_clinical_vars) > 0) {
  stop(
    "Missing clinical variables: ",
    paste(missing_clinical_vars, collapse = ", ")
  )
}

# ----------------------------------------------------------
# 2. Validate sample IDs and sample order
# ----------------------------------------------------------

if (anyNA(clinical_WGCNA_traits$Sample) ||
    anyDuplicated(clinical_WGCNA_traits$Sample)) {
  stop("Clinical Sample IDs contain missing or duplicated values.")
}

if (is.null(rownames(MEs_biological)) ||
    anyNA(rownames(MEs_biological)) ||
    anyDuplicated(rownames(MEs_biological))) {
  stop("MEs_biological row names are missing or duplicated.")
}

if (!identical(
  rownames(MEs_biological),
  clinical_WGCNA_traits$Sample
)) {
  stop(
    "Sample order mismatch between MEs_biological ",
    "and clinical_WGCNA_traits."
  )
}

if (nrow(MEs_biological) != nrow(clinical_WGCNA_traits)) {
  stop("Clinical and eigengene matrices have different sample sizes.")
}

# ----------------------------------------------------------
# 3. Select clinical variables
# ----------------------------------------------------------

clinical_cox <- clinical_WGCNA_traits[
  ,
  required_clinical_vars,
  drop = FALSE
]

# ----------------------------------------------------------
# 4. Validate Age and survival variables
# ----------------------------------------------------------

if (!is.numeric(clinical_cox$Age_trait)) {
  stop("Age_trait must be numeric. Do not use factor integer codes.")
}

if (!is.numeric(clinical_cox$OS_time)) {
  stop("OS_time must be numeric.")
}

if (!is.numeric(clinical_cox$OS_event)) {
  stop("OS_event must be numeric and coded as 0/1.")
}

invalid_event_values <- unique(
  clinical_cox$OS_event[
    !is.na(clinical_cox$OS_event) &
      !clinical_cox$OS_event %in% c(0, 1)
  ]
)

if (length(invalid_event_values) > 0) {
  stop(
    "Unexpected OS_event values: ",
    paste(invalid_event_values, collapse = ", ")
  )
}

if (any(
  !is.na(clinical_cox$Age_trait) &
  !is.finite(clinical_cox$Age_trait)
)) {
  stop("Age_trait contains non-finite values.")
}

# ----------------------------------------------------------
# 5. Validate and standardize categorical variables
#
# These operations preserve the Step 8A clinical coding.
# They do not reclassify unknown values.
# ----------------------------------------------------------

# Gender
# ----------------------------------------------------------
# Gender: accept only 0/1 or Female/Male coding
# ----------------------------------------------------------

gender_values <- as.character(
  clinical_cox$Gender_trait
)

observed_gender <- unique(
  gender_values[!is.na(gender_values)]
)

if (all(observed_gender %in% c("0", "1"))) {
  
  clinical_cox$Gender_trait <- factor(
    gender_values,
    levels = c("0", "1"),
    labels = c("Female", "Male")
  )
  
} else if (all(observed_gender %in% c("Female", "Male"))) {
  
  clinical_cox$Gender_trait <- factor(
    gender_values,
    levels = c("Female", "Male")
  )
  
} else {
  
  stop(
    "Unexpected Gender_trait coding: ",
    paste(observed_gender, collapse = ", ")
  )
}

# ----------------------------------------------------------
# Validate categorical clinical labels before factor conversion
# ----------------------------------------------------------

allowed_categories <- list(
  
  Stage_trait = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  ),
  
  T_trait = c(
    "T1", "T2", "T3", "T4"
  ),
  
  N_trait = c(
    "N0", "N1", "N2", "N3"
  ),
  
  M_trait = c(
    "M0", "M1"
  )
)

for (variable_name in names(allowed_categories)) {
  
  observed_values <- unique(
    as.character(
      clinical_cox[[variable_name]]
    )
  )
  
  observed_values <- observed_values[
    !is.na(observed_values)
  ]
  
  unexpected_values <- setdiff(
    observed_values,
    allowed_categories[[variable_name]]
  )
  
  if (length(unexpected_values) > 0) {
    
    stop(
      "Unexpected categories in ",
      variable_name,
      ": ",
      paste(unexpected_values, collapse = ", ")
    )
  }
}

# Stage
clinical_cox$Stage_trait <- factor(
  as.character(clinical_cox$Stage_trait),
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )
)

# T stage
clinical_cox$T_trait <- factor(
  as.character(clinical_cox$T_trait),
  levels = c("T1", "T2", "T3", "T4")
)

# N stage
clinical_cox$N_trait <- factor(
  as.character(clinical_cox$N_trait),
  levels = c("N0", "N1", "N2", "N3")
)

# M stage
clinical_cox$M_trait <- factor(
  as.character(clinical_cox$M_trait),
  levels = c("M0", "M1")
)

# ----------------------------------------------------------
# 6. Confirm categorical coding
# ----------------------------------------------------------

expected_factor_levels <- list(
  Gender_trait = c("Female", "Male"),
  Stage_trait = c(
    "Stage I", "Stage II", "Stage III", "Stage IV"
  ),
  T_trait = c("T1", "T2", "T3", "T4"),
  N_trait = c("N0", "N1", "N2", "N3"),
  M_trait = c("M0", "M1")
)

for (variable_name in names(expected_factor_levels)) {
  
  observed_levels <- levels(
    clinical_cox[[variable_name]]
  )
  
  if (!identical(
    observed_levels,
    expected_factor_levels[[variable_name]]
  )) {
    stop(
      "Factor level definition failed for ",
      variable_name
    )
  }
}

# ----------------------------------------------------------
# 7. Check module selection based on Step 8C and Step 8C-1
# ----------------------------------------------------------

if (!all(
  c("Module", "FDR", "Model_status") %in%
  colnames(Cox_results)
)) {
  stop("Cox_results does not contain the required columns.")
}

if (!all(
  c("Module", "PH_P", "PH_status") %in%
  colnames(ph_summary)
)) {
  stop("ph_summary does not contain the required columns.")
}

selected_modules <- Cox_results$Module[
  !is.na(Cox_results$FDR) &
    Cox_results$FDR < 0.05 &
    Cox_results$Model_status == "OK"
]

selected_modules <- unique(selected_modules)

if (length(selected_modules) == 0) {
  stop(
    "No FDR-significant and successfully fitted modules ",
    "were identified in Step 8C."
  )
}

missing_ph_results <- setdiff(
  selected_modules,
  ph_summary$Module
)

if (length(missing_ph_results) > 0) {
  stop(
    "PH assessment is missing for selected modules: ",
    paste(missing_ph_results, collapse = ", ")
  )
}

if (!all(selected_modules %in% colnames(MEs_biological))) {
  stop("Some selected modules are missing from MEs_biological.")
}

# ----------------------------------------------------------
# 8. Construct aligned Cox design data
# ----------------------------------------------------------

module_data <- as.data.frame(
  MEs_biological[
    ,
    selected_modules,
    drop = FALSE
  ],
  check.names = FALSE
)

cox_design <- cbind(
  clinical_cox,
  module_data
)

if (!identical(
  cox_design$Sample,
  rownames(MEs_biological)
)) {
  stop("Sample order changed during Cox design construction.")
}

# ----------------------------------------------------------
# 9. QC summary
# ----------------------------------------------------------

cat("\n========== Step 8C-2A QC ==========\n")

cat("Clinical samples:", nrow(clinical_cox), "\n")
cat("Biological modules:", ncol(MEs_biological), "\n")
cat("Selected modules:", length(selected_modules), "\n")

cat("\nSelected modules:\n")
print(selected_modules)

cat("\nGender distribution:\n")
print(table(cox_design$Gender_trait, useNA = "ifany"))

cat("\nStage distribution:\n")
print(table(cox_design$Stage_trait, useNA = "ifany"))

cat("\nT-stage distribution:\n")
print(table(cox_design$T_trait, useNA = "ifany"))

cat("\nN-stage distribution:\n")
print(table(cox_design$N_trait, useNA = "ifany"))

cat("\nM-stage distribution:\n")
print(table(cox_design$M_trait, useNA = "ifany"))

cat("\n===================================\n")


############################################################
# Step 8C-2B
# Clinical covariate overlap assessment
#
# Descriptive QC only.
# No statistical model is fitted.
############################################################

cat("\n========== Clinical Covariate Overlap ==========\n")

clinical_cross_tables <- list(
  Stage_vs_T = table(
    cox_design$Stage_trait,
    cox_design$T_trait,
    useNA = "ifany"
  ),
  
  Stage_vs_N = table(
    cox_design$Stage_trait,
    cox_design$N_trait,
    useNA = "ifany"
  ),
  
  Stage_vs_M = table(
    cox_design$Stage_trait,
    cox_design$M_trait,
    useNA = "ifany"
  ),
  
  T_vs_N = table(
    cox_design$T_trait,
    cox_design$N_trait,
    useNA = "ifany"
  ),
  
  T_vs_M = table(
    cox_design$T_trait,
    cox_design$M_trait,
    useNA = "ifany"
  ),
  
  N_vs_M = table(
    cox_design$N_trait,
    cox_design$M_trait,
    useNA = "ifany"
  )
)

for (table_name in names(clinical_cross_tables)) {
  
  cat("\n---", table_name, "---\n")
  
  print(clinical_cross_tables[[table_name]])
}

# ----------------------------------------------------------
# Save descriptive overlap audit
# ----------------------------------------------------------

saveRDS(
  clinical_cross_tables,
  file = file.path(
    result_dir,
    "Workflow-04-2M_WGCNA_clinical_covariate_overlap.rds"
  )
)

cat("\nClinical covariate overlap assessment completed.\n")

############################################################
# Step 8C-2C
# Age and Gender descriptive QC
############################################################

cat("\n========== Age and Gender QC ==========\n")

cat("\nAge summary:\n")
print(summary(cox_design$Age_trait))

cat("\nGender distribution:\n")
print(table(cox_design$Gender_trait, useNA = "ifany"))

cat("\nAge by Gender:\n")

age_by_gender <- split(
  cox_design$Age_trait,
  cox_design$Gender_trait
)

age_by_gender_summary <- lapply(
  age_by_gender,
  summary
)

print(age_by_gender_summary)

cat("\n========================================\n")

############################################################
# Step 8C-2D
# Define clinical adjustment frameworks
#
# Model A: Age + Gender + Stage
# Model B: Age + Gender + T + N + M
############################################################

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

primary_model_label <- paste(
  primary_covariates,
  collapse = " + "
)

sensitivity_model_label <- paste(
  sensitivity_covariates,
  collapse = " + "
)

cat("\nPrimary clinical adjustment:\n")
cat(primary_model_label, "\n")

cat("\nSensitivity clinical adjustment:\n")
cat(sensitivity_model_label, "\n")


############################################################
# Step 8C-2E
# Classify modules according to PH assessment
############################################################

module_ph_status <- ph_summary[
  match(selected_modules, ph_summary$Module),
  ,
  drop = FALSE
]

if (!identical(
  module_ph_status$Module,
  selected_modules
)) {
  stop("Module matching failed during PH classification.")
}

module_ph_status$PH_class <- ifelse(
  is.na(module_ph_status$PH_P),
  "PH assessment unavailable",
  ifelse(
    module_ph_status$PH_P < 0.05,
    "Potential PH violation",
    "No evidence of PH violation"
  )
)

module_ph_status$Recommended_model <- ifelse(
  module_ph_status$PH_class ==
    "No evidence of PH violation",
  "Conventional Cox candidate",
  ifelse(
    module_ph_status$PH_class ==
      "Potential PH violation",
    "Time-dependent / sensitivity candidate",
    "Further assessment required"
  )
)

conventional_modules <- module_ph_status$Module[
  module_ph_status$PH_class ==
    "No evidence of PH violation"
]

time_dependent_modules <- module_ph_status$Module[
  module_ph_status$PH_class ==
    "Potential PH violation"
]

ph_unavailable_modules <- module_ph_status$Module[
  module_ph_status$PH_class ==
    "PH assessment unavailable"
]

cat("\n========== Module PH classification ==========\n")

print(
  module_ph_status[
    ,
    c(
      "Module",
      "PH_P",
      "PH_class",
      "Recommended_model"
    )
  ],
  row.names = FALSE
)

cat("\nConventional Cox candidates:\n")
print(conventional_modules)

cat("\nTime-dependent / sensitivity candidates:\n")
print(time_dependent_modules)

cat("\nModules requiring further PH assessment:\n")
print(ph_unavailable_modules)

cat("\n==============================================\n")


############################################################
# Step 8C-2F
# Save model design
#
# No multivariable Cox model is fitted here.
############################################################

model_design <- list(
  
  analysis = "Tumor-only WGCNA module survival analysis",
  
  input_sample_n = nrow(cox_design),
  
  selected_modules = selected_modules,
  
  conventional_modules = conventional_modules,
  
  time_dependent_modules = time_dependent_modules,
  
  ph_unavailable_modules = ph_unavailable_modules,
  
  module_ph_status = module_ph_status,
  
  primary_clinical_model = primary_model_label,
  
  sensitivity_clinical_model = sensitivity_model_label,
  
  survival_eligibility = paste(
    "OS_time finite and >= 0;",
    "OS_event coded as 0 or 1;",
    "OS_time = 0 retained."
  ),
  
  primary_model_note = paste(
    "Conventional multivariable Cox models are",
    "considered for modules without evidence of",
    "PH violation."
  ),
  
  time_dependent_note = paste(
    "Potential PH violations are classified as",
    "time-dependent or sensitivity candidates.",
    "No time-dependent Cox model is fitted in this step."
  ),
  
  multiple_testing_note = paste(
    "Module selection is based on Step 8C.",
    "No additional BH correction is performed",
    "for covariate coefficients in the multivariable models."
  ),
  
  software = list(
    R = as.character(getRversion()),
    survival = as.character(
      utils::packageVersion("survival")
    ),
    here = as.character(
      utils::packageVersion("here")
    )
  ),
  
  random_seed = 20261106,
  
  session_info = utils::sessionInfo(),
  
  note = paste(
    "This checkpoint records the model design only.",
    "No final multivariable Cox model is fitted here."
  )
)

saveRDS(
  model_design,
  file = file.path(
    result_dir,
    "Workflow-04-2M_WGCNA_OS_multivariable_model_design.rds"
  )
)

write.table(
  module_ph_status,
  file = file.path(
    result_dir,
    "Workflow-04-2M_WGCNA_OS_multivariable_model_design.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2 model design saved.\n")

############################################################
# Step 8C-2G
# Define module-specific complete-case cohorts
#
# No Cox model is fitted at this step.
############################################################

# ----------------------------------------------------------
# 1. Define survival eligibility
# ----------------------------------------------------------

survival_eligible <- (
  !is.na(cox_design$OS_time) &
    is.finite(cox_design$OS_time) &
    cox_design$OS_time >= 0 &
    !is.na(cox_design$OS_event) &
    cox_design$OS_event %in% c(0, 1)
)

cat("\nSurvival-eligible samples:",
    sum(survival_eligible), "\n")

cat("Eligible events:",
    sum(cox_design$OS_event[survival_eligible] == 1),
    "\n")

cat("Eligible censored:",
    sum(cox_design$OS_event[survival_eligible] == 0),
    "\n")

cat("Eligible OS.time = 0:",
    sum(cox_design$OS_time[survival_eligible] == 0),
    "\n")


# ----------------------------------------------------------
# 2. Define complete-case cohort function
# ----------------------------------------------------------

get_complete_cohort <- function(
    data,
    module,
    covariates,
    survival_eligible
) {
  
  required_vars <- c(
    "Sample",
    "OS_time",
    "OS_event",
    module,
    covariates
  )
  
  missing_vars <- setdiff(
    required_vars,
    colnames(data)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      "Missing variables for ",
      module,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  if (length(survival_eligible) != nrow(data)) {
    stop("survival_eligible has an incorrect length.")
  }
  
  complete_idx <- (
    survival_eligible &
      complete.cases(
        data[, c(module, covariates), drop = FALSE]
      ) &
      is.finite(data[[module]])
  )
  
  cohort <- data[
    complete_idx,
    required_vars,
    drop = FALSE
  ]
  
  # Remove unused factor levels in the model cohort.
  factor_vars <- covariates[
    vapply(
      cohort[, covariates, drop = FALSE],
      is.factor,
      logical(1)
    )
  ]
  
  for (variable_name in factor_vars) {
    cohort[[variable_name]] <- droplevels(
      cohort[[variable_name]]
    )
  }
  
  rownames(cohort) <- NULL
  
  cohort
}


# ----------------------------------------------------------
# 3. Define module-specific primary cohorts
# ----------------------------------------------------------

primary_cohorts <- lapply(
  selected_modules,
  function(module) {
    
    get_complete_cohort(
      data = cox_design,
      module = module,
      covariates = primary_covariates,
      survival_eligible = survival_eligible
    )
  }
)

names(primary_cohorts) <- selected_modules


# ----------------------------------------------------------
# 4. Define module-specific sensitivity cohorts
# ----------------------------------------------------------

sensitivity_cohorts <- lapply(
  selected_modules,
  function(module) {
    
    get_complete_cohort(
      data = cox_design,
      module = module,
      covariates = sensitivity_covariates,
      survival_eligible = survival_eligible
    )
  }
)

names(sensitivity_cohorts) <- selected_modules


############################################################
# Step 8C-2H
# Final QC of module-specific complete-case cohorts
############################################################

# ----------------------------------------------------------
# 1. Summarize one cohort
# ----------------------------------------------------------

summarize_cohort <- function(
    data,
    module,
    model_type
) {
  
  result <- data.frame(
    Module = module,
    Model_type = model_type,
    N = nrow(data),
    Events = sum(data$OS_event == 1),
    Censored = sum(data$OS_event == 0),
    ZeroTime = sum(data$OS_time == 0),
    stringsAsFactors = FALSE
  )
  
  if (model_type == "Primary") {
    
    result$Age_N <- sum(!is.na(data$Age_trait))
    result$Gender_N <- sum(!is.na(data$Gender_trait))
    result$Stage_N <- sum(!is.na(data$Stage_trait))
    
  } else {
    
    result$Age_N <- sum(!is.na(data$Age_trait))
    result$Gender_N <- sum(!is.na(data$Gender_trait))
    result$T_N <- sum(!is.na(data$T_trait))
    result$N_N <- sum(!is.na(data$N_trait))
    result$M_N <- sum(!is.na(data$M_trait))
  }
  
  result
}


# ----------------------------------------------------------
# 2. Primary cohort QC
# ----------------------------------------------------------

primary_qc <- do.call(
  rbind,
  lapply(
    selected_modules,
    function(module) {
      
      summarize_cohort(
        data = primary_cohorts[[module]],
        module = module,
        model_type = "Primary"
      )
    }
  )
)

rownames(primary_qc) <- NULL


# ----------------------------------------------------------
# 3. Sensitivity cohort QC
# ----------------------------------------------------------

sensitivity_qc <- do.call(
  rbind,
  lapply(
    selected_modules,
    function(module) {
      
      summarize_cohort(
        data = sensitivity_cohorts[[module]],
        module = module,
        model_type = "Sensitivity"
      )
    }
  )
)

rownames(sensitivity_qc) <- NULL


# ----------------------------------------------------------
# 4. Display QC results
# ----------------------------------------------------------

cat("\n========== Primary cohort QC ==========\n")

print(primary_qc, row.names = FALSE)

cat("\n========== Sensitivity cohort QC ==========\n")

print(sensitivity_qc, row.names = FALSE)


# ----------------------------------------------------------
# 5. Detailed categorical distribution by module
# ----------------------------------------------------------

primary_factor_distributions <- lapply(
  selected_modules,
  function(module) {
    
    dat <- primary_cohorts[[module]]
    
    list(
      Module = module,
      Gender = table(dat$Gender_trait),
      Stage = table(dat$Stage_trait)
    )
  }
)

names(primary_factor_distributions) <- selected_modules


sensitivity_factor_distributions <- lapply(
  selected_modules,
  function(module) {
    
    dat <- sensitivity_cohorts[[module]]
    
    list(
      Module = module,
      Gender = table(dat$Gender_trait),
      T = table(dat$T_trait),
      N = table(dat$N_trait),
      M = table(dat$M_trait)
    )
  }
)

names(sensitivity_factor_distributions) <- selected_modules


# ----------------------------------------------------------
# 6. Save cohort QC
# ----------------------------------------------------------

saveRDS(
  list(
    survival_eligible = survival_eligible,
    primary_qc = primary_qc,
    sensitivity_qc = sensitivity_qc,
    primary_cohorts = primary_cohorts,
    sensitivity_cohorts = sensitivity_cohorts,
    primary_factor_distributions =
      primary_factor_distributions,
    sensitivity_factor_distributions =
      sensitivity_factor_distributions
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds"
  )
)

write.table(
  primary_qc,
  file = file.path(
    result_dir,
    "Workflow-04-2N_WGCNA_OS_primary_cohort_QC.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sensitivity_qc,
  file = file.path(
    result_dir,
    "Workflow-04-2N_WGCNA_OS_sensitivity_cohort_QC.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2H cohort QC completed.\n")


############################################################
# Step 8C-2I
# Primary conventional multivariable Cox models
#
# Model:
# Module eigengene + Age + Gender + Stage
#
# Only modules without evidence of PH violation
# are considered for conventional Cox modelling.
############################################################

# ----------------------------------------------------------
# 1. Check conventional modules
# ----------------------------------------------------------

if (length(conventional_modules) == 0) {
  
  warning(
    "No modules were classified as conventional Cox candidates. ",
    "No primary multivariable Cox models will be fitted."
  )
  
}

missing_primary_cohorts <- setdiff(
  conventional_modules,
  names(primary_cohorts)
)

if (length(missing_primary_cohorts) > 0) {
  stop(
    "Primary cohorts are missing for: ",
    paste(missing_primary_cohorts, collapse = ", ")
  )
}


# ----------------------------------------------------------
# 2. Define model fitting function
# ----------------------------------------------------------

fit_primary_cox <- function(module, data) {
  
  if (nrow(data) == 0) {
    return(list(
      model = NULL,
      status = "Not fitted: empty complete-case cohort",
      warnings = character(0)
    ))
  }
  
  # Rename module eigengene to a stable formula term.
  model_data <- data.frame(
    OS_time = data$OS_time,
    OS_event = data$OS_event,
    ME_value = data[[module]],
    Age_trait = data$Age_trait,
    Gender_trait = droplevels(data$Gender_trait),
    Stage_trait = droplevels(data$Stage_trait)
  )
  
  # Check sufficient observations and events.
  if (sum(model_data$OS_event == 1) < 1) {
    return(list(
      model = NULL,
      status = "Not fitted: no death events",
      warnings = character(0)
    ))
  }
  
  # Check variation in continuous variables.
  if (length(unique(model_data$ME_value)) < 2 ||
      length(unique(model_data$Age_trait)) < 2) {
    return(list(
      model = NULL,
      status = "Not fitted: insufficient variation in continuous covariates",
      warnings = character(0)
    ))
  }
  
  # Check categorical variables have at least two observed levels.
  if (nlevels(model_data$Gender_trait) < 2 ||
      nlevels(model_data$Stage_trait) < 2) {
    return(list(
      model = NULL,
      status = "Not fitted: insufficient observed categorical levels",
      warnings = character(0)
    ))
  }
  
  # Explicit model formula.
  model_formula <- survival::Surv(
    OS_time,
    OS_event
  ) ~ ME_value + Age_trait + Gender_trait + Stage_trait
  
  warning_messages <- character(0)
  
  model <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        model_formula,
        data = model_data,
        ties = "efron",
        x = TRUE,
        y = TRUE,
        model = TRUE
      ),
      warning = function(w) {
        warning_messages <<- c(
          warning_messages,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  
  if (inherits(model, "error")) {
    return(list(
      model = NULL,
      status = paste(
        "Model error:",
        conditionMessage(model)
      ),
      warnings = warning_messages
    ))
  }
  
  if (anyNA(stats::coef(model))) {
    return(list(
      model = model,
      status = "Fitted with non-estimable coefficient(s)",
      warnings = warning_messages
    ))
  }
  
  list(
    model = model,
    status = "OK",
    warnings = warning_messages
  )
}


# ----------------------------------------------------------
# 3. Fit primary models
# ----------------------------------------------------------

primary_model_objects <- lapply(
  conventional_modules,
  function(module) {
    
    fit_primary_cox(
      module = module,
      data = primary_cohorts[[module]]
    )
  }
)

names(primary_model_objects) <- conventional_modules

primary_models <- lapply(
  primary_model_objects,
  function(x) x$model
)


# ----------------------------------------------------------
# 4. Extract coefficient-level results
# ----------------------------------------------------------

extract_primary_cox_results <- function(
    model,
    module,
    model_status
) {
  
  if (is.null(model)) {
    return(NULL)
  }
  
  coefficient_table <- summary(model)$coefficients
  
  confidence_table <- summary(model)$conf.int
  
  coefficient_names <- rownames(coefficient_table)
  
  result <- data.frame(
    Module = module,
    Variable = coefficient_names,
    Beta = coefficient_table[, "coef"],
    HR = confidence_table[, "exp(coef)"],
    Lower95CI = confidence_table[, "lower .95"],
    Upper95CI = confidence_table[, "upper .95"],
    P_value = coefficient_table[, "Pr(>|z|)"],
    Model_status = model_status,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  result
}


primary_result_list <- lapply(
  conventional_modules,
  function(module) {
    
    extract_primary_cox_results(
      model = primary_models[[module]],
      module = module,
      model_status = primary_model_objects[[module]]$status
    )
  }
)

primary_result_list <- Filter(
  Negate(is.null),
  primary_result_list
)

if (length(primary_result_list) > 0) {
  
  primary_results <- do.call(
    rbind,
    primary_result_list
  )
  
  rownames(primary_results) <- NULL
  
} else {
  
  primary_results <- data.frame(
    Module = character(),
    Variable = character(),
    Beta = numeric(),
    HR = numeric(),
    Lower95CI = numeric(),
    Upper95CI = numeric(),
    P_value = numeric(),
    Model_status = character(),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------
# 5. Model-level QC
# ----------------------------------------------------------

primary_model_qc <- do.call(
  rbind,
  lapply(
    conventional_modules,
    function(module) {
      
      model <- primary_models[[module]]
      
      dat <- primary_cohorts[[module]]
      
      model_status <- primary_model_objects[[module]]$status
      
      if (is.null(model)) {
        
        return(data.frame(
          Module = module,
          N = nrow(dat),
          Events = sum(dat$OS_event == 1),
          Censored = sum(dat$OS_event == 0),
          ZeroTime = sum(dat$OS_time == 0),
          Concordance = NA_real_,
          Concordance_SE = NA_real_,
          Model_status = model_status,
          Model_warning = paste(
            primary_model_objects[[module]]$warnings,
            collapse = " | "
          ),
          stringsAsFactors = FALSE
        ))
      }
      
      model_summary <- summary(model)
      
      data.frame(
        Module = module,
        N = model$n,
        Events = model$nevent,
        Censored = model$n - model$nevent,
        ZeroTime = sum(dat$OS_time == 0),
        Concordance = unname(
          model_summary$concordance[1]
        ),
        Concordance_SE = unname(
          model_summary$concordance[2]
        ),
        Model_status = model_status,
        Model_warning = paste(
          primary_model_objects[[module]]$warnings,
          collapse = " | "
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(primary_model_qc) <- NULL


# ----------------------------------------------------------
# 6. Display results
# ----------------------------------------------------------

cat("\n========== Primary multivariable Cox results ==========\n")

print(primary_results, row.names = FALSE)

cat("\n========== Primary model QC ==========\n")

print(primary_model_qc, row.names = FALSE)


# ----------------------------------------------------------
# 7. Save primary Cox results
# ----------------------------------------------------------

saveRDS(
  list(
    models = primary_models,
    model_objects = primary_model_objects,
    results = primary_results,
    model_qc = primary_model_qc,
    conventional_modules = conventional_modules,
    covariates = primary_covariates
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2O_WGCNA_OS_primary_conventional_Cox.rds"
  )
)

write.table(
  primary_results,
  file = file.path(
    result_dir,
    "Workflow-04-2O_WGCNA_OS_primary_conventional_Cox.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_model_qc,
  file = file.path(
    result_dir,
    "Workflow-04-2O_WGCNA_OS_primary_conventional_Cox_QC.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2I completed.\n")


############################################################
# Step 8C-2J
# PH assumption assessment for primary multivariable Cox
#
# Schoenfeld residual-based tests
############################################################

# ----------------------------------------------------------
# 1. Identify successfully fitted models
# ----------------------------------------------------------

successful_primary_modules <- names(primary_models)[
  vapply(
    primary_models,
    function(x) !is.null(x),
    logical(1)
  )
]

cat("\nModels available for multivariable PH assessment:\n")
print(successful_primary_modules)


# ----------------------------------------------------------
# 2. Run cox.zph for each primary model
# ----------------------------------------------------------

primary_zph_objects <- list()
primary_ph_result_list <- list()

for (module in successful_primary_modules) {
  
  model <- primary_models[[module]]
  
  zph_object <- tryCatch(
    survival::cox.zph(
      model,
      transform = "km"
    ),
    error = function(e) e
  )
  
  if (inherits(zph_object, "error")) {
    
    primary_ph_result_list[[module]] <- data.frame(
      Module = module,
      Variable = NA_character_,
      ChiSquare = NA_real_,
      DF = NA_real_,
      P_value = NA_real_,
      PH_status = paste(
        "cox.zph error:",
        conditionMessage(zph_object)
      ),
      stringsAsFactors = FALSE
    )
    
    primary_zph_objects[[module]] <- NULL
    
    next
  }
  
  primary_zph_objects[[module]] <- zph_object
  
  zph_table <- as.data.frame(
    zph_object$table
  )
  
  zph_table$Variable <- rownames(zph_table)
  rownames(zph_table) <- NULL
  
  primary_ph_result_list[[module]] <- data.frame(
    Module = module,
    Variable = zph_table$Variable,
    ChiSquare = zph_table$chisq,
    DF = zph_table$df,
    P_value = zph_table$p,
    PH_status = "OK",
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------
# 3. Combine PH results
# ----------------------------------------------------------

if (length(primary_ph_result_list) > 0) {
  
  primary_ph_results <- do.call(
    rbind,
    primary_ph_result_list
  )
  
  rownames(primary_ph_results) <- NULL
  
} else {
  
  primary_ph_results <- data.frame(
    Module = character(),
    Variable = character(),
    ChiSquare = numeric(),
    DF = numeric(),
    P_value = numeric(),
    PH_status = character(),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------
# 4. Add interpretation
# ----------------------------------------------------------

primary_ph_results$PH_interpretation <- ifelse(
  is.na(primary_ph_results$P_value),
  "Not assessed / unavailable",
  ifelse(
    primary_ph_results$P_value < 0.05,
    "Potential PH violation",
    "No evidence of PH violation"
  )
)


# ----------------------------------------------------------
# 5. Display PH results
# ----------------------------------------------------------

cat("\n========== Multivariable PH results ==========\n")

print(
  primary_ph_results,
  row.names = FALSE
)


# ----------------------------------------------------------
# 6. Save PH results
# ----------------------------------------------------------

saveRDS(
  list(
    zph_objects = primary_zph_objects,
    results = primary_ph_results
  ),
  file = file.path(
    result_dir,
    "Workflow-04-2P_WGCNA_OS_multivariable_PH.rds"
  )
)

write.table(
  primary_ph_results,
  file = file.path(
    result_dir,
    "Workflow-04-2P_WGCNA_OS_multivariable_PH.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nStep 8C-2J completed.\n")


############################################################
# Step 8C-2K
# Schoenfeld residual diagnostic plots
#
# Generate plots for all successfully fitted primary models.
############################################################

ph_plot_file <- file.path(
  result_dir,
  "Workflow-04-2Q_WGCNA_OS_multivariable_PH_diagnostic_plots.pdf"
)

pdf(
  file = ph_plot_file,
  width = 8,
  height = 7,
  onefile = TRUE
)

for (module in names(primary_zph_objects)) {
  
  zph_object <- primary_zph_objects[[module]]
  
  if (is.null(zph_object)) {
    next
  }
  
  zph_terms <- rownames(zph_object$table)
  
  # Exclude GLOBAL from variable-specific plots.
  plot_terms <- setdiff(
    zph_terms,
    "GLOBAL"
  )
  
  for (term in plot_terms) {
    
    plot(
      zph_object,
      var = term,
      main = paste(
        module,
        "-",
        term
      )
    )
    
    abline(
      h = 0,
      lty = 2
    )
  }
}

dev.off()

cat(
  "\nMultivariable PH diagnostic plots saved:\n",
  ph_plot_file,
  "\n"
)

############################################################
# Step 8C-2L
# Reproducibility checkpoint
#
# Explicitly save key analysis objects.
# No complete workspace image is saved.
############################################################

# ----------------------------------------------------------
# 1. Define checkpoint object list
# ----------------------------------------------------------

checkpoint_objects <- list(
  
  # Clinical data
  clinical_aligned = clinical_aligned,
  clinical_WGCNA_traits = clinical_WGCNA_traits,
  clinical_cox = clinical_cox,
  
  # WGCNA data
  datExpr_MAD = datExpr_MAD,
  gene_MAD = gene_MAD,
  net = net,
  moduleColors = moduleColors,
  moduleSizes = moduleSizes,
  MEs = MEs,
  MEs_biological = MEs_biological,
  
  # Univariate survival results
  Cox_results = Cox_results,
  Cox_sig = Cox_sig,
  
  # Univariate PH results
  ph_summary = ph_summary,
  ph_results = ph_results,
  
  # Multivariable model design
  cox_design = cox_design,
  selected_modules = selected_modules,
  module_ph_status = module_ph_status,
  model_design = model_design,
  
  # Complete-case cohorts
  primary_cohorts = primary_cohorts,
  sensitivity_cohorts = sensitivity_cohorts,
  primary_qc = primary_qc,
  sensitivity_qc = sensitivity_qc,
  
  # Primary multivariable Cox
  primary_models = primary_models,
  primary_model_objects = primary_model_objects,
  primary_results = primary_results,
  primary_model_qc = primary_model_qc,
  
  # Multivariable PH assessment
  primary_zph_objects = primary_zph_objects,
  primary_ph_results = primary_ph_results
)

# ----------------------------------------------------------
# 2. Save checkpoint
# ----------------------------------------------------------

checkpoint_file <- file.path(
  result_dir,
  "Workflow-04-2R_WGCNA_OS_multivariable_analysis_checkpoint.rds"
)

saveRDS(
  checkpoint_objects,
  file = checkpoint_file
)

# ----------------------------------------------------------
# 3. Save software environment
# ----------------------------------------------------------

session_info_file <- file.path(
  result_dir,
  "Workflow-04-2R_WGCNA_OS_multivariable_sessionInfo.txt"
)

writeLines(
  capture.output(sessionInfo()),
  con = session_info_file
)

# ----------------------------------------------------------
# 4. Save checkpoint object inventory
# ----------------------------------------------------------

checkpoint_inventory <- data.frame(
  Object = names(checkpoint_objects),
  Class = vapply(
    checkpoint_objects,
    function(x) paste(class(x), collapse = "; "),
    character(1)
  ),
  stringsAsFactors = FALSE
)

write.table(
  checkpoint_inventory,
  file = file.path(
    result_dir,
    "Workflow-04-2R_WGCNA_OS_multivariable_checkpoint_inventory.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ----------------------------------------------------------
# 5. Verify output files
# ----------------------------------------------------------

checkpoint_files <- c(
  checkpoint_file,
  session_info_file,
  file.path(
    result_dir,
    "Workflow-04-2R_WGCNA_OS_multivariable_checkpoint_inventory.txt"
  )
)

checkpoint_verification <- data.frame(
  File = basename(checkpoint_files),
  Exists = file.exists(checkpoint_files),
  Size_bytes = ifelse(
    file.exists(checkpoint_files),
    file.info(checkpoint_files)$size,
    NA_real_
  ),
  stringsAsFactors = FALSE
)

print(
  checkpoint_verification,
  row.names = FALSE
)

if (!all(checkpoint_verification$Exists)) {
  stop("One or more checkpoint files were not created.")
}

cat("\nStep 8C-2 reproducibility checkpoint completed.\n")


