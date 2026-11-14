
# ============================================================
# GitHub1_Paper1_differential_expression_GSE32863.R
#
# Project:
# Integrated transcriptomic and single-cell characterization
# of AOC3, MS4A7, and IL7R across the TME in LUAD
#
# Purpose:
# Independent expression assessment of AOC3, MS4A7, and IL7R
# in the GSE32863 lung adenocarcinoma cohort.
#
# Manuscript correspondence:
# Methods: Section 2.2, Differential Expression Analysis
# Figure: Figure 2D
#
# Analysis design:
# - 116 samples from 59 subjects
# - 58 lung tumor samples
# - 58 adjacent-normal samples
# - 57 complete tumor-normal pairs
# - One tumor-only and one normal-only sample
# - Paired analysis restricted to the 57 complete pairs
# - One probe per gene selected by highest mean expression
#   across all 116 samples
# - Two-sided paired Wilcoxon signed-rank tests
#
# R version:
# 4.5.1 or compatible R 4.5.x
#
# ============================================================


# ============================================================
# 00. Configuration and project paths
# ============================================================

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. ",
    "Please install it before running this script."
  )
}

project_root <- here::here()

input_dir <- file.path(
  project_root,
  "data",
  "input_data"
)

manifest_dir <- file.path(
  project_root,
  "data",
  "sample_manifests"
)

figure_dir <- file.path(
  project_root,
  "results",
  "figures"
)

table_dir <- file.path(
  project_root,
  "results",
  "tables"
)

stat_dir <- file.path(
  project_root,
  "results",
  "statistical_outputs"
)

environment_dir <- file.path(
  project_root,
  "environment"
)

output_dirs <- c(
  manifest_dir,
  figure_dir,
  table_dir,
  stat_dir,
  environment_dir
)

for (dir_path in output_dirs) {
  dir.create(
    dir_path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 01. Input files and package dependencies
# ============================================================

series_matrix_file <- file.path(
  input_dir,
  "GSE32863_series_matrix.txt.gz"
)

platform_annotation_file <- file.path(
  input_dir,
  "GPL6884-11607.txt"
)

required_packages <- c(
  "here",
  "GEOquery",
  "Biobase",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "The following required packages are missing: ",
    paste(missing_packages, collapse = ", "),
    ". Please install them before running this script."
  )
}

if (!file.exists(series_matrix_file)) {
  stop(
    "GSE32863 series matrix file not found: ",
    series_matrix_file
  )
}

if (!file.exists(platform_annotation_file)) {
  stop(
    "GPL6884 platform annotation file not found: ",
    platform_annotation_file
  )
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(ggplot2)
})


# ============================================================
# 02. Analysis configuration
# ============================================================

target_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expected_sample_count <- 116L
expected_tumor_count <- 58L
expected_normal_count <- 58L
expected_subject_count <- 59L
expected_pair_count <- 57L
expected_singleton_count <- 2L

group_levels <- c(
  "Normal",
  "Tumor"
)

figure_width <- 11
figure_height <- 4.2
figure_dpi <- 600


# ============================================================
# 03. Load GSE32863 expression data and sample metadata
# ============================================================

message("Loading GSE32863 series matrix...")

gse_object <- GEOquery::getGEO(
  filename = series_matrix_file
)

if (inherits(gse_object, "ExpressionSet")) {
  
  eset <- gse_object
  
} else if (is.list(gse_object)) {
  
  eset_index <- which(
    vapply(
      gse_object,
      inherits,
      logical(1),
      what = "ExpressionSet"
    )
  )
  
  if (length(eset_index) != 1L) {
    stop(
      "The GEO series matrix did not return exactly one ",
      "ExpressionSet object."
    )
  }
  
  eset <- gse_object[[eset_index]]
  
} else {
  
  stop(
    "Unexpected object returned by GEOquery::getGEO()."
  )
}

expr <- Biobase::exprs(eset)

pdata <- as.data.frame(
  Biobase::pData(eset),
  stringsAsFactors = FALSE
)

if (is.null(rownames(expr))) {
  stop("The expression matrix has no probe identifiers.")
}

if (is.null(colnames(expr))) {
  stop("The expression matrix has no sample identifiers.")
}

if (!is.numeric(expr)) {
  stop("The expression matrix is not numeric.")
}

if (anyDuplicated(rownames(expr)) > 0L) {
  stop("Duplicated probe identifiers detected.")
}

if (anyDuplicated(colnames(expr)) > 0L) {
  stop("Duplicated sample identifiers detected.")
}

if (ncol(expr) != expected_sample_count) {
  stop(
    "Unexpected sample count. Expected ",
    expected_sample_count,
    ", observed ",
    ncol(expr),
    "."
  )
}

if (nrow(expr) != 48803L) {
  warning(
    "The number of expression probes differs from the ",
    "previously audited matrix (48,803). Verify the ",
    "series matrix input file."
  )
}


# ============================================================
# 04. Load GPL6884 platform annotation
# ============================================================

message("Loading GPL6884 platform annotation...")

anno <- utils::read.delim(
  platform_annotation_file,
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

required_annotation_columns <- c(
  "ID",
  "Symbol"
)

if (!all(required_annotation_columns %in% colnames(anno))) {
  stop(
    "The platform annotation does not contain the required ",
    "columns: ID and Symbol."
  )
}

anno$ID <- as.character(anno$ID)
anno$Symbol <- trimws(as.character(anno$Symbol))

if (anyNA(anno$ID) || any(anno$ID == "")) {
  stop("Missing or empty probe IDs in platform annotation.")
}

if (anyDuplicated(anno$ID) > 0L) {
  stop(
    "Duplicated probe IDs were detected in the platform ",
    "annotation. Resolve annotation ambiguity before ",
    "continuing."
  )
}


# ============================================================
# 05. Construct and align the complete sample manifest
# ============================================================

metadata_columns <- c(
  "geo_accession",
  "subject:ch1",
  "tissue:ch1",
  "type:ch1",
  "title",
  "source_name_ch1"
)

missing_metadata_columns <- setdiff(
  metadata_columns,
  colnames(pdata)
)

if (length(missing_metadata_columns) > 0L) {
  stop(
    "Missing required sample metadata columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

sample_manifest <- data.frame(
  GEO_ID = as.character(pdata$geo_accession),
  Subject = as.character(pdata[["subject:ch1"]]),
  Tissue = as.character(pdata[["tissue:ch1"]]),
  Sample_Type = as.character(pdata[["type:ch1"]]),
  Title = as.character(pdata$title),
  Source_Name = as.character(pdata$source_name_ch1),
  stringsAsFactors = FALSE
)

if (anyNA(sample_manifest$GEO_ID) ||
    any(sample_manifest$GEO_ID == "")) {
  stop("Missing GEO accession identifiers in sample metadata.")
}

if (anyDuplicated(sample_manifest$GEO_ID) > 0L) {
  stop("Duplicated GEO accessions in sample metadata.")
}

# Explicitly align sample metadata to expression matrix columns.
sample_index <- match(
  colnames(expr),
  sample_manifest$GEO_ID
)

if (anyNA(sample_index)) {
  stop(
    "Some expression matrix sample IDs are absent ",
    "from the sample metadata."
  )
}

sample_manifest <- sample_manifest[
  sample_index,
  ,
  drop = FALSE
]

if (!identical(
  colnames(expr),
  sample_manifest$GEO_ID
)) {
  stop(
    "Expression matrix and sample metadata are not ",
    "in identical sample order after alignment."
  )
}

if (anyNA(sample_manifest$Subject) ||
    any(sample_manifest$Subject == "")) {
  stop("Missing subject identifiers in sample metadata.")
}

if (nrow(sample_manifest) != expected_sample_count) {
  stop("Unexpected number of sample metadata records.")
}


# ============================================================
# 06. Define tumor and adjacent-normal sample groups
# ============================================================

sample_manifest$Group <- NA_character_

sample_manifest$Group[
  sample_manifest$Tissue == "Lung tumor"
] <- "Tumor"

sample_manifest$Group[
  sample_manifest$Tissue == "Normal lung"
] <- "Normal"

if (anyNA(sample_manifest$Group)) {
  stop(
    "Unexpected tissue categories were found. ",
    "Please inspect the tissue metadata."
  )
}

sample_manifest$Group <- factor(
  sample_manifest$Group,
  levels = group_levels
)

group_counts <- table(sample_manifest$Group)

if (as.integer(group_counts["Tumor"]) != expected_tumor_count ||
    as.integer(group_counts["Normal"]) != expected_normal_count) {
  stop(
    "Unexpected tumor/normal sample counts. ",
    "Expected 58 tumor and 58 normal samples."
  )
}

if (length(unique(sample_manifest$Subject)) !=
    expected_subject_count) {
  stop(
    "Unexpected number of unique subjects. ",
    "Expected 59 subjects."
  )
}


# ============================================================
# 07. Identify complete tumor-normal pairs by subject
# ============================================================

tumor_manifest <- sample_manifest[
  sample_manifest$Group == "Tumor",
  ,
  drop = FALSE
]

normal_manifest <- sample_manifest[
  sample_manifest$Group == "Normal",
  ,
  drop = FALSE
]

tumor_subject_counts <- table(tumor_manifest$Subject)
normal_subject_counts <- table(normal_manifest$Subject)

if (any(tumor_subject_counts > 1L) ||
    any(normal_subject_counts > 1L)) {
  stop(
    "At least one subject has multiple samples within ",
    "the same tissue group. Pairing is ambiguous."
  )
}

paired_subjects <- intersect(
  tumor_manifest$Subject,
  normal_manifest$Subject
)

paired_manifest <- data.frame(
  Subject = paired_subjects,
  GEO_ID_Tumor = tumor_manifest$GEO_ID[
    match(paired_subjects, tumor_manifest$Subject)
  ],
  GEO_ID_Normal = normal_manifest$GEO_ID[
    match(paired_subjects, normal_manifest$Subject)
  ],
  stringsAsFactors = FALSE
)

if (nrow(paired_manifest) != expected_pair_count) {
  stop(
    "Unexpected number of complete pairs. Expected ",
    expected_pair_count,
    ", observed ",
    nrow(paired_manifest),
    "."
  )
}

if (anyDuplicated(paired_manifest$Subject) > 0L ||
    anyDuplicated(paired_manifest$GEO_ID_Tumor) > 0L ||
    anyDuplicated(paired_manifest$GEO_ID_Normal) > 0L) {
  stop("Duplicate subjects or sample IDs in paired manifest.")
}

paired_sample_ids <- c(
  paired_manifest$GEO_ID_Tumor,
  paired_manifest$GEO_ID_Normal
)

sample_manifest$Included_in_Paired_Analysis <-
  sample_manifest$GEO_ID %in% paired_sample_ids

singleton_manifest <- sample_manifest[
  !sample_manifest$Included_in_Paired_Analysis,
  ,
  drop = FALSE
]

if (nrow(singleton_manifest) != expected_singleton_count) {
  stop(
    "Unexpected number of singleton samples. Expected ",
    expected_singleton_count,
    "."
  )
}

# Save complete, paired, and singleton sample manifests.
utils::write.csv(
  sample_manifest,
  file.path(
    manifest_dir,
    "GSE32863_complete_sample_manifest.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  paired_manifest,
  file.path(
    manifest_dir,
    "GSE32863_57_complete_pairs.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  singleton_manifest,
  file.path(
    manifest_dir,
    "GSE32863_singleton_samples_excluded_from_paired_analysis.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 08. Align the platform annotation to expression probes
# ============================================================

annotation_index <- match(
  rownames(expr),
  anno$ID
)

matched_probe <- !is.na(annotation_index)

annotation_match_audit <- data.frame(
  Metric = c(
    "Expression probes",
    "Annotation records",
    "Matched expression probes",
    "Unmatched expression probes",
    "Annotation records not represented in expression matrix"
  ),
  Count = c(
    nrow(expr),
    nrow(anno),
    sum(matched_probe),
    sum(!matched_probe),
    sum(!anno$ID %in% rownames(expr))
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  annotation_match_audit,
  file.path(
    stat_dir,
    "GSE32863_probe_annotation_match_audit.csv"
  ),
  row.names = FALSE
)

if (sum(matched_probe) == 0L) {
  stop("No expression probes could be matched to GPL6884.")
}

expr_annotated <- expr[
  matched_probe,
  ,
  drop = FALSE
]

anno_aligned <- anno[
  annotation_index[matched_probe],
  ,
  drop = FALSE
]

if (!identical(
  rownames(expr_annotated),
  anno_aligned$ID
)) {
  stop(
    "Matched expression probes and annotation probes ",
    "are not in identical order."
  )
}


# ============================================================
# 09. Select one probe per gene using the highest mean expression
# ============================================================

# Historical Figure 2D probe-selection strategy:
# retain the probe with the highest mean expression across
# all 116 GSE32863 samples for each target gene.

target_probe_index <- which(
  !is.na(anno_aligned$Symbol) &
    anno_aligned$Symbol %in% target_genes &
    anno_aligned$Symbol != ""
)

if (length(target_probe_index) == 0L) {
  stop("No probes were annotated for the target genes.")
}

target_probe_expr <- expr_annotated[
  target_probe_index,
  ,
  drop = FALSE
]

target_probe_anno <- anno_aligned[
  target_probe_index,
  ,
  drop = FALSE
]

target_probe_anno$MeanExpr <- rowMeans(
  target_probe_expr,
  na.rm = TRUE
)

target_probe_anno <- target_probe_anno[
  is.finite(target_probe_anno$MeanExpr),
  ,
  drop = FALSE
]

if (nrow(target_probe_anno) == 0L) {
  stop("No finite expression values were found for target probes.")
}

# Deterministic ordering:
# descending mean expression, followed by probe ID.
target_probe_anno <- target_probe_anno[
  order(
    target_probe_anno$Symbol,
    -target_probe_anno$MeanExpr,
    target_probe_anno$ID
  ),
  ,
  drop = FALSE
]

# The first probe for each gene has the highest mean expression.
selected_probe_index <- match(
  target_genes,
  target_probe_anno$Symbol
)

if (anyNA(selected_probe_index)) {
  missing_genes <- target_genes[
    is.na(selected_probe_index)
  ]
  
  stop(
    "No valid probe was found for: ",
    paste(missing_genes, collapse = ", ")
  )
}

selected_probes <- target_probe_anno[
  selected_probe_index,
  ,
  drop = FALSE
]

if (anyDuplicated(selected_probes$Symbol) > 0L) {
  stop("More than one probe was selected for a target gene.")
}

if (!setequal(selected_probes$Symbol, target_genes)) {
  stop("Selected probe symbols do not match the target genes.")
}

# Save all candidate probes and the selected probe per gene.
utils::write.csv(
  target_probe_anno,
  file.path(
    table_dir,
    "GSE32863_target_gene_candidate_probes.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  selected_probes,
  file.path(
    table_dir,
    "GSE32863_selected_probe_per_gene.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 10. Extract selected-probe expression values
# ============================================================

selected_expr <- expr[
  selected_probes$ID,
  ,
  drop = FALSE
]

rownames(selected_expr) <- selected_probes$Symbol

selected_expr <- selected_expr[
  target_genes,
  ,
  drop = FALSE
]

if (!identical(
  rownames(selected_expr),
  target_genes
)) {
  stop("Selected expression matrix has unexpected gene order.")
}

if (any(!is.finite(selected_expr))) {
  stop(
    "Non-finite expression values detected for selected probes."
  )
}

expression_wide <- data.frame(
  GEO_ID = colnames(selected_expr),
  Subject = sample_manifest$Subject,
  Group = as.character(sample_manifest$Group),
  t(selected_expr),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

utils::write.csv(
  expression_wide,
  file.path(
    table_dir,
    "GSE32863_selected_gene_expression_all_samples.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. Construct the paired expression dataset
# ============================================================

paired_expression_list <- vector(
  "list",
  length(target_genes)
)

names(paired_expression_list) <- target_genes

for (gene in target_genes) {
  
  tumor_values <- selected_expr[
    gene,
    paired_manifest$GEO_ID_Tumor
  ]
  
  normal_values <- selected_expr[
    gene,
    paired_manifest$GEO_ID_Normal
  ]
  
  paired_expression_list[[gene]] <- data.frame(
    Subject = paired_manifest$Subject,
    Normal = as.numeric(normal_values),
    Tumor = as.numeric(tumor_values),
    stringsAsFactors = FALSE
  )
}

paired_expression_long <- do.call(
  rbind,
  lapply(
    target_genes,
    function(gene) {
      
      gene_data <- paired_expression_list[[gene]]
      
      rbind(
        data.frame(
          Gene = gene,
          Subject = gene_data$Subject,
          Group = "Normal",
          Expression = gene_data$Normal,
          stringsAsFactors = FALSE
        ),
        data.frame(
          Gene = gene,
          Subject = gene_data$Subject,
          Group = "Tumor",
          Expression = gene_data$Tumor,
          stringsAsFactors = FALSE
        )
      )
    }
  )
)

paired_expression_long$Group <- factor(
  paired_expression_long$Group,
  levels = group_levels
)

utils::write.csv(
  paired_expression_long,
  file.path(
    table_dir,
    "GSE32863_paired_expression_long_format.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. Paired two-sided Wilcoxon signed-rank tests
# ============================================================

wilcoxon_results <- lapply(
  target_genes,
  function(gene) {
    
    gene_data <- paired_expression_list[[gene]]
    
    test_result <- stats::wilcox.test(
      x = gene_data$Tumor,
      y = gene_data$Normal,
      paired = TRUE,
      alternative = "two.sided",
      exact = FALSE,
      correct = TRUE
    )
    
    data.frame(
      Gene = gene,
      N_Pairs = nrow(gene_data),
      Median_Normal = median(gene_data$Normal),
      Median_Tumor = median(gene_data$Tumor),
      Median_Paired_Difference = median(
        gene_data$Tumor - gene_data$Normal
      ),
      Wilcoxon_Statistic = unname(test_result$statistic),
      P_Value = test_result$p.value,
      stringsAsFactors = FALSE
    )
  }
)

wilcoxon_results <- do.call(
  rbind,
  wilcoxon_results
)

# BH-adjusted P values are provided as an additional
# multiple-testing summary across the three gene-level tests.
wilcoxon_results$P_Value_BH <- stats::p.adjust(
  wilcoxon_results$P_Value,
  method = "BH"
)

utils::write.csv(
  wilcoxon_results,
  file.path(
    stat_dir,
    "GSE32863_paired_Wilcoxon_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13. Prepare statistical annotations for Figure 2D
# ============================================================

format_p_value <- function(p_value) {
  
  if (is.na(p_value)) {
    return("P = NA")
  }
  
  if (p_value < 0.001) {
    return("P < 0.001")
  }
  
  paste0(
    "P = ",
    format.pval(
      p_value,
      digits = 3,
      eps = 0.001
    )
  )
}

annotation_data <- do.call(
  rbind,
  lapply(
    target_genes,
    function(gene) {
      
      gene_values <- paired_expression_long$Expression[
        paired_expression_long$Gene == gene
      ]
      
      gene_range <- range(
        gene_values,
        na.rm = TRUE
      )
      
      y_span <- diff(gene_range)
      
      if (!is.finite(y_span) || y_span == 0) {
        y_span <- max(abs(gene_range), 1)
      }
      
      p_value <- wilcoxon_results$P_Value[
        wilcoxon_results$Gene == gene
      ]
      
      data.frame(
        Gene = gene,
        x = 1.5,
        y = max(gene_values, na.rm = TRUE) + 0.16 * y_span,
        label = format_p_value(p_value),
        stringsAsFactors = FALSE
      )
    }
  )
)


# ============================================================
# 14. Generate Figure 2D
# ============================================================

plot_data <- paired_expression_long

plot_data$Group <- factor(
  plot_data$Group,
  levels = group_levels
)

figure2D <- ggplot(
  plot_data,
  aes(
    x = Group,
    y = Expression,
    fill = Group
  )
) +
  
  geom_violin(
    trim = FALSE,
    alpha = 0.65,
    colour = NA,
    width = 0.85
  ) +
  
  # Each line connects the matched samples from the same subject.
  geom_line(
    aes(
      group = Subject
    ),
    colour = "grey55",
    alpha = 0.35,
    linewidth = 0.35
  ) +
  
  geom_boxplot(
    width = 0.16,
    fill = "white",
    colour = "black",
    linewidth = 0.45,
    outlier.shape = NA
  ) +
  
  geom_segment(
    data = annotation_data,
    aes(
      x = 1,
      xend = 2,
      y = y,
      yend = y
    ),
    inherit.aes = FALSE,
    linewidth = 0.45
  ) +
  
  geom_text(
    data = annotation_data,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3.5,
    vjust = -0.5
  ) +
  
  facet_wrap(
    ~ Gene,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Normal" = "#4DBBD5",
      "Tumor" = "#E64B35"
    )
  ) +
  
  scale_y_continuous(
    name = "Expression",
    expand = expansion(
      mult = c(0.05, 0.30)
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  labs(
    x = NULL,
    y = "Expression"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 11
    ),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    plot.margin = margin(
      t = 15,
      r = 12,
      b = 10,
      l = 10
    )
  )

print(figure2D)


# ============================================================
# 15. Save Figure 2D
# ============================================================

figure_pdf <- file.path(
  figure_dir,
  "Figure2D_GSE32863.pdf"
)

figure_tiff <- file.path(
  figure_dir,
  "Figure2D_GSE32863.tiff"
)

ggplot2::ggsave(
  filename = figure_pdf,
  plot = figure2D,
  width = figure_width,
  height = figure_height,
  device = if (capabilities("cairo")) {
    grDevices::cairo_pdf
  } else {
    grDevices::pdf
  },
  limitsize = FALSE
)

ggplot2::ggsave(
  filename = figure_tiff,
  plot = figure2D,
  width = figure_width,
  height = figure_height,
  dpi = figure_dpi,
  compression = "lzw",
  limitsize = FALSE
)


# ============================================================
# 16. Save analysis summary and software environment
# ============================================================

analysis_summary <- data.frame(
  Item = c(
    "Dataset",
    "Platform",
    "Total samples",
    "Tumor samples",
    "Adjacent-normal samples",
    "Unique subjects",
    "Complete tumor-normal pairs",
    "Singleton samples excluded from paired testing",
    "Target genes",
    "Probe selection",
    "Statistical test",
    "Alternative hypothesis",
    "Multiple-testing summary",
    "Manuscript figure"
  ),
  Value = c(
    "GSE32863",
    "GPL6884",
    ncol(expr),
    sum(sample_manifest$Group == "Tumor"),
    sum(sample_manifest$Group == "Normal"),
    length(unique(sample_manifest$Subject)),
    nrow(paired_manifest),
    nrow(singleton_manifest),
    paste(target_genes, collapse = ", "),
    paste(
      "One probe per gene with the highest mean expression",
      "across all 116 samples"
    ),
    "Paired Wilcoxon signed-rank test",
    "Two-sided",
    "BH adjustment across the three gene-level tests",
    "Figure 2D"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_summary,
  file.path(
    stat_dir,
    "GSE32863_analysis_summary.csv"
  ),
  row.names = FALSE
)

# Record the R and package environment used for this run.
session_info_file <- file.path(
  environment_dir,
  "GSE32863_sessionInfo.txt"
)

capture.output(
  sessionInfo(),
  file = session_info_file
)

package_versions <- data.frame(
  Package = required_packages,
  Version = vapply(
    required_packages,
    function(pkg) {
      as.character(utils::packageVersion(pkg))
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  package_versions,
  file.path(
    stat_dir,
    "GSE32863_package_versions.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. Final execution report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("GSE32863 Figure 2D analysis completed.\n")
cat("============================================================\n")

cat(
  "Expression matrix:",
  nrow(expr),
  "probes x",
  ncol(expr),
  "samples\n"
)

cat(
  "Tumor samples:",
  sum(sample_manifest$Group == "Tumor"),
  "\n"
)

cat(
  "Adjacent-normal samples:",
  sum(sample_manifest$Group == "Normal"),
  "\n"
)

cat(
  "Unique subjects:",
  length(unique(sample_manifest$Subject)),
  "\n"
)

cat(
  "Complete tumor-normal pairs:",
  nrow(paired_manifest),
  "\n"
)

cat(
  "Singleton samples excluded from paired testing:",
  nrow(singleton_manifest),
  "\n"
)

cat("\nSelected probes:\n")

print(
  selected_probes[
    ,
    c("Symbol", "ID", "MeanExpr"),
    drop = FALSE
  ],
  row.names = FALSE
)

cat("\nPaired Wilcoxon results:\n")

print(
  wilcoxon_results,
  row.names = FALSE
)

cat("\nFigure files:\n")
cat(figure_pdf, "\n")
cat(figure_tiff, "\n")

cat("\nAnalysis outputs were saved under:\n")
cat(file.path(project_root, "results"), "\n")

cat("\n============================================================\n")
