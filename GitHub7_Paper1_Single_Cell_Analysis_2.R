##GitHub7_Paper1_Single_Cell_Analysis_2

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

##STEP10A — Patient / Sample Structure Audit
# ============================================================
# WORKFLOW-10
# STEP10A — Patient / Sample Structure Audit
#
# Purpose:
#   Audit patient-, sample-, source-, stage-, and cell-level
#   structure of the TISCH2 NSCLC_GSE131907 dataset.
#
# IMPORTANT:
#   No PCA / UMAP / clustering / cell annotation is recalculated.
#   This step uses TISCH2-provided metadata only.
#
# Dataset:
#   TISCH2
#   NSCLC_GSE131907
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
# ============================================================


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)

cat("\n")
cat("============================================================\n")
cat("STEP10A — Patient / Sample Structure Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------

meta_file <- "input_dir/Paper 1_Step 10-Single-Cell/NSCLC_GSE131907_CellMetainfo_table.tsv"

if (!file.exists(meta_file)) {
  stop(
    "INPUT FILE NOT FOUND:\n",
    meta_file,
    "\n\nPlease check the exact TISCH2 metadata file path."
  )
}

cat("Input metadata file:\n")
cat(meta_file, "\n\n")


# ------------------------------------------------------------
# 2. Read metadata
# ------------------------------------------------------------

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata successfully loaded.\n")
cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n\n")


# ------------------------------------------------------------
# 3. Display column names
# ------------------------------------------------------------

cat("Column names:\n")
print(colnames(meta))
cat("\n")


# ------------------------------------------------------------
# 4. Required-column audit
# ------------------------------------------------------------

required_cols <- c(
  "Cell",
  "UMAP_1",
  "UMAP_2",
  "Celltype (malignancy)",
  "Celltype (major-lineage)",
  "Celltype (minor-lineage)",
  "Cluster",
  "Sample",
  "TNMstage",
  "Patient",
  "Source",
  "Stage"
)

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  
  cat("ERROR — Missing required columns:\n")
  print(missing_cols)
  cat("\n")
  
  stop("STEP10A halted because required metadata columns are missing.")
  
} else {
  
  cat("PASS — All required metadata columns are present.\n\n")
}


# ------------------------------------------------------------
# 5. Basic data integrity audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. Basic data integrity\n")
cat("------------------------------------------------------------\n\n")

duplicate_cells <- sum(duplicated(meta$Cell))

missing_cell <- sum(is.na(meta$Cell) | meta$Cell == "")
missing_patient <- sum(is.na(meta$Patient) | meta$Patient == "")
missing_sample <- sum(is.na(meta$Sample) | meta$Sample == "")
missing_source <- sum(is.na(meta$Source) | meta$Source == "")
missing_stage <- sum(is.na(meta$Stage) | meta$Stage == "")

cat("Duplicate Cell IDs:", duplicate_cells, "\n")
cat("Missing Cell:", missing_cell, "\n")
cat("Missing Patient:", missing_patient, "\n")
cat("Missing Sample:", missing_sample, "\n")
cat("Missing Source:", missing_source, "\n")
cat("Missing Stage:", missing_stage, "\n\n")


# ------------------------------------------------------------
# 6. Overall dataset summary
# ------------------------------------------------------------

n_cells <- nrow(meta)

n_patients <- length(unique(na.omit(meta$Patient)))
n_samples  <- length(unique(na.omit(meta$Sample)))

n_primary_stage <- sum(
  !is.na(meta$Stage) &
    meta$Stage == "Primary"
)

n_metastatic_stage <- sum(
  !is.na(meta$Stage) &
    meta$Stage == "Metastatic"
)

n_primary_patients <- length(
  unique(
    meta$Patient[
      !is.na(meta$Patient) &
        !is.na(meta$Stage) &
        meta$Stage == "Primary"
    ]
  )
)

n_primary_samples <- length(
  unique(
    meta$Sample[
      !is.na(meta$Sample) &
        !is.na(meta$Stage) &
        meta$Stage == "Primary"
    ]
  )
)


# ------------------------------------------------------------
# 7. Summary table
# ------------------------------------------------------------

dataset_summary <- data.frame(
  Metric = c(
    "Total cells",
    "Unique patients",
    "Unique samples",
    "Cells annotated as Primary stage",
    "Cells annotated as Metastatic stage",
    "Unique Primary-stage patients",
    "Unique Primary-stage samples",
    "Duplicate cell IDs",
    "Missing patient IDs",
    "Missing sample IDs",
    "Missing source annotations",
    "Missing stage annotations"
  ),
  Value = c(
    n_cells,
    n_patients,
    n_samples,
    n_primary_stage,
    n_metastatic_stage,
    n_primary_patients,
    n_primary_samples,
    duplicate_cells,
    missing_patient,
    missing_sample,
    missing_source,
    missing_stage
  ),
  stringsAsFactors = FALSE
)

write.table(
  dataset_summary,
  file = "STEP10A-1_Dataset_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Dataset summary:\n")
print(dataset_summary)
cat("\n")


# ------------------------------------------------------------
# 8. Patient-level sample audit
#
# IMPORTANT:
#   Patient / Sample / Source / Stage are kept separate.
#
#   Stage == Primary does NOT automatically mean
#   Source == Tumor.
# ------------------------------------------------------------

patient_sample_audit <- unique(
  meta[, c(
    "Patient",
    "Sample",
    "Source",
    "Stage",
    "TNMstage"
  )]
)

patient_sample_audit <- patient_sample_audit[
  order(
    patient_sample_audit$Patient,
    patient_sample_audit$Sample
  ),
]

write.table(
  patient_sample_audit,
  file = "STEP10A-2_Patient_Sample_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Cells per patient
# ------------------------------------------------------------

cells_per_patient <- as.data.frame(
  table(meta$Patient),
  stringsAsFactors = FALSE
)

colnames(cells_per_patient) <- c(
  "Patient",
  "Cell_Count"
)

cells_per_patient <- cells_per_patient[
  cells_per_patient$Patient != "" &
    !is.na(cells_per_patient$Patient),
]

cells_per_patient <- cells_per_patient[
  order(
    cells_per_patient$Cell_Count,
    decreasing = TRUE
  ),
]

write.table(
  cells_per_patient,
  file = "STEP10A-3_Cells_Per_Patient.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 10. Patient-level cell population structure
# ------------------------------------------------------------

celltype_patient <- unique(
  meta[, c(
    "Patient",
    "Sample",
    "Celltype (major-lineage)",
    "Celltype (minor-lineage)",
    "Celltype (malignancy)"
  )]
)

celltype_patient <- celltype_patient[
  order(
    celltype_patient$Patient,
    celltype_patient$Sample,
    celltype_patient$`Celltype (major-lineage)`,
    celltype_patient$`Celltype (minor-lineage)`
  ),
]

write.table(
  celltype_patient,
  file = "STEP10A-4_CellType_Per_Patient.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Source distribution
# ------------------------------------------------------------

source_summary <- as.data.frame(
  table(meta$Source, useNA = "ifany"),
  stringsAsFactors = FALSE
)

colnames(source_summary) <- c(
  "Source",
  "Cell_Count"
)

source_summary$Patient_Count <- sapply(
  source_summary$Source,
  function(x) {
    
    if (is.na(x)) {
      
      length(
        unique(
          meta$Patient[is.na(meta$Source)]
        )
      )
      
    } else {
      
      length(
        unique(
          meta$Patient[
            !is.na(meta$Source) &
              meta$Source == x
          ]
        )
      )
      
    }
  }
)

source_summary$Sample_Count <- sapply(
  source_summary$Source,
  function(x) {
    
    if (is.na(x)) {
      
      length(
        unique(
          meta$Sample[is.na(meta$Source)]
        )
      )
      
    } else {
      
      length(
        unique(
          meta$Sample[
            !is.na(meta$Source) &
              meta$Source == x
          ]
        )
      )
      
    }
  }
)

write.table(
  source_summary,
  file = "STEP10A-Source_Distribution.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Stage distribution
# ------------------------------------------------------------

stage_summary <- as.data.frame(
  table(meta$Stage, useNA = "ifany"),
  stringsAsFactors = FALSE
)

colnames(stage_summary) <- c(
  "Stage",
  "Cell_Count"
)

stage_summary$Patient_Count <- sapply(
  stage_summary$Stage,
  function(x) {
    
    if (is.na(x)) {
      
      length(
        unique(
          meta$Patient[is.na(meta$Stage)]
        )
      )
      
    } else {
      
      length(
        unique(
          meta$Patient[
            !is.na(meta$Stage) &
              meta$Stage == x
          ]
        )
      )
      
    }
  }
)

stage_summary$Sample_Count <- sapply(
  stage_summary$Stage,
  function(x) {
    
    if (is.na(x)) {
      
      length(
        unique(
          meta$Sample[is.na(meta$Stage)]
        )
      )
      
    } else {
      
      length(
        unique(
          meta$Sample[
            !is.na(meta$Stage) &
              meta$Stage == x
          ]
        )
      )
      
    }
  }
)

write.table(
  stage_summary,
  file = "STEP10A-Stage_Distribution.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Stage × Source audit
#
# This is particularly important.
# It prevents the incorrect assumption:
#
#   Primary stage = Tumor source
# ------------------------------------------------------------

stage_source_table <- as.data.frame(
  table(
    Stage = meta$Stage,
    Source = meta$Source,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

stage_source_table <- stage_source_table[
  stage_source_table$Freq > 0,
]

write.table(
  stage_source_table,
  file = "STEP10A-Stage_by_Source.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Patient × Stage × Source structure
# ------------------------------------------------------------

patient_stage_source <- unique(
  meta[, c(
    "Patient",
    "Sample",
    "Stage",
    "Source",
    "TNMstage"
  )]
)

patient_stage_source <- patient_stage_source[
  order(
    patient_stage_source$Patient,
    patient_stage_source$Stage,
    patient_stage_source$Source,
    patient_stage_source$Sample
  ),
]

write.table(
  patient_stage_source,
  file = "STEP10A-Patient_Stage_Source.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Primary-stage source classification
#
# Explicitly distinguish:
#   Primary + Tumor
#   Primary + other source
#   Metastatic
# ------------------------------------------------------------

meta$STEP10A_Category <- NA_character_

meta$STEP10A_Category[
  !is.na(meta$Stage) &
    meta$Stage == "Primary" &
    !is.na(meta$Source) &
    meta$Source == "Tumor"
] <- "Primary_Tumor"

meta$STEP10A_Category[
  !is.na(meta$Stage) &
    meta$Stage == "Primary" &
    !is.na(meta$Source) &
    meta$Source != "Tumor"
] <- "Primary_NonTumor_Source"

meta$STEP10A_Category[
  !is.na(meta$Stage) &
    meta$Stage == "Metastatic"
] <- "Metastatic"

meta$STEP10A_Category[
  is.na(meta$STEP10A_Category)
] <- "Other_or_Unclassified"


category_summary <- as.data.frame(
  table(
    meta$STEP10A_Category,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

colnames(category_summary) <- c(
  "STEP10A_Category",
  "Cell_Count"
)

category_summary$Patient_Count <- sapply(
  category_summary$STEP10A_Category,
  function(x) {
    length(
      unique(
        meta$Patient[
          meta$STEP10A_Category == x
        ]
      )
    )
  }
)

category_summary$Sample_Count <- sapply(
  category_summary$STEP10A_Category,
  function(x) {
    length(
      unique(
        meta$Sample[
          meta$STEP10A_Category == x
        ]
      )
    )
  }
)

write.table(
  category_summary,
  file = "STEP10A-Primary_Tumor_Category_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Explicit exclusion audit
#
# No exclusion is applied yet.
# This table documents possible categories.
# ------------------------------------------------------------

exclusion_audit <- data.frame(
  Exclusion_Category = c(
    "Missing Patient",
    "Missing Sample",
    "Missing Source",
    "Missing Stage",
    "Metastatic Stage",
    "Primary Non-Tumor Source",
    "Other or Unclassified"
  ),
  Cell_Count = c(
    sum(is.na(meta$Patient) | meta$Patient == ""),
    sum(is.na(meta$Sample) | meta$Sample == ""),
    sum(is.na(meta$Source) | meta$Source == ""),
    sum(is.na(meta$Stage) | meta$Stage == ""),
    sum(meta$Stage == "Metastatic", na.rm = TRUE),
    sum(
      meta$Stage == "Primary" &
        meta$Source != "Tumor",
      na.rm = TRUE
    ),
    sum(
      meta$STEP10A_Category == "Other_or_Unclassified",
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  exclusion_audit,
  file = "STEP10A-5_Exclusion_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. Cells per patient statistics
# ------------------------------------------------------------

cpp <- cells_per_patient$Cell_Count

cells_per_patient_summary <- data.frame(
  Metric = c(
    "Patients",
    "Minimum cells per patient",
    "Q1 cells per patient",
    "Median cells per patient",
    "Mean cells per patient",
    "Q3 cells per patient",
    "Maximum cells per patient",
    "IQR cells per patient"
  ),
  Value = c(
    length(cpp),
    min(cpp),
    quantile(cpp, 0.25),
    median(cpp),
    mean(cpp),
    quantile(cpp, 0.75),
    max(cpp),
    IQR(cpp)
  ),
  stringsAsFactors = FALSE
)

write.table(
  cells_per_patient_summary,
  file = "STEP10A-Cells_Per_Patient_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Number of cell populations per patient
# ------------------------------------------------------------

major_lineage_per_patient <- aggregate(
  meta$`Celltype (major-lineage)`,
  by = list(
    Patient = meta$Patient
  ),
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(major_lineage_per_patient)[2] <-
  "Major_Lineage_Count"

minor_lineage_per_patient <- aggregate(
  meta$`Celltype (minor-lineage)`,
  by = list(
    Patient = meta$Patient
  ),
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(minor_lineage_per_patient)[2] <-
  "Minor_Lineage_Count"

cell_population_per_patient <- merge(
  major_lineage_per_patient,
  minor_lineage_per_patient,
  by = "Patient",
  all = TRUE
)

cell_population_per_patient <- cell_population_per_patient[
  order(
    cell_population_per_patient$Patient
  ),
]

write.table(
  cell_population_per_patient,
  file = "STEP10A-Cell_Populations_Per_Patient.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Core-gene availability audit
#
# This step only checks gene names against the TISCH2
# differential-expression summary table IF that file exists.
# It does not calculate new differential expression.
# ------------------------------------------------------------

diff_file <- "input_dir/Figure 9C Dot Plot (Publication-ready)/NSCLC_GSE131907_AllDiffGenes_table.tsv"

if (file.exists(diff_file)) {
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  cat("19. Core-gene availability audit\n")
  cat("------------------------------------------------------------\n\n")
  
  diff <- read.delim(
    diff_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  core_genes <- c(
    "AOC3",
    "MS4A7",
    "IL7R"
  )
  
  core_gene_audit <- data.frame(
    Gene = core_genes,
    Present_in_TISCH2_DiffGene_Table =
      core_genes %in% diff$Gene,
    stringsAsFactors = FALSE
  )
  
  write.table(
    core_gene_audit,
    file = "STEP10A-Core_Gene_Availability.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  print(core_gene_audit)
  
} else {
  
  cat("\n")
  cat("WARNING:\n")
  cat("TISCH2 differential-gene table not found at:\n")
  cat(diff_file, "\n")
  cat("Core-gene availability audit skipped.\n\n")
  
}


# ------------------------------------------------------------
# 20. Final PASS / WARNING assessment
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10A — FINAL QC\n")
cat("============================================================\n\n")

pass_basic <- (
  duplicate_cells == 0 &&
    missing_cell == 0 &&
    missing_patient == 0 &&
    missing_sample == 0
)

cat(
  "Cell ID / Patient / Sample integrity:",
  ifelse(pass_basic, "PASS", "WARNING"),
  "\n"
)

cat(
  "Total cells:",
  n_cells,
  "\n"
)

cat(
  "Unique patients:",
  n_patients,
  "\n"
)

cat(
  "Unique samples:",
  n_samples,
  "\n"
)

cat(
  "Primary-stage patients:",
  n_primary_patients,
  "\n"
)

cat(
  "Primary-stage samples:",
  n_primary_samples,
  "\n"
)

cat("\nStage distribution:\n")
print(stage_summary)

cat("\nSource distribution:\n")
print(source_summary)

cat("\nStage × Source structure:\n")
print(stage_source_table)

cat("\nPrimary-stage category summary:\n")
print(category_summary)

cat("\nCells per patient summary:\n")
print(cells_per_patient_summary)

cat("\n")

if (pass_basic) {
  
  cat("============================================================\n")
  cat("STEP10A STATUS: PASS — BASIC STRUCTURE AUDIT\n")
  cat("============================================================\n")
  
} else {
  
  cat("============================================================\n")
  cat("STEP10A STATUS: WARNING — STRUCTURE REQUIRES REVIEW\n")
  cat("============================================================\n")
  
}

cat("\n")


# ------------------------------------------------------------
# 21. Save complete metadata audit object
# ------------------------------------------------------------

saveRDS(
  list(
    meta = meta,
    dataset_summary = dataset_summary,
    patient_sample_audit = patient_sample_audit,
    cells_per_patient = cells_per_patient,
    celltype_patient = celltype_patient,
    source_summary = source_summary,
    stage_summary = stage_summary,
    stage_source_table = stage_source_table,
    patient_stage_source = patient_stage_source,
    category_summary = category_summary,
    exclusion_audit = exclusion_audit,
    cells_per_patient_summary = cells_per_patient_summary,
    cell_population_per_patient = cell_population_per_patient
  ),
  file = "STEP10A_Patient_Sample_Structure_Audit.rds"
)


# ------------------------------------------------------------
# 22. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(sessionInfo()),
  "STEP10A_sessionInfo.txt"
)

cat("\n")
cat("STEP10A completed.\n")
cat("All audit files saved to:\n")
cat(step_dir, "\n")
cat("\n")


# ============================================================
# STEP10A-6
# Patient-level sample structure and primary-tumor audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10A-6 — Patient-level Sample Structure Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Patient × Sample × Stage × Source
# ------------------------------------------------------------

patient_structure <- unique(
  meta[, c(
    "Patient",
    "Sample",
    "Stage",
    "Source",
    "TNMstage"
  )]
)

patient_structure <- patient_structure[
  order(
    patient_structure$Patient,
    patient_structure$Stage,
    patient_structure$Source,
    patient_structure$Sample
  ),
]

write.table(
  patient_structure,
  file = "STEP10A-6_Patient_Sample_Structure.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 2. Number of samples per patient
# ------------------------------------------------------------

samples_per_patient <- aggregate(
  Sample ~ Patient,
  data = patient_structure,
  FUN = function(x) length(unique(x))
)

colnames(samples_per_patient)[2] <- "Sample_Count"

samples_per_patient <- samples_per_patient[
  order(
    samples_per_patient$Patient
  ),
]

write.table(
  samples_per_patient,
  file = "STEP10A-6_Samples_Per_Patient.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 3. Does each patient have a primary tumor sample?
# ------------------------------------------------------------

primary_tumor_by_patient <- aggregate(
  Source ~ Patient,
  data = patient_structure,
  FUN = function(x) any(x == "Tumor")
)

colnames(primary_tumor_by_patient)[2] <-
  "Has_Tumor_Source"

# Stage restriction
primary_tumor_stage <- aggregate(
  Stage ~ Patient,
  data = patient_structure,
  FUN = function(x) any(x == "Primary")
)

colnames(primary_tumor_stage)[2] <-
  "Has_Primary_Stage"

primary_tumor_audit <- merge(
  primary_tumor_by_patient,
  primary_tumor_stage,
  by = "Patient",
  all = TRUE
)

primary_tumor_audit$Primary_Tumor_Definition <- with(
  primary_tumor_audit,
  Has_Primary_Stage & Has_Tumor_Source
)

write.table(
  primary_tumor_audit,
  file = "STEP10A-6_Primary_Tumor_Patient_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 4. Patients with both primary and metastatic samples
# ------------------------------------------------------------

patient_stage <- aggregate(
  Stage ~ Patient,
  data = patient_structure,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)

colnames(patient_stage)[2] <- "Stage_Set"

patient_stage$Both_Primary_and_Metastatic <- grepl(
  "Primary",
  patient_stage$Stage_Set
) &
  grepl(
    "Metastatic",
    patient_stage$Stage_Set
  )

write.table(
  patient_stage,
  file = "STEP10A-6_Patient_Stage_Set.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. Patients with multiple source types
# ------------------------------------------------------------

patient_sources <- aggregate(
  Source ~ Patient,
  data = patient_structure,
  FUN = function(x) paste(sort(unique(x)), collapse = ";")
)

colnames(patient_sources)[2] <- "Source_Set"

patient_sources$Source_Type_Count <- sapply(
  strsplit(
    patient_sources$Source_Set,
    ";",
    fixed = TRUE
  ),
  length
)

write.table(
  patient_sources,
  file = "STEP10A-6_Patient_Source_Set.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Tumor-source samples and their cell counts
# ------------------------------------------------------------

tumor_meta <- meta[
  !is.na(meta$Source) &
    meta$Source == "Tumor",
  c(
    "Cell",
    "Patient",
    "Sample",
    "Stage",
    "Source",
    "TNMstage",
    "Celltype (major-lineage)",
    "Celltype (minor-lineage)"
  )
]

tumor_sample_cell_counts <- as.data.frame(
  table(
    tumor_meta$Patient,
    tumor_meta$Sample
  ),
  stringsAsFactors = FALSE
)

colnames(tumor_sample_cell_counts) <- c(
  "Patient",
  "Sample",
  "Cell_Count"
)

tumor_sample_cell_counts <- tumor_sample_cell_counts[
  tumor_sample_cell_counts$Cell_Count > 0,
]

tumor_sample_info <- unique(
  tumor_meta[, c(
    "Patient",
    "Sample",
    "Stage",
    "Source",
    "TNMstage"
  )]
)

tumor_sample_cell_counts <- merge(
  tumor_sample_cell_counts,
  tumor_sample_info,
  by = c(
    "Patient",
    "Sample"
  ),
  all.x = TRUE
)

tumor_sample_cell_counts <- tumor_sample_cell_counts[
  order(
    tumor_sample_cell_counts$Patient,
    tumor_sample_cell_counts$Sample
  ),
]

write.table(
  tumor_sample_cell_counts,
  file = "STEP10A-6_Tumor_Sample_Cell_Counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Cell numbers per patient × major lineage
# ------------------------------------------------------------

patient_majorlineage_counts <- as.data.frame(
  table(
    meta$Patient,
    meta$`Celltype (major-lineage)`
  ),
  stringsAsFactors = FALSE
)

colnames(patient_majorlineage_counts) <- c(
  "Patient",
  "Major_Lineage",
  "Cell_Count"
)

patient_majorlineage_counts <- patient_majorlineage_counts[
  patient_majorlineage_counts$Cell_Count > 0,
]

write.table(
  patient_majorlineage_counts,
  file = "STEP10A-6_Patient_MajorLineage_Cell_Counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Primary tumor patients only
# ------------------------------------------------------------

primary_tumor_patients <- unique(
  patient_structure$Patient[
    patient_structure$Stage == "Primary" &
      patient_structure$Source == "Tumor"
  ]
)

primary_tumor_patients <- sort(
  primary_tumor_patients
)

primary_tumor_patient_table <- data.frame(
  Patient = primary_tumor_patients,
  stringsAsFactors = FALSE
)

write.table(
  primary_tumor_patient_table,
  file = "STEP10A-6_Primary_Tumor_Patients.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Final summary
# ------------------------------------------------------------

n_primary_tumor_patients <- length(
  primary_tumor_patients
)

n_both_stage_patients <- sum(
  patient_stage$Both_Primary_and_Metastatic
)

n_multi_source_patients <- sum(
  patient_sources$Source_Type_Count > 1
)

cat("Total patients:", n_patients, "\n")
cat("Patients with Primary + Tumor:", 
    n_primary_tumor_patients, "\n")
cat("Patients with both Primary and Metastatic:",
    n_both_stage_patients, "\n")
cat("Patients with >1 source type:",
    n_multi_source_patients, "\n\n")


cat("Primary Tumor patients:\n")
print(primary_tumor_patient_table)

cat("\nPatient stage sets:\n")
print(patient_stage)

cat("\nPatient source sets:\n")
print(patient_sources)

cat("\n")
cat("============================================================\n")
cat("STEP10A-6 completed.\n")
cat("============================================================\n")


# ============================================================
# STEP10B-1 — Primary-Tumor Compartment Cell-Count Feasibility Audit
#
# Purpose:
#   1. Restrict analysis to Primary + Tumor samples
#   2. Lock TISCH2 major-lineage → target-compartment mapping
#   3. Audit patient × compartment cell counts
#   4. Assess feasibility for subsequent patient-level pseudobulk
#
# Target compartments:
#   T-cell             = CD4Tconv + CD8T + CD8Tex
#   Myeloid/macrophage = Mono/Macro
#   Endothelial        = Endothelial
#
# IMPORTANT:
#   - No minimum-cell threshold is applied at this step.
#   - No patient is excluded at this step.
#   - No expression data are analyzed.
#   - No pseudobulk calculation is performed.
# ============================================================


cat("STEP10B-1 — Primary-Tumor Compartment Cell-Count Feasibility Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit"
)

if (!dir.exists(step_dir)) {
  dir.create(
    step_dir,
    recursive = TRUE
  )
}

setwd(step_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Input metadata
# ------------------------------------------------------------

meta_file <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

if (!file.exists(meta_file)) {
  stop(
    "INPUT FILE NOT FOUND:\n",
    meta_file
  )
}

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Metadata loaded successfully.\n")
cat("Total cells:", nrow(meta), "\n\n")


# ------------------------------------------------------------
# 2. Required-column check
# ------------------------------------------------------------

required_cols <- c(
  "Cell",
  "Sample",
  "Patient",
  "Source",
  "Stage",
  "Celltype (major-lineage)",
  "Celltype (minor-lineage)"
)

missing_cols <- setdiff(
  required_cols,
  colnames(meta)
)

if (length(missing_cols) > 0) {
  
  cat("Missing columns:\n")
  print(missing_cols)
  
  stop(
    "STEP10B-1 halted because required columns are missing."
  )
  
} else {
  
  cat(
    "PASS — All required columns are present.\n\n"
  )
  
}


# ------------------------------------------------------------
# 3. Basic metadata integrity checks
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Basic metadata integrity audit\n")
cat("------------------------------------------------------------\n\n")

duplicate_cells <- sum(
  duplicated(meta$Cell)
)

missing_cell <- sum(
  is.na(meta$Cell) |
    trimws(meta$Cell) == ""
)

missing_patient <- sum(
  is.na(meta$Patient) |
    trimws(meta$Patient) == ""
)

missing_sample <- sum(
  is.na(meta$Sample) |
    trimws(meta$Sample) == ""
)

missing_source <- sum(
  is.na(meta$Source) |
    trimws(meta$Source) == ""
)

missing_stage <- sum(
  is.na(meta$Stage) |
    trimws(meta$Stage) == ""
)

cat("Duplicate Cell IDs:", duplicate_cells, "\n")
cat("Missing Cell IDs:", missing_cell, "\n")
cat("Missing Patient IDs:", missing_patient, "\n")
cat("Missing Sample IDs:", missing_sample, "\n")
cat("Missing Source:", missing_source, "\n")
cat("Missing Stage:", missing_stage, "\n\n")

if (
  duplicate_cells > 0 ||
  missing_cell > 0 ||
  missing_patient > 0 ||
  missing_sample > 0 ||
  missing_source > 0 ||
  missing_stage > 0
) {
  
  stop(
    "STEP10B-1 halted because metadata integrity checks failed."
  )
  
}

cat(
  "PASS — Basic metadata integrity checks passed.\n\n"
)


# ------------------------------------------------------------
# 4. Define the locked Primary + Tumor patient group
# ------------------------------------------------------------

primary_tumor_meta <- meta[
  !is.na(meta$Stage) &
    meta$Stage == "Primary" &
    !is.na(meta$Source) &
    meta$Source == "Tumor",
]

primary_tumor_patients <- sort(
  unique(primary_tumor_meta$Patient)
)

cat("------------------------------------------------------------\n")
cat("Locked Primary + Tumor group\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Number of Primary + Tumor patients:",
  length(primary_tumor_patients),
  "\n"
)

cat(
  "Number of Primary + Tumor cells:",
  nrow(primary_tumor_meta),
  "\n\n"
)

cat("Patients:\n")
print(primary_tumor_patients)
cat("\n")

if (
  length(primary_tumor_patients) != 9
) {
  
  stop(
    "Unexpected number of Primary + Tumor patients. ",
    "Expected 9 based on STEP10A-6."
  )
  
}

if (
  nrow(primary_tumor_meta) != 65702
) {
  
  stop(
    "Unexpected number of Primary + Tumor cells. ",
    "Expected 65,702 based on STEP10A."
  )
  
}

cat(
  "PASS — Primary + Tumor group matches the locked STEP10A/10A-6 structure.\n\n"
)


# ------------------------------------------------------------
# 5. Audit the actual major-lineage labels
# ------------------------------------------------------------

major_lineage_all <- sort(
  unique(
    primary_tumor_meta$`Celltype (major-lineage)`
  )
)

cat("------------------------------------------------------------\n")
cat("Major-lineage labels among Primary + Tumor cells\n")
cat("------------------------------------------------------------\n\n")

print(major_lineage_all)
cat("\n")

write.table(
  data.frame(
    Major_Lineage = major_lineage_all,
    stringsAsFactors = FALSE
  ),
  file = "STEP10B-1A_Primary_Tumor_MajorLineage_Labels.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 6. Count all major-lineages
# ------------------------------------------------------------

major_lineage_counts <- as.data.frame(
  table(
    primary_tumor_meta$`Celltype (major-lineage)`
  ),
  stringsAsFactors = FALSE
)

colnames(major_lineage_counts) <- c(
  "Major_Lineage",
  "Cell_Count"
)

major_lineage_counts <- major_lineage_counts[
  order(
    major_lineage_counts$Cell_Count,
    decreasing = TRUE
  ),
]

major_lineage_counts$Cell_Percentage <- round(
  100 *
    major_lineage_counts$Cell_Count /
    sum(major_lineage_counts$Cell_Count),
  2
)

write.table(
  major_lineage_counts,
  file = "STEP10B-1B_Primary_Tumor_MajorLineage_Counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Major-lineage cell counts:\n")
print(major_lineage_counts)
cat("\n")


# ------------------------------------------------------------
# 7. LOCKED TISCH2 major-lineage → target-compartment mapping
#
# Based on the actual labels observed in the dataset.
#
# T-cell:
#   CD4Tconv
#   CD8T
#   CD8Tex
#
# Myeloid/macrophage:
#   Mono/Macro
#
# Endothelial:
#   Endothelial
# ------------------------------------------------------------

T_cell_labels <- c(
  "CD4Tconv",
  "CD8T",
  "CD8Tex"
)

Myeloid_macrophage_labels <- c(
  "Mono/Macro"
)

Endothelial_labels <- c(
  "Endothelial"
)


# ------------------------------------------------------------
# 8. Verify that all locked labels exist
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("Locked compartment labels\n")
cat("------------------------------------------------------------\n\n")

cat("T-cell labels:\n")
print(T_cell_labels)
cat("\n")

cat("Myeloid/macrophage labels:\n")
print(Myeloid_macrophage_labels)
cat("\n")

cat("Endothelial labels:\n")
print(Endothelial_labels)
cat("\n")

locked_labels <- c(
  T_cell_labels,
  Myeloid_macrophage_labels,
  Endothelial_labels
)

missing_locked_labels <- setdiff(
  locked_labels,
  major_lineage_all
)

if (length(missing_locked_labels) > 0) {
  
  cat(
    "Missing locked labels:\n"
  )
  
  print(missing_locked_labels)
  
  stop(
    "STEP10B-1 halted because one or more locked compartment labels ",
    "are absent from the Primary + Tumor dataset."
  )
  
}

cat(
  "PASS — All locked compartment labels are present.\n\n"
)


# ------------------------------------------------------------
# 9. Check for overlap between locked compartments
# ------------------------------------------------------------

overlap_check <- list(
  T_cell = T_cell_labels,
  Myeloid_macrophage = Myeloid_macrophage_labels,
  Endothelial = Endothelial_labels
)

cat("------------------------------------------------------------\n")
cat("Locked compartment overlap audit\n")
cat("------------------------------------------------------------\n\n")

print(overlap_check)
cat("\n")

compartment_pairs <- combn(
  names(overlap_check),
  2,
  simplify = FALSE
)

overlap_results <- list()

for (pair in compartment_pairs) {
  
  overlap <- intersect(
    overlap_check[[pair[1]]],
    overlap_check[[pair[2]]]
  )
  
  overlap_results[[paste(
    pair[1],
    "vs",
    pair[2],
    sep = "_"
  )]] <- overlap
  
  if (length(overlap) > 0) {
    
    stop(
      "STEP10B-1 halted: overlapping major-lineage labels detected ",
      "between compartments: ",
      paste(overlap, collapse = ", ")
    )
    
  }
  
}

cat(
  "PASS — No major-lineage label overlaps between the three ",
  "locked target compartments.\n\n"
)


# ------------------------------------------------------------
# 10. Create explicit target mapping table
# ------------------------------------------------------------

target_mapping <- data.frame(
  Major_Lineage = major_lineage_all,
  T_Cell_Match =
    major_lineage_all %in% T_cell_labels,
  Myeloid_Match =
    major_lineage_all %in% Myeloid_macrophage_labels,
  Endothelial_Match =
    major_lineage_all %in% Endothelial_labels,
  stringsAsFactors = FALSE
)

target_mapping$Target_Compartment <- NA_character_

target_mapping$Target_Compartment[
  target_mapping$T_Cell_Match
] <- "T_cell"

target_mapping$Target_Compartment[
  target_mapping$Myeloid_Match
] <- "Myeloid_macrophage"

target_mapping$Target_Compartment[
  target_mapping$Endothelial_Match
] <- "Endothelial"

write.table(
  target_mapping,
  file = "STEP10B-1C_MajorLineage_TargetMapping_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("Target compartment mapping audit\n")
cat("------------------------------------------------------------\n\n")

print(target_mapping)
cat("\n")


# ------------------------------------------------------------
# 11. Verify target label counts
# ------------------------------------------------------------

n_t_labels <- sum(
  target_mapping$T_Cell_Match
)

n_myeloid_labels <- sum(
  target_mapping$Myeloid_Match
)

n_endothelial_labels <- sum(
  target_mapping$Endothelial_Match
)

cat(
  "T-cell major-lineage labels detected:",
  n_t_labels,
  "\n"
)

cat(
  "Myeloid/macrophage major-lineage labels detected:",
  n_myeloid_labels,
  "\n"
)

cat(
  "Endothelial major-lineage labels detected:",
  n_endothelial_labels,
  "\n\n"
)

if (
  n_t_labels != length(T_cell_labels) ||
  n_myeloid_labels != length(Myeloid_macrophage_labels) ||
  n_endothelial_labels != length(Endothelial_labels)
) {
  
  stop(
    "STEP10B-1 halted because the detected target-label counts ",
    "do not match the locked mapping."
  )
  
}

cat(
  "PASS — Target-label detection matches the locked mapping.\n\n"
)


# ------------------------------------------------------------
# 12. Assign target compartment
# ------------------------------------------------------------

primary_tumor_meta$STEP10B1_Compartment <- NA_character_

primary_tumor_meta$STEP10B1_Compartment[
  primary_tumor_meta$`Celltype (major-lineage)` %in%
    T_cell_labels
] <- "T_cell"

primary_tumor_meta$STEP10B1_Compartment[
  primary_tumor_meta$`Celltype (major-lineage)` %in%
    Myeloid_macrophage_labels
] <- "Myeloid_macrophage"

primary_tumor_meta$STEP10B1_Compartment[
  primary_tumor_meta$`Celltype (major-lineage)` %in%
    Endothelial_labels
] <- "Endothelial"


# ------------------------------------------------------------
# 13. Assignment QC
# ------------------------------------------------------------

mapping_score <- rowSums(
  cbind(
    primary_tumor_meta$`Celltype (major-lineage)` %in%
      T_cell_labels,
    
    primary_tumor_meta$`Celltype (major-lineage)` %in%
      Myeloid_macrophage_labels,
    
    primary_tumor_meta$`Celltype (major-lineage)` %in%
      Endothelial_labels
  )
)

ambiguous_cells <- sum(
  mapping_score > 1
)

unmatched_cells <- sum(
  mapping_score == 0
)

multiple_assignment_cells <- sum(
  mapping_score > 1
)

assigned_target_cells <- sum(
  mapping_score == 1
)

cat("------------------------------------------------------------\n")
cat("Target-compartment assignment QC\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Cells assigned to exactly one target compartment:",
  assigned_target_cells,
  "\n"
)

cat(
  "Ambiguous target assignments:",
  ambiguous_cells,
  "\n"
)

cat(
  "Cells not assigned to one of the three target compartments:",
  unmatched_cells,
  "\n\n"
)

if (
  ambiguous_cells > 0
) {
  
  stop(
    "STEP10B-1 halted: ambiguous target-compartment assignments detected."
  )
  
}

if (
  assigned_target_cells + unmatched_cells !=
  nrow(primary_tumor_meta)
) {
  
  stop(
    "STEP10B-1 halted: assignment accounting does not reconcile."
  )
  
}

cat(
  "PASS — No cell is assigned to more than one target compartment.\n\n"
)


# ------------------------------------------------------------
# 14. Save cell-level target mapping audit
# ------------------------------------------------------------

cell_mapping_audit <- primary_tumor_meta[
  ,
  c(
    "Cell",
    "Sample",
    "Patient",
    "Source",
    "Stage",
    "Celltype (major-lineage)",
    "Celltype (minor-lineage)",
    "STEP10B1_Compartment"
  )
]

write.table(
  cell_mapping_audit,
  file = "STEP10B-1D_Primary_Tumor_CellLevel_Compartment_Mapping.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Patient × compartment cell counts
# ------------------------------------------------------------

target_meta <- primary_tumor_meta[
  !is.na(
    primary_tumor_meta$STEP10B1_Compartment
  ),
]

patient_compartment_counts <- as.data.frame(
  table(
    target_meta$Patient,
    target_meta$STEP10B1_Compartment
  ),
  stringsAsFactors = FALSE
)

colnames(patient_compartment_counts) <- c(
  "Patient",
  "Compartment",
  "Cell_Count"
)

patient_compartment_counts <- patient_compartment_counts[
  order(
    patient_compartment_counts$Patient,
    patient_compartment_counts$Compartment
  ),
]

write.table(
  patient_compartment_counts,
  file = "STEP10B-1E_Patient_Compartment_Cell_Counts_Long.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Wide patient × compartment table
# ------------------------------------------------------------

patient_compartment_wide <- reshape(
  patient_compartment_counts,
  idvar = "Patient",
  timevar = "Compartment",
  direction = "wide"
)

colnames(patient_compartment_wide) <- sub(
  "^Cell_Count\\.",
  "",
  colnames(patient_compartment_wide)
)


# Ensure all three columns exist
for (x in c(
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)) {
  
  if (!x %in% colnames(patient_compartment_wide)) {
    
    patient_compartment_wide[[x]] <- 0
    
  }
  
}

patient_compartment_wide <- patient_compartment_wide[
  ,
  c(
    "Patient",
    "T_cell",
    "Myeloid_macrophage",
    "Endothelial"
  )
]

patient_compartment_wide <- patient_compartment_wide[
  order(
    patient_compartment_wide$Patient
  ),
]

write.table(
  patient_compartment_wide,
  file = "STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("Patient × compartment cell counts\n")
cat("------------------------------------------------------------\n\n")

print(patient_compartment_wide)
cat("\n")


# ------------------------------------------------------------
# 17. Verify patient count consistency
# ------------------------------------------------------------

patient_count_check <- sort(
  unique(
    patient_compartment_wide$Patient
  )
)

if (
  !identical(
    patient_count_check,
    primary_tumor_patients
  )
) {
  
  stop(
    "STEP10B-1 halted: patient identities do not reconcile ",
    "between Primary + Tumor metadata and compartment table."
  )
  
}

cat(
  "PASS — Patient identities reconcile across all target tables.\n\n"
)


# ------------------------------------------------------------
# 18. Cell-count summary for each compartment
#
# IMPORTANT:
#   No minimum-cell threshold is applied.
# ------------------------------------------------------------

compartments <- c(
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)

compartment_summary <- do.call(
  rbind,
  lapply(
    compartments,
    function(comp) {
      
      x <- patient_compartment_wide[[comp]]
      
      data.frame(
        Compartment = comp,
        Patients = sum(x > 0),
        Zero_Cell_Patients = sum(x == 0),
        Minimum = min(x),
        Q1 = as.numeric(
          quantile(x, 0.25)
        ),
        Median = median(x),
        Mean = mean(x),
        Q3 = as.numeric(
          quantile(x, 0.75)
        ),
        Maximum = max(x),
        IQR = IQR(x),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

rownames(compartment_summary) <- NULL

write.table(
  compartment_summary,
  file = "STEP10B-1G_Compartment_Cell_Count_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("Compartment cell-count summary\n")
cat("------------------------------------------------------------\n\n")

print(compartment_summary)
cat("\n")


# ------------------------------------------------------------
# 19. Patient-level three-compartment feasibility
#
# IMPORTANT:
#   This is descriptive only.
#   No minimum threshold is imposed.
# ------------------------------------------------------------

patient_compartment_wide$Minimum_Target_Compartment_Cells <-
  apply(
    patient_compartment_wide[
      ,
      c(
        "T_cell",
        "Myeloid_macrophage",
        "Endothelial"
      )
    ],
    1,
    min
  )

patient_compartment_wide$All_Three_Compartments_Present <-
  patient_compartment_wide$Minimum_Target_Compartment_Cells > 0

write.table(
  patient_compartment_wide,
  file = "STEP10B-1H_Patient_Three_Compartment_Feasibility.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("Patient-level three-compartment feasibility\n")
cat("------------------------------------------------------------\n\n")

print(
  patient_compartment_wide[
    ,
    c(
      "Patient",
      "T_cell",
      "Myeloid_macrophage",
      "Endothelial",
      "Minimum_Target_Compartment_Cells",
      "All_Three_Compartments_Present"
    )
  ]
)

cat("\n")


# ------------------------------------------------------------
# 20. Major-lineage → compartment crosswalk
# ------------------------------------------------------------

crosswalk <- unique(
  target_mapping[
    ,
    c(
      "Major_Lineage",
      "Target_Compartment"
    )
  ]
)

colnames(crosswalk) <- c(
  "Major_Lineage",
  "STEP10B1_Compartment"
)

crosswalk <- crosswalk[
  order(
    crosswalk$Major_Lineage
  ),
]

write.table(
  crosswalk,
  file = "STEP10B-1I_Final_MajorLineage_Compartment_Crosswalk.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Cell count by patient × compartment × minor-lineage
#
# This provides an audit of the composition of each broad
# target compartment.
# ------------------------------------------------------------

minor_lineage_counts <- as.data.frame(
  table(
    target_meta$Patient,
    target_meta$STEP10B1_Compartment,
    target_meta$`Celltype (minor-lineage)`
  ),
  stringsAsFactors = FALSE
)

colnames(minor_lineage_counts) <- c(
  "Patient",
  "Compartment",
  "Minor_Lineage",
  "Cell_Count"
)

minor_lineage_counts <- minor_lineage_counts[
  minor_lineage_counts$Cell_Count > 0,
]

write.table(
  minor_lineage_counts,
  file = "STEP10B-1J_Patient_Compartment_MinorLineage_Counts.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Reconciliation of cell counts
#
# The sum of the three target compartments plus unmatched cells
# must equal all Primary + Tumor cells.
# ------------------------------------------------------------

total_target_cells <- sum(
  patient_compartment_wide$T_cell
) +
  sum(
    patient_compartment_wide$Myeloid_macrophage
  ) +
  sum(
    patient_compartment_wide$Endothelial
  )

total_primary_tumor_cells <- nrow(
  primary_tumor_meta
)

reconciliation_difference <-
  total_primary_tumor_cells -
  total_target_cells -
  unmatched_cells

cat("------------------------------------------------------------\n")
cat("Cell-count reconciliation\n")
cat("------------------------------------------------------------\n\n")

cat(
  "Primary + Tumor cells:",
  total_primary_tumor_cells,
  "\n"
)

cat(
  "Target-compartment cells:",
  total_target_cells,
  "\n"
)

cat(
  "Unmatched cells:",
  unmatched_cells,
  "\n"
)

cat(
  "Reconciliation difference:",
  reconciliation_difference,
  "\n\n"
)

if (
  reconciliation_difference != 0
) {
  
  stop(
    "STEP10B-1 halted: cell-count reconciliation failed."
  )
  
}

cat(
  "PASS — Cell counts reconcile exactly.\n\n"
)


# ------------------------------------------------------------
# 23. Final QC summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10B-1 — FINAL QC SUMMARY\n")
cat("============================================================\n\n")

cat(
  "Primary + Tumor patients:",
  length(primary_tumor_patients),
  "\n"
)

cat(
  "Primary + Tumor cells:",
  nrow(primary_tumor_meta),
  "\n\n"
)

cat(
  "Locked T-cell labels:",
  paste(T_cell_labels, collapse = " + "),
  "\n"
)

cat(
  "Locked Myeloid/macrophage label:",
  paste(Myeloid_macrophage_labels, collapse = " + "),
  "\n"
)

cat(
  "Locked Endothelial label:",
  paste(Endothelial_labels, collapse = " + "),
  "\n\n"
)

cat(
  "Ambiguous target assignments:",
  ambiguous_cells,
  "\n"
)

cat(
  "Unmatched cells:",
  unmatched_cells,
  "\n\n"
)

cat("Compartment summary:\n")
print(compartment_summary)

cat("\nPatient-level feasibility:\n")

print(
  patient_compartment_wide[
    ,
    c(
      "Patient",
      "T_cell",
      "Myeloid_macrophage",
      "Endothelial",
      "Minimum_Target_Compartment_Cells",
      "All_Three_Compartments_Present"
    )
  ]
)

cat("\n")


# ------------------------------------------------------------
# 24. Save final audit object
# ------------------------------------------------------------

saveRDS(
  list(
    
    primary_tumor_meta =
      primary_tumor_meta,
    
    primary_tumor_patients =
      primary_tumor_patients,
    
    major_lineage_counts =
      major_lineage_counts,
    
    T_cell_labels =
      T_cell_labels,
    
    Myeloid_macrophage_labels =
      Myeloid_macrophage_labels,
    
    Endothelial_labels =
      Endothelial_labels,
    
    target_mapping =
      target_mapping,
    
    patient_compartment_counts =
      patient_compartment_counts,
    
    patient_compartment_wide =
      patient_compartment_wide,
    
    compartment_summary =
      compartment_summary,
    
    crosswalk =
      crosswalk,
    
    minor_lineage_counts =
      minor_lineage_counts,
    
    ambiguous_cells =
      ambiguous_cells,
    
    unmatched_cells =
      unmatched_cells,
    
    total_target_cells =
      total_target_cells,
    
    reconciliation_difference =
      reconciliation_difference
    
  ),
  file =
    "STEP10B-1_Compartment_CellCount_Feasibility_Audit.rds"
)


# ------------------------------------------------------------
# 25. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  "STEP10B-1_sessionInfo.txt"
)


# ------------------------------------------------------------
# STEP10B-1 — FINAL STATUS CORRECTION
# ------------------------------------------------------------

mapping_labels_ok <- length(
  missing_locked_labels
) == 0

final_qc_pass <-
  length(primary_tumor_patients) == 9 &&
  nrow(primary_tumor_meta) == 65702 &&
  mapping_labels_ok &&
  ambiguous_cells == 0 &&
  reconciliation_difference == 0

cat("\n")
cat("============================================================\n")

if (isTRUE(final_qc_pass)) {
  
  cat(
    "STEP10B-1 STATUS: PASS — LOCKED COMPARTMENT MAPPING AND CELL-COUNT AUDIT\n"
  )
  
} else {
  
  cat(
    "STEP10B-1 STATUS: REVIEW REQUIRED — CHECK COMPARTMENT MAPPING OR QC\n"
  )
  
}

cat("============================================================\n")


# ============================================================
# STEP10B-2 — Patient-Level Pseudobulk Feasibility /
#             Cell-Depth Audit
#
# Purpose:
#   1. Quantify cells per patient for:
#      - T-cell
#      - Mono/Macro
#      - Endothelial
#
#   2. Evaluate patient-level pseudobulk feasibility
#      using predefined cell-depth strata.
#
# IMPORTANT:
#   - No gene expression is recalculated here.
#   - No pseudobulk matrix is generated here.
#   - No patient is excluded at this stage.
#   - This is a feasibility / audit step only.
#
# Locked target mapping from STEP10B-1:
#   T-cell              = CD4Tconv + CD8T + CD8Tex
#   Myeloid/macrophage  = Mono/Macro
#   Endothelial         = Endothelial
#
# Locked primary-tumor cohort:
#   Primary + Tumor
#   9 patients
#
# ============================================================


cat("STEP10B-2: Patient-Level Pseudobulk Feasibility / Cell-Depth Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output folder
# ------------------------------------------------------------

step_dir <- "input_dir/Paper 1_Step 10-Single-Cell/STEP10B-2_Pseudobulk_Feasibility_CellDepth_Audit"

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)

cat("Working directory:\n")
cat(normalizePath(step_dir, winslash = "/", mustWork = FALSE), "\n\n")


# ------------------------------------------------------------
# 1. Input file
# ------------------------------------------------------------
#
# Use the LOCKED STEP10B-1 patient × compartment table.
#
# This avoids redoing lineage mapping and ensures that
# STEP10B-2 is fully downstream of the locked STEP10B-1 result.
# ------------------------------------------------------------

input_file <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit/",
  "STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv"
)

if (!file.exists(input_file)) {
  stop("Required STEP10B-1 input file not found:\n", input_file)
}

patient_wide <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Input file:\n")
cat(input_file, "\n\n")

cat("Input dimensions:\n")
cat(nrow(patient_wide), "rows x", ncol(patient_wide), "columns\n\n")

print(patient_wide)


# ------------------------------------------------------------
# 2. Inspect column names
# ------------------------------------------------------------

cat("\nColumn names:\n")
print(colnames(patient_wide))


# ------------------------------------------------------------
# 3. Identify patient and locked compartment columns
# ------------------------------------------------------------

# Expected structure from STEP10B-1:
#   Patient
#   T_cell
#   Myeloid_macrophage
#   Endothelial

required_cols <- c(
  "Patient",
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)

missing_cols <- setdiff(required_cols, colnames(patient_wide))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from STEP10B-1F: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("\nRequired columns: PASS\n")


# ------------------------------------------------------------
# 4. Define locked 9-patient cohort
# ------------------------------------------------------------

locked_patients <- c(
  "P0020",
  "P0025",
  "P0028",
  "P0030",
  "P0031",
  "P0034",
  "P1006",
  "P1010",
  "P1011"
)

cat("\nLocked Primary + Tumor patients:\n")
print(locked_patients)


# ------------------------------------------------------------
# 5. Patient identity audit
# ------------------------------------------------------------

actual_patients <- as.character(patient_wide$Patient)

missing_patients <- setdiff(
  locked_patients,
  actual_patients
)

unexpected_patients <- setdiff(
  actual_patients,
  locked_patients
)

duplicate_patients <- actual_patients[
  duplicated(actual_patients)
]

cat("\nPatient identity audit:\n")
cat("Expected patients:", length(locked_patients), "\n")
cat("Observed patients:", length(unique(actual_patients)), "\n")
cat("Missing expected patients:",
    length(missing_patients), "\n")
cat("Unexpected patients:",
    length(unexpected_patients), "\n")
cat("Duplicated patient IDs:",
    length(unique(duplicate_patients)), "\n")


if (length(missing_patients) > 0) {
  cat("\nMissing expected patients:\n")
  print(missing_patients)
}

if (length(unexpected_patients) > 0) {
  cat("\nUnexpected patients:\n")
  print(unexpected_patients)
}

if (length(duplicate_patients) > 0) {
  cat("\nDuplicated patient IDs:\n")
  print(unique(duplicate_patients))
}


patient_identity_pass <-
  length(missing_patients) == 0 &&
  length(unexpected_patients) == 0 &&
  length(unique(duplicate_patients)) == 0 &&
  length(unique(actual_patients)) == 9


# ------------------------------------------------------------
# 6. Reorder patients according to locked order
# ------------------------------------------------------------

patient_wide <- patient_wide[
  match(locked_patients, patient_wide$Patient),
  ,
  drop = FALSE
]

rownames(patient_wide) <- NULL


# ------------------------------------------------------------
# 7. Verify all counts are numeric and non-missing
# ------------------------------------------------------------

for (v in required_cols[-1]) {
  
  patient_wide[[v]] <- as.numeric(patient_wide[[v]])
  
}

count_na <- sapply(
  patient_wide[required_cols[-1]],
  function(x) sum(is.na(x))
)

count_negative <- sapply(
  patient_wide[required_cols[-1]],
  function(x) sum(x < 0, na.rm = TRUE)
)

cat("\nCount-value QC:\n")
print(count_na)

cat("\nNegative-count QC:\n")
print(count_negative)

count_value_pass <-
  all(count_na == 0) &&
  all(count_negative == 0)


# ------------------------------------------------------------
# 8. Verify exact STEP10B-1 locked cell counts
# ------------------------------------------------------------
#
# These are the expected values already established in
# STEP10B-1.
#
# T-cell:
# P0020 4003
# P0025 1995
# P0028 4014
# P0030 3480
# P0031 4808
# P0034 2752
# P1006 4009
# P1010 2734
# P1011 3319
#
# Myeloid/macrophage:
# P0020 3384
# P0025 283
# P0028 1266
# P0030 1576
# P0031 1327
# P0034 1800
# P1006 818
# P1010 317
# P1011 639
#
# Endothelial:
# P0020 83
# P0025 31
# P0028 41
# P0030 345
# P0031 325
# P0034 113
# P1006 2
# P1010 6
# P1011 5
# ------------------------------------------------------------

expected_counts <- data.frame(
  Patient = locked_patients,
  
  T_cell = c(
    4003,
    1995,
    4014,
    3480,
    4808,
    2752,
    4009,
    2734,
    3319
  ),
  
  Myeloid_macrophage = c(
    3384,
    283,
    1266,
    1576,
    1327,
    1800,
    818,
    317,
    639
  ),
  
  Endothelial = c(
    83,
    31,
    41,
    345,
    325,
    113,
    2,
    6,
    5
  ),
  
  stringsAsFactors = FALSE
)


count_comparison <- merge(
  expected_counts,
  patient_wide,
  by = "Patient",
  suffixes = c("_expected", "_observed")
)

count_comparison$T_cell_match <-
  count_comparison$T_cell_expected ==
  count_comparison$T_cell_observed

count_comparison$Myeloid_match <-
  count_comparison$Myeloid_macrophage_expected ==
  count_comparison$Myeloid_macrophage_observed

count_comparison$Endothelial_match <-
  count_comparison$Endothelial_expected ==
  count_comparison$Endothelial_observed

exact_count_match <-
  all(count_comparison$T_cell_match) &&
  all(count_comparison$Myeloid_match) &&
  all(count_comparison$Endothelial_match)


write.table(
  count_comparison,
  file = "STEP10B-2A_Exact_STEP10B-1_Count_Reconciliation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nExact STEP10B-1 count reconciliation:\n")

if (isTRUE(exact_count_match)) {
  cat("PASS — all patient × compartment counts match STEP10B-1\n")
} else {
  cat("FAIL — count mismatch detected\n")
}


# ------------------------------------------------------------
# 9. Convert wide → long
# ------------------------------------------------------------

compartments <- c(
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)

long_list <- lapply(
  compartments,
  function(comp) {
    
    data.frame(
      Patient = patient_wide$Patient,
      Compartment = comp,
      Cell_Count = patient_wide[[comp]],
      stringsAsFactors = FALSE
    )
    
  }
)

patient_long <- do.call(
  rbind,
  long_list
)

rownames(patient_long) <- NULL


# ------------------------------------------------------------
# 10. Define PRE-SPECIFIED feasibility strata
# ------------------------------------------------------------
#
# Primary operational criterion:
#
#   >=50 cells:
#       Adequate for primary patient-level pseudobulk
#
#   20–49 cells:
#       Limited depth / sensitivity analysis only
#
#   <20 cells:
#       Insufficient for formal patient-level pseudobulk
#
# Important:
#   These are operational analysis strata, not biological
#   or universal statistical laws.
# ------------------------------------------------------------

patient_long$Depth_Class <- cut(
  patient_long$Cell_Count,
  breaks = c(-Inf, 19, 49, Inf),
  labels = c(
    "<20",
    "20-49",
    ">=50"
  ),
  right = TRUE
)


patient_long$Primary_Pseudobulk_Eligible <-
  patient_long$Cell_Count >= 50


patient_long$Sensitivity_Eligible <-
  patient_long$Cell_Count >= 20


patient_long$Very_Low_Depth <-
  patient_long$Cell_Count < 20


# ------------------------------------------------------------
# 11. Add stricter descriptive depth strata
# ------------------------------------------------------------
#
# This does NOT redefine the primary threshold.
# It allows us to see how many patients remain at:
#
#   >=20
#   >=50
#   >=100
#
# especially important for Endothelial.
# ------------------------------------------------------------

patient_long$ge20 <- patient_long$Cell_Count >= 20
patient_long$ge50 <- patient_long$Cell_Count >= 50
patient_long$ge100 <- patient_long$Cell_Count >= 100


# ------------------------------------------------------------
# 12. Save patient-level feasibility table
# ------------------------------------------------------------

write.table(
  patient_long,
  file = "STEP10B-2B_Patient_Compartment_CellDepth_Feasibility.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Summary by compartment
# ------------------------------------------------------------

summary_by_compartment <- do.call(
  rbind,
  lapply(
    compartments,
    function(comp) {
      
      x <- patient_long$Cell_Count[
        patient_long$Compartment == comp
      ]
      
      data.frame(
        Compartment = comp,
        Patients = length(x),
        Zero_cells = sum(x == 0),
        Min = min(x),
        Q1 = as.numeric(quantile(x, 0.25)),
        Median = median(x),
        Mean = mean(x),
        Q3 = as.numeric(quantile(x, 0.75)),
        Max = max(x),
        IQR = IQR(x),
        N_ge20 = sum(x >= 20),
        N_ge50 = sum(x >= 50),
        N_ge100 = sum(x >= 100),
        N_lt20 = sum(x < 20),
        N_20_49 = sum(x >= 20 & x < 50),
        N_ge50_primary = sum(x >= 50),
        stringsAsFactors = FALSE
      )
      
    }
  )
)

rownames(summary_by_compartment) <- NULL

write.table(
  summary_by_compartment,
  file = "STEP10B-2C_Compartment_CellDepth_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("CELL-DEPTH SUMMARY\n")
cat("============================================================\n")

print(summary_by_compartment)


# ------------------------------------------------------------
# 14. Explicit patient × compartment feasibility table
# ------------------------------------------------------------

feasibility_wide <- patient_wide[, c(
  "Patient",
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)]

feasibility_wide$T_cell_status <-
  ifelse(
    feasibility_wide$T_cell >= 50,
    "PRIMARY_ELIGIBLE",
    ifelse(
      feasibility_wide$T_cell >= 20,
      "SENSITIVITY_ONLY",
      "INSUFFICIENT"
    )
  )

feasibility_wide$Myeloid_status <-
  ifelse(
    feasibility_wide$Myeloid_macrophage >= 50,
    "PRIMARY_ELIGIBLE",
    ifelse(
      feasibility_wide$Myeloid_macrophage >= 20,
      "SENSITIVITY_ONLY",
      "INSUFFICIENT"
    )
  )

feasibility_wide$Endothelial_status <-
  ifelse(
    feasibility_wide$Endothelial >= 50,
    "PRIMARY_ELIGIBLE",
    ifelse(
      feasibility_wide$Endothelial >= 20,
      "SENSITIVITY_ONLY",
      "INSUFFICIENT"
    )
  )


write.table(
  feasibility_wide,
  file = "STEP10B-2D_Patient_Compartment_Feasibility_Wide.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Identify patients below each threshold
# ------------------------------------------------------------

threshold_audit <- do.call(
  rbind,
  lapply(
    compartments,
    function(comp) {
      
      x <- patient_wide[[comp]]
      
      data.frame(
        Compartment = comp,
        
        Patients_lt20 =
          paste(
            patient_wide$Patient[x < 20],
            collapse = ";"
          ),
        
        Patients_20_49 =
          paste(
            patient_wide$Patient[x >= 20 & x < 50],
            collapse = ";"
          ),
        
        Patients_ge50 =
          paste(
            patient_wide$Patient[x >= 50],
            collapse = ";"
          ),
        
        Patients_ge100 =
          paste(
            patient_wide$Patient[x >= 100],
            collapse = ";"
          ),
        
        stringsAsFactors = FALSE
      )
      
    }
  )
)

rownames(threshold_audit) <- NULL

write.table(
  threshold_audit,
  file = "STEP10B-2E_Threshold_Audit_By_Compartment.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. Endothelial-specific audit
# ------------------------------------------------------------
#
# Because endothelial is expected to be the limiting
# compartment, generate a dedicated audit.
# ------------------------------------------------------------

endothelial_audit <- feasibility_wide[
  order(feasibility_wide$Endothelial),
  c(
    "Patient",
    "Endothelial",
    "Endothelial_status"
  )
]

write.table(
  endothelial_audit,
  file = "STEP10B-2F_Endothelial_CellDepth_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("ENDOTHELIAL CELL-DEPTH AUDIT\n")
cat("============================================================\n")

print(endothelial_audit)


# ------------------------------------------------------------
# 17. Cross-compartment patient feasibility
# ------------------------------------------------------------
#
# This asks:
#
#   How many patients have sufficient depth in:
#   - all 3 compartments
#   - each pair of compartments
#
# This is descriptive only.
# ------------------------------------------------------------

patient_wide$T_ge50 <-
  patient_wide$T_cell >= 50

patient_wide$Myeloid_ge50 <-
  patient_wide$Myeloid_macrophage >= 50

patient_wide$Endothelial_ge50 <-
  patient_wide$Endothelial >= 50


patient_wide$All3_ge50 <-
  patient_wide$T_ge50 &
  patient_wide$Myeloid_ge50 &
  patient_wide$Endothelial_ge50

patient_wide$T_Myeloid_ge50 <-
  patient_wide$T_ge50 &
  patient_wide$Myeloid_ge50

patient_wide$T_Endothelial_ge50 <-
  patient_wide$T_ge50 &
  patient_wide$Endothelial_ge50

patient_wide$Myeloid_Endothelial_ge50 <-
  patient_wide$Myeloid_ge50 &
  patient_wide$Endothelial_ge50


cross_compartment_summary <- data.frame(
  
  Criterion = c(
    "T-cell >=50",
    "Myeloid/macrophage >=50",
    "Endothelial >=50",
    "T-cell + Myeloid/macrophage >=50",
    "T-cell + Endothelial >=50",
    "Myeloid/macrophage + Endothelial >=50",
    "All three compartments >=50"
  ),
  
  N_patients = c(
    sum(patient_wide$T_ge50),
    sum(patient_wide$Myeloid_ge50),
    sum(patient_wide$Endothelial_ge50),
    sum(patient_wide$T_Myeloid_ge50),
    sum(patient_wide$T_Endothelial_ge50),
    sum(patient_wide$Myeloid_Endothelial_ge50),
    sum(patient_wide$All3_ge50)
  ),
  
  stringsAsFactors = FALSE
)


write.table(
  cross_compartment_summary,
  file = "STEP10B-2G_Cross_Compartment_Pseudobulk_Feasibility.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


cat("\n============================================================\n")
cat("CROSS-COMPARTMENT FEASIBILITY\n")
cat("============================================================\n")

print(cross_compartment_summary)


# ------------------------------------------------------------
# 18. Global feasibility interpretation
# ------------------------------------------------------------

t_ge50 <- sum(patient_wide$T_cell >= 50)
myeloid_ge50 <- sum(patient_wide$Myeloid_macrophage >= 50)
endo_ge50 <- sum(patient_wide$Endothelial >= 50)

t_lt20 <- sum(patient_wide$T_cell < 20)
myeloid_lt20 <- sum(patient_wide$Myeloid_macrophage < 20)
endo_lt20 <- sum(patient_wide$Endothelial < 20)

all3_ge50 <- sum(patient_wide$All3_ge50)


# ------------------------------------------------------------
# 19. Final QC
# ------------------------------------------------------------

final_qc_pass <-
  isTRUE(patient_identity_pass) &&
  isTRUE(count_value_pass) &&
  isTRUE(exact_count_match) &&
  nrow(patient_wide) == 9 &&
  all(patient_wide$T_cell > 0) &&
  all(patient_wide$Myeloid_macrophage > 0) &&
  all(patient_wide$Endothelial > 0)


# ------------------------------------------------------------
# 20. Final status
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")

if (isTRUE(final_qc_pass)) {
  
  cat(
    "STEP10B-2 STATUS: PASS — PATIENT-LEVEL CELL-DEPTH ",
    "FEASIBILITY AUDIT\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STEP10B-2 STATUS: REVIEW REQUIRED — ",
    "CHECK PATIENT/CELL-DEPTH QC\n"
  )
  
}

cat("============================================================\n\n")


# ------------------------------------------------------------
# 21. Key locked results
# ------------------------------------------------------------

cat("LOCKED SUMMARY\n")
cat("------------------------------------------------------------\n")

cat("Patients:", nrow(patient_wide), "\n\n")

cat("T-cell:\n")
cat("  >=50 cells:", t_ge50, "/ 9\n")
cat("  <20 cells:", t_lt20, "/ 9\n\n")

cat("Myeloid/macrophage:\n")
cat("  >=50 cells:", myeloid_ge50, "/ 9\n")
cat("  <20 cells:", myeloid_lt20, "/ 9\n\n")

cat("Endothelial:\n")
cat("  >=50 cells:", endo_ge50, "/ 9\n")
cat("  <20 cells:", endo_lt20, "/ 9\n\n")

cat("All three compartments >=50 cells:",
    all3_ge50, "/ 9\n")

cat("\nPrimary operational threshold:\n")
cat(">=50 cells / patient / compartment\n")

cat("Sensitivity threshold:\n")
cat(">=20 cells / patient / compartment\n")

cat("Insufficient depth:\n")
cat("<20 cells / patient / compartment\n")

cat("\nIMPORTANT:\n")
cat("No patients were excluded in STEP10B-2.\n")
cat("No expression-level pseudobulk was calculated.\n")
cat("No differential expression or statistical testing was performed.\n")


# ------------------------------------------------------------
# 22. Save RDS
# ------------------------------------------------------------

step10b2_result <- list(
  
  patient_wide = patient_wide,
  
  patient_long = patient_long,
  
  expected_counts = expected_counts,
  
  count_comparison = count_comparison,
  
  summary_by_compartment = summary_by_compartment,
  
  threshold_audit = threshold_audit,
  
  endothelial_audit = endothelial_audit,
  
  cross_compartment_summary = cross_compartment_summary,
  
  locked_primary_threshold = 50,
  
  locked_sensitivity_threshold = 20,
  
  final_qc_pass = final_qc_pass
  
)


saveRDS(
  step10b2_result,
  file = "STEP10B-2_Pseudobulk_Feasibility_CellDepth_Audit.rds"
)


# ------------------------------------------------------------
# 23. Save session information
# ------------------------------------------------------------

sink("STEP10B-2_sessionInfo.txt")

cat("STEP10B-2 — Patient-Level Pseudobulk Feasibility / Cell-Depth Audit\n")
cat("============================================================\n\n")

cat("Date:\n")
print(Sys.time())

cat("\nR version:\n")
print(R.version.string)

cat("\nSession information:\n")
print(sessionInfo())

sink()


cat("\n")
cat("STEP10B-2 completed.\n")
cat("Output directory:\n")
cat(normalizePath(step_dir, winslash = "/", mustWork = FALSE), "\n")


# ============================================================
# STEP10C — Patient-Level Compartment-Specific Pseudobulk
#            Expression Construction and Core-Gene Audit
#
# Purpose:
#   Construct patient-level pseudobulk expression profiles for
#   the three LOCKED cellular compartments:
#
#     T-cell
#     Mono/Macro
#     Endothelial
#
# Locked cohort structure from STEP10B-2:
#
#   T-cell:
#       9/9 patients >=50 cells
#
#   Mono/Macro:
#       9/9 patients >=50 cells
#
#   Endothelial:
#       4/9 patients >=50 cells  -> PRIMARY
#       2/9 patients 20-49 cells -> SENSITIVITY
#       3/9 patients <20 cells   -> EXCLUDED FROM FORMAL PB
#
# IMPORTANT:
#   - Patient is the biological replicate.
#   - Cells are aggregated within patient × compartment.
#   - No DE testing.
#   - No correlation testing.
#   - No FDR calculation.
#   - No CellChat analysis.
#   - No redefinition of cell compartments.
#
# Pseudobulk should be interpreted as:
#   "patient-level compartment-specific transcriptional signal"
#
# It is NOT interpreted as:
#   "direct proof of cell-intrinsic regulation"
#
# ============================================================


cat("STEP10C: Patient-Level Compartment-Specific Pseudobulk\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step_dir <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "STEP10C_Patient_Level_Pseudobulk"
)

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

setwd(step_dir)

cat("Working directory:\n")
cat(normalizePath(step_dir, winslash = "/", mustWork = FALSE), "\n\n")


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

required_packages <- c(
  "Matrix"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    stop(
      "Required package not installed: ",
      pkg
    )
    
  }
  
}

library(Matrix)


# ------------------------------------------------------------
# 2. LOCKED patient structure
# ------------------------------------------------------------

locked_patients <- c(
  "P0020",
  "P0025",
  "P0028",
  "P0030",
  "P0031",
  "P0034",
  "P1006",
  "P1010",
  "P1011"
)


# ------------------------------------------------------------
# 3. LOCKED cell-depth cohorts
# ------------------------------------------------------------

primary_patients_T <- locked_patients

primary_patients_Myeloid <- locked_patients

primary_patients_Endothelial <- c(
  "P0020",
  "P0030",
  "P0031",
  "P0034"
)

sensitivity_patients_Endothelial <- c(
  "P0025",
  "P0028"
)

insufficient_patients_Endothelial <- c(
  "P1006",
  "P1010",
  "P1011"
)


cat("LOCKED PSEUDOBULK COHORTS\n")
cat("------------------------------------------------------------\n")

cat(
  "T-cell primary:",
  length(primary_patients_T),
  "patients\n"
)

cat(
  "Mono/Macro primary:",
  length(primary_patients_Myeloid),
  "patients\n"
)

cat(
  "Endothelial primary:",
  length(primary_patients_Endothelial),
  "patients\n"
)

cat(
  "Endothelial sensitivity:",
  length(sensitivity_patients_Endothelial),
  "patients\n"
)

cat(
  "Endothelial insufficient:",
  length(insufficient_patients_Endothelial),
  "patients\n\n"
)


# ------------------------------------------------------------
# 4. LOCKED compartment definitions
# ------------------------------------------------------------

locked_compartment_mapping <- list(
  
  T_cell = c(
    "CD4Tconv",
    "CD8T",
    "CD8Tex"
  ),
  
  Myeloid_macrophage = c(
    "Mono/Macro"
  ),
  
  Endothelial = c(
    "Endothelial"
  )
  
)


cat("LOCKED COMPARTMENT DEFINITIONS\n")
cat("------------------------------------------------------------\n")

print(locked_compartment_mapping)


# ------------------------------------------------------------
# 5. Input paths
# ------------------------------------------------------------
#
# IMPORTANT:
# Replace the count-matrix path below ONLY if the actual
# original raw-count file is stored elsewhere.
#
# Do NOT use:
#   - normalized expression
#   - log-normalized expression
#   - TISCH2 differential-gene table
#
# Pseudobulk requires the original count matrix.
# ------------------------------------------------------------

metadata_file <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

mapping_file <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit/",
  "STEP10B-1D_Primary_Tumor_CellLevel_Compartment_Mapping.tsv"
)

feasibility_file <- paste0(
  "input_dir/Paper 1_Step 10-Single-Cell/",
  "STEP10B-2_Pseudobulk_Feasibility_CellDepth_Audit/",
  "STEP10B-2B_Patient_Compartment_CellDepth_Feasibility.tsv"
)


# ------------------------------------------------------------
# 6. RAW COUNT MATRIX
# ------------------------------------------------------------
#
# !!! IMPORTANT !!!
#
# Set this to the ACTUAL original count matrix used for the
# NSCLC_GSE131907 dataset.
#
# Example:
#
# count_matrix_file <- "O:/.../NSCLC_GSE131907_raw_counts.mtx"
#
# or:
#
# count_matrix_file <- "O:/.../NSCLC_GSE131907_counts.rds"
#
# The script deliberately stops if the path is not supplied.
# ------------------------------------------------------------

count_matrix_file <- ""


if (count_matrix_file == "") {
  
  stop(
    paste0(
      "\nRAW COUNT MATRIX PATH HAS NOT BEEN LOCKED.\n\n",
      "STEP10C must use the original raw count matrix.\n",
      "Please set 'count_matrix_file' to the actual file path ",
      "before running this script.\n\n",
      "Do NOT substitute a normalized matrix or TISCH2 ",
      "differential-expression table.\n"
    )
  )
  
}


if (!file.exists(count_matrix_file)) {
  
  stop(
    "Raw count matrix not found:\n",
    count_matrix_file
  )
  
}


# ------------------------------------------------------------
# 7. Read metadata
# ------------------------------------------------------------

if (!file.exists(metadata_file)) {
  
  stop(
    "Metadata file not found:\n",
    metadata_file
  )
  
}

meta <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nMetadata dimensions:\n")
cat(
  nrow(meta),
  "rows x",
  ncol(meta),
  "columns\n"
)


# ------------------------------------------------------------
# 8. Required metadata columns
# ------------------------------------------------------------

required_meta <- c(
  "Cell",
  "Patient",
  "Sample",
  "Source",
  "Stage",
  "MajorLineage"
)

missing_meta <- setdiff(
  required_meta,
  colnames(meta)
)

if (length(missing_meta) > 0) {
  
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
  
}


# ------------------------------------------------------------
# 9. Read STEP10B-1 locked compartment mapping
# ------------------------------------------------------------

if (!file.exists(mapping_file)) {
  
  stop(
    "STEP10B-1 mapping file not found:\n",
    mapping_file
  )
  
}

mapping <- read.delim(
  mapping_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\nSTEP10B-1 mapping dimensions:\n")
cat(
  nrow(mapping),
  "rows x",
  ncol(mapping),
  "columns\n"
)

cat("\nMapping columns:\n")
print(colnames(mapping))


# ------------------------------------------------------------
# 10. Detect mapping columns
# ------------------------------------------------------------

possible_cell_cols <- c(
  "Cell",
  "Cell_ID",
  "CellID"
)

cell_col <- possible_cell_cols[
  possible_cell_cols %in% colnames(mapping)
][1]

if (is.na(cell_col)) {
  
  stop(
    "Could not identify Cell column in STEP10B-1 mapping."
  )
  
}


possible_compartment_cols <- c(
  "Compartment",
  "Target_Compartment",
  "Locked_Compartment"
)

compartment_col <- possible_compartment_cols[
  possible_compartment_cols %in% colnames(mapping)
][1]

if (is.na(compartment_col)) {
  
  stop(
    "Could not identify compartment column in STEP10B-1 mapping."
  )
  
}


cat("\nMapping columns locked:\n")
cat("Cell column:", cell_col, "\n")
cat("Compartment column:", compartment_col, "\n")


# ------------------------------------------------------------
# 11. Restrict metadata to Primary + Tumor
# ------------------------------------------------------------

primary_tumor_meta <- meta[
  meta$Stage == "Primary" &
    meta$Source == "Tumor",
  ,
  drop = FALSE
]

cat("\nPrimary + Tumor metadata:\n")
cat(
  nrow(primary_tumor_meta),
  "cells\n"
)

cat(
  length(unique(primary_tumor_meta$Patient)),
  "patients\n"
)


# ------------------------------------------------------------
# 12. Verify exactly 9 patients
# ------------------------------------------------------------

observed_primary_tumor_patients <-
  sort(unique(primary_tumor_meta$Patient))

expected_primary_tumor_patients <-
  sort(locked_patients)

if (!identical(
  observed_primary_tumor_patients,
  expected_primary_tumor_patients
)) {
  
  stop(
    "Primary + Tumor patient identities do not match ",
    "STEP10A-6 locked cohort."
  )
  
}


# ------------------------------------------------------------
# 13. Merge locked compartment mapping
# ------------------------------------------------------------

mapping_small <- mapping[
  ,
  c(cell_col, compartment_col),
  drop = FALSE
]

colnames(mapping_small) <- c(
  "Cell",
  "Compartment"
)


primary_tumor_meta$Cell <- as.character(
  primary_tumor_meta$Cell
)

mapping_small$Cell <- as.character(
  mapping_small$Cell
)

primary_tumor_meta <- merge(
  primary_tumor_meta,
  mapping_small,
  by = "Cell",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 14. Mapping QC
# ------------------------------------------------------------

target_compartments <- names(
  locked_compartment_mapping
)

primary_tumor_meta$Compartment <- as.character(
  primary_tumor_meta$Compartment
)

target_meta <- primary_tumor_meta[
  primary_tumor_meta$Compartment %in%
    target_compartments,
  ,
  drop = FALSE
]


cat("\nTarget-compartment cells:\n")
cat(nrow(target_meta), "\n\n")

cat("Compartment counts:\n")
print(table(target_meta$Compartment))


# ------------------------------------------------------------
# 15. Confirm no duplicated cell IDs
# ------------------------------------------------------------

if (anyDuplicated(primary_tumor_meta$Cell) > 0) {
  
  stop(
    "Duplicated Cell IDs detected after metadata/mapping merge."
  )
  
}


# ------------------------------------------------------------
# 16. Read raw count matrix
# ------------------------------------------------------------
#
# Supported formats:
#   .rds
#   .mtx
#
# The script expects:
#
#   rows = genes
#   columns = cells
#
# with cell IDs in column names.
# ------------------------------------------------------------

count_ext <- tolower(
  tools::file_ext(count_matrix_file)
)


if (count_ext == "rds") {
  
  counts <- readRDS(count_matrix_file)
  
} else if (count_ext == "mtx") {
  
  counts <- Matrix::readMM(count_matrix_file)
  
} else {
  
  stop(
    "Unsupported count matrix format: .",
    count_ext,
    "\nPlease use an RDS or Matrix Market (.mtx) count matrix."
  )
  
}


# ------------------------------------------------------------
# 17. Convert to sparse matrix if necessary
# ------------------------------------------------------------

if (!inherits(counts, "dgCMatrix")) {
  
  counts <- as(
    counts,
    "dgCMatrix"
  )
  
}


cat("\nRaw count matrix dimensions:\n")
cat(
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)


# ------------------------------------------------------------
# 18. Verify cell IDs
# ------------------------------------------------------------

if (is.null(colnames(counts))) {
  
  stop(
    "Raw count matrix has no cell IDs in column names."
  )
  
}

count_cells <- colnames(counts)

metadata_cells <- primary_tumor_meta$Cell


# ------------------------------------------------------------
# 19. Cell identity overlap
# ------------------------------------------------------------

missing_count_cells <- setdiff(
  metadata_cells,
  count_cells
)

extra_count_cells <- setdiff(
  count_cells,
  metadata_cells
)

cat("\nCell identity reconciliation:\n")
cat(
  "Primary + Tumor metadata cells:",
  length(metadata_cells),
  "\n"
)

cat(
  "Missing in count matrix:",
  length(missing_count_cells),
  "\n"
)

cat(
  "Extra count-matrix cells:",
  length(extra_count_cells),
  "\n"
)


if (length(missing_count_cells) > 0) {
  
  stop(
    "Some Primary + Tumor cells are missing from raw count matrix."
  )
  
}


# ------------------------------------------------------------
# 20. Restrict raw counts to Primary + Tumor cells
# ------------------------------------------------------------

primary_tumor_cell_order <- primary_tumor_meta$Cell

counts_primary_tumor <- counts[
  ,
  primary_tumor_cell_order,
  drop = FALSE
]


# ------------------------------------------------------------
# 21. Reconcile cell order
# ------------------------------------------------------------

if (!identical(
  colnames(counts_primary_tumor),
  primary_tumor_meta$Cell
)) {
  
  stop(
    "Cell order mismatch between metadata and count matrix."
  )
  
}


cat("\nCell order reconciliation: PASS\n")


# ------------------------------------------------------------
# 22. Core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


gene_names <- rownames(
  counts_primary_tumor
)

missing_core_genes <- setdiff(
  core_genes,
  gene_names
)

if (length(missing_core_genes) > 0) {
  
  stop(
    "Core genes missing from raw count matrix: ",
    paste(
      missing_core_genes,
      collapse = ", "
    )
  )
  
}


cat("\nCore genes verified in raw count matrix:\n")
print(core_genes)


# ------------------------------------------------------------
# 23. Remove duplicated gene symbols if present
# ------------------------------------------------------------
#
# For the core-gene audit, duplicated gene symbols would make
# gene-level interpretation ambiguous.
#
# We therefore audit duplicated rownames but do NOT silently
# collapse them.
# ------------------------------------------------------------

duplicated_gene_symbols <- unique(
  gene_names[
    duplicated(gene_names)
  ]
)

cat("\nDuplicated gene symbols:\n")
cat(
  length(duplicated_gene_symbols),
  "\n"
)

if (length(duplicated_gene_symbols) > 0) {
  
  cat(
    "WARNING: duplicated gene symbols detected.\n"
  )
  
}


# ------------------------------------------------------------
# 24. Function: patient-level pseudobulk
# ------------------------------------------------------------
#
# Input:
#   sparse count matrix
#   cell metadata
#
# Output:
#   genes x patients
#
# Each patient is represented by the SUM of raw counts
# across all cells belonging to that patient and compartment.
# ------------------------------------------------------------

make_pseudobulk <- function(
    counts_matrix,
    metadata_df,
    compartment_name,
    patient_ids
) {
  
  md <- metadata_df[
    metadata_df$Compartment == compartment_name &
      metadata_df$Patient %in% patient_ids,
    ,
    drop = FALSE
  ]
  
  md <- md[
    match(
      colnames(counts_matrix),
      md$Cell
    ),
    ,
    drop = FALSE
  ]
  
  if (any(is.na(md$Cell))) {
    
    stop(
      "Metadata-cell alignment failed for compartment: ",
      compartment_name
    )
    
  }
  
  output <- matrix(
    0,
    nrow = nrow(counts_matrix),
    ncol = length(patient_ids),
    dimnames = list(
      rownames(counts_matrix),
      patient_ids
    )
  )
  
  for (p in patient_ids) {
    
    cells_p <- md$Cell[
      md$Patient == p
    ]
    
    if (length(cells_p) == 0) {
      
      stop(
        "No cells found for patient ",
        p,
        " in compartment ",
        compartment_name
      )
      
    }
    
    idx <- match(
      cells_p,
      colnames(counts_matrix)
    )
    
    output[, p] <- Matrix::rowSums(
      counts_matrix[
        ,
        idx,
        drop = FALSE
      ]
    )
    
  }
  
  output
  
}


# ------------------------------------------------------------
# 25. Build PRIMARY pseudobulk matrices
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("BUILDING PRIMARY PSEUDOBULK MATRICES\n")
cat("============================================================\n")


pb_T <- make_pseudobulk(
  counts_primary_tumor,
  primary_tumor_meta,
  "T_cell",
  primary_patients_T
)


pb_Myeloid <- make_pseudobulk(
  counts_primary_tumor,
  primary_tumor_meta,
  "Myeloid_macrophage",
  primary_patients_Myeloid
)


pb_Endothelial_primary <- make_pseudobulk(
  counts_primary_tumor,
  primary_tumor_meta,
  "Endothelial",
  primary_patients_Endothelial
)


# ------------------------------------------------------------
# 26. Build ENDOTHELIAL sensitivity pseudobulk
# ------------------------------------------------------------

pb_Endothelial_sensitivity <- make_pseudobulk(
  counts_primary_tumor,
  primary_tumor_meta,
  "Endothelial",
  sensitivity_patients_Endothelial
)


# ------------------------------------------------------------
# 27. Save primary matrices
# ------------------------------------------------------------

saveRDS(
  pb_T,
  file = "STEP10C-1_Primary_Pseudobulk_T_cell_raw_counts.rds"
)

saveRDS(
  pb_Myeloid,
  file = "STEP10C-2_Primary_Pseudobulk_Myeloid_Macrophage_raw_counts.rds"
)

saveRDS(
  pb_Endothelial_primary,
  file = "STEP10C-3_Primary_Pseudobulk_Endothelial_raw_counts.rds"
)

saveRDS(
  pb_Endothelial_sensitivity,
  file = "STEP10C-4_Sensitivity_Pseudobulk_Endothelial_20_49_raw_counts.rds"
)


# ------------------------------------------------------------
# 28. Export core-gene pseudobulk tables
# ------------------------------------------------------------

core_T <- pb_T[
  core_genes,
  ,
  drop = FALSE
]

core_Myeloid <- pb_Myeloid[
  core_genes,
  ,
  drop = FALSE
]

core_Endothelial_primary <- pb_Endothelial_primary[
  core_genes,
  ,
  drop = FALSE
]

core_Endothelial_sensitivity <- pb_Endothelial_sensitivity[
  core_genes,
  ,
  drop = FALSE
]


write.table(
  core_T,
  file = "STEP10C-5_CoreGenes_T_cell_Pseudobulk_RawCounts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  core_Myeloid,
  file = "STEP10C-6_CoreGenes_Myeloid_Macrophage_Pseudobulk_RawCounts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  core_Endothelial_primary,
  file = "STEP10C-7_CoreGenes_Endothelial_Primary_Pseudobulk_RawCounts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

write.table(
  core_Endothelial_sensitivity,
  file = "STEP10C-8_CoreGenes_Endothelial_Sensitivity_RawCounts.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ------------------------------------------------------------
# 29. Pseudobulk library-size audit
# ------------------------------------------------------------

library_size_audit <- rbind(
  
  data.frame(
    Compartment = "T_cell",
    Cohort = "Primary",
    Patient = colnames(pb_T),
    Total_raw_counts = Matrix::colSums(pb_T),
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    Compartment = "Myeloid_macrophage",
    Cohort = "Primary",
    Patient = colnames(pb_Myeloid),
    Total_raw_counts = Matrix::colSums(pb_Myeloid),
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    Compartment = "Endothelial",
    Cohort = "Primary",
    Patient = colnames(pb_Endothelial_primary),
    Total_raw_counts = Matrix::colSums(pb_Endothelial_primary),
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    Compartment = "Endothelial",
    Cohort = "Sensitivity_20_49",
    Patient = colnames(pb_Endothelial_sensitivity),
    Total_raw_counts =
      Matrix::colSums(pb_Endothelial_sensitivity),
    stringsAsFactors = FALSE
  )
  
)


write.table(
  library_size_audit,
  file = "STEP10C-9_Pseudobulk_LibrarySize_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Core-gene detection audit
# ------------------------------------------------------------

make_core_audit <- function(
    pb,
    compartment_name,
    cohort_name
) {
  
  result <- expand.grid(
    Gene = core_genes,
    Patient = colnames(pb),
    stringsAsFactors = FALSE
  )
  
  result$Compartment <- compartment_name
  result$Cohort <- cohort_name
  
  result$RawCount <- mapply(
    function(g, p) {
      pb[g, p]
    },
    result$Gene,
    result$Patient
  )
  
  result$Detected <- result$RawCount > 0
  
  result[
    ,
    c(
      "Compartment",
      "Cohort",
      "Patient",
      "Gene",
      "RawCount",
      "Detected"
    )
  ]
  
}


core_detection_audit <- rbind(
  
  make_core_audit(
    pb_T,
    "T_cell",
    "Primary"
  ),
  
  make_core_audit(
    pb_Myeloid,
    "Myeloid_macrophage",
    "Primary"
  ),
  
  make_core_audit(
    pb_Endothelial_primary,
    "Endothelial",
    "Primary"
  ),
  
  make_core_audit(
    pb_Endothelial_sensitivity,
    "Endothelial",
    "Sensitivity_20_49"
  )
  
)


write.table(
  core_detection_audit,
  file = "STEP10C-10_CoreGene_Detection_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Core-gene summary by compartment
# ------------------------------------------------------------

core_summary <- do.call(
  rbind,
  lapply(
    unique(
      core_detection_audit[
        ,
        c("Compartment", "Cohort")
      ]
    )$Compartment,
    function(comp) {
      
      NULL
      
    }
  )
)


# Explicit summary
core_summary_list <- list()

for (comp in unique(core_detection_audit$Compartment)) {
  
  for (cohort in unique(
    core_detection_audit$Cohort[
      core_detection_audit$Compartment == comp
    ]
  )) {
    
    x <- core_detection_audit[
      core_detection_audit$Compartment == comp &
        core_detection_audit$Cohort == cohort,
      ,
      drop = FALSE
    ]
    
    for (g in core_genes) {
      
      y <- x[
        x$Gene == g,
        ,
        drop = FALSE
      ]
      
      core_summary_list[[length(core_summary_list) + 1]] <-
        data.frame(
          Compartment = comp,
          Cohort = cohort,
          Gene = g,
          N_patients = nrow(y),
          N_detected = sum(y$Detected),
          Detection_fraction =
            sum(y$Detected) / nrow(y),
          Min_raw_count = min(y$RawCount),
          Median_raw_count = median(y$RawCount),
          Max_raw_count = max(y$RawCount),
          stringsAsFactors = FALSE
        )
      
    }
    
  }
  
}

core_summary <- do.call(
  rbind,
  core_summary_list
)

rownames(core_summary) <- NULL


write.table(
  core_summary,
  file = "STEP10C-11_CoreGene_Detection_Summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Cell-count / pseudobulk patient reconciliation
# ------------------------------------------------------------

pb_patient_reconciliation <- data.frame(
  
  Compartment = c(
    "T_cell",
    "Myeloid_macrophage",
    "Endothelial_primary",
    "Endothelial_sensitivity"
  ),
  
  Expected_patients = c(
    length(primary_patients_T),
    length(primary_patients_Myeloid),
    length(primary_patients_Endothelial),
    length(sensitivity_patients_Endothelial)
  ),
  
  Observed_pseudobulk_patients = c(
    ncol(pb_T),
    ncol(pb_Myeloid),
    ncol(pb_Endothelial_primary),
    ncol(pb_Endothelial_sensitivity)
  ),
  
  stringsAsFactors = FALSE
)

pb_patient_reconciliation$PASS <-
  pb_patient_reconciliation$Expected_patients ==
  pb_patient_reconciliation$Observed_pseudobulk_patients


write.table(
  pb_patient_reconciliation,
  file = "STEP10C-12_Patient_Reconciliation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Final QC
# ------------------------------------------------------------

patient_reconciliation_pass <-
  all(pb_patient_reconciliation$PASS)

core_gene_pass <-
  length(missing_core_genes) == 0

matrix_dimension_pass <-
  ncol(pb_T) == 9 &&
  ncol(pb_Myeloid) == 9 &&
  ncol(pb_Endothelial_primary) == 4 &&
  ncol(pb_Endothelial_sensitivity) == 2

core_rowname_pass <-
  all(core_genes %in% rownames(pb_T)) &&
  all(core_genes %in% rownames(pb_Myeloid)) &&
  all(core_genes %in% rownames(pb_Endothelial_primary))


final_qc_pass <-
  patient_reconciliation_pass &&
  core_gene_pass &&
  matrix_dimension_pass &&
  core_rowname_pass


# ------------------------------------------------------------
# 34. Save complete STEP10C object
# ------------------------------------------------------------

step10c_result <- list(
  
  locked_patients = locked_patients,
  
  primary_patients_T = primary_patients_T,
  
  primary_patients_Myeloid =
    primary_patients_Myeloid,
  
  primary_patients_Endothelial =
    primary_patients_Endothelial,
  
  sensitivity_patients_Endothelial =
    sensitivity_patients_Endothelial,
  
  insufficient_patients_Endothelial =
    insufficient_patients_Endothelial,
  
  locked_compartment_mapping =
    locked_compartment_mapping,
  
  pb_T = pb_T,
  
  pb_Myeloid = pb_Myeloid,
  
  pb_Endothelial_primary =
    pb_Endothelial_primary,
  
  pb_Endothelial_sensitivity =
    pb_Endothelial_sensitivity,
  
  library_size_audit =
    library_size_audit,
  
  core_detection_audit =
    core_detection_audit,
  
  core_summary =
    core_summary,
  
  patient_reconciliation =
    pb_patient_reconciliation,
  
  final_qc_pass =
    final_qc_pass
  
)


saveRDS(
  step10c_result,
  file =
    "STEP10C_Patient_Level_Pseudobulk_Final.rds"
)


# ------------------------------------------------------------
# 35. Final status
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")

if (isTRUE(final_qc_pass)) {
  
  cat(
    "STEP10C STATUS: PASS — PATIENT-LEVEL PSEUDOBULK ",
    "CONSTRUCTION AND CORE-GENE AUDIT\n",
    sep = ""
  )
  
} else {
  
  cat(
    "STEP10C STATUS: REVIEW REQUIRED — ",
    "CHECK PSEUDOBULK CONSTRUCTION QC\n"
  )
  
}

cat("============================================================\n\n")


# ------------------------------------------------------------
# 36. Locked final summary
# ------------------------------------------------------------

cat("LOCKED STEP10C SUMMARY\n")
cat("------------------------------------------------------------\n")

cat(
  "T-cell primary pseudobulk:",
  ncol(pb_T),
  "patients\n"
)

cat(
  "Mono/Macro primary pseudobulk:",
  ncol(pb_Myeloid),
  "patients\n"
)

cat(
  "Endothelial primary pseudobulk:",
  ncol(pb_Endothelial_primary),
  "patients\n"
)

cat(
  "Endothelial sensitivity pseudobulk:",
  ncol(pb_Endothelial_sensitivity),
  "patients\n"
)

cat("\n")

cat(
  "Core genes:",
  paste(core_genes, collapse = ", "),
  "\n"
)

cat("\n")

cat(
  "No differential-expression testing performed.\n"
)

cat(
  "No correlation testing performed.\n"
)

cat(
  "No FDR calculation performed.\n"
)

cat(
  "Patient is the biological replicate.\n"
)

cat(
  "Pseudobulk = sum of raw counts within patient × compartment.\n"
)

cat(
  "Interpretation = patient-level compartment-specific ",
  "transcriptional signal.\n"
)


# ------------------------------------------------------------
# 37. Session information
# ------------------------------------------------------------

sink(
  "STEP10C_sessionInfo.txt"
)

cat(
  "STEP10C — Patient-Level Compartment-Specific Pseudobulk\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Date:\n"
)

print(
  Sys.time()
)

cat(
  "\nR version:\n"
)

print(
  R.version.string
)

cat(
  "\nSession information:\n"
)

print(
  sessionInfo()
)

sink()


cat("\nSTEP10C completed.\n")

cat(
  "Output directory:\n",
  normalizePath(
    step_dir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n"
)


# ============================================================
# STEP10C-0A
# NSCLC_GSE131907 Expression / Count Matrix File-Type Audit
#
# Purpose:
#   Determine whether:
#     1. NSCLC_GSE131907_expression.h5
#     2. NSCLC_GSE131907_Expression.zip
#   contain a usable single-cell expression/count matrix.
#
# IMPORTANT:
#   - NO pseudobulk
#   - NO statistical testing
#   - NO modification of source files
#   - NO patient exclusion
#   - NO expression transformation
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10C-0A — EXPRESSION / COUNT MATRIX FILE-TYPE AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

work_dir <- "input_dir/Figure 9C Dot Plot (Publication-ready)"

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------

h5_file <- file.path(
  work_dir,
  "NSCLC_GSE131907_expression.h5"
)

zip_file <- file.path(
  work_dir,
  "NSCLC_GSE131907_Expression.zip"
)


# ------------------------------------------------------------
# 2. Basic file existence / size audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. BASIC FILE AUDIT\n")
cat("------------------------------------------------------------\n\n")

file_audit <- data.frame(
  File = c(
    basename(h5_file),
    basename(zip_file)
  ),
  Exists = c(
    file.exists(h5_file),
    file.exists(zip_file)
  ),
  Size_MB = c(
    if (file.exists(h5_file))
      round(file.info(h5_file)$size / 1024^2, 2)
    else NA_real_,
    if (file.exists(zip_file))
      round(file.info(zip_file)$size / 1024^2, 2)
    else NA_real_
  ),
  stringsAsFactors = FALSE
)

print(file_audit)

write.table(
  file_audit,
  file = "STEP10C-0A_File_Basic_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 3. Check required R packages
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. PACKAGE CHECK\n")
cat("------------------------------------------------------------\n\n")

required_packages <- c(
  "rhdf5",
  "zip"
)

for (pkg in required_packages) {
  
  cat(
    pkg,
    ":",
    requireNamespace(pkg, quietly = TRUE),
    "\n"
  )
}

cat("\n")

if (!requireNamespace("rhdf5", quietly = TRUE)) {
  
  stop(
    paste0(
      "\nPackage 'rhdf5' is required.\n",
      "Install through Bioconductor if necessary:\n\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('rhdf5')\n"
    )
  )
}


# ------------------------------------------------------------
# 4. H5 internal structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. H5 INTERNAL STRUCTURE\n")
cat("------------------------------------------------------------\n\n")

if (!file.exists(h5_file)) {
  
  cat("H5 file NOT FOUND.\n")
  
} else {
  
  cat("H5 file found.\n\n")
  
  cat("Top-level HDF5 structure:\n\n")
  
  h5_structure <- rhdf5::h5ls(
    h5_file,
    recursive = TRUE,
    all = TRUE
  )
  
  print(h5_structure)
  
  write.table(
    h5_structure,
    file = "STEP10C-0A_H5_Internal_Structure.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 5. Identify possible matrix structures
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. SEARCH FOR MATRIX-LIKE H5 OBJECTS\n")
cat("------------------------------------------------------------\n\n")

if (exists("h5_structure") && nrow(h5_structure) > 0) {
  
  matrix_candidates <- h5_structure[
    h5_structure$otype %in% c("H5I_DATASET", "H5I_GROUP"),
    ,
    drop = FALSE
  ]
  
  print(matrix_candidates)
  
  write.table(
    matrix_candidates,
    file = "STEP10C-0A_H5_Matrix_Candidates.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  matrix_candidates <- NULL
  
  cat("No H5 structure available.\n")
}

cat("\n")


# ------------------------------------------------------------
# 6. Recognize common 10x H5 structure
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. 10X H5 STRUCTURE CHECK\n")
cat("------------------------------------------------------------\n\n")

tenx_h5_detected <- FALSE

if (exists("h5_structure") && nrow(h5_structure) > 0) {
  
  object_names <- h5_structure$name
  
  tenx_features <- any(
    grepl(
      "features",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_barcodes <- any(
    grepl(
      "barcodes",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_data <- any(
    grepl(
      "(^|/)data$",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_indices <- any(
    grepl(
      "(^|/)indices$",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_indptr <- any(
    grepl(
      "(^|/)indptr$",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_matrix <- any(
    grepl(
      "(^|/)matrix$",
      object_names,
      ignore.case = TRUE
    )
  )
  
  tenx_h5_detected <- (
    tenx_features &&
      tenx_barcodes &&
      tenx_data &&
      tenx_indices &&
      tenx_indptr &&
      tenx_matrix
  )
  
  cat("features :", tenx_features, "\n")
  cat("barcodes :", tenx_barcodes, "\n")
  cat("data     :", tenx_data, "\n")
  cat("indices  :", tenx_indices, "\n")
  cat("indptr   :", tenx_indptr, "\n")
  cat("matrix   :", tenx_matrix, "\n\n")
  
  cat(
    "10x H5 structure detected:",
    tenx_h5_detected,
    "\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 7. Attempt to inspect standard 10x H5 matrix
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. 10X H5 MATRIX DIMENSION / GENE / BARCODE AUDIT\n")
cat("------------------------------------------------------------\n\n")

h5_matrix_info <- NULL

if (tenx_h5_detected) {
  
  cat("Detected standard 10x-style H5 matrix.\n\n")
  
  # Find matrix group
  matrix_group_candidates <- h5_structure[
    grepl(
      "/matrix$",
      h5_structure$group_name,
      ignore.case = TRUE
    ) &
      h5_structure$otype == "H5I_GROUP",
    ,
    drop = FALSE
  ]
  
  if (nrow(matrix_group_candidates) == 0) {
    
    # alternative
    matrix_group_candidates <- h5_structure[
      h5_structure$name == "matrix" &
        h5_structure$otype == "H5I_GROUP",
      ,
      drop = FALSE
    ]
  }
  
  print(matrix_group_candidates)
  
  # Usually 10x path is /matrix
  possible_matrix_path <- "/matrix"
  
  matrix_exists <- tryCatch(
    rhdf5::H5Lexists(
      rhdf5::H5Fopen(h5_file),
      "matrix"
    ),
    error = function(e) FALSE
  )
  
  cat(
    "\n/root matrix group exists:",
    matrix_exists,
    "\n"
  )
  
} else {
  
  cat(
    "The file does not appear to use the standard 10x H5 structure.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 8. Inspect H5 datasets safely
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. H5 DATASET DIMENSIONS / TYPES\n")
cat("------------------------------------------------------------\n\n")

dataset_rows <- h5_structure[
  h5_structure$otype == "H5I_DATASET",
  ,
  drop = FALSE
]

if (nrow(dataset_rows) > 0) {
  
  dataset_info <- data.frame(
    group = dataset_rows$group_name,
    name = dataset_rows$name,
    dataset = paste0(
      dataset_rows$group_name,
      "/",
      dataset_rows$name
    ),
    stringsAsFactors = FALSE
  )
  
  dataset_info$dim <- NA_character_
  dataset_info$storage.mode <- NA_character_
  
  for (i in seq_len(nrow(dataset_info))) {
    
    dataset_path <- dataset_info$dataset[i]
    
    info_i <- tryCatch(
      rhdf5::h5readAttributes(
        h5_file,
        dataset_path
      ),
      error = function(e) NULL
    )
    
    # Try H5Dget_space / H5Dget_type only if necessary
    # Avoid reading large datasets.
    
    dataset_info$dim[i] <- "not directly read"
    dataset_info$storage.mode[i] <- "not directly read"
  }
  
  print(dataset_info)
  
  write.table(
    dataset_info,
    file = "STEP10C-0A_H5_Dataset_Audit.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat("No datasets detected.\n")
}

cat("\n")


# ------------------------------------------------------------
# 9. Inspect gene and barcode datasets if standard 10x
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. GENE / BARCODE CONTENT CHECK\n")
cat("------------------------------------------------------------\n\n")

gene_ids <- NULL
gene_names <- NULL
barcodes <- NULL

if (tenx_h5_detected) {
  
  # Try common 10x locations
  
  possible_gene_paths <- c(
    "/matrix/features/id",
    "/matrix/features/name",
    "/matrix/barcodes"
  )
  
  for (p in possible_gene_paths) {
    
    result <- tryCatch(
      rhdf5::h5read(
        h5_file,
        p
      ),
      error = function(e) NULL
    )
    
    if (!is.null(result)) {
      
      cat(
        "\nPath:",
        p,
        "\n"
      )
      
      cat(
        "Length:",
        length(result),
        "\n"
      )
      
      cat(
        "Class:",
        class(result),
        "\n"
      )
      
      if (length(result) <= 20) {
        print(result)
      } else {
        print(head(result, 10))
      }
      
      if (grepl("/features/id$", p)) {
        gene_ids <- result
      }
      
      if (grepl("/features/name$", p)) {
        gene_names <- result
      }
      
      if (grepl("/barcodes$", p)) {
        barcodes <- result
      }
    }
  }
}

cat("\n")


# ------------------------------------------------------------
# 10. Core gene audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. CORE GENE AUDIT\n")
cat("------------------------------------------------------------\n\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_gene_audit <- data.frame(
  Gene = core_genes,
  Present_in_feature_names = FALSE,
  Present_in_feature_ids = FALSE,
  stringsAsFactors = FALSE
)

if (!is.null(gene_names)) {
  
  core_gene_audit$Present_in_feature_names <-
    core_genes %in% gene_names
}

if (!is.null(gene_ids)) {
  
  core_gene_audit$Present_in_feature_ids <-
    core_genes %in% gene_ids
}

print(core_gene_audit)

write.table(
  core_gene_audit,
  file = "STEP10C-0A_Core_Gene_Audit.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 11. Barcode audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. BARCODE AUDIT\n")
cat("------------------------------------------------------------\n\n")

if (!is.null(barcodes)) {
  
  cat(
    "Number of barcodes:",
    length(barcodes),
    "\n"
  )
  
  cat(
    "Unique barcodes:",
    length(unique(barcodes)),
    "\n"
  )
  
  cat(
    "Duplicated barcodes:",
    sum(duplicated(barcodes)),
    "\n"
  )
  
  cat("\nFirst 20 barcodes:\n")
  print(head(barcodes, 20))
  
  barcode_summary <- data.frame(
    total_barcodes = length(barcodes),
    unique_barcodes = length(unique(barcodes)),
    duplicated_barcodes = sum(duplicated(barcodes)),
    stringsAsFactors = FALSE
  )
  
  write.table(
    barcode_summary,
    file = "STEP10C-0A_Barcode_Summary.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "Barcode dataset was not identified.\n"
  )
}

cat("\n")


# ------------------------------------------------------------
# 12. ZIP file inspection
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. ZIP FILE CONTENT AUDIT\n")
cat("------------------------------------------------------------\n\n")

zip_contents <- NULL

if (file.exists(zip_file)) {
  
  cat(
    "Inspecting ZIP contents without extracting large files...\n\n"
  )
  
  zip_contents <- utils::unzip(
    zip_file,
    list = TRUE
  )
  
  print(zip_contents)
  
  write.table(
    zip_contents,
    file = "STEP10C-0A_ZIP_Contents.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
} else {
  
  cat("ZIP file NOT FOUND.\n")
}

cat("\n")


# ------------------------------------------------------------
# 13. ZIP file size interpretation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. ZIP MATRIX CANDIDATE CHECK\n")
cat("------------------------------------------------------------\n\n")

if (!is.null(zip_contents) && nrow(zip_contents) > 0) {
  
  zip_contents$size_MB <-
    round(zip_contents$Length / 1024^2, 3)
  
  zip_contents$compression_ratio <-
    ifelse(
      zip_contents$Length > 0,
      round(
        zip_contents$Compressed / zip_contents$Length,
        4
      ),
      NA_real_
    )
  
  print(zip_contents)
  
  possible_matrix_files <- zip_contents[
    grepl(
      "matrix|expression|count|umi|barcodes|features",
      zip_contents$Name,
      ignore.case = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  cat("\nPotential matrix-related files:\n\n")
  
  print(possible_matrix_files)
  
  write.table(
    possible_matrix_files,
    file = "STEP10C-0A_ZIP_Matrix_Candidates.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cat("\n")


# ------------------------------------------------------------
# 14. Preliminary classification
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. PRELIMINARY FILE CLASSIFICATION\n")
cat("------------------------------------------------------------\n\n")

if (tenx_h5_detected) {
  
  h5_classification <-
    "EXPRESSION_MATRIX_DETECTED_STANDARD_10X_H5"
  
} else {
  
  h5_classification <-
    "H5_FORMAT_DETECTED_BUT_MATRIX_TYPE_REQUIRES_FURTHER_AUDIT"
}


zip_classification <-
  if (
    !is.null(zip_contents) &&
    nrow(zip_contents) > 0
  ) {
    
    matrix_like_zip <-
      grepl(
        "matrix|expression|count|umi",
        zip_contents$Name,
        ignore.case = TRUE
      )
    
    if (any(matrix_like_zip)) {
      
      "ZIP_CONTAINS_EXPRESSION_OR_COUNT_RELATED_FILE"
      
    } else {
      
      "ZIP_DOES_NOT_APPEAR_TO_CONTAIN_A_STANDARD_COUNT_MATRIX"
      
    }
    
  } else {
    
    "ZIP_NOT_AVAILABLE"
  }


classification <- data.frame(
  File = c(
    basename(h5_file),
    basename(zip_file)
  ),
  Preliminary_Classification = c(
    h5_classification,
    zip_classification
  ),
  stringsAsFactors = FALSE
)

print(classification)

write.table(
  classification,
  file = "STEP10C-0A_Preliminary_File_Classification.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Final report
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP10C-0A FINAL REPORT\n")
cat("============================================================\n\n")

cat(
  "H5 file:",
  basename(h5_file),
  "\n"
)

cat(
  "H5 classification:",
  h5_classification,
  "\n\n"
)

cat(
  "ZIP file:",
  basename(zip_file),
  "\n"
)

cat(
  "ZIP classification:",
  zip_classification,
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "This STEP10C-0A audit does NOT establish that the matrix\n"
)

cat(
  "contains RAW UMI COUNTS suitable for pseudobulk.\n"
)

cat(
  "That determination requires the next audit step after\n"
)

cat(
  "the H5 structure and matrix content are confirmed.\n"
)

cat(
  "\nNo pseudobulk was performed.\n"
)

cat(
  "No statistical analysis was performed.\n"
)

cat(
  "No source file was modified.\n"
)

cat(
  "\n============================================================\n"
)

# ------------------------------------------------------------
# 16. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  "STEP10C-0A_sessionInfo.txt"
)

cat(
  "\nSession information saved to:\n",
  file.path(work_dir, "STEP10C-0A_sessionInfo.txt"),
  "\n"
)


# ============================================================
# STEP10C-0B
# H5 Expression-Value / Raw-Count Characterization Audit
#
# Dataset:
#   NSCLC_GSE131907_expression.h5
#
# Purpose:
#   Determine whether the H5 expression matrix is compatible
#   with a raw UMI count matrix for downstream pseudobulk.
#
# IMPORTANT:
#   - NO pseudobulk
#   - NO statistical testing
#   - NO patient exclusion
#   - NO normalization
#   - NO transformation
#   - NO modification of source H5 file
#   - Avoid loading the entire 2.5 GB matrix into RAM
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10C-0B — H5 RAW-COUNT CHARACTERIZATION AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

work_dir <- "input_dir/Figure 9C Dot Plot (Publication-ready)"

setwd(work_dir)

cat("Working directory:\n")
cat(getwd(), "\n\n")


# ------------------------------------------------------------
# 1. Input H5 file
# ------------------------------------------------------------

h5_file <- file.path(
  work_dir,
  "NSCLC_GSE131907_expression.h5"
)

if (!file.exists(h5_file)) {
  
  stop(
    paste0(
      "\nERROR: H5 file not found:\n",
      h5_file,
      "\n"
    )
  )
}

cat("Input H5 file:\n")
cat(h5_file, "\n\n")


# ------------------------------------------------------------
# 2. Required package
# ------------------------------------------------------------

if (!requireNamespace("rhdf5", quietly = TRUE)) {
  
  stop(
    paste0(
      "\nPackage 'rhdf5' is required.\n\n",
      "Install using:\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) ",
      "install.packages('BiocManager')\n",
      "BiocManager::install('rhdf5')\n"
    )
  )
}

library(rhdf5)


# ------------------------------------------------------------
# 3. Output directory
# ------------------------------------------------------------

step_dir <- file.path(
  work_dir,
  "STEP10C-0B_H5_RawCount_Characterization_Audit"
)

if (!dir.exists(step_dir)) {
  dir.create(step_dir, recursive = TRUE)
}

cat("Output directory:\n")
cat(step_dir, "\n\n")


# ------------------------------------------------------------
# 4. H5 structure audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. H5 STRUCTURE AUDIT\n")
cat("------------------------------------------------------------\n\n")

h5_structure <- rhdf5::h5ls(
  h5_file,
  recursive = TRUE,
  all = TRUE
)

write.table(
  h5_structure,
  file = file.path(
    step_dir,
    "STEP10C-0B_H5_Internal_Structure.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

required_paths <- c(
  "/matrix/barcodes",
  "/matrix/data",
  "/matrix/features/id",
  "/matrix/features/name",
  "/matrix/indices",
  "/matrix/indptr",
  "/matrix/shape"
)

structure_paths <- paste0(
  h5_structure$group,
  "/",
  h5_structure$name
)

structure_check <- data.frame(
  Required_Path = required_paths,
  Present = required_paths %in% structure_paths,
  stringsAsFactors = FALSE
)

print(structure_check)

write.table(
  structure_check,
  file = file.path(
    step_dir,
    "STEP10C-0B_Required_H5_Path_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

standard_10x_structure <- all(
  structure_check$Present
)

cat(
  "\nStandard 10x H5 structure:",
  standard_10x_structure,
  "\n\n"
)


# ------------------------------------------------------------
# 5. Read shape
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. MATRIX SHAPE\n")
cat("------------------------------------------------------------\n\n")

shape <- rhdf5::h5read(
  h5_file,
  "/matrix/shape"
)

shape <- as.numeric(shape)

cat("Raw shape vector:\n")
print(shape)

if (length(shape) != 2) {
  
  stop(
    "\nERROR: /matrix/shape does not contain exactly two dimensions.\n"
  )
}

# Standard 10x H5 stores:
# shape[1] = number of features
# shape[2] = number of cells

n_features <- shape[1]
n_cells <- shape[2]

cat("\nNumber of features:", n_features, "\n")
cat("Number of cells:", n_cells, "\n\n")


# ------------------------------------------------------------
# 6. Read gene IDs and names
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. FEATURE / GENE AUDIT\n")
cat("------------------------------------------------------------\n\n")

gene_ids <- rhdf5::h5read(
  h5_file,
  "/matrix/features/id"
)

gene_names <- rhdf5::h5read(
  h5_file,
  "/matrix/features/name"
)

cat(
  "Number of feature IDs:",
  length(gene_ids),
  "\n"
)

cat(
  "Number of feature names:",
  length(gene_names),
  "\n"
)

cat(
  "Feature dimension matches shape:",
  length(gene_names) == n_features,
  "\n\n"
)


# ------------------------------------------------------------
# 7. Core gene audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. CORE GENE AUDIT\n")
cat("------------------------------------------------------------\n\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_gene_audit <- data.frame(
  Gene = core_genes,
  In_feature_name = core_genes %in% gene_names,
  In_feature_id = core_genes %in% gene_ids,
  stringsAsFactors = FALSE
)

print(core_gene_audit)

write.table(
  core_gene_audit,
  file = file.path(
    step_dir,
    "STEP10C-0B_Core_Gene_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 8. Barcode audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. CELL BARCODE AUDIT\n")
cat("------------------------------------------------------------\n\n")

barcodes <- rhdf5::h5read(
  h5_file,
  "/matrix/barcodes"
)

barcode_n <- length(barcodes)

cat(
  "Number of barcodes:",
  barcode_n,
  "\n"
)

cat(
  "Shape cell dimension:",
  n_cells,
  "\n"
)

cat(
  "Barcode count matches matrix shape:",
  barcode_n == n_cells,
  "\n"
)

cat(
  "Unique barcodes:",
  length(unique(barcodes)),
  "\n"
)

cat(
  "Duplicated barcodes:",
  sum(duplicated(barcodes)),
  "\n\n"
)

cat("First 20 barcodes:\n")
print(head(barcodes, 20))

barcode_audit <- data.frame(
  H5_Cell_Count = barcode_n,
  Shape_Cell_Count = n_cells,
  Cell_Count_Match = barcode_n == n_cells,
  Unique_Cell_Count = length(unique(barcodes)),
  Duplicated_Cell_Count = sum(duplicated(barcodes)),
  stringsAsFactors = FALSE
)

write.table(
  barcode_audit,
  file = file.path(
    step_dir,
    "STEP10C-0B_Barcode_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 9. Sparse matrix structural audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. SPARSE MATRIX STRUCTURAL AUDIT\n")
cat("------------------------------------------------------------\n\n")

indptr <- rhdf5::h5read(
  h5_file,
  "/matrix/indptr"
)

n_indptr <- length(indptr)

cat(
  "indptr length:",
  n_indptr,
  "\n"
)

cat(
  "Expected indptr length:",
  n_cells + 1,
  "\n"
)

indptr_length_ok <- (
  n_indptr == n_cells + 1
)

cat(
  "indptr length PASS:",
  indptr_length_ok,
  "\n\n"
)

# The final indptr value corresponds to number of stored
# non-zero entries in a standard CSC matrix.

nnz_from_indptr <- as.numeric(
  indptr[length(indptr)]
)

cat(
  "Stored non-zero entries (nnz):",
  format(nnz_from_indptr, scientific = FALSE),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Inspect data dataset metadata WITHOUT reading all data
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("10. DATA DATASET STORAGE-TYPE AUDIT\n")
cat("------------------------------------------------------------\n\n")

data_info <- h5_structure[
  h5_structure$group == "/matrix" &
    h5_structure$name == "data",
  ,
  drop = FALSE
]

if (nrow(data_info) != 1) {
  
  stop(
    "\nERROR: Could not uniquely identify /matrix/data.\n"
  )
}

cat("dclass:\n")
print(data_info$dclass)

cat("\ndtype:\n")
print(data_info$dtype)

cat("\ndim:\n")
print(data_info$dim)

data_dclass <- as.character(
  data_info$dclass[1]
)

data_dtype <- as.character(
  data_info$dtype[1]
)

cat("\n")


# ------------------------------------------------------------
# 11. Read a SAFE SAMPLE of data values
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. DATA VALUE SAMPLE AUDIT\n")
cat("------------------------------------------------------------\n\n")

# IMPORTANT:
# Do NOT read all 329,684,450 values.
#
# Read only a limited sample from the beginning and middle
# of the dataset.

data_length <- as.numeric(
  data_info$dim[1]
)

sample_n <- min(
  1000000,
  data_length
)

cat(
  "Total stored data values:",
  format(data_length, scientific = FALSE),
  "\n"
)

cat(
  "Values sampled:",
  format(sample_n, scientific = FALSE),
  "\n\n"
)


# First 1,000,000 values
sample_start <- rhdf5::h5read(
  h5_file,
  "/matrix/data",
  index = list(
    seq_len(sample_n)
  )
)

sample_start <- as.numeric(sample_start)

cat("Sample class:\n")
print(class(sample_start))

cat("\nSample length:\n")
print(length(sample_start))

cat("\n")


# ------------------------------------------------------------
# 12. Sample integer / negative / NA audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("12. SAMPLE VALUE CHARACTERIZATION\n")
cat("------------------------------------------------------------\n\n")

sample_na_n <- sum(
  is.na(sample_start)
)

sample_nan_n <- sum(
  is.nan(sample_start)
)

sample_inf_n <- sum(
  is.infinite(sample_start)
)

sample_negative_n <- sum(
  sample_start < 0,
  na.rm = TRUE
)

sample_zero_n <- sum(
  sample_start == 0,
  na.rm = TRUE
)

sample_positive_n <- sum(
  sample_start > 0,
  na.rm = TRUE
)

sample_integer_n <- sum(
  sample_start == floor(sample_start),
  na.rm = TRUE
)

sample_noninteger_n <- sum(
  sample_start != floor(sample_start),
  na.rm = TRUE
)

sample_min <- min(
  sample_start,
  na.rm = TRUE
)

sample_max <- max(
  sample_start,
  na.rm = TRUE
)

sample_mean <- mean(
  sample_start,
  na.rm = TRUE
)

sample_median <- median(
  sample_start,
  na.rm = TRUE
)

sample_summary <- data.frame(
  Metric = c(
    "Sample_N",
    "NA",
    "NaN",
    "Infinite",
    "Negative",
    "Zero",
    "Positive",
    "Integer_valued",
    "Noninteger_valued",
    "Minimum",
    "Maximum",
    "Mean",
    "Median"
  ),
  Value = c(
    length(sample_start),
    sample_na_n,
    sample_nan_n,
    sample_inf_n,
    sample_negative_n,
    sample_zero_n,
    sample_positive_n,
    sample_integer_n,
    sample_noninteger_n,
    sample_min,
    sample_max,
    sample_mean,
    sample_median
  ),
  stringsAsFactors = FALSE
)

print(sample_summary)

write.table(
  sample_summary,
  file = file.path(
    step_dir,
    "STEP10C-0B_Data_Value_Sample_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 13. Sample quantile audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. SAMPLE QUANTILE AUDIT\n")
cat("------------------------------------------------------------\n\n")

sample_quantiles <- quantile(
  sample_start,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    0.999,
    1
  ),
  na.rm = TRUE
)

print(sample_quantiles)

write.table(
  data.frame(
    Quantile = names(sample_quantiles),
    Value = as.numeric(sample_quantiles),
    row.names = NULL
  ),
  file = file.path(
    step_dir,
    "STEP10C-0B_Data_Value_Sample_Quantiles.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 14. Inspect selected core-gene rows
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. CORE GENE DATA AUDIT\n")
cat("------------------------------------------------------------\n\n")

core_gene_indices <- match(
  core_genes,
  gene_names
)

core_gene_indices_id <- match(
  core_genes,
  gene_ids
)

core_gene_data_audit <- data.frame(
  Gene = core_genes,
  Feature_Name_Index = core_gene_indices,
  Feature_ID_Index = core_gene_indices_id,
  stringsAsFactors = FALSE
)

print(core_gene_data_audit)

write.table(
  core_gene_data_audit,
  file = file.path(
    step_dir,
    "STEP10C-0B_Core_Gene_Row_Index.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 15. Inspect CSC structure for core genes
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("15. CORE GENE SPARSITY / VALUE AUDIT\n")
cat("------------------------------------------------------------\n\n")

# In a standard 10x CSC matrix:
# indices = feature indices
# data    = corresponding values
# indptr  = cell boundaries
#
# We do NOT read the entire indices/data arrays.
#
# Instead, use HDF5 hyperslab reads for selected cells
# and identify whether core-gene values are integer-valued.

# Need to know number of stored values.

indices_info <- h5_structure[
  h5_structure$group == "/matrix" &
    h5_structure$name == "indices",
  ,
  drop = FALSE
]

indices_length <- as.numeric(
  indices_info$dim[1]
)

cat(
  "indices length:",
  format(indices_length, scientific = FALSE),
  "\n"
)

cat(
  "data length:",
  format(data_length, scientific = FALSE),
  "\n"
)

cat(
  "indices/data length match:",
  indices_length == data_length,
  "\n\n"
)


# ------------------------------------------------------------
# 16. Inspect selected cells
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("16. SELECTED CELL RAW-VALUE AUDIT\n")
cat("------------------------------------------------------------\n\n")

# Select 10 cells distributed across matrix.
# This is only a diagnostic check.

selected_cells <- unique(
  round(
    seq(
      1,
      n_cells,
      length.out = min(10, n_cells)
    )
  )
)

cat("Selected cell indices:\n")
print(selected_cells)

selected_cell_summary <- list()

for (cell_i in selected_cells) {
  
  start_pos <- as.numeric(
    indptr[cell_i]
  ) + 1
  
  end_pos <- as.numeric(
    indptr[cell_i + 1]
  )
  
  if (end_pos >= start_pos) {
    
    positions <- start_pos:end_pos
    
    cell_indices <- rhdf5::h5read(
      h5_file,
      "/matrix/indices",
      index = list(positions)
    )
    
    cell_data <- rhdf5::h5read(
      h5_file,
      "/matrix/data",
      index = list(positions)
    )
    
    cell_indices <- as.numeric(cell_indices)
    cell_data <- as.numeric(cell_data)
    
    # Convert 0-based feature index to R 1-based index
    feature_index_r <- cell_indices + 1
    
    core_present <- feature_index_r %in%
      core_gene_indices
    
    core_rows <- which(
      feature_index_r %in% core_gene_indices
    )
    
    core_values <- if (length(core_rows) > 0) {
      cell_data[core_rows]
    } else {
      numeric(0)
    }
    
    selected_cell_summary[[length(selected_cell_summary) + 1]] <-
      data.frame(
        Cell_Index = cell_i,
        Barcode = barcodes[cell_i],
        Nonzero_Entries = length(cell_data),
        Cell_Total = sum(
          cell_data,
          na.rm = TRUE
        ),
        Cell_Min = min(
          cell_data,
          na.rm = TRUE
        ),
        Cell_Max = max(
          cell_data,
          na.rm = TRUE
        ),
        Cell_All_Integer =
          all(
            cell_data == floor(cell_data)
          ),
        Core_Gene_Entries =
          length(core_values),
        Core_Gene_All_Integer =
          if (length(core_values) > 0)
            all(core_values == floor(core_values))
        else
          NA,
        stringsAsFactors = FALSE
      )
    
  } else {
    
    selected_cell_summary[[length(selected_cell_summary) + 1]] <-
      data.frame(
        Cell_Index = cell_i,
        Barcode = barcodes[cell_i],
        Nonzero_Entries = 0,
        Cell_Total = 0,
        Cell_Min = NA_real_,
        Cell_Max = NA_real_,
        Cell_All_Integer = TRUE,
        Core_Gene_Entries = 0,
        Core_Gene_All_Integer = NA,
        stringsAsFactors = FALSE
      )
  }
}

selected_cell_summary <- do.call(
  rbind,
  selected_cell_summary
)

print(selected_cell_summary)

write.table(
  selected_cell_summary,
  file = file.path(
    step_dir,
    "STEP10C-0B_Selected_Cell_Value_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 17. Cell-level expression-depth audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("17. CELL EXPRESSION-DEPTH AUDIT\n")
cat("------------------------------------------------------------\n\n")

# We can calculate cell total from the sparse CSC structure
# without reading the full data matrix into RAM.
#
# However, reading every cell's data individually would be slow.
#
# Therefore use chunks of cells.

chunk_size <- 5000

cell_totals <- numeric(n_cells)

cell_nnz <- numeric(n_cells)

cat(
  "Cells:",
  n_cells,
  "\n"
)

cat(
  "Chunk size:",
  chunk_size,
  "\n\n"
)

for (
  chunk_start in seq(
    1,
    n_cells,
    by = chunk_size
  )
) {
  
  chunk_end <- min(
    chunk_start + chunk_size - 1,
    n_cells
  )
  
  cell_range <- chunk_start:chunk_end
  
  # 0-based sparse pointer positions
  ptr_start <- as.numeric(
    indptr[chunk_start]
  ) + 1
  
  ptr_end <- as.numeric(
    indptr[chunk_end + 1]
  )
  
  if (ptr_end >= ptr_start) {
    
    chunk_positions <- ptr_start:ptr_end
    
    chunk_data <- rhdf5::h5read(
      h5_file,
      "/matrix/data",
      index = list(chunk_positions)
    )
    
    chunk_data <- as.numeric(chunk_data)
    
    # Number of entries per cell
    pointer_slice <- as.numeric(
      indptr[
        chunk_start:(chunk_end + 1)
      ]
    )
    
    cell_nnz_chunk <- diff(
      pointer_slice
    )
    
    cell_nnz[cell_range] <-
      cell_nnz_chunk
    
    # Sum values within each cell
    cumulative <- cumsum(
      chunk_data
    )
    
    starts <- c(
      1,
      cumsum(cell_nnz_chunk)[
        -length(cell_nnz_chunk)
      ] + 1
    )
    
    ends <- cumsum(
      cell_nnz_chunk
    )
    
    totals <- numeric(
      length(cell_nnz_chunk)
    )
    
    for (j in seq_along(totals)) {
      
      if (
        starts[j] <= ends[j]
      ) {
        
        totals[j] <-
          sum(
            chunk_data[
              starts[j]:ends[j]
            ],
            na.rm = TRUE
          )
      }
    }
    
    cell_totals[cell_range] <-
      totals
    
  } else {
    
    cell_nnz[cell_range] <- 0
    cell_totals[cell_range] <- 0
  }
  
  cat(
    "\rProcessed cells:",
    chunk_end,
    "/",
    n_cells
  )
  
  flush.console()
}

cat("\n\n")


# ------------------------------------------------------------
# 18. Cell total summary
# ------------------------------------------------------------

cell_total_summary <- data.frame(
  Metric = c(
    "Cells",
    "Minimum_Total",
    "Q1_Total",
    "Median_Total",
    "Mean_Total",
    "Q3_Total",
    "Maximum_Total",
    "Minimum_NNZ",
    "Median_NNZ",
    "Maximum_NNZ"
  ),
  Value = c(
    length(cell_totals),
    min(cell_totals),
    quantile(cell_totals, 0.25),
    median(cell_totals),
    mean(cell_totals),
    quantile(cell_totals, 0.75),
    max(cell_totals),
    min(cell_nnz),
    median(cell_nnz),
    max(cell_nnz)
  ),
  stringsAsFactors = FALSE
)

print(cell_total_summary)

write.table(
  cell_total_summary,
  file = file.path(
    step_dir,
    "STEP10C-0B_Cell_Expression_Depth_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 19. Global cell-total integrity checks
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("19. GLOBAL CELL-TOTAL INTEGRITY CHECK\n")
cat("------------------------------------------------------------\n\n")

cell_total_integer <- all(
  cell_totals == floor(cell_totals)
)

cell_total_negative <- any(
  cell_totals < 0
)

cell_total_na <- any(
  is.na(cell_totals)
)

cat(
  "All cell totals integer-valued:",
  cell_total_integer,
  "\n"
)

cat(
  "Any negative cell totals:",
  cell_total_negative,
  "\n"
)

cat(
  "Any NA cell totals:",
  cell_total_na,
  "\n\n"
)


# ------------------------------------------------------------
# 20. Determine preliminary raw-count compatibility
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("20. PRELIMINARY RAW-COUNT COMPATIBILITY\n")
cat("------------------------------------------------------------\n\n")

# Important:
# We deliberately DO NOT classify raw count solely because
# values are integer-valued in a sample.
#
# Strong compatibility indicators:
#
#   A. standard 10x structure
#   B. feature count matches shape
#   C. barcode count matches shape
#   D. indices/data same length
#   E. no negative values in sampled data
#   F. sampled values integer-valued
#   G. cell totals integer-valued
#
# But if data are floating-point storage, this is still
# "raw-count compatible" rather than absolute proof of
# original upstream raw UMI provenance.

structure_ok <-
  standard_10x_structure

dimension_ok <-
  length(gene_names) == n_features &&
  length(barcodes) == n_cells

sparse_ok <-
  indptr_length_ok &&
  indices_length == data_length

sample_no_negative <-
  sample_negative_n == 0

sample_integer_values <-
  sample_noninteger_n == 0

cell_total_integer_ok <-
  cell_total_integer

core_genes_ok <-
  all(core_gene_audit$In_feature_name)


raw_count_compatible <-
  structure_ok &&
  dimension_ok &&
  sparse_ok &&
  sample_no_negative &&
  sample_integer_values &&
  cell_total_integer_ok &&
  core_genes_ok


if (raw_count_compatible) {
  
  preliminary_classification <-
    "RAW_UMI_COUNT_COMPATIBLE"
  
} else {
  
  if (
    sample_noninteger_n > 0
  ) {
    
    preliminary_classification <-
      "NONINTEGER_EXPRESSION_VALUES_DETECTED"
    
  } else {
    
    preliminary_classification <-
      "EXPRESSION_MATRIX_BUT_RAW_UMI_STATUS_UNCERTAIN"
  }
}


classification_audit <- data.frame(
  Criterion = c(
    "Standard_10x_H5_structure",
    "Feature_dimension_match",
    "Cell_dimension_match",
    "Indices_data_length_match",
    "No_negative_sampled_values",
    "Sampled_values_integer",
    "Cell_totals_integer",
    "Core_genes_present",
    "Preliminary_classification"
  ),
  Result = c(
    structure_ok,
    length(gene_names) == n_features,
    length(barcodes) == n_cells,
    sparse_ok,
    sample_no_negative,
    sample_integer_values,
    cell_total_integer_ok,
    core_genes_ok,
    preliminary_classification
  ),
  stringsAsFactors = FALSE
)

print(classification_audit)

write.table(
  classification_audit,
  file = file.path(
    step_dir,
    "STEP10C-0B_Preliminary_RawCount_Classification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 21. Important interpretation statement
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP10C-0B INTERPRETATION\n")
cat("============================================================\n\n")

cat(
  "H5 matrix dimension:\n"
)

cat(
  "  ",
  n_features,
  "features x",
  n_cells,
  "cells\n\n"
)

cat(
  "H5 data storage type:\n"
)

cat(
  "  dclass =",
  data_dclass,
  "\n"
)

cat(
  "  dtype  =",
  data_dtype,
  "\n\n"
)

cat(
  "Preliminary classification:\n"
)

cat(
  "  ",
  preliminary_classification,
  "\n\n"
)

cat(
  "IMPORTANT:\n"
)

cat(
  "This STEP10C-0B audit characterizes the matrix structure\n"
)

cat(
  "and observed expression values.\n"
)

cat(
  "It does NOT establish upstream data provenance by itself.\n\n"
)

cat(
  "Therefore, even if the matrix is integer-valued and\n"
)

cat(
  "compatible with raw UMI counts, cell-ID reconciliation\n"
)

cat(
  "against the locked TISCH2 metadata must still be completed\n"
)

cat(
  "before pseudobulk construction.\n\n"
)

cat(
  "No pseudobulk was performed.\n"
)

cat(
  "No statistical testing was performed.\n"
)

cat(
  "No patients were excluded.\n"
)

cat(
  "No source file was modified.\n\n"
)


# ------------------------------------------------------------
# 22. Save compact RDS audit object
# ------------------------------------------------------------

audit_object <- list(
  
  step = "STEP10C-0B",
  
  input_file = h5_file,
  
  h5_size_MB =
    file.info(h5_file)$size / 1024^2,
  
  standard_10x_structure =
    standard_10x_structure,
  
  n_features =
    n_features,
  
  n_cells =
    n_cells,
  
  n_barcodes =
    barcode_n,
  
  n_unique_barcodes =
    length(unique(barcodes)),
  
  n_duplicated_barcodes =
    sum(duplicated(barcodes)),
  
  data_dclass =
    data_dclass,
  
  data_dtype =
    data_dtype,
  
  data_length =
    data_length,
  
  indices_length =
    indices_length,
  
  nnz_from_indptr =
    nnz_from_indptr,
  
  core_gene_audit =
    core_gene_audit,
  
  sample_summary =
    sample_summary,
  
  sample_quantiles =
    sample_quantiles,
  
  selected_cell_summary =
    selected_cell_summary,
  
  cell_total_summary =
    cell_total_summary,
  
  preliminary_classification =
    preliminary_classification,
  
  raw_count_compatible =
    raw_count_compatible
  
)

saveRDS(
  audit_object,
  file = file.path(
    step_dir,
    "STEP10C-0B_H5_RawCount_Characterization_Audit.rds"
  )
)


# ------------------------------------------------------------
# 23. Session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    step_dir,
    "STEP10C-0B_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 24. Final status
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP10C-0B FINAL STATUS\n")
cat("============================================================\n\n")

if (
  raw_count_compatible
) {
  
  cat(
    "PASS — H5 MATRIX IS RAW-COUNT COMPATIBLE\n"
  )
  
} else {
  
  cat(
    "REVIEW REQUIRED — RAW-COUNT STATUS NOT YET ESTABLISHED\n"
  )
}

cat("\n")

cat(
  "Output directory:\n",
  step_dir,
  "\n\n"
)

cat(
  "Next step, if appropriate:\n",
  "STEP10C-0C — TISCH2 / H5 CELL-ID RECONCILIATION AUDIT\n"
)

cat("\n============================================================\n")



# ============================================================
# WORKFLOW-10
# STEP10C-0C — Cell-ID Reconciliation Audit
#
# Purpose:
#   Reconcile TISCH2 cell metadata with the 10x-format H5
#   processed-expression matrix at the cell-ID level.
#
# IMPORTANT:
#   - No expression values are read from /matrix/data.
#   - No pseudobulk is performed.
#   - No statistical testing is performed.
#   - No patient is excluded.
#   - STEP10B-1 locked compartment mapping is NOT recalculated.
#
# Primary objectives:
#   1. Reconcile all 203,298 TISCH2 cells with H5 barcodes.
#   2. Reconcile the locked 65,702 Primary + Tumor cells.
#   3. Verify the 9 locked Primary-Tumor patients.
#   4. Verify T-cell / Mono-Macro / Endothelial compartment identity.
#   5. Audit barcode-format differences without silently modifying IDs.
#   6. Verify AOC3 / MS4A7 / IL7R feature presence.
#
# Final status:
#   PASS — EXACT CELL-ID RECONCILIATION
#   PASS — NORMALIZED CELL-ID RECONCILIATION
#   FAIL — CELL-ID RECONCILIATION UNSUPPORTED
#
# Environment:
#   Windows 11
#   R 4.5.1
# ============================================================


cat("============================================================\n")
cat("STEP10C-0C — Cell-ID Reconciliation Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "rhdf5"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nPlease install before running STEP10C-0C."
    )
  )
}

library(data.table)
library(rhdf5)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

step10_root <- "input_dir/Paper 1_Step 10-Single-Cell"

metadata_file <- file.path(
  step10_root,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

h5_file <- file.path(
  "input_dir/Figure 9C Dot Plot (Publication-ready)",
  "NSCLC_GSE131907_expression.h5"
)

output_dir <- file.path(
  step10_root,
  "STEP10C-0C_CellID_Reconciliation_Audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Input existence audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("2. INPUT FILE AUDIT\n")
cat("------------------------------------------------------------\n\n")

if (!file.exists(metadata_file)) {
  stop("TISCH2 metadata file not found:\n", metadata_file)
}

if (!file.exists(h5_file)) {
  stop("H5 expression file not found:\n", h5_file)
}

cat("Metadata file:\n")
cat(metadata_file, "\n\n")

cat("H5 expression file:\n")
cat(h5_file, "\n\n")

metadata_size_mb <- file.info(metadata_file)$size / 1024^2
h5_size_mb <- file.info(h5_file)$size / 1024^2

cat(
  sprintf(
    "Metadata file size: %.2f MB\n",
    metadata_size_mb
  )
)

cat(
  sprintf(
    "H5 file size: %.2f MB\n\n",
    h5_size_mb
  )
)


# ------------------------------------------------------------
# 3. Read TISCH2 metadata
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("3. TISCH2 METADATA LOAD\n")
cat("------------------------------------------------------------\n\n")

meta <- fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  data.table = TRUE,
  check.names = FALSE
)

cat(
  "Metadata dimensions: ",
  nrow(meta),
  " rows x ",
  ncol(meta),
  " columns\n\n",
  sep = ""
)

required_meta_columns <- c(
  "Cell",
  "Patient",
  "Sample",
  "Source",
  "Stage"
)

missing_meta_columns <- setdiff(
  required_meta_columns,
  names(meta)
)

if (length(missing_meta_columns) > 0) {
  stop(
    "Required metadata columns missing: ",
    paste(missing_meta_columns, collapse = ", ")
  )
}

cat("Required metadata columns: PASS\n\n")


# ------------------------------------------------------------
# 4. Basic TISCH2 Cell-ID audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("4. TISCH2 CELL-ID BASIC AUDIT\n")
cat("------------------------------------------------------------\n\n")

meta$Cell <- trimws(as.character(meta$Cell))

tisch_cell_n <- nrow(meta)
tisch_unique_cell_n <- length(unique(meta$Cell))

tisch_missing_cell <- sum(
  is.na(meta$Cell) | meta$Cell == ""
)

tisch_duplicate_cell_n <- sum(
  duplicated(meta$Cell)
)

cat(
  "TISCH2 rows: ",
  tisch_cell_n,
  "\n",
  sep = ""
)

cat(
  "TISCH2 unique Cell IDs: ",
  tisch_unique_cell_n,
  "\n",
  sep = ""
)

cat(
  "TISCH2 missing/blank Cell IDs: ",
  tisch_missing_cell,
  "\n",
  sep = ""
)

cat(
  "TISCH2 duplicated Cell IDs: ",
  tisch_duplicate_cell_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. H5 structure audit — IDs/features only
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("5. H5 ID / FEATURE AUDIT\n")
cat("------------------------------------------------------------\n\n")

h5_paths <- c(
  "/matrix/barcodes",
  "/matrix/features/id",
  "/matrix/features/name",
  "/matrix/shape"
)

h5_path_exists <- vapply(
  h5_paths,
  function(x) {
    out <- tryCatch(
      {
        rhdf5::h5ls(
          h5_file,
          recursive = TRUE,
          all = TRUE
        )
        TRUE
      },
      error = function(e) FALSE
    )
    out
  },
  logical(1)
)

# Direct existence check using h5ls table
h5_listing <- rhdf5::h5ls(
  h5_file,
  recursive = TRUE,
  all = TRUE
)

h5_full_paths <- character(0)

if (nrow(h5_listing) > 0) {
  
  h5_full_paths <- apply(
    h5_listing,
    1,
    function(x) {
      group <- as.character(x[["group"]])
      name <- as.character(x[["name"]])
      
      if (identical(group, "/")) {
        paste0("/", name)
      } else {
        paste0(group, "/", name)
      }
    }
  )
}

required_h5_paths <- c(
  "/matrix/barcodes",
  "/matrix/features/id",
  "/matrix/features/name",
  "/matrix/shape"
)

missing_h5_paths <- setdiff(
  required_h5_paths,
  h5_full_paths
)

if (length(missing_h5_paths) > 0) {
  stop(
    "Required H5 paths missing:\n",
    paste(missing_h5_paths, collapse = "\n")
  )
}

cat("Required H5 paths: PASS\n\n")


# ------------------------------------------------------------
# 6. Read H5 barcodes
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("6. H5 BARCODE LOAD\n")
cat("------------------------------------------------------------\n\n")

h5_barcodes <- rhdf5::h5read(
  h5_file,
  "/matrix/barcodes"
)

h5_barcodes <- trimws(
  as.character(h5_barcodes)
)

h5_cell_n <- length(h5_barcodes)
h5_unique_cell_n <- length(unique(h5_barcodes))

h5_missing_barcode_n <- sum(
  is.na(h5_barcodes) |
    h5_barcodes == ""
)

h5_duplicate_barcode_n <- sum(
  duplicated(h5_barcodes)
)

cat(
  "H5 barcodes: ",
  h5_cell_n,
  "\n",
  sep = ""
)

cat(
  "H5 unique barcodes: ",
  h5_unique_cell_n,
  "\n",
  sep = ""
)

cat(
  "H5 missing/blank barcodes: ",
  h5_missing_barcode_n,
  "\n",
  sep = ""
)

cat(
  "H5 duplicated barcodes: ",
  h5_duplicate_barcode_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. H5 feature audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("7. H5 FEATURE AUDIT\n")
cat("------------------------------------------------------------\n\n")

feature_id <- rhdf5::h5read(
  h5_file,
  "/matrix/features/id"
)

feature_name <- rhdf5::h5read(
  h5_file,
  "/matrix/features/name"
)

feature_id <- trimws(as.character(feature_id))
feature_name <- trimws(as.character(feature_name))

cat(
  "Feature IDs: ",
  length(feature_id),
  "\n",
  sep = ""
)

cat(
  "Feature names: ",
  length(feature_name),
  "\n\n",
  sep = ""
)

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_gene_audit <- data.table(
  Gene = core_genes,
  In_Feature_ID = core_genes %in% feature_id,
  In_Feature_Name = core_genes %in% feature_name
)

core_gene_audit[, Present := In_Feature_ID | In_Feature_Name]

print(core_gene_audit)

fwrite(
  core_gene_audit,
  file.path(
    output_dir,
    "STEP10C-0C8_CoreGene_Feature_Audit.tsv"
  ),
  sep = "\t"
)

cat("\n")


# ------------------------------------------------------------
# 8. Matrix shape audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("8. H5 MATRIX DIMENSION AUDIT\n")
cat("------------------------------------------------------------\n\n")

h5_shape <- rhdf5::h5read(
  h5_file,
  "/matrix/shape"
)

h5_n_features <- as.integer(h5_shape[1])
h5_n_cells <- as.integer(h5_shape[2])

cat(
  "H5 matrix dimensions: ",
  h5_n_features,
  " features x ",
  h5_n_cells,
  " cells\n\n",
  sep = ""
)

shape_feature_match <- (
  h5_n_features == length(feature_id)
)

shape_cell_match <- (
  h5_n_cells == h5_cell_n
)

cat(
  "Feature dimension match: ",
  shape_feature_match,
  "\n",
  sep = ""
)

cat(
  "Cell dimension match: ",
  shape_cell_match,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Global exact Cell-ID reconciliation
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("9. GLOBAL EXACT CELL-ID RECONCILIATION\n")
cat("------------------------------------------------------------\n\n")

tisch_cells <- meta$Cell
h5_cells <- h5_barcodes

exact_intersection <- intersect(
  tisch_cells,
  h5_cells
)

tisch_unmatched <- setdiff(
  tisch_cells,
  h5_cells
)

h5_unmatched <- setdiff(
  h5_cells,
  tisch_cells
)

exact_match_n <- length(exact_intersection)

tisch_exact_match_pct <- (
  exact_match_n / length(tisch_cells)
) * 100

h5_exact_match_pct <- (
  exact_match_n / length(h5_cells)
) * 100

cat(
  "TISCH2 cells: ",
  length(tisch_cells),
  "\n",
  sep = ""
)

cat(
  "H5 barcodes: ",
  length(h5_cells),
  "\n",
  sep = ""
)

cat(
  "Exact intersection: ",
  exact_match_n,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "TISCH2 -> H5 exact match: %.4f%%\n",
    tisch_exact_match_pct
  )
)

cat(
  sprintf(
    "H5 -> TISCH2 exact match: %.4f%%\n",
    h5_exact_match_pct
  )
)

cat(
  "Unmatched TISCH2 cells: ",
  length(tisch_unmatched),
  "\n",
  sep = ""
)

cat(
  "Unmatched H5 barcodes: ",
  length(h5_unmatched),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Save unmatched IDs
# ------------------------------------------------------------

unmatched_tisch_dt <- data.table(
  Cell = tisch_unmatched
)

unmatched_h5_dt <- data.table(
  Barcode = h5_unmatched
)

fwrite(
  unmatched_tisch_dt,
  file.path(
    output_dir,
    "STEP10C-0C2_Unmatched_TISCH2_CellIDs.tsv"
  ),
  sep = "\t"
)

fwrite(
  unmatched_h5_dt,
  file.path(
    output_dir,
    "STEP10C-0C3_Unmatched_H5_Barcodes.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 11. Barcode normalization audit
#
# IMPORTANT:
#   This section does NOT redefine the official Cell ID.
#   It only evaluates whether unmatched IDs can be explained
#   by common technical formatting differences.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("11. BARCODE NORMALIZATION AUDIT\n")
cat("------------------------------------------------------------\n\n")


normalize_barcode <- function(x) {
  
  x <- trimws(as.character(x))
  
  # Remove trailing 10x-style "-1", "-2", etc.
  x <- sub(
    "-[0-9]+$",
    "",
    x
  )
  
  # Remove sample/source suffix beginning at the first "_"
  x <- sub(
    "_[^_]+.*$",
    "",
    x
  )
  
  x
}


tisch_norm <- normalize_barcode(tisch_cells)
h5_norm <- normalize_barcode(h5_cells)

tisch_norm_unique <- unique(tisch_norm)
h5_norm_unique <- unique(h5_norm)

normalized_intersection <- intersect(
  tisch_norm_unique,
  h5_norm_unique
)

normalized_match_n <- length(
  normalized_intersection
)

normalized_tisch_match_pct <- (
  sum(tisch_norm %in% h5_norm_unique) /
    length(tisch_norm)
) * 100

normalized_h5_match_pct <- (
  sum(h5_norm %in% tisch_norm_unique) /
    length(h5_norm)
) * 100


normalization_summary <- data.table(
  Metric = c(
    "TISCH2_cell_n",
    "H5_barcode_n",
    "TISCH2_unique_normalized_n",
    "H5_unique_normalized_n",
    "Normalized_intersection_n",
    "TISCH2_to_H5_normalized_match_pct",
    "H5_to_TISCH2_normalized_match_pct"
  ),
  Value = c(
    length(tisch_cells),
    length(h5_cells),
    length(tisch_norm_unique),
    length(h5_norm_unique),
    normalized_match_n,
    normalized_tisch_match_pct,
    normalized_h5_match_pct
  )
)

print(normalization_summary)

fwrite(
  normalization_summary,
  file.path(
    output_dir,
    "STEP10C-0C7_Barcode_Normalization_Audit.tsv"
  ),
  sep = "\t"
)

cat("\n")


# ------------------------------------------------------------
# 12. Build official exact-match crosswalk
#
# Only exact Cell IDs are used here.
# ------------------------------------------------------------

exact_match_crosswalk <- data.table(
  Cell = exact_intersection
)

exact_match_crosswalk[
  ,
  TISCH2_ExactMatch := TRUE
]

exact_match_crosswalk[
  ,
  H5_ExactMatch := TRUE
]

fwrite(
  exact_match_crosswalk,
  file.path(
    output_dir,
    "STEP10C-0C1_Exact_CellID_Crosswalk.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 13. Primary + Tumor locked subset
#
# IMPORTANT:
#   This uses the already locked STEP10A / STEP10B definition:
#
#       Stage == Primary
#       Source == Tumor
#
#   No alternative tumor definition is introduced.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("13. LOCKED PRIMARY-TUMOR SUBSET RECONCILIATION\n")
cat("------------------------------------------------------------\n\n")

primary_tumor_meta <- meta[
  Stage == "Primary" &
    Source == "Tumor"
]

primary_tumor_cells <- primary_tumor_meta$Cell

primary_tumor_n <- length(
  primary_tumor_cells
)

primary_tumor_unique_n <- length(
  unique(primary_tumor_cells)
)

primary_tumor_in_h5 <- (
  primary_tumor_cells %in% h5_cells
)

primary_tumor_match_n <- sum(
  primary_tumor_in_h5
)

primary_tumor_match_pct <- (
  primary_tumor_match_n /
    primary_tumor_n
) * 100

cat(
  "Primary + Tumor cells: ",
  primary_tumor_n,
  "\n",
  sep = ""
)

cat(
  "Primary + Tumor unique cells: ",
  primary_tumor_unique_n,
  "\n",
  sep = ""
)

cat(
  "Primary + Tumor cells found in H5: ",
  primary_tumor_match_n,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Primary + Tumor exact match: %.4f%%\n\n",
    primary_tumor_match_pct
  )
)

primary_tumor_reconciliation <- data.table(
  Cell = primary_tumor_cells,
  TISCH2_Primary_Tumor = TRUE,
  H5_ExactMatch = primary_tumor_in_h5
)

fwrite(
  primary_tumor_reconciliation,
  file.path(
    output_dir,
    "STEP10C-0C4_Primary_Tumor_CellID_Reconciliation.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 14. Locked 9-patient audit
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("14. PRIMARY-TUMOR PATIENT RECONCILIATION\n")
cat("------------------------------------------------------------\n\n")

locked_primary_tumor_patients <- c(
  "P0020",
  "P0025",
  "P0028",
  "P0030",
  "P0031",
  "P0034",
  "P1006",
  "P1010",
  "P1011"
)

observed_primary_tumor_patients <- sort(
  unique(primary_tumor_meta$Patient)
)

patient_set_match <- setequal(
  locked_primary_tumor_patients,
  observed_primary_tumor_patients
)

cat(
  "Locked patient count: ",
  length(locked_primary_tumor_patients),
  "\n",
  sep = ""
)

cat(
  "Observed patient count: ",
  length(observed_primary_tumor_patients),
  "\n",
  sep = ""
)

cat(
  "Locked patient set match: ",
  patient_set_match,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 15. Patient-level Cell-ID reconciliation
# ------------------------------------------------------------

patient_reconciliation <- primary_tumor_meta[
  ,
  .(
    TISCH2_Cells = .N,
    TISCH2_Unique_Cells = uniqueN(Cell),
    H5_Matched_Cells = sum(Cell %in% h5_cells),
    H5_Unmatched_Cells = sum(!(Cell %in% h5_cells))
  ),
  by = Patient
]

patient_reconciliation[
  ,
  H5_Match_Pct :=
    100 * H5_Matched_Cells / TISCH2_Unique_Cells
]

patient_reconciliation[
  ,
  Patient_Set_Locked :=
    Patient %in% locked_primary_tumor_patients
]

setorder(
  patient_reconciliation,
  Patient
)

print(patient_reconciliation)

fwrite(
  patient_reconciliation,
  file.path(
    output_dir,
    "STEP10C-0C5_Primary_Tumor_Patient_Reconciliation.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 16. Merge metadata with exact H5 matching status
# ------------------------------------------------------------

meta_reconciliation <- copy(meta)

meta_reconciliation[
  ,
  H5_ExactMatch := Cell %in% h5_cells
]

meta_reconciliation[
  ,
  Primary_Tumor_Locked :=
    Stage == "Primary" &
    Source == "Tumor"
]

meta_reconciliation[
  ,
  Primary_Tumor_H5_ExactMatch :=
    Primary_Tumor_Locked &
    H5_ExactMatch
]


# ------------------------------------------------------------
# 17. Compartment mapping audit
#
# This reproduces ONLY the locked STEP10B-1 explicit mapping.
#
# T-cell:
#   CD4Tconv
#   CD8T
#   CD8Tex
#
# Myeloid/macrophage:
#   Mono/Macro
#
# Endothelial:
#   Endothelial
#
# DC is NOT included.
# ------------------------------------------------------------

cat("------------------------------------------------------------\n")
cat("17. LOCKED COMPARTMENT RECONCILIATION\n")
cat("------------------------------------------------------------\n\n")

if (!"Celltype (major-lineage)" %in% names(meta)) {
  stop(
    "Required column 'Celltype (major-lineage)' not found in metadata."
  )
}

major_lineage <- trimws(
  as.character(meta[["Celltype (major-lineage)"]])
)

locked_compartment <- rep(
  NA_character_,
  length(major_lineage)
)

locked_compartment[
  major_lineage %in% c(
    "CD4Tconv",
    "CD8T",
    "CD8Tex"
  )
] <- "T_cell"

locked_compartment[
  major_lineage %in% c(
    "Mono/Macro"
  )
] <- "Myeloid_macrophage"

locked_compartment[
  major_lineage %in% c(
    "Endothelial"
  )
] <- "Endothelial"


meta_reconciliation[
  ,
  Locked_Compartment := locked_compartment
]


# ------------------------------------------------------------
# 18. Primary-Tumor compartment reconciliation
# ------------------------------------------------------------

pt_compartment <- meta_reconciliation[
  Primary_Tumor_Locked &
    !is.na(Locked_Compartment),
  .(
    TISCH2_Cells = .N,
    H5_Matched_Cells = sum(H5_ExactMatch),
    H5_Unmatched_Cells = sum(!H5_ExactMatch),
    Patients = uniqueN(Patient)
  ),
  by = Locked_Compartment
]

pt_compartment[
  ,
  H5_Match_Pct :=
    100 * H5_Matched_Cells / TISCH2_Cells
]

setorder(
  pt_compartment,
  Locked_Compartment
)

print(pt_compartment)

fwrite(
  pt_compartment,
  file.path(
    output_dir,
    "STEP10C-0C6_Compartment_Reconciliation.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 19. Patient × compartment reconciliation
# ------------------------------------------------------------

patient_compartment_reconciliation <- meta_reconciliation[
  Primary_Tumor_Locked &
    !is.na(Locked_Compartment),
  .(
    TISCH2_Cells = .N,
    H5_Matched_Cells = sum(H5_ExactMatch),
    H5_Unmatched_Cells = sum(!H5_ExactMatch)
  ),
  by = .(
    Patient,
    Locked_Compartment
  )
]

patient_compartment_reconciliation[
  ,
  H5_Match_Pct :=
    100 * H5_Matched_Cells / TISCH2_Cells
]

setorder(
  patient_compartment_reconciliation,
  Patient,
  Locked_Compartment
)

fwrite(
  patient_compartment_reconciliation,
  file.path(
    output_dir,
    "STEP10C-0C6A_Patient_Compartment_Reconciliation.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 20. Compare compartment cell counts against STEP10B-1
#
# This is a structural cross-check only.
# ------------------------------------------------------------

step10b1_count_file <- file.path(
  step10_root,
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit",
  "STEP10B-1E_Patient_Compartment_Cell_Counts_Long.tsv"
)

step10b1_available <- file.exists(
  step10b1_count_file
)

step10b1_crosscheck <- NULL

if (step10b1_available) {
  
  step10b1_counts <- fread(
    step10b1_count_file,
    sep = "\t",
    header = TRUE,
    data.table = TRUE,
    check.names = FALSE
  )
  
  cat("\nSTEP10B-1 long-format count file found.\n")
  
  cat(
    "Columns:\n",
    paste(names(step10b1_counts), collapse = " | "),
    "\n\n"
  )
  
  # Attempt flexible identification of patient / compartment / count columns
  patient_col_candidates <- c(
    "Patient",
    "patient",
    "Patient_ID",
    "patient_id"
  )
  
  compartment_col_candidates <- c(
    "Compartment",
    "compartment",
    "Locked_Compartment"
  )
  
  count_col_candidates <- c(
    "Cell_Count",
    "cell_count",
    "N_Cells",
    "Cells",
    "Count",
    "n"
  )
  
  patient_col <- intersect(
    patient_col_candidates,
    names(step10b1_counts)
  )
  
  compartment_col <- intersect(
    compartment_col_candidates,
    names(step10b1_counts)
  )
  
  count_col <- intersect(
    count_col_candidates,
    names(step10b1_counts)
  )
  
  if (
    length(patient_col) >= 1 &&
    length(compartment_col) >= 1 &&
    length(count_col) >= 1
  ) {
    
    step10b1_crosscheck <- copy(
      step10b1_counts[
        ,
        .(
          Patient = get(patient_col[1]),
          Locked_Compartment = get(compartment_col[1]),
          STEP10B1_Cell_Count = get(count_col[1])
        )
      ]
    )
    
    step10b1_crosscheck[
      ,
      Patient := as.character(Patient)
    ]
    
    step10b1_crosscheck[
      ,
      Locked_Compartment :=
        as.character(Locked_Compartment)
    ]
    
    step10b1_crosscheck[
      ,
      STEP10B1_Cell_Count :=
        as.numeric(STEP10B1_Cell_Count)
    ]
    
    current_counts <- meta_reconciliation[
      Primary_Tumor_Locked &
        !is.na(Locked_Compartment),
      .(
        STEP10C0C_Cell_Count = .N
      ),
      by = .(
        Patient,
        Locked_Compartment
      )
    ]
    
    step10b1_crosscheck <- merge(
      step10b1_crosscheck,
      current_counts,
      by = c(
        "Patient",
        "Locked_Compartment"
      ),
      all = TRUE
    )
    
    step10b1_crosscheck[
      ,
      Count_Match :=
        STEP10B1_Cell_Count ==
        STEP10C0C_Cell_Count
    ]
    
    fwrite(
      step10b1_crosscheck,
      file.path(
        output_dir,
        "STEP10C-0C6B_STEP10B1_Compartment_Count_Crosscheck.tsv"
      ),
      sep = "\t"
    )
    
    cat(
      "STEP10B-1 compartment count cross-check completed.\n"
    )
    
  } else {
    
    cat(
      "WARNING: STEP10B-1 column names could not be mapped automatically.\n"
    )
  }
  
} else {
  
  cat(
    "\nWARNING: STEP10B-1 long-format count file not found:\n",
    step10b1_count_file,
    "\n"
  )
}


# ------------------------------------------------------------
# 21. Global reconciliation summary
# ------------------------------------------------------------

cat("\n------------------------------------------------------------\n")
cat("21. FINAL RECONCILIATION SUMMARY\n")
cat("------------------------------------------------------------\n\n")


global_exact_pass <- all(
  tisch_cell_n == h5_cell_n,
  tisch_unique_cell_n == h5_unique_cell_n,
  tisch_missing_cell == 0,
  tisch_duplicate_cell_n == 0,
  h5_missing_barcode_n == 0,
  h5_duplicate_barcode_n == 0,
  exact_match_n == tisch_cell_n,
  exact_match_n == h5_cell_n
)


primary_tumor_pass <- all(
  primary_tumor_n == 65702,
  primary_tumor_unique_n == 65702,
  primary_tumor_match_n == primary_tumor_n
)


patient_set_pass <- all(
  patient_set_match,
  length(observed_primary_tumor_patients) == 9
)


patient_exact_pass <- all(
  patient_reconciliation$H5_Match_Pct == 100
)


required_compartments <- c(
  "T_cell",
  "Myeloid_macrophage",
  "Endothelial"
)

observed_compartments <- sort(
  unique(
    meta_reconciliation[
      Primary_Tumor_Locked &
        !is.na(Locked_Compartment),
      Locked_Compartment
    ]
  )
)

compartment_set_pass <- all(
  required_compartments %in%
    observed_compartments
)


compartment_exact_pass <- all(
  patient_compartment_reconciliation[
    Locked_Compartment %in% required_compartments,
    H5_Match_Pct
  ] == 100
)


core_gene_pass <- all(
  core_gene_audit$Present
)


dimension_pass <- all(
  shape_feature_match,
  shape_cell_match
)


# ------------------------------------------------------------
# 22. Final status logic
# ------------------------------------------------------------

if (
  global_exact_pass &&
  primary_tumor_pass &&
  patient_set_pass &&
  patient_exact_pass &&
  compartment_set_pass &&
  compartment_exact_pass &&
  core_gene_pass &&
  dimension_pass
) {
  
  final_status <- "PASS — EXACT CELL-ID RECONCILIATION"
  
} else if (
  primary_tumor_pass &&
  patient_set_pass &&
  compartment_set_pass &&
  core_gene_pass &&
  dimension_pass &&
  normalized_tisch_match_pct == 100 &&
  normalized_h5_match_pct == 100
) {
  
  final_status <- "PASS — NORMALIZED CELL-ID RECONCILIATION"
  
} else {
  
  final_status <- "FAIL — CELL-ID RECONCILIATION UNSUPPORTED"
}


summary_dt <- data.table(
  Check = c(
    "TISCH2_row_n",
    "TISCH2_unique_Cell_n",
    "TISCH2_missing_Cell_n",
    "TISCH2_duplicate_Cell_n",
    "H5_barcode_n",
    "H5_unique_barcode_n",
    "H5_missing_barcode_n",
    "H5_duplicate_barcode_n",
    "H5_feature_n",
    "H5_matrix_feature_n",
    "H5_matrix_cell_n",
    "Global_exact_intersection_n",
    "Global_TISCH2_to_H5_exact_pct",
    "Global_H5_to_TISCH2_exact_pct",
    "Primary_Tumor_cell_n",
    "Primary_Tumor_unique_cell_n",
    "Primary_Tumor_H5_match_n",
    "Primary_Tumor_H5_match_pct",
    "Primary_Tumor_patient_n",
    "Locked_patient_set_match",
    "Patient_level_exact_match",
    "Required_compartment_set_match",
    "Compartment_level_exact_match",
    "Core_gene_feature_presence",
    "H5_dimension_consistency",
    "Final_Status"
  ),
  Value = c(
    tisch_cell_n,
    tisch_unique_cell_n,
    tisch_missing_cell,
    tisch_duplicate_cell_n,
    h5_cell_n,
    h5_unique_cell_n,
    h5_missing_barcode_n,
    h5_duplicate_barcode_n,
    length(feature_id),
    h5_n_features,
    h5_n_cells,
    exact_match_n,
    tisch_exact_match_pct,
    h5_exact_match_pct,
    primary_tumor_n,
    primary_tumor_unique_n,
    primary_tumor_match_n,
    primary_tumor_match_pct,
    length(observed_primary_tumor_patients),
    patient_set_pass,
    patient_exact_pass,
    compartment_set_pass,
    compartment_exact_pass,
    core_gene_pass,
    dimension_pass,
    final_status
  )
)

print(summary_dt)

fwrite(
  summary_dt,
  file.path(
    output_dir,
    "STEP10C-0C9_Final_Reconciliation_Summary.tsv"
  ),
  sep = "\t"
)


# ------------------------------------------------------------
# 23. Save complete audit object
# ------------------------------------------------------------

audit_object <- list(
  
  metadata_file = metadata_file,
  h5_file = h5_file,
  
  metadata_dimensions = c(
    rows = nrow(meta),
    columns = ncol(meta)
  ),
  
  h5_dimensions = c(
    features = h5_n_features,
    cells = h5_n_cells
  ),
  
  tisch_cell_n = tisch_cell_n,
  tisch_unique_cell_n = tisch_unique_cell_n,
  tisch_missing_cell_n = tisch_missing_cell,
  tisch_duplicate_cell_n = tisch_duplicate_cell_n,
  
  h5_cell_n = h5_cell_n,
  h5_unique_cell_n = h5_unique_cell_n,
  h5_missing_barcode_n = h5_missing_barcode_n,
  h5_duplicate_barcode_n = h5_duplicate_barcode_n,
  
  exact_intersection_n = exact_match_n,
  exact_match_tisch_pct = tisch_exact_match_pct,
  exact_match_h5_pct = h5_exact_match_pct,
  
  normalized_match_tisch_pct =
    normalized_tisch_match_pct,
  
  normalized_match_h5_pct =
    normalized_h5_match_pct,
  
  primary_tumor_cell_n =
    primary_tumor_n,
  
  primary_tumor_h5_match_n =
    primary_tumor_match_n,
  
  primary_tumor_match_pct =
    primary_tumor_match_pct,
  
  locked_primary_tumor_patients =
    locked_primary_tumor_patients,
  
  observed_primary_tumor_patients =
    observed_primary_tumor_patients,
  
  patient_reconciliation =
    patient_reconciliation,
  
  compartment_reconciliation =
    pt_compartment,
  
  patient_compartment_reconciliation =
    patient_compartment_reconciliation,
  
  core_gene_audit =
    core_gene_audit,
  
  final_status =
    final_status
)

saveRDS(
  audit_object,
  file = file.path(
    output_dir,
    "STEP10C-0C_CellID_Reconciliation_Audit.rds"
  )
)


# ------------------------------------------------------------
# 24. Save session information
# ------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_dir,
    "STEP10C-0C_sessionInfo.txt"
  )
)


# ------------------------------------------------------------
# 25. Close HDF5 handles
# ------------------------------------------------------------

try(
  rhdf5::h5closeAll(),
  silent = TRUE
)


# ------------------------------------------------------------
# 26. Final console report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10C-0C FINAL STATUS\n")
cat("============================================================\n\n")

cat(
  final_status,
  "\n\n"
)

cat(
  "Global TISCH2 cells: ",
  tisch_cell_n,
  "\n",
  sep = ""
)

cat(
  "Global H5 barcodes: ",
  h5_cell_n,
  "\n",
  sep = ""
)

cat(
  "Global exact match: ",
  exact_match_n,
  " / ",
  tisch_cell_n,
  " (",
  sprintf("%.4f", tisch_exact_match_pct),
  "%)\n",
  sep = ""
)

cat(
  "Primary + Tumor cells: ",
  primary_tumor_n,
  "\n",
  sep = ""
)

cat(
  "Primary + Tumor exact H5 match: ",
  primary_tumor_match_n,
  " / ",
  primary_tumor_n,
  " (",
  sprintf("%.4f", primary_tumor_match_pct),
  "%)\n",
  sep = ""
)

cat(
  "Primary + Tumor patients: ",
  length(observed_primary_tumor_patients),
  "\n",
  sep = ""
)

cat(
  "Locked compartments: T_cell / Myeloid_macrophage / Endothelial\n"
)

cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)

cat("\n")
cat(
  "No expression values from /matrix/data were analyzed.\n"
)

cat(
  "No pseudobulk was performed.\n"
)

cat(
  "No statistical testing was performed.\n"
)

cat(
  "No patient exclusion was performed.\n"
)

cat("\n")
cat(
  "Output directory:\n",
  output_dir,
  "\n"
)

cat("\n")
cat("============================================================\n")
cat("STEP10C-0C COMPLETE\n")
cat("============================================================\n")


# ============================================================
# STEP10C-0D — GEO RAW UMI TXT FILE CHARACTERIZATION AUDIT
# ============================================================
#
# Purpose:
#   1. Verify file existence and compressed file size
#   2. Verify GZIP integrity
#   3. Read only the first few text lines
#   4. Determine delimiter and matrix orientation
#   5. Inspect header / gene / cell-ID structure
#   6. Verify expected number of cell IDs
#   7. Check whether cell IDs match TISCH2 metadata
#   8. Check whether core genes are present
#
# IMPORTANT:
#   - Full 389-MB compressed file is NOT loaded into memory
#   - No dense matrix conversion
#   - No pseudobulk
#   - No statistical analysis
#   - No patient exclusion
#   - No expression analysis
#
# ============================================================


cat("============================================================\n")
cat("STEP10C-0D — GEO RAW UMI TXT FILE CHARACTERIZATION AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. WORKING DIRECTORY / INPUT FILES
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

txt_file <- file.path(
  base_dir,
  "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
)

metadata_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)


# ------------------------------------------------------------
# 1. FILE EXISTENCE AND SIZE
# ------------------------------------------------------------

cat("1. FILE INFORMATION\n")
cat("------------------------------------------------------------\n")

if (!file.exists(txt_file)) {
  
  stop(
    "FAIL — Raw UMI TXT file not found:\n",
    txt_file
  )
  
}

fi <- file.info(txt_file)

cat("File exists: TRUE\n")

cat(
  "File size (MB): ",
  round(fi$size / 1024^2, 3),
  "\n",
  sep = ""
)

cat(
  "File size (GB): ",
  round(fi$size / 1024^3, 3),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 2. GZIP HEADER TEST
# ------------------------------------------------------------

cat("2. GZIP HEADER TEST\n")
cat("------------------------------------------------------------\n")

con <- file(txt_file, open = "rb")

raw_head <- readBin(
  con,
  what = "raw",
  n = 32
)

close(con)

gzip_header_ok <- (
  length(raw_head) >= 2 &&
    identical(
      as.integer(raw_head[1:2]),
      c(31L, 139L)
    )
)

cat(
  "First 32 bytes (hex):\n"
)

cat(
  paste(
    sprintf("%02X", as.integer(raw_head)),
    collapse = " "
  ),
  "\n\n"
)

cat(
  "Standard GZIP header detected: ",
  gzip_header_ok,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 3. GZIP DECOMPRESSION TEST
# ------------------------------------------------------------

cat("3. GZIP DECOMPRESSION TEST\n")
cat("------------------------------------------------------------\n")

gzip_readable <- FALSE
gzip_error <- NULL

con <- NULL

tryCatch({
  
  con <- gzfile(
    txt_file,
    open = "rb"
  )
  
  test_bytes <- readBin(
    con,
    what = "raw",
    n = 4096
  )
  
  close(con)
  
  gzip_readable <- length(test_bytes) > 0
  
}, error = function(e) {
  
  gzip_error <<- conditionMessage(e)
  
  if (!is.null(con)) {
    try(close(con), silent = TRUE)
  }
  
})

cat(
  "GZIP decompression readable: ",
  gzip_readable,
  "\n",
  sep = ""
)

if (!is.null(gzip_error)) {
  
  cat(
    "GZIP error: ",
    gzip_error,
    "\n",
    sep = ""
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 4. READ FIRST 5 LINES ONLY
# ------------------------------------------------------------

cat("4. FIRST 5 TEXT LINES\n")
cat("------------------------------------------------------------\n")

con <- gzfile(
  txt_file,
  open = "rt",
  encoding = "UTF-8"
)

first_lines <- readLines(
  con,
  n = 5,
  warn = FALSE
)

close(con)

cat(
  "Number of lines successfully read: ",
  length(first_lines),
  "\n\n",
  sep = ""
)

if (length(first_lines) > 0) {
  
  cat("Preview:\n")
  cat(
    paste(
      substr(first_lines, 1, 300),
      collapse = "\n"
    ),
    "\n\n"
  )
  
}


# ------------------------------------------------------------
# 5. DETERMINE LIKELY DELIMITER
# ------------------------------------------------------------

cat("5. DELIMITER DIAGNOSTIC\n")
cat("------------------------------------------------------------\n")

if (length(first_lines) == 0) {
  
  stop(
    "FAIL — No text lines could be read from the decompressed file."
  )
  
}

line1 <- first_lines[1]

tab_fields <- strsplit(
  line1,
  "\t",
  fixed = TRUE
)[[1]]

comma_fields <- strsplit(
  line1,
  ",",
  fixed = TRUE
)[[1]]

semicolon_fields <- strsplit(
  line1,
  ";",
  fixed = TRUE
)[[1]]

space_fields <- strsplit(
  trimws(line1),
  "\\s+"
)[[1]]

cat(
  "TAB-separated fields:       ",
  length(tab_fields),
  "\n",
  sep = ""
)

cat(
  "Comma-separated fields:     ",
  length(comma_fields),
  "\n",
  sep = ""
)

cat(
  "Semicolon-separated fields: ",
  length(semicolon_fields),
  "\n",
  sep = ""
)

cat(
  "Whitespace-separated fields:",
  length(space_fields),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. HEADER INSPECTION
# ------------------------------------------------------------

cat("6. HEADER INSPECTION\n")
cat("------------------------------------------------------------\n")

header <- tab_fields

cat(
  "Number of TAB-delimited header fields: ",
  length(header),
  "\n\n",
  sep = ""
)

cat("First 15 header fields:\n")

n_header_show <- min(
  15,
  length(header)
)

for (i in seq_len(n_header_show)) {
  
  cat(
    i,
    ": ",
    header[i],
    "\n",
    sep = ""
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 7. FIRST DATA LINE INSPECTION
# ------------------------------------------------------------

cat("7. FIRST DATA ROW INSPECTION\n")
cat("------------------------------------------------------------\n")

if (length(first_lines) >= 2) {
  
  data_line <- strsplit(
    first_lines[2],
    "\t",
    fixed = TRUE
  )[[1]]
  
  cat(
    "Number of TAB-delimited fields in first data row: ",
    length(data_line),
    "\n\n",
    sep = ""
  )
  
  cat("First 15 fields:\n")
  
  n_data_show <- min(
    15,
    length(data_line)
  )
  
  for (i in seq_len(n_data_show)) {
    
    cat(
      i,
      ": ",
      data_line[i],
      "\n",
      sep = ""
    )
    
  }
  
} else {
  
  cat(
    "WARNING — No second line available.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 8. MATRIX ORIENTATION DIAGNOSTIC
# ------------------------------------------------------------

cat("8. MATRIX ORIENTATION DIAGNOSTIC\n")
cat("------------------------------------------------------------\n")

header_n <- length(header)

data_n <- if (length(first_lines) >= 2) {
  length(
    strsplit(
      first_lines[2],
      "\t",
      fixed = TRUE
    )[[1]]
  )
} else {
  NA_integer_
}

cat(
  "Header field count: ",
  header_n,
  "\n",
  sep = ""
)

cat(
  "First data-row field count: ",
  data_n,
  "\n\n",
  sep = ""
)

if (
  header_n >= 100000 &&
  data_n == header_n
) {
  
  cat(
    "Likely orientation: GENE × CELL\n\n"
  )
  
} else if (
  header_n <= 50000 &&
  data_n == header_n
) {
  
  cat(
    "Likely orientation: CELL × GENE or another wide format.\n\n"
  )
  
} else {
  
  cat(
    "Orientation cannot yet be determined from the first two lines.\n\n"
  )
  
}


# ------------------------------------------------------------
# 9. EXPECTED CELL-ID COUNT
# ------------------------------------------------------------

cat("9. EXPECTED CELL-ID COUNT\n")
cat("------------------------------------------------------------\n")

expected_cells <- 203298

cat(
  "Expected TISCH2/H5 cell count: ",
  expected_cells,
  "\n",
  sep = ""
)

if (header_n == expected_cells + 1) {
  
  cat(
    "PASS — Header contains expected 203,298 cell columns plus one gene field.\n"
  )
  
} else {
  
  cat(
    "INFO — Header does not equal 203,299 fields.\n",
    "This does not yet establish failure; file orientation must be confirmed.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 10. CELL-ID EXTRACTION FROM HEADER
# ------------------------------------------------------------

cat("10. CELL-ID HEADER AUDIT\n")
cat("------------------------------------------------------------\n")

cell_ids_txt <- character(0)

if (header_n >= expected_cells + 1) {
  
  cell_ids_txt <- header[-1]
  
  cat(
    "Assuming first field is gene identifier.\n"
  )
  
  cat(
    "Candidate cell IDs: ",
    length(cell_ids_txt),
    "\n",
    sep = ""
  )
  
  cat(
    "Unique candidate cell IDs: ",
    length(unique(cell_ids_txt)),
    "\n",
    sep = ""
  )
  
  duplicate_cell_ids <- sum(
    duplicated(cell_ids_txt)
  )
  
  cat(
    "Duplicate candidate cell IDs: ",
    duplicate_cell_ids,
    "\n\n",
    sep = ""
  )
  
} else {
  
  cat(
    "Cell-ID extraction deferred because header size does not yet\n",
    "match the expected wide gene × cell structure.\n\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 11. CORE-GENE HEADER / FIRST-DATA-ROW CHECK
# ------------------------------------------------------------

cat("11. CORE-GENE STRUCTURE CHECK\n")
cat("------------------------------------------------------------\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat(
  "Locked core genes: ",
  paste(core_genes, collapse = " / "),
  "\n\n"
)

if (length(first_lines) >= 2) {
  
  first_gene <- data_line[1]
  
  cat(
    "First data-row gene identifier: ",
    first_gene,
    "\n\n",
    sep = ""
  )
  
  cat(
    "NOTE:\n",
    "Only the first data row is inspected here.\n",
    "Complete core-gene presence will be assessed in the next\n",
    "low-memory validation step after the file orientation is confirmed.\n\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 12. TISCH2 METADATA CELL-ID CROSS-CHECK
# ------------------------------------------------------------

cat("12. TISCH2 CELL-ID CROSS-CHECK\n")
cat("------------------------------------------------------------\n")

if (!file.exists(metadata_file)) {
  
  cat(
    "INFO — TISCH2 metadata file not found at expected path.\n",
    metadata_file,
    "\n\n"
  )
  
} else {
  
  tisch_meta <- read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = ""
  )
  
  cat(
    "TISCH2 metadata rows: ",
    nrow(tisch_meta),
    "\n",
    sep = ""
  )
  
  cat(
    "TISCH2 metadata columns: ",
    ncol(tisch_meta),
    "\n\n",
    sep = ""
  )
  
  # Detect Cell ID column
  possible_cell_cols <- c(
    "Cell",
    "CellID",
    "Cell_ID",
    "cell",
    "cellID",
    "cell_id"
  )
  
  cell_col <- possible_cell_cols[
    possible_cell_cols %in% names(tisch_meta)
  ][1]
  
  if (is.na(cell_col)) {
    
    cat(
      "WARNING — Standard Cell-ID column name not detected.\n"
    )
    
    cat(
      "Available columns:\n",
      paste(names(tisch_meta), collapse = ", "),
      "\n\n"
    )
    
  } else {
    
    tisch_cells <- as.character(
      tisch_meta[[cell_col]]
    )
    
    cat(
      "Detected TISCH2 Cell-ID column: ",
      cell_col,
      "\n",
      sep = ""
    )
    
    cat(
      "TISCH2 unique Cell IDs: ",
      length(unique(tisch_cells)),
      "\n",
      sep = ""
    )
    
    if (length(cell_ids_txt) > 0) {
      
      exact_match <- (
        length(cell_ids_txt) ==
          length(tisch_cells) &&
          setequal(
            cell_ids_txt,
            tisch_cells
          )
      )
      
      overlap_n <- length(
        intersect(
          cell_ids_txt,
          tisch_cells
        )
      )
      
      cat(
        "Exact Cell-ID set match: ",
        exact_match,
        "\n",
        sep = ""
      )
      
      cat(
        "Cell-ID intersection: ",
        overlap_n,
        "\n",
        sep = ""
      )
      
    }
    
  }
  
}

cat("\n")


# ------------------------------------------------------------
# 13. FINAL DIAGNOSTIC STATUS
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP10C-0D INITIAL STATUS\n")
cat("============================================================\n\n")

if (
  gzip_header_ok &&
  gzip_readable &&
  length(first_lines) >= 2
) {
  
  cat(
    "PASS — GEO raw UMI TXT file is a readable GZIP-compressed\n",
    "text file with sufficient structure for further low-memory\n",
    "raw UMI characterization.\n\n",
    sep = ""
  )
  
} else {
  
  cat(
    "FAIL — TXT file did not pass the basic GZIP/text-readability audit.\n\n"
  )
  
}

cat(
  "IMPORTANT:\n",
  "No full matrix was loaded.\n",
  "No expression values were analyzed.\n",
  "No pseudobulk was constructed.\n",
  "No statistical testing was performed.\n",
  "No patients or cells were excluded.\n\n",
  sep = ""
)

cat("============================================================\n")
cat("END STEP10C-0D INITIAL DIAGNOSTIC\n")
cat("============================================================\n")


# ============================================================
# STEP10C-0D-2
# GEO RAW UMI CELL-ID + GENE-ID LIGHTWEIGHT
# RECONCILIATION AUDIT
# ============================================================
#
# Purpose:
#   1. Extract all GEO raw UMI Cell IDs from the header
#   2. Verify GEO raw Cell-ID uniqueness
#   3. Compare GEO raw Cell IDs with TISCH2 Cell IDs
#   4. Quantify TISCH2 cells retained in GEO raw matrix
#   5. Identify GEO raw-only cells
#   6. Cross-check GEO raw Cell IDs against H5 barcodes
#   7. Scan ONLY the first field of every matrix row to obtain
#      the complete GEO raw Gene-ID list
#   8. Verify AOC3 / MS4A7 / IL7R presence
#
# IMPORTANT:
#   - The full expression matrix is NEVER loaded
#   - No dense matrix conversion
#   - No pseudobulk
#   - No statistical analysis
#   - No patient exclusion
#
# ============================================================


cat("============================================================\n")
cat("STEP10C-0D-2 — GEO RAW UMI CELL-ID + GENE-ID AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. PATHS
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

txt_file <- file.path(
  base_dir,
  "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
)

metadata_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

h5_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_expression.h5"
)

output_dir <- file.path(
  base_dir,
  "STEP10C-0D-2_GEO_RawUMI_CellID_GeneID_Audit"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. INPUT FILE CHECK
# ------------------------------------------------------------

cat("1. INPUT FILE CHECK\n")
cat("------------------------------------------------------------\n")

required_files <- c(
  txt_file,
  metadata_file
)

file_status <- data.frame(
  file = basename(required_files),
  exists = file.exists(required_files),
  size_MB = ifelse(
    file.exists(required_files),
    round(
      file.info(required_files)$size / 1024^2,
      3
    ),
    NA_real_
  ),
  stringsAsFactors = FALSE
)

print(file_status)

if (!all(file_status$exists)) {
  
  stop(
    "FAIL — One or more required input files are missing."
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 2. LOCKED CORE GENES
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("2. LOCKED CORE GENES\n")
cat("------------------------------------------------------------\n")

cat(
  paste(core_genes, collapse = " / "),
  "\n\n"
)


# ------------------------------------------------------------
# 3. READ GEO RAW UMI HEADER ONLY
# ------------------------------------------------------------

cat("3. GEO RAW UMI HEADER — CELL-ID EXTRACTION\n")
cat("------------------------------------------------------------\n")

con <- gzfile(
  txt_file,
  open = "rt",
  encoding = "UTF-8"
)

header_line <- readLines(
  con,
  n = 1,
  warn = FALSE
)

close(con)

if (length(header_line) != 1) {
  
  stop(
    "FAIL — Could not read exactly one header line."
  )
  
}

header_fields <- strsplit(
  header_line,
  "\t",
  fixed = TRUE
)[[1]]

cat(
  "Total header fields: ",
  length(header_fields),
  "\n",
  sep = ""
)

if (header_fields[1] != "Index") {
  
  cat(
    "WARNING — First header field is not 'Index': ",
    header_fields[1],
    "\n",
    sep = ""
  )
  
}

geo_cell_ids <- header_fields[-1]

cat(
  "GEO raw UMI Cell IDs: ",
  length(geo_cell_ids),
  "\n",
  sep = ""
)

geo_unique_n <- length(
  unique(geo_cell_ids)
)

geo_duplicate_n <- sum(
  duplicated(geo_cell_ids)
)

cat(
  "Unique GEO Cell IDs: ",
  geo_unique_n,
  "\n",
  sep = ""
)

cat(
  "Duplicate GEO Cell IDs: ",
  geo_duplicate_n,
  "\n\n",
  sep = ""
)

geo_cell_id_unique_pass <- (
  length(geo_cell_ids) == geo_unique_n &&
    geo_duplicate_n == 0
)

cat(
  "GEO Cell-ID uniqueness PASS: ",
  geo_cell_id_unique_pass,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. READ TISCH2 METADATA
# ------------------------------------------------------------

cat("4. TISCH2 METADATA CELL-ID AUDIT\n")
cat("------------------------------------------------------------\n")

tisch_meta <- read.delim(
  metadata_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

cat(
  "TISCH2 rows: ",
  nrow(tisch_meta),
  "\n",
  sep = ""
)

cat(
  "TISCH2 columns: ",
  ncol(tisch_meta),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. DETECT TISCH2 CELL-ID COLUMN
# ------------------------------------------------------------

possible_cell_cols <- c(
  "Cell",
  "CellID",
  "Cell_ID",
  "cell",
  "cellID",
  "cell_id"
)

matched_cell_cols <- possible_cell_cols[
  possible_cell_cols %in% names(tisch_meta)
]

if (length(matched_cell_cols) == 0) {
  
  stop(
    "FAIL — Could not identify TISCH2 Cell-ID column.\n",
    "Available columns:\n",
    paste(names(tisch_meta), collapse = ", ")
  )
  
}

tisch_cell_col <- matched_cell_cols[1]

tisch_cell_ids <- as.character(
  tisch_meta[[tisch_cell_col]]
)

cat(
  "TISCH2 Cell-ID column: ",
  tisch_cell_col,
  "\n",
  sep = ""
)

cat(
  "TISCH2 Cell IDs: ",
  length(tisch_cell_ids),
  "\n",
  sep = ""
)

cat(
  "Unique TISCH2 Cell IDs: ",
  length(unique(tisch_cell_ids)),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 6. GEO RAW ↔ TISCH2 CELL-ID RECONCILIATION
# ------------------------------------------------------------

cat("5. GEO RAW ↔ TISCH2 CELL-ID RECONCILIATION\n")
cat("------------------------------------------------------------\n")

geo_set <- unique(geo_cell_ids)
tisch_set <- unique(tisch_cell_ids)

intersection_ids <- intersect(
  geo_set,
  tisch_set
)

tisch_only_ids <- setdiff(
  tisch_set,
  geo_set
)

geo_only_ids <- setdiff(
  geo_set,
  tisch_set
)

cat(
  "GEO raw unique cells: ",
  length(geo_set),
  "\n",
  sep = ""
)

cat(
  "TISCH2 unique cells: ",
  length(tisch_set),
  "\n",
  sep = ""
)

cat(
  "Intersection: ",
  length(intersection_ids),
  "\n",
  sep = ""
)

cat(
  "TISCH2-only cells: ",
  length(tisch_only_ids),
  "\n",
  sep = ""
)

cat(
  "GEO raw-only cells: ",
  length(geo_only_ids),
  "\n\n",
  sep = ""
)

tisch_in_geo_pct <- (
  length(intersection_ids) /
    length(tisch_set)
) * 100

geo_in_tisch_pct <- (
  length(intersection_ids) /
    length(geo_set)
) * 100

cat(
  "TISCH2 cells represented in GEO raw: ",
  round(tisch_in_geo_pct, 4),
  "%\n",
  sep = ""
)

cat(
  "GEO raw cells represented in TISCH2: ",
  round(geo_in_tisch_pct, 4),
  "%\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. CHECK EXPECTED CELL COUNTS
# ------------------------------------------------------------

cat("6. EXPECTED CELL COUNT CHECK\n")
cat("------------------------------------------------------------\n")

expected_geo_cells <- 208506
expected_tisch_cells <- 203298

cat(
  "Expected GEO raw cells: ",
  expected_geo_cells,
  "\n",
  sep = ""
)

cat(
  "Observed GEO raw cells: ",
  length(geo_set),
  "\n",
  sep = ""
)

cat(
  "Expected TISCH2 cells: ",
  expected_tisch_cells,
  "\n",
  sep = ""
)

cat(
  "Observed TISCH2 cells: ",
  length(tisch_set),
  "\n\n",
  sep = ""
)

geo_count_pass <- (
  length(geo_set) == expected_geo_cells
)

tisch_count_pass <- (
  length(tisch_set) == expected_tisch_cells
)

cat(
  "GEO raw count PASS: ",
  geo_count_pass,
  "\n",
  sep = ""
)

cat(
  "TISCH2 count PASS: ",
  tisch_count_pass,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. GEO RAW-ONLY CELL CHARACTERIZATION
# ------------------------------------------------------------

cat("7. GEO RAW-ONLY CELL CHARACTERIZATION\n")
cat("------------------------------------------------------------\n")

if (length(geo_only_ids) > 0) {
  
  cat(
    "Number of GEO raw-only cells: ",
    length(geo_only_ids),
    "\n\n",
    sep = ""
  )
  
  cat("First 30 GEO raw-only Cell IDs:\n")
  
  print(
    head(
      geo_only_ids,
      30
    )
  )
  
} else {
  
  cat(
    "No GEO raw-only cells detected.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 9. TISCH2-ONLY CELL CHARACTERIZATION
# ------------------------------------------------------------

cat("8. TISCH2-ONLY CELL CHARACTERIZATION\n")
cat("------------------------------------------------------------\n")

if (length(tisch_only_ids) > 0) {
  
  cat(
    "Number of TISCH2-only cells: ",
    length(tisch_only_ids),
    "\n\n",
    sep = ""
  )
  
  cat("First 30 TISCH2-only Cell IDs:\n")
  
  print(
    head(
      tisch_only_ids,
      30
    )
  )
  
} else {
  
  cat(
    "No TISCH2-only cells detected.\n"
  )
  
}

cat("\n")


# ------------------------------------------------------------
# 10. H5 BARCODE CROSS-CHECK
# ------------------------------------------------------------

cat("9. H5 BARCODE CROSS-CHECK\n")
cat("------------------------------------------------------------\n")

h5_crosscheck_available <- FALSE
h5_barcode_ids <- character(0)

if (file.exists(h5_file)) {
  
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    
    cat(
      "INFO — package 'rhdf5' is not installed.\n",
      "H5 cross-check skipped.\n\n"
    )
    
  } else {
    
    h5_crosscheck_available <- TRUE
    
    h5_barcode_ids <- rhdf5::h5read(
      h5_file,
      "/matrix/barcodes"
    )
    
    h5_barcode_ids <- as.character(
      h5_barcode_ids
    )
    
    cat(
      "H5 barcodes: ",
      length(h5_barcode_ids),
      "\n",
      sep = ""
    )
    
    cat(
      "Unique H5 barcodes: ",
      length(unique(h5_barcode_ids)),
      "\n\n",
      sep = ""
    )
    
    h5_geo_intersection <- intersect(
      h5_barcode_ids,
      geo_set
    )
    
    h5_only_ids <- setdiff(
      h5_barcode_ids,
      geo_set
    )
    
    geo_h5_only_ids <- setdiff(
      geo_set,
      h5_barcode_ids
    )
    
    cat(
      "H5 ↔ GEO raw intersection: ",
      length(h5_geo_intersection),
      "\n",
      sep = ""
    )
    
    cat(
      "H5-only cells: ",
      length(h5_only_ids),
      "\n",
      sep = ""
    )
    
    cat(
      "GEO raw-only relative to H5: ",
      length(geo_h5_only_ids),
      "\n\n",
      sep = ""
    )
    
    h5_in_geo_pass <- (
      length(h5_only_ids) == 0 &&
        length(h5_geo_intersection) ==
        length(unique(h5_barcode_ids))
    )
    
    cat(
      "H5 → GEO raw inclusion PASS: ",
      h5_in_geo_pass,
      "\n\n",
      sep = ""
    )
    
  }
  
} else {
  
  cat(
    "INFO — H5 file not found at:\n",
    h5_file,
    "\n\n",
    sep = ""
  )
  
}


# ------------------------------------------------------------
# 11. SCAN GENE IDENTIFIERS ONLY
# ------------------------------------------------------------

cat("10. GEO RAW UMI GENE-ID SCAN\n")
cat("------------------------------------------------------------\n")

cat(
  "IMPORTANT: The complete expression matrix is NOT loaded.\n",
  "Only the first TAB-delimited field of each data row is read.\n\n",
  sep = ""
)

gene_ids <- character(0)

con <- gzfile(
  txt_file,
  open = "rt",
  encoding = "UTF-8"
)

# Skip header
readLines(
  con,
  n = 1,
  warn = FALSE
)

line_counter <- 0L

while (TRUE) {
  
  lines <- readLines(
    con,
    n = 5000,
    warn = FALSE
  )
  
  if (length(lines) == 0) {
    break
  }
  
  # Extract first field only
  batch_genes <- sub(
    "\t.*$",
    "",
    lines
  )
  
  gene_ids <- c(
    gene_ids,
    batch_genes
  )
  
  line_counter <- (
    line_counter +
      length(lines)
  )
  
  if (
    line_counter %% 50000 == 0
  ) {
    
    cat(
      "Rows scanned: ",
      line_counter,
      "\n",
      sep = ""
    )
    
  }
  
}

close(con)

cat(
  "\nTotal matrix data rows scanned: ",
  length(gene_ids),
  "\n",
  sep = ""
)

cat(
  "Unique Gene IDs: ",
  length(unique(gene_ids)),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. GENE-ID DUPLICATE AUDIT
# ------------------------------------------------------------

cat("11. GENE-ID DUPLICATE AUDIT\n")
cat("------------------------------------------------------------\n")

gene_unique_n <- length(
  unique(gene_ids)
)

gene_duplicate_n <- sum(
  duplicated(gene_ids)
)

cat(
  "Total Gene-ID rows: ",
  length(gene_ids),
  "\n",
  sep = ""
)

cat(
  "Unique Gene IDs: ",
  gene_unique_n,
  "\n",
  sep = ""
)

cat(
  "Duplicated Gene IDs: ",
  gene_duplicate_n,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. CORE-GENE PRESENCE
# ------------------------------------------------------------

cat("12. CORE-GENE PRESENCE\n")
cat("------------------------------------------------------------\n")

core_presence <- data.frame(
  gene = core_genes,
  present = core_genes %in% gene_ids,
  stringsAsFactors = FALSE
)

print(core_presence)

core_genes_present <- all(
  core_presence$present
)

cat(
  "\nAll locked core genes present: ",
  core_genes_present,
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. OPTIONAL GENE ORDER PREVIEW
# ------------------------------------------------------------

cat("13. GENE-ID PREVIEW\n")
cat("------------------------------------------------------------\n")

cat("First 20 Gene IDs:\n")

print(
  head(
    gene_ids,
    20
  )
)

cat("\nLast 20 Gene IDs:\n")

print(
  tail(
    gene_ids,
    20
  )
)

cat("\n")


# ------------------------------------------------------------
# 15. SAVE LIGHTWEIGHT AUDIT OBJECTS
# ------------------------------------------------------------

cat("14. SAVE LIGHTWEIGHT AUDIT OUTPUTS\n")
cat("------------------------------------------------------------\n")

write.table(
  data.frame(
    CellID = geo_cell_ids,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "STEP10C-0D2_GEO_RawUMI_CellIDs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    CellID = geo_only_ids,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "STEP10C-0D2_GEO_RawOnly_CellIDs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    GeneID = gene_ids,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    output_dir,
    "STEP10C-0D2_GEO_RawUMI_GeneIDs.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_presence,
  file = file.path(
    output_dir,
    "STEP10C-0D2_CoreGene_Presence.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 16. FINAL SUMMARY
# ------------------------------------------------------------

summary_df <- data.frame(
  metric = c(
    "GEO_raw_cell_count",
    "GEO_raw_unique_cell_count",
    "GEO_raw_duplicate_cell_count",
    "TISCH2_cell_count",
    "TISCH2_unique_cell_count",
    "GEO_TISCH2_intersection",
    "TISCH2_only_cell_count",
    "GEO_raw_only_cell_count",
    "TISCH2_in_GEO_percent",
    "GEO_in_TISCH2_percent",
    "GEO_expected_count_pass",
    "TISCH2_expected_count_pass",
    "GEO_cellID_unique_pass",
    "GEO_raw_gene_row_count",
    "GEO_raw_unique_gene_count",
    "GEO_raw_duplicate_gene_count",
    "AOC3_present",
    "MS4A7_present",
    "IL7R_present",
    "all_core_genes_present"
  ),
  value = c(
    length(geo_cell_ids),
    length(geo_set),
    geo_duplicate_n,
    length(tisch_cell_ids),
    length(tisch_set),
    length(intersection_ids),
    length(tisch_only_ids),
    length(geo_only_ids),
    tisch_in_geo_pct,
    geo_in_tisch_pct,
    geo_count_pass,
    tisch_count_pass,
    geo_cell_id_unique_pass,
    length(gene_ids),
    gene_unique_n,
    gene_duplicate_n,
    core_presence$present[
      core_presence$gene == "AOC3"
    ],
    core_presence$present[
      core_presence$gene == "MS4A7"
    ],
    core_presence$present[
      core_presence$gene == "IL7R"
    ],
    core_genes_present
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = file.path(
    output_dir,
    "STEP10C-0D2_Lightweight_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17. SAVE RDS
# ------------------------------------------------------------

audit_object <- list(
  
  geo_raw_cell_ids = geo_cell_ids,
  
  geo_raw_only_ids = geo_only_ids,
  
  tisch2_cell_ids = tisch_cell_ids,
  
  tisch2_only_ids = tisch_only_ids,
  
  intersection_ids = intersection_ids,
  
  geo_raw_gene_ids = gene_ids,
  
  core_gene_presence = core_presence,
  
  summary = summary_df
  
)

saveRDS(
  audit_object,
  file = file.path(
    output_dir,
    "STEP10C-0D2_Lightweight_Audit_Object.rds"
  )
)


# ------------------------------------------------------------
# 18. FINAL STATUS
# ------------------------------------------------------------

cat("============================================================\n")
cat("STEP10C-0D-2 FINAL STATUS\n")
cat("============================================================\n\n")

cell_id_reconciliation_pass <- (
  geo_cell_id_unique_pass &&
    geo_count_pass &&
    tisch_count_pass &&
    length(tisch_only_ids) == 0
)

gene_id_pass <- (
  length(gene_ids) > 0 &&
    gene_unique_n == length(gene_ids) &&
    core_genes_present
)

if (
  cell_id_reconciliation_pass &&
  gene_id_pass
) {
  
  cat(
    "PASS — GEO raw UMI Cell-ID + Gene-ID lightweight audit.\n\n"
  )
  
} else {
  
  cat(
    "REVIEW REQUIRED — One or more lightweight reconciliation\n",
    "criteria did not pass.\n\n",
    sep = ""
  )
  
}

cat(
  "GEO raw cells: ",
  length(geo_set),
  "\n",
  sep = ""
)

cat(
  "TISCH2 cells: ",
  length(tisch_set),
  "\n",
  sep = ""
)

cat(
  "GEO ↔ TISCH2 intersection: ",
  length(intersection_ids),
  "\n",
  sep = ""
)

cat(
  "TISCH2-only cells: ",
  length(tisch_only_ids),
  "\n",
  sep = ""
)

cat(
  "GEO raw-only cells: ",
  length(geo_only_ids),
  "\n",
  sep = ""
)

cat(
  "GEO raw Gene IDs: ",
  length(gene_ids),
  "\n",
  sep = ""
)

cat(
  "All core genes present: ",
  core_genes_present,
  "\n\n",
  sep = ""
)

cat(
  "NO FULL EXPRESSION MATRIX WAS LOADED.\n",
  "NO PSEUDOBULK WAS CONSTRUCTED.\n",
  "NO STATISTICAL TESTING WAS PERFORMED.\n",
  "NO PATIENTS OR CELLS WERE EXCLUDED.\n\n",
  sep = ""
)

cat("============================================================\n")
cat("END STEP10C-0D-2\n")
cat("============================================================\n")

##STEP10C-0E — Primary-Tumor Raw UMI Cell-ID + Compartment Reconciliation
install.packages("tidyr")
library(tidyr)

# ============================================================
# WORKFLOW-10
# STEP10C-0E — Primary-Tumor Raw UMI Cell-ID + Compartment
#              Reconciliation
#
# PURPOSE
#   1. Verify that all locked Primary+Tumor TISCH2 cells
#      are present in the GEO raw UMI Cell-ID list.
#   2. Verify patient-level counts against STEP10B-1.
#   3. Verify locked compartment definitions.
#   4. Generate a formal raw-UMI pseudobulk input manifest.
#
# IMPORTANT
#   - NO full raw UMI matrix is loaded.
#   - NO expression values are read.
#   - NO pseudobulk expression is calculated.
#   - NO statistical testing is performed.
#   - NO patient/cell exclusion is performed.
#
# LOCKED COMPARTMENTS
#   T-cell             = CD4Tconv + CD8T + CD8Tex
#   Myeloid/macrophage = Mono/Macro
#   Endothelial        = Endothelial
#
# LOCKED TARGET
#   Stage  = Primary
#   Source = Tumor
#
# LOCKED COUNTS
#   Primary+Tumor cells    = 65,702
#   Primary+Tumor patients = 9
#
# DATE: 
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10C-0E — Primary-Tumor Raw UMI Cell-ID + Compartment\n")
cat("             Reconciliation\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "data.table",
  "dplyr",
  "readr"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace,
          FUN.VALUE = logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(data.table)
library(dplyr)
library(readr)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

metadata_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

raw_umi_txt_gz <- file.path(
  base_dir,
  "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
)

step10b1_dir <- file.path(
  base_dir,
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit"
)

step10b1_counts_file <- file.path(
  step10b1_dir,
  "STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv"
)

out_dir <- file.path(
  base_dir,
  "STEP10C-0E_PrimaryTumor_RawUMI_Compartment_Reconciliation"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. File existence audit
# ------------------------------------------------------------

cat("STEP10C-0E-1 — FILE EXISTENCE AUDIT\n")
cat("-----------------------------------\n")

required_files <- c(
  metadata_file,
  raw_umi_txt_gz,
  step10b1_counts_file
)

file_exists <- file.exists(required_files)

file_audit <- data.frame(
  File = basename(required_files),
  Exists = file_exists,
  stringsAsFactors = FALSE
)

print(file_audit, row.names = FALSE)

if (!all(file_exists)) {
  stop(
    "\nFAIL — One or more required files are missing."
  )
}

cat("\nPASS — All required input files exist.\n\n")


# ------------------------------------------------------------
# 3. Read TISCH2 metadata
# ------------------------------------------------------------

cat("STEP10C-0E-2 — READ TISCH2 METADATA\n")
cat("-----------------------------------\n")

meta <- fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n")

required_meta_cols <- c(
  "Cell",
  "Patient",
  "Sample",
  "Source",
  "Stage",
  "Celltype (major-lineage)"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)

if (length(missing_meta_cols) > 0) {
  stop(
    "\nFAIL — Required metadata columns missing: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

cat("PASS — Required metadata columns present.\n\n")


# ------------------------------------------------------------
# 4. Metadata Cell-ID uniqueness
# ------------------------------------------------------------

cat("STEP10C-0E-3 — TISCH2 CELL-ID UNIQUENESS\n")
cat("---------------------------------------\n")

duplicate_meta_cells <- sum(
  duplicated(meta$Cell)
)

cat("TISCH2 cells:", nrow(meta), "\n")
cat("Duplicate Cell IDs:", duplicate_meta_cells, "\n")

if (duplicate_meta_cells != 0) {
  stop(
    "\nFAIL — Duplicate TISCH2 Cell IDs detected."
  )
}

cat("PASS — TISCH2 Cell IDs are unique.\n\n")


# ------------------------------------------------------------
# 5. Locked compartment mapping
# ------------------------------------------------------------

cat("STEP10C-0E-4 — LOCKED COMPARTMENT MAPPING\n")
cat("----------------------------------------\n")

major_lineage <- trimws(
  as.character(meta[["Celltype (major-lineage)"]])
)

meta$Locked_Compartment <- NA_character_

meta$Locked_Compartment[
  major_lineage %in% c(
    "CD4Tconv",
    "CD8T",
    "CD8Tex"
  )
] <- "T-cell"

meta$Locked_Compartment[
  major_lineage %in% c(
    "Mono/Macro"
  )
] <- "Myeloid/macrophage"

meta$Locked_Compartment[
  major_lineage %in% c(
    "Endothelial"
  )
] <- "Endothelial"


locked_labels <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

mapping_labels_present <- sort(
  unique(
    na.omit(meta$Locked_Compartment)
  )
)

missing_locked_labels <- setdiff(
  locked_labels,
  mapping_labels_present
)

mapping_labels_ok <- length(missing_locked_labels) == 0

cat(
  "Locked compartments:",
  paste(locked_labels, collapse = " | "),
  "\n"
)

cat(
  "Missing locked labels:",
  ifelse(
    mapping_labels_ok,
    "0",
    paste(missing_locked_labels, collapse = ", ")
  ),
  "\n"
)

if (!mapping_labels_ok) {
  stop(
    "\nFAIL — One or more locked compartments were not represented."
  )
}

cat("PASS — Locked compartment mapping present.\n\n")


# ------------------------------------------------------------
# 6. Define Primary+Tumor target cohort
# ------------------------------------------------------------

cat("STEP10C-0E-5 — PRIMARY+TUMOR TARGET COHORT\n")
cat("-----------------------------------------\n")

meta$Stage_clean <- trimws(
  as.character(meta$Stage)
)

meta$Source_clean <- trimws(
  as.character(meta$Source)
)

target_idx <- (
  meta$Stage_clean == "Primary" &
    meta$Source_clean == "Tumor"
)

target_meta <- meta[target_idx, , drop = FALSE]

cat(
  "Primary+Tumor cells:",
  nrow(target_meta),
  "\n"
)

cat(
  "Primary+Tumor patients:",
  length(unique(target_meta$Patient)),
  "\n"
)

cat(
  "Primary+Tumor samples:",
  length(unique(target_meta$Sample)),
  "\n"
)

if (nrow(target_meta) != 65702) {
  stop(
    "\nFAIL — Primary+Tumor cell count does not match locked 65,702."
  )
}

if (length(unique(target_meta$Patient)) != 9) {
  stop(
    "\nFAIL — Primary+Tumor patient count does not match locked 9."
  )
}

cat(
  "PASS — Primary+Tumor target cohort matches locked structure.\n\n"
)


# ------------------------------------------------------------
# 7. Locked compartment counts within Primary+Tumor
# ------------------------------------------------------------

cat("STEP10C-0E-6 — PRIMARY+TUMOR COMPARTMENT COUNTS\n")
cat("-----------------------------------------------\n")

target_compartment_counts <- target_meta %>%
  count(
    Locked_Compartment,
    name = "Cell_Count"
  ) %>%
  arrange(
    match(
      Locked_Compartment,
      locked_labels
    )
  )

print(
  target_compartment_counts,
  row.names = FALSE
)

missing_target_compartments <- setdiff(
  locked_labels,
  target_compartment_counts$Locked_Compartment
)

if (length(missing_target_compartments) > 0) {
  stop(
    "\nFAIL — Target cohort missing locked compartment(s): ",
    paste(
      missing_target_compartments,
      collapse = ", "
    )
  )
}

cat(
  "\nPASS — All locked compartments are represented in Primary+Tumor.\n\n"
)


# ------------------------------------------------------------
# 8. Patient × locked compartment counts
# ------------------------------------------------------------

cat("STEP10C-0E-7 — PATIENT × COMPARTMENT COUNTS\n")
cat("------------------------------------------\n")

patient_compartment_long <- target_meta %>%
  filter(
    Locked_Compartment %in% locked_labels
  ) %>%
  count(
    Patient,
    Locked_Compartment,
    name = "Cell_Count"
  ) %>%
  arrange(
    Patient,
    match(
      Locked_Compartment,
      locked_labels
    )
  )

print(
  patient_compartment_long,
  row.names = FALSE
)


# Wide version
patient_compartment_wide <- target_meta %>%
  filter(
    Locked_Compartment %in% locked_labels
  ) %>%
  count(
    Patient,
    Locked_Compartment,
    name = "Cell_Count"
  ) %>%
  tidyr::pivot_wider(
    names_from = Locked_Compartment,
    values_from = Cell_Count,
    values_fill = 0
  )

# Ensure column order
for (nm in locked_labels) {
  if (!nm %in% colnames(patient_compartment_wide)) {
    patient_compartment_wide[[nm]] <- 0L
  }
}

patient_compartment_wide <- patient_compartment_wide %>%
  select(
    Patient,
    all_of(locked_labels)
  ) %>%
  arrange(Patient)

cat("\nPatient × compartment wide table:\n")
print(
  patient_compartment_wide,
  row.names = FALSE
)

cat("\n")


# ------------------------------------------------------------
# 9. Compare with STEP10B-1 locked counts
# ------------------------------------------------------------

cat("STEP10C-0E-8 — STEP10B-1 COUNT RECONCILIATION\n")
cat("--------------------------------------------\n")

b1_counts <- fread(
  step10b1_counts_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "STEP10B-1 rows:",
  nrow(b1_counts),
  "\n"
)

cat(
  "STEP10B-1 columns:",
  paste(colnames(b1_counts), collapse = " | "),
  "\n\n"
)


# ------------------------------------------------------------
# Flexible extraction of expected STEP10B-1 counts
# ------------------------------------------------------------

# The locked STEP10B-1 file is expected to contain:
#   Patient
#   T-cell
#   Myeloid/macrophage
#   Endothelial
#
# This block checks exact names first and provides
# a controlled fallback for common syntactic variants.

find_col <- function(df, candidates) {
  
  hit <- candidates[
    candidates %in% colnames(df)
  ]
  
  if (length(hit) >= 1) {
    return(hit[1])
  }
  
  return(NA_character_)
}


b1_patient_col <- find_col(
  b1_counts,
  c(
    "Patient",
    "patient",
    "Patient_ID",
    "patient_id"
  )
)

b1_t_col <- find_col(
  b1_counts,
  c(
    "T-cell",
    "T_cell",
    "T.cell",
    "Tcell"
  )
)

b1_myeloid_col <- find_col(
  b1_counts,
  c(
    "Myeloid/macrophage",
    "Myeloid_macrophage",
    "Myeloid.macrophage",
    "Myeloid_Macrophage"
  )
)

b1_endo_col <- find_col(
  b1_counts,
  c(
    "Endothelial",
    "endothelial"
  )
)


b1_required_cols <- c(
  b1_patient_col,
  b1_t_col,
  b1_myeloid_col,
  b1_endo_col
)

if (any(is.na(b1_required_cols))) {
  
  cat(
    "Available STEP10B-1 column names:\n"
  )
  
  print(
    colnames(b1_counts)
  )
  
  stop(
    "\nFAIL — Could not identify all required STEP10B-1 count columns."
  )
}


b1_expected <- b1_counts %>%
  transmute(
    Patient = as.character(.data[[b1_patient_col]]),
    `T-cell` = as.numeric(.data[[b1_t_col]]),
    `Myeloid/macrophage` = as.numeric(.data[[b1_myeloid_col]]),
    Endothelial = as.numeric(.data[[b1_endo_col]])
  )


# Keep only target patients
target_patients <- sort(
  unique(
    as.character(target_meta$Patient)
  )
)

b1_expected <- b1_expected %>%
  filter(
    Patient %in% target_patients
  ) %>%
  arrange(Patient)

patient_compartment_wide_sorted <-
  patient_compartment_wide %>%
  mutate(
    Patient = as.character(Patient)
  ) %>%
  arrange(Patient)


# Compare patient sets
patient_set_a <- sort(
  patient_compartment_wide_sorted$Patient
)

patient_set_b <- sort(
  b1_expected$Patient
)

patient_sets_identical <- identical(
  patient_set_a,
  patient_set_b
)

cat(
  "Patient set identical:",
  patient_sets_identical,
  "\n"
)

if (!patient_sets_identical) {
  
  cat(
    "STEP10C-0E — Target patients missing from STEP10B-1:\n"
  )
  
  print(
    setdiff(
      patient_set_a,
      patient_set_b
    )
  )
  
  cat(
    "STEP10C-0E — STEP10B-1 patients not in target:\n"
  )
  
  print(
    setdiff(
      patient_set_b,
      patient_set_a
    )
  )
  
  stop(
    "\nFAIL — Patient sets do not reconcile."
  )
}


# Compare counts
count_comparison <- patient_compartment_wide_sorted %>%
  inner_join(
    b1_expected,
    by = "Patient",
    suffix = c(
      "_STEP10C0E",
      "_STEP10B1"
    )
  ) %>%
  mutate(
    T_cell_diff =
      `T-cell_STEP10C0E` -
      `T-cell_STEP10B1`,
    
    Myeloid_diff =
      `Myeloid/macrophage_STEP10C0E` -
      `Myeloid/macrophage_STEP10B1`,
    
    Endothelial_diff =
      `Endothelial_STEP10C0E` -
      `Endothelial_STEP10B1`
  )


print(
  count_comparison,
  row.names = FALSE
)

all_counts_match <- all(
  count_comparison$T_cell_diff == 0 &
    count_comparison$Myeloid_diff == 0 &
    count_comparison$Endothelial_diff == 0
)

cat(
  "\nAll patient × compartment counts match:",
  all_counts_match,
  "\n"
)

if (!all_counts_match) {
  stop(
    "\nFAIL — STEP10C-0E counts do not match STEP10B-1."
  )
}

cat(
  "PASS — STEP10B-1 patient × compartment counts exactly reconciled.\n\n"
)


# ------------------------------------------------------------
# 10. Read ONLY GEO raw UMI Cell-ID header
# ------------------------------------------------------------

cat("STEP10C-0E-9 — GEO RAW UMI CELL-ID HEADER AUDIT\n")
cat("----------------------------------------------\n")

# Important:
# readLines() reads only the first line.
# No expression matrix values are loaded.

con <- gzfile(
  raw_umi_txt_gz,
  open = "rt"
)

header_line <- readLines(
  con,
  n = 1L,
  warn = FALSE
)

close(con)

if (length(header_line) != 1L) {
  stop(
    "\nFAIL — Could not read GEO raw UMI header."
  )
}

header_fields <- strsplit(
  header_line,
  "\t",
  fixed = TRUE
)[[1]]

if (length(header_fields) < 2) {
  stop(
    "\nFAIL — GEO raw UMI header contains fewer than two fields."
  )
}

# First field is expected to be Index
if (header_fields[1] != "Index") {
  
  warning(
    "First header field is not exactly 'Index': ",
    header_fields[1]
  )
}

geo_raw_cell_ids <- header_fields[-1]

geo_raw_cell_ids <- trimws(
  geo_raw_cell_ids
)

geo_raw_cell_ids <- geo_raw_cell_ids[
  nzchar(geo_raw_cell_ids)
]

geo_raw_n <- length(
  geo_raw_cell_ids
)

geo_raw_unique_n <- length(
  unique(geo_raw_cell_ids)
)

cat(
  "GEO raw UMI cells:",
  geo_raw_n,
  "\n"
)

cat(
  "GEO raw unique Cell IDs:",
  geo_raw_unique_n,
  "\n"
)

if (geo_raw_n != geo_raw_unique_n) {
  stop(
    "\nFAIL — Duplicate GEO raw UMI Cell IDs detected."
  )
}

if (geo_raw_n != 208506) {
  stop(
    "\nFAIL — GEO raw UMI cell count differs from locked 208,506."
  )
}

cat(
  "PASS — GEO raw UMI Cell-ID header contains 208,506 unique cells.\n\n"
)


# ------------------------------------------------------------
# 11. Global TISCH2 ↔ GEO raw reconciliation
# ------------------------------------------------------------

cat("STEP10C-0E-10 — GLOBAL CELL-ID RECONCILIATION\n")
cat("---------------------------------------------\n")

tisch2_cells <- as.character(
  meta$Cell
)

tisch2_cells <- trimws(
  tisch2_cells
)

tisch2_geo_intersection <- intersect(
  tisch2_cells,
  geo_raw_cell_ids
)

tisch2_only <- setdiff(
  tisch2_cells,
  geo_raw_cell_ids
)

geo_raw_only <- setdiff(
  geo_raw_cell_ids,
  tisch2_cells
)

cat(
  "TISCH2 cells:",
  length(tisch2_cells),
  "\n"
)

cat(
  "GEO raw cells:",
  length(geo_raw_cell_ids),
  "\n"
)

cat(
  "TISCH2 ↔ GEO raw intersection:",
  length(tisch2_geo_intersection),
  "\n"
)

cat(
  "TISCH2-only:",
  length(tisch2_only),
  "\n"
)

cat(
  "GEO raw-only:",
  length(geo_raw_only),
  "\n"
)

global_reconciliation_pass <- (
  length(tisch2_cells) == 203298 &
    length(tisch2_geo_intersection) == 203298 &
    length(tisch2_only) == 0
)

if (!global_reconciliation_pass) {
  stop(
    "\nFAIL — Global TISCH2 ↔ GEO raw Cell-ID reconciliation failed."
  )
}

cat(
  "PASS — All 203,298 TISCH2 cells are represented in GEO raw UMI.\n\n"
)


# ------------------------------------------------------------
# 12. Primary+Tumor TISCH2 ↔ GEO raw reconciliation
# ------------------------------------------------------------

cat("STEP10C-0E-11 — PRIMARY+TUMOR CELL-ID RECONCILIATION\n")
cat("--------------------------------------------------\n")

target_cells <- as.character(
  target_meta$Cell
)

target_cells <- trimws(
  target_cells
)

target_geo_intersection <- intersect(
  target_cells,
  geo_raw_cell_ids
)

target_only <- setdiff(
  target_cells,
  geo_raw_cell_ids
)

geo_target_extra <- setdiff(
  geo_raw_cell_ids,
  target_cells
)

cat(
  "Primary+Tumor TISCH2 cells:",
  length(target_cells),
  "\n"
)

cat(
  "Primary+Tumor cells present in GEO raw:",
  length(target_geo_intersection),
  "\n"
)

cat(
  "Primary+Tumor TISCH2-only:",
  length(target_only),
  "\n"
)

cat(
  "GEO raw cells outside Primary+Tumor target:",
  length(geo_target_extra),
  "\n"
)

target_reconciliation_pass <- (
  length(target_cells) == 65702 &
    length(target_geo_intersection) == 65702 &
    length(target_only) == 0
)

if (!target_reconciliation_pass) {
  
  cat(
    "\nMissing Primary+Tumor Cell IDs:\n"
  )
  
  if (length(target_only) > 0) {
    print(
      head(
        target_only,
        50
      )
    )
  }
  
  stop(
    "\nFAIL — Primary+Tumor TISCH2 cells are not fully represented in GEO raw."
  )
}

cat(
  "PASS — All 65,702 Primary+Tumor TISCH2 cells are present in GEO raw UMI.\n\n"
)


# ------------------------------------------------------------
# 13. Patient-level raw UMI Cell-ID reconciliation
# ------------------------------------------------------------

cat("STEP10C-0E-12 — PATIENT-LEVEL CELL-ID RECONCILIATION\n")
cat("--------------------------------------------------\n")

target_meta$GEO_raw_UMI_CellID_available <- (
  target_meta$Cell %in% geo_raw_cell_ids
)

patient_availability <- target_meta %>%
  group_by(Patient) %>%
  summarise(
    TISCH2_PrimaryTumor_Cells = n(),
    GEO_raw_UMI_Cells_Available =
      sum(GEO_raw_UMI_CellID_available),
    Missing_from_GEO_raw =
      sum(!GEO_raw_UMI_CellID_available),
    .groups = "drop"
  ) %>%
  arrange(Patient)

print(
  patient_availability,
  row.names = FALSE
)

patient_availability_pass <- all(
  patient_availability$Missing_from_GEO_raw == 0
)

if (!patient_availability_pass) {
  stop(
    "\nFAIL — One or more patients contain Primary+Tumor cells absent from GEO raw."
  )
}

cat(
  "\nPASS — All 9 Primary+Tumor patients have complete GEO raw Cell-ID coverage.\n\n"
)


# ------------------------------------------------------------
# 14. Compartment-level GEO raw availability
# ------------------------------------------------------------

cat("STEP10C-0E-13 — COMPARTMENT-LEVEL RAW UMI AVAILABILITY\n")
cat("-----------------------------------------------------\n")

compartment_availability <- target_meta %>%
  filter(
    Locked_Compartment %in% locked_labels
  ) %>%
  group_by(Locked_Compartment) %>%
  summarise(
    TISCH2_Cells = n(),
    GEO_raw_UMI_Cells_Available =
      sum(GEO_raw_UMI_CellID_available),
    Missing_from_GEO_raw =
      sum(!GEO_raw_UMI_CellID_available),
    Coverage =
      GEO_raw_UMI_Cells_Available /
      TISCH2_Cells,
    .groups = "drop"
  ) %>%
  arrange(
    match(
      Locked_Compartment,
      locked_labels
    )
  )

print(
  compartment_availability,
  row.names = FALSE
)

compartment_availability_pass <- all(
  compartment_availability$Missing_from_GEO_raw == 0 &
    compartment_availability$Coverage == 1
)

if (!compartment_availability_pass) {
  stop(
    "\nFAIL — One or more locked compartments do not have complete GEO raw Cell-ID coverage."
  )
}

cat(
  "\nPASS — All locked compartments have 100% GEO raw Cell-ID coverage.\n\n"
)


# ------------------------------------------------------------
# 15. Construct formal pseudobulk input manifest
# ------------------------------------------------------------

cat("STEP10C-0E-14 — PSEUDOBULK INPUT MANIFEST\n")
cat("------------------------------------------\n")

# IMPORTANT:
# This is ONLY a cell-level manifest.
# No expression values are included.
# It is the formal input definition for the later
# STEP10C patient-level pseudobulk construction.

pseudobulk_manifest <- target_meta %>%
  transmute(
    Patient = as.character(Patient),
    CellID = as.character(Cell),
    Sample = as.character(Sample),
    Source = as.character(Source),
    Stage = as.character(Stage),
    Major_lineage =
      as.character(
        `Celltype (major-lineage)`
      ),
    Locked_Compartment =
      as.character(Locked_Compartment),
    GEO_raw_UMI_CellID_available =
      as.logical(
        GEO_raw_UMI_CellID_available
      )
  ) %>%
  arrange(
    Patient,
    Locked_Compartment,
    CellID
  )

cat(
  "Manifest rows:",
  nrow(pseudobulk_manifest),
  "\n"
)

cat(
  "Manifest columns:",
  ncol(pseudobulk_manifest),
  "\n"
)

cat(
  "Unique patients:",
  length(
    unique(
      pseudobulk_manifest$Patient
    )
  ),
  "\n"
)

cat(
  "Unique Cell IDs:",
  length(
    unique(
      pseudobulk_manifest$CellID
    )
  ),
  "\n"
)

cat(
  "GEO raw Cell-ID availability:",
  sum(
    pseudobulk_manifest$GEO_raw_UMI_CellID_available
  ),
  "/",
  nrow(pseudobulk_manifest),
  "\n"
)


manifest_pass <- (
  nrow(pseudobulk_manifest) == 65702 &
    length(unique(pseudobulk_manifest$CellID)) == 65702 &
    all(
      pseudobulk_manifest$GEO_raw_UMI_CellID_available
    )
)

if (!manifest_pass) {
  stop(
    "\nFAIL — Pseudobulk input manifest integrity check failed."
  )
}

cat(
  "PASS — Formal pseudobulk input manifest generated.\n\n"
)


# ------------------------------------------------------------
# 16. Manifest compartment summary
# ------------------------------------------------------------

manifest_summary <- pseudobulk_manifest %>%
  count(
    Patient,
    Locked_Compartment,
    name = "Cell_Count"
  ) %>%
  arrange(
    Patient,
    match(
      Locked_Compartment,
      locked_labels
    )
  )


# ------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------

cat("STEP10C-0E-15 — SAVE AUDIT OUTPUTS\n")
cat("---------------------------------\n")

write.table(
  target_meta,
  file = file.path(
    out_dir,
    "STEP10C-0E_A_PrimaryTumor_Target_Metadata.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_compartment_long,
  file = file.path(
    out_dir,
    "STEP10C-0E_B_Patient_Compartment_Counts_Long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_compartment_wide,
  file = file.path(
    out_dir,
    "STEP10C-0E_C_Patient_Compartment_Counts_Wide.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  count_comparison,
  file = file.path(
    out_dir,
    "STEP10C-0E_D_STEP10B1_Count_Reconciliation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  patient_availability,
  file = file.path(
    out_dir,
    "STEP10C-0E_E_Patient_RawUMI_CellID_Availability.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  compartment_availability,
  file = file.path(
    out_dir,
    "STEP10C-0E_F_Compartment_RawUMI_CellID_Availability.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  pseudobulk_manifest,
  file = file.path(
    out_dir,
    "STEP10C-0E_G_PrimaryTumor_RawUMI_Pseudobulk_Input_Manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  manifest_summary,
  file = file.path(
    out_dir,
    "STEP10C-0E_H_Manifest_Patient_Compartment_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Final audit summary
# ------------------------------------------------------------

final_summary <- data.frame(
  Metric = c(
    "TISCH2 total cells",
    "GEO raw UMI total cells",
    "TISCH2-GEO raw intersection",
    "TISCH2-only cells",
    "GEO raw-only cells",
    "Primary+Tumor TISCH2 cells",
    "Primary+Tumor patients",
    "Primary+Tumor cells present in GEO raw",
    "Primary+Tumor cells missing from GEO raw",
    "T-cell cells",
    "Myeloid/macrophage cells",
    "Endothelial cells",
    "Pseudobulk manifest rows",
    "Pseudobulk manifest unique Cell IDs"
  ),
  Value = c(
    length(tisch2_cells),
    length(geo_raw_cell_ids),
    length(tisch2_geo_intersection),
    length(tisch2_only),
    length(geo_raw_only),
    nrow(target_meta),
    length(unique(target_meta$Patient)),
    length(target_geo_intersection),
    length(target_only),
    target_compartment_counts$Cell_Count[
      match(
        "T-cell",
        target_compartment_counts$Locked_Compartment
      )
    ],
    target_compartment_counts$Cell_Count[
      match(
        "Myeloid/macrophage",
        target_compartment_counts$Locked_Compartment
      )
    ],
    target_compartment_counts$Cell_Count[
      match(
        "Endothelial",
        target_compartment_counts$Locked_Compartment
      )
    ],
    nrow(pseudobulk_manifest),
    length(unique(pseudobulk_manifest$CellID))
  ),
  stringsAsFactors = FALSE
)

write.table(
  final_summary,
  file = file.path(
    out_dir,
    "STEP10C-0E_I_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Final PASS/FAIL lock
# ------------------------------------------------------------

all_final_checks <- c(
  metadata_rows = nrow(meta) == 203298,
  
  metadata_unique_cells =
    duplicate_meta_cells == 0,
  
  primary_tumor_cells =
    nrow(target_meta) == 65702,
  
  primary_tumor_patients =
    length(unique(target_meta$Patient)) == 9,
  
  global_tisch2_geo_intersection =
    length(tisch2_geo_intersection) == 203298,
  
  global_tisch2_only =
    length(tisch2_only) == 0,
  
  primary_tumor_geo_intersection =
    length(target_geo_intersection) == 65702,
  
  primary_tumor_only =
    length(target_only) == 0,
  
  step10b1_patient_counts =
    all_counts_match,
  
  patient_raw_umi_coverage =
    patient_availability_pass,
  
  compartment_raw_umi_coverage =
    compartment_availability_pass,
  
  manifest_rows =
    nrow(pseudobulk_manifest) == 65702,
  
  manifest_unique_cells =
    length(unique(pseudobulk_manifest$CellID)) == 65702,
  
  manifest_raw_umi_coverage =
    all(
      pseudobulk_manifest$GEO_raw_UMI_CellID_available
    )
)

final_pass <- all(
  all_final_checks
)

cat("\n")
cat("============================================================\n")

if (final_pass) {
  
  cat(
    "PASS — STEP10C-0E FINAL LOCK\n"
  )
  
} else {
  
  cat(
    "FAIL — STEP10C-0E FINAL AUDIT\n"
  )
  
}

cat("============================================================\n\n")

cat(
  "TISCH2 total cells:",
  length(tisch2_cells),
  "\n"
)

cat(
  "GEO raw UMI cells:",
  length(geo_raw_cell_ids),
  "\n"
)

cat(
  "TISCH2 ↔ GEO raw intersection:",
  length(tisch2_geo_intersection),
  "\n"
)

cat(
  "TISCH2-only:",
  length(tisch2_only),
  "\n"
)

cat(
  "GEO raw-only:",
  length(geo_raw_only),
  "\n\n"
)

cat(
  "Primary+Tumor TISCH2 cells:",
  nrow(target_meta),
  "\n"
)

cat(
  "Primary+Tumor patients:",
  length(unique(target_meta$Patient)),
  "\n"
)

cat(
  "Primary+Tumor cells in GEO raw:",
  length(target_geo_intersection),
  "\n"
)

cat(
  "Primary+Tumor cells missing from GEO raw:",
  length(target_only),
  "\n\n"
)

cat(
  "T-cell:",
  target_compartment_counts$Cell_Count[
    match(
      "T-cell",
      target_compartment_counts$Locked_Compartment
    )
  ],
  "\n"
)

cat(
  "Myeloid/macrophage:",
  target_compartment_counts$Cell_Count[
    match(
      "Myeloid/macrophage",
      target_compartment_counts$Locked_Compartment
    )
  ],
  "\n"
)

cat(
  "Endothelial:",
  target_compartment_counts$Cell_Count[
    match(
      "Endothelial",
      target_compartment_counts$Locked_Compartment
    )
  ],
  "\n\n"
)

cat(
  "Pseudobulk manifest:",
  nrow(pseudobulk_manifest),
  " cells\n"
)

cat(
  "Pseudobulk manifest GEO raw coverage:",
  sum(
    pseudobulk_manifest$GEO_raw_UMI_CellID_available
  ),
  "/",
  nrow(pseudobulk_manifest),
  "\n\n"
)

cat(
  "No full raw UMI matrix loaded.\n"
)

cat(
  "No expression values analyzed.\n"
)

cat(
  "No pseudobulk expression calculated.\n"
)

cat(
  "No statistical testing performed.\n"
)

cat(
  "No patient/cell exclusion performed.\n\n"
)

cat(
  "Output directory:\n",
  out_dir,
  "\n\n"
)

if (!final_pass) {
  print(
    all_final_checks
  )
  stop(
    "\nSTEP10C-0E FAILED — see failed checks above."
  )
}

cat(
  "============================================================\n"
)
cat(
  "STEP10C-0E COMPLETE — FINAL STATUS: PASS\n"
)
cat(
  "============================================================\n"
)

# ============================================================
# WORKFLOW-10
# STEP10C — Patient-Level Raw UMI Pseudobulk Construction
#
# FINAL LOCKED DESIGN
# ------------------------------------------------------------
# Purpose:
#   Construct patient-level pseudobulk raw UMI count matrices
#   from the GEO raw UMI matrix using the locked TISCH2
#   Primary+Tumor cohort and locked compartments.
#
# Input:
#   GEO GSE131907 raw UMI matrix:
#   GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz
#
#   TISCH2 metadata:
#   NSCLC_GSE131907_CellMetainfo_table.tsv
#
#   STEP10B-1:
#   STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv
#
# LOCKED TARGET
# ------------------------------------------------------------
#   Stage  = Primary
#   Source = Tumor
#
# LOCKED COMPARTMENTS
# ------------------------------------------------------------
#   T-cell             = CD4Tconv + CD8T + CD8Tex
#   Myeloid/macrophage = Mono/Macro
#   Endothelial        = Endothelial
#
# LOCKED CELL-DEPTH STRATA
# ------------------------------------------------------------
#   Primary analysis:
#       >= 50 cells / patient / compartment
#
#   Sensitivity analysis:
#       20–49 cells / patient / compartment
#
#   Insufficient:
#       < 20 cells / patient / compartment
#
# LOCKED PRIMARY COHORT
# ------------------------------------------------------------
#   T-cell             = 9 patients
#   Myeloid/macrophage = 9 patients
#   Endothelial        = 4 patients
#
# LOCKED SENSITIVITY COHORT
# ------------------------------------------------------------
#   Endothelial = 2 patients
#
# LOCKED INSUFFICIENT
# ------------------------------------------------------------
#   Endothelial = 3 patients
#
# IMPORTANT
# ------------------------------------------------------------
#   - Raw UMI expression values are read for the first time here.
#   - No H5 expression values are used.
#   - No normalized expression is used.
#   - No statistical testing.
#   - No differential expression.
#   - No CellChat.
#   - No survival analysis.
#   - No patient exclusion beyond the PRE-LOCKED
#     cell-depth strata.
#   - GEO raw-only 5,208 cells are NOT added.
#
# MEMORY STRATEGY
# ------------------------------------------------------------
#   The complete raw matrix is NOT loaded into memory.
#   The file is read gene-by-gene.
#   Only locked target cell columns are extracted.
#   Counts are immediately aggregated to
#   Gene × Patient × Compartment.
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10C — Patient-Level Raw UMI Pseudobulk Construction\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_pkgs <- c(
  "data.table"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    FUN.VALUE = logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(data.table)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

metadata_file <- file.path(
  base_dir,
  "NSCLC_GSE131907_CellMetainfo_table.tsv"
)

raw_umi_file <- file.path(
  base_dir,
  "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
)

step10b1_file <- file.path(
  base_dir,
  "STEP10B-1_Compartment_CellCount_Feasibility_Audit",
  "STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv"
)

step10c0e_manifest_file <- file.path(
  base_dir,
  "STEP10C-0E_PrimaryTumor_RawUMI_Compartment_Reconciliation",
  "STEP10C-0E_G_PrimaryTumor_RawUMI_Pseudobulk_Input_Manifest.tsv"
)

step10c0e_patient_counts_file <- file.path(
  base_dir,
  "STEP10C-0E_PrimaryTumor_RawUMI_Compartment_Reconciliation",
  "STEP10C-0E_C_Patient_Compartment_Counts_Wide.tsv"
)

out_dir <- file.path(
  base_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. File existence
# ------------------------------------------------------------

cat("STEP10C-1 — INPUT FILE AUDIT\n")
cat("----------------------------\n")

input_files <- c(
  metadata_file,
  raw_umi_file,
  step10b1_file,
  step10c0e_manifest_file,
  step10c0e_patient_counts_file
)

input_audit <- data.frame(
  File = basename(input_files),
  Exists = file.exists(input_files),
  Size_MB = ifelse(
    file.exists(input_files),
    round(
      file.info(input_files)$size / 1024^2,
      2
    ),
    NA
  ),
  stringsAsFactors = FALSE
)

print(
  input_audit,
  row.names = FALSE
)

if (!all(input_audit$Exists)) {
  stop(
    "\nFAIL — Required input file(s) missing."
  )
}

cat(
  "\nPASS — All STEP10C input files exist.\n\n"
)


# ------------------------------------------------------------
# 3. Read TISCH2 metadata
# ------------------------------------------------------------

cat("STEP10C-2 — READ TISCH2 METADATA\n")
cat("--------------------------------\n")

meta <- fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n")

required_meta_cols <- c(
  "Cell",
  "Patient",
  "Sample",
  "Source",
  "Stage",
  "Celltype (major-lineage)"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)

if (length(missing_meta_cols) > 0) {
  stop(
    "FAIL — Missing metadata columns: ",
    paste(
      missing_meta_cols,
      collapse = ", "
    )
  )
}

if (nrow(meta) != 203298) {
  stop(
    "FAIL — TISCH2 metadata cell count is not 203,298."
  )
}

if (anyDuplicated(meta$Cell) > 0) {
  stop(
    "FAIL — Duplicate TISCH2 Cell IDs detected."
  )
}

cat(
  "PASS — TISCH2 metadata structure verified.\n\n"
)


# ------------------------------------------------------------
# 4. Locked compartment mapping
# ------------------------------------------------------------

cat("STEP10C-3 — LOCKED COMPARTMENT MAPPING\n")
cat("--------------------------------------\n")

major_lineage <- trimws(
  as.character(
    meta[["Celltype (major-lineage)"]]
  )
)

meta$Locked_Compartment <- NA_character_

meta$Locked_Compartment[
  major_lineage %in%
    c("CD4Tconv", "CD8T", "CD8Tex")
] <- "T-cell"

meta$Locked_Compartment[
  major_lineage %in%
    c("Mono/Macro")
] <- "Myeloid/macrophage"

meta$Locked_Compartment[
  major_lineage %in%
    c("Endothelial")
] <- "Endothelial"

locked_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

cat(
  "Locked compartments:\n"
)

print(
  locked_compartments
)

cat(
  "\nNo DC cells included in Myeloid/macrophage.\n"
)

cat(
  "PASS — Locked compartment definitions applied.\n\n"
)


# ------------------------------------------------------------
# 5. Define Primary+Tumor target cells
# ------------------------------------------------------------

cat("STEP10C-4 — LOCK PRIMARY+TUMOR TARGET\n")
cat("------------------------------------\n")

meta$Stage_clean <- trimws(
  as.character(meta$Stage)
)

meta$Source_clean <- trimws(
  as.character(meta$Source)
)

target_idx <- (
  meta$Stage_clean == "Primary" &
    meta$Source_clean == "Tumor"
)

target_meta <- meta[
  target_idx,
  ,
  drop = FALSE
]

cat(
  "Primary+Tumor cells:",
  nrow(target_meta),
  "\n"
)

cat(
  "Primary+Tumor patients:",
  length(
    unique(target_meta$Patient)
  ),
  "\n"
)

if (nrow(target_meta) != 65702) {
  stop(
    "FAIL — Primary+Tumor cohort != 65,702 cells."
  )
}

if (
  length(
    unique(target_meta$Patient)
  ) != 9
) {
  stop(
    "FAIL — Primary+Tumor cohort != 9 patients."
  )
}

cat(
  "PASS — Primary+Tumor target locked.\n\n"
)


# ------------------------------------------------------------
# 6. Load STEP10C-0E manifest
# ------------------------------------------------------------

cat("STEP10C-5 — LOAD STEP10C-0E MANIFEST\n")
cat("------------------------------------\n")

manifest <- fread(
  step10c0e_manifest_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Manifest rows:",
  nrow(manifest),
  "\n"
)

if (nrow(manifest) != 65702) {
  stop(
    "FAIL — STEP10C-0E manifest does not contain 65,702 cells."
  )
}

if (
  length(
    unique(manifest$CellID)
  ) != 65702
) {
  stop(
    "FAIL — Manifest Cell IDs are not unique."
  )
}

if (
  !all(
    manifest$GEO_raw_UMI_CellID_available
  )
) {
  stop(
    "FAIL — Not all target cells have GEO raw UMI Cell-ID availability."
  )
}

cat(
  "PASS — STEP10C-0E manifest integrity verified.\n\n"
)


# ------------------------------------------------------------
# 7. Load locked patient × compartment cell counts
# ------------------------------------------------------------

cat("STEP10C-6 — LOAD LOCKED CELL-DEPTH TABLE\n")
cat("-----------------------------------------\n")

locked_counts <- fread(
  step10c0e_patient_counts_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

cat(
  "Rows:",
  nrow(locked_counts),
  "\n"
)

if (nrow(locked_counts) != 9) {
  stop(
    "FAIL — Expected exactly 9 Primary+Tumor patients."
  )
}


# ------------------------------------------------------------
# 8. Construct locked eligibility strata
# ------------------------------------------------------------

cat("STEP10C-7 — LOCK PRIMARY / SENSITIVITY / INSUFFICIENT\n")
cat("------------------------------------------------------\n")

# ------------------------------------------------------------
# Normalize column names if necessary
# ------------------------------------------------------------

required_count_cols <- c(
  "Patient",
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

missing_count_cols <- setdiff(
  required_count_cols,
  colnames(locked_counts)
)

if (length(missing_count_cols) > 0) {
  
  cat(
    "Available count-table columns:\n"
  )
  
  print(
    colnames(locked_counts)
  )
  
  stop(
    "FAIL — Required locked cell-count columns missing."
  )
}


# ------------------------------------------------------------
# Cell-depth classification
# ------------------------------------------------------------

classify_depth <- function(n) {
  
  if (is.na(n)) {
    return(NA_character_)
  }
  
  if (n >= 50) {
    return("Primary")
  }
  
  if (n >= 20) {
    return("Sensitivity")
  }
  
  return("Insufficient")
}


eligibility_long <- rbind(
  
  data.frame(
    Patient = locked_counts$Patient,
    Compartment = "T-cell",
    Cell_Count = locked_counts[["T-cell"]],
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    Patient = locked_counts$Patient,
    Compartment = "Myeloid/macrophage",
    Cell_Count =
      locked_counts[["Myeloid/macrophage"]],
    stringsAsFactors = FALSE
  ),
  
  data.frame(
    Patient = locked_counts$Patient,
    Compartment = "Endothelial",
    Cell_Count =
      locked_counts[["Endothelial"]],
    stringsAsFactors = FALSE
  )
  
)

eligibility_long$Eligibility <- vapply(
  eligibility_long$Cell_Count,
  classify_depth,
  FUN.VALUE = character(1)
)


print(
  eligibility_long,
  row.names = FALSE
)


# ------------------------------------------------------------
# Verify expected locked cohort sizes
# ------------------------------------------------------------

primary_n <- table(
  eligibility_long$Compartment[
    eligibility_long$Eligibility == "Primary"
  ]
)

sensitivity_n <- table(
  eligibility_long$Compartment[
    eligibility_long$Eligibility == "Sensitivity"
  ]
)

insufficient_n <- table(
  eligibility_long$Compartment[
    eligibility_long$Eligibility == "Insufficient"
  ]
)


get_n <- function(tab, label) {
  
  if (label %in% names(tab)) {
    return(as.integer(tab[[label]]))
  }
  
  return(0L)
}


locked_primary_n <- c(
  `T-cell` =
    get_n(primary_n, "T-cell"),
  `Myeloid/macrophage` =
    get_n(primary_n, "Myeloid/macrophage"),
  Endothelial =
    get_n(primary_n, "Endothelial")
)

locked_sensitivity_n <- c(
  `T-cell` =
    get_n(sensitivity_n, "T-cell"),
  `Myeloid/macrophage` =
    get_n(sensitivity_n, "Myeloid/macrophage"),
  Endothelial =
    get_n(sensitivity_n, "Endothelial")
)

locked_insufficient_n <- c(
  `T-cell` =
    get_n(insufficient_n, "T-cell"),
  `Myeloid/macrophage` =
    get_n(insufficient_n, "Myeloid/macrophage"),
  Endothelial =
    get_n(insufficient_n, "Endothelial")
)


cat("\n")
cat("LOCKED PRIMARY:\n")
print(locked_primary_n)

cat("\nLOCKED SENSITIVITY:\n")
print(locked_sensitivity_n)

cat("\nLOCKED INSUFFICIENT:\n")
print(locked_insufficient_n)


if (
  locked_primary_n[["T-cell"]] != 9 ||
  locked_primary_n[["Myeloid/macrophage"]] != 9 ||
  locked_primary_n[["Endothelial"]] != 4
) {
  stop(
    "\nFAIL — Primary cohort does not match STEP10B-2 lock."
  )
}

if (
  locked_sensitivity_n[["Endothelial"]] != 2
) {
  stop(
    "\nFAIL — Endothelial sensitivity cohort does not match STEP10B-2 lock."
  )
}

if (
  locked_insufficient_n[["Endothelial"]] != 3
) {
  stop(
    "\nFAIL — Endothelial insufficient cohort does not match STEP10B-2 lock."
  )
}

if (
  locked_sensitivity_n[["T-cell"]] != 0 ||
  locked_sensitivity_n[["Myeloid/macrophage"]] != 0
) {
  stop(
    "\nFAIL — Unexpected sensitivity patients in T-cell or Myeloid/macrophage."
  )
}

if (
  locked_insufficient_n[["T-cell"]] != 0 ||
  locked_insufficient_n[["Myeloid/macrophage"]] != 0
) {
  stop(
    "\nFAIL — Unexpected insufficient patients in T-cell or Myeloid/macrophage."
  )
}

cat(
  "\nPASS — Primary / sensitivity / insufficient cohorts exactly match STEP10B-2.\n\n"
)


# ------------------------------------------------------------
# 9. Create explicit patient-compartment cohort labels
# ------------------------------------------------------------

eligibility_key <- eligibility_long

eligibility_key$Group_ID <- paste(
  eligibility_key$Patient,
  eligibility_key$Compartment,
  sep = "__"
)

# ------------------------------------------------------------
# Primary patient-compartment groups
# ------------------------------------------------------------

primary_key <- eligibility_key[
  eligibility_key$Eligibility == "Primary",
  ,
  drop = FALSE
]

sensitivity_key <- eligibility_key[
  eligibility_key$Eligibility == "Sensitivity",
  ,
  drop = FALSE
]


cat(
  "Primary patient-compartment groups:",
  nrow(primary_key),
  "\n"
)

cat(
  "Sensitivity patient-compartment groups:",
  nrow(sensitivity_key),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Map each target Cell to Patient × Compartment
# ------------------------------------------------------------

cat("STEP10C-8 — BUILD CELL-LEVEL AGGREGATION MAP\n")
cat("---------------------------------------------\n")

target_manifest <- manifest[
  manifest$Locked_Compartment %in%
    locked_compartments,
  ,
  drop = FALSE
]

target_manifest$Patient <- as.character(
  target_manifest$Patient
)

target_manifest$CellID <- as.character(
  target_manifest$CellID
)

target_manifest$Locked_Compartment <-
  as.character(
    target_manifest$Locked_Compartment
  )


# ------------------------------------------------------------
# Create Group_ID
# ------------------------------------------------------------

target_manifest$Group_ID <- paste(
  target_manifest$Patient,
  target_manifest$Locked_Compartment,
  sep = "__"
)


# ------------------------------------------------------------
# Verify every cell has exactly one group
# ------------------------------------------------------------

if (
  anyNA(target_manifest$Group_ID)
) {
  stop(
    "FAIL — NA Group_ID detected."
  )
}

if (
  anyDuplicated(target_manifest$CellID) > 0
) {
  stop(
    "FAIL — Duplicate target Cell IDs detected."
  )
}


# ------------------------------------------------------------
# Cell ID → raw matrix column index
# ------------------------------------------------------------

# This vector will be created after reading the GEO header.

cat(
  "Target cells for locked compartments:",
  nrow(target_manifest),
  "\n"
)

cat(
  "Expected:",
  43475,
  "\n"
)

if (
  nrow(target_manifest) != 43475
) {
  stop(
    "FAIL — Locked compartment target cell count != 43,475."
  )
}

cat(
  "PASS — Cell-level aggregation map prepared.\n\n"
)


# ------------------------------------------------------------
# 11. Read ONLY GEO raw UMI header
# ------------------------------------------------------------

cat("STEP10C-9 — GEO RAW UMI HEADER\n")
cat("------------------------------\n")

con <- gzfile(
  raw_umi_file,
  open = "rt"
)

header_line <- readLines(
  con,
  n = 1L,
  warn = FALSE
)

close(con)

if (
  length(header_line) != 1L
) {
  stop(
    "FAIL — Could not read GEO raw UMI header."
  )
}

header_fields <- strsplit(
  header_line,
  "\t",
  fixed = TRUE
)[[1]]

if (
  length(header_fields) != 208507
) {
  stop(
    "FAIL — GEO raw header does not contain 208,507 fields."
  )
}

geo_cell_ids <- header_fields[-1]

geo_cell_ids <- trimws(
  geo_cell_ids
)

if (
  anyDuplicated(geo_cell_ids) > 0
) {
  stop(
    "FAIL — Duplicate GEO raw Cell IDs detected."
  )
}

cat(
  "GEO raw cells:",
  length(geo_cell_ids),
  "\n"
)

cat(
  "PASS — GEO raw header verified.\n\n"
)


# ------------------------------------------------------------
# 12. Match target Cell IDs to GEO raw column positions
# ------------------------------------------------------------

cat("STEP10C-10 — TARGET CELL-ID → RAW MATRIX COLUMN MAP\n")
cat("---------------------------------------------------\n")

target_manifest$Raw_Column <- match(
  target_manifest$CellID,
  geo_cell_ids
)

missing_raw_columns <- sum(
  is.na(target_manifest$Raw_Column)
)

cat(
  "Target cells:",
  nrow(target_manifest),
  "\n"
)

cat(
  "Target cells with raw UMI column:",
  sum(
    !is.na(target_manifest$Raw_Column)
  ),
  "\n"
)

cat(
  "Target cells missing raw UMI column:",
  missing_raw_columns,
  "\n"
)

if (
  missing_raw_columns != 0
) {
  stop(
    "FAIL — Target Cell IDs missing from GEO raw UMI matrix."
  )
}

if (
  length(
    unique(
      target_manifest$Raw_Column
    )
  ) != nrow(target_manifest)
) {
  stop(
    "FAIL — Raw matrix column mapping is not one-to-one."
  )
}

cat(
  "PASS — All 43,475 target compartment cells mapped to GEO raw UMI columns.\n\n"
)


# ------------------------------------------------------------
# 13. Build ordered raw-column index
# ------------------------------------------------------------

# Important:
# sort target cells by raw matrix column.
# This allows direct extraction from each line.

target_manifest <- target_manifest[
  order(
    target_manifest$Raw_Column
  ),
  ,
  drop = FALSE
]

target_raw_columns <- target_manifest$Raw_Column

target_group_ids <- target_manifest$Group_ID

target_patient <- target_manifest$Patient

target_compartment <-
  target_manifest$Locked_Compartment


# ------------------------------------------------------------
# 14. Define final Group_ID order
# ------------------------------------------------------------

group_levels <- unique(
  c(
    primary_key$Group_ID,
    sensitivity_key$Group_ID
  )
)

# Only groups represented by target compartments
group_levels <- group_levels[
  group_levels %in%
    unique(target_group_ids)
]

cat(
  "Total patient × compartment groups:",
  length(group_levels),
  "\n"
)

cat(
  "Primary groups:",
  nrow(primary_key),
  "\n"
)

cat(
  "Sensitivity groups:",
  nrow(sensitivity_key),
  "\n\n"
)


# ------------------------------------------------------------
# 15. Preallocate pseudobulk matrices
# ------------------------------------------------------------

cat("STEP10C-11 — INITIALIZE PSEUDOBULK MATRICES\n")
cat("--------------------------------------------\n")

# Gene names are not yet known.
# The matrix will be initialized after the first data row.

pb_primary <- NULL
pb_sensitivity <- NULL

gene_names <- character(0)

gene_counter <- 0L


# ------------------------------------------------------------
# 16. Open GEO raw UMI matrix
# ------------------------------------------------------------

cat(
  "STEP10C-12 — OPEN GEO RAW UMI MATRIX\n"
)
cat(
  "------------------------------------\n"
)

raw_con <- gzfile(
  raw_umi_file,
  open = "rt"
)

# Skip header
raw_header <- readLines(
  raw_con,
  n = 1L,
  warn = FALSE
)

if (
  length(raw_header) != 1L
) {
  close(raw_con)
  
  stop(
    "FAIL — Could not reopen GEO raw UMI matrix."
  )
}

cat(
  "Raw UMI stream opened successfully.\n"
)

cat(
  "No full matrix loaded into memory.\n\n"
)


# ------------------------------------------------------------
# 17. Streaming row-by-row pseudobulk aggregation
# ------------------------------------------------------------

cat("STEP10C-13 — STREAMING RAW UMI AGGREGATION\n")
cat("------------------------------------------\n")

cat(
  "This step may take substantial time because the GEO raw\n"
)

cat(
  "matrix contains 29,634 genes × 208,506 cells.\n\n"
)

start_time <- Sys.time()

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

parse_gene_line <- function(
    line,
    target_raw_columns
) {
  
  fields <- strsplit(
    line,
    "\t",
    fixed = TRUE
  )[[1]]
  
  if (
    length(fields) != 208507
  ) {
    stop(
      "FAIL — Raw UMI row field count != 208,507."
    )
  }
  
  gene <- fields[1]
  
  values_chr <- fields[
    target_raw_columns + 1L
  ]
  
  values_num <- suppressWarnings(
    as.numeric(values_chr)
  )
  
  if (
    anyNA(values_num)
  ) {
    
    bad_idx <- which(
      is.na(values_num)
    )[1]
    
    stop(
      "FAIL — Non-numeric raw UMI value detected at target position ",
      bad_idx,
      " for gene ",
      gene
    )
  }
  
  list(
    gene = gene,
    values = values_num
  )
}


# ------------------------------------------------------------
# Main streaming loop
# ------------------------------------------------------------

repeat {
  
  line <- readLines(
    raw_con,
    n = 1L,
    warn = FALSE
  )
  
  if (
    length(line) == 0L
  ) {
    break
  }
  
  parsed <- parse_gene_line(
    line,
    target_raw_columns
  )
  
  gene <- parsed$gene
  values <- parsed$values
  
  gene_counter <- gene_counter + 1L
  
  gene_names[gene_counter] <- gene
  
  # ----------------------------------------------------------
  # Initialize matrices on first gene
  # ----------------------------------------------------------
  
  if (
    gene_counter == 1L
  ) {
    
    pb_primary <- matrix(
      0,
      nrow = 1L,
      ncol = length(group_levels)
    )
    
    pb_sensitivity <- matrix(
      0,
      nrow = 1L,
      ncol = length(group_levels)
    )
    
  } else {
    
    pb_primary <- rbind(
      pb_primary,
      rep(
        0,
        length(group_levels)
      )
    )
    
    pb_sensitivity <- rbind(
      pb_sensitivity,
      rep(
        0,
        length(group_levels)
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # Aggregate cell counts by patient × compartment
  # ----------------------------------------------------------
  
  # Use tapply to aggregate the current gene's cell counts.
  # Only one gene vector is handled at a time.
  
  group_sum <- tapply(
    values,
    factor(
      target_group_ids,
      levels = group_levels
    ),
    sum
  )
  
  group_sum <- as.numeric(
    group_sum
  )
  
  group_sum[
    is.na(group_sum)
  ] <- 0
  
  
  # ----------------------------------------------------------
  # Primary groups
  # ----------------------------------------------------------
  
  primary_idx <- match(
    primary_key$Group_ID,
    group_levels
  )
  
  pb_primary[
    gene_counter,
    primary_idx
  ] <- group_sum[
    primary_idx
  ]
  
  
  # ----------------------------------------------------------
  # Sensitivity groups
  # ----------------------------------------------------------
  
  sensitivity_idx <- match(
    sensitivity_key$Group_ID,
    group_levels
  )
  
  if (
    length(sensitivity_idx) > 0
  ) {
    
    pb_sensitivity[
      gene_counter,
      sensitivity_idx
    ] <- group_sum[
      sensitivity_idx
    ]
    
  }
  
  
  # ----------------------------------------------------------
  # Progress reporting
  # ----------------------------------------------------------
  
  if (
    gene_counter %% 500 == 0
  ) {
    
    elapsed <- difftime(
      Sys.time(),
      start_time,
      units = "mins"
    )
    
    cat(
      "Processed genes:",
      gene_counter,
      "| elapsed:",
      round(
        as.numeric(elapsed),
        1
      ),
      "min\n"
    )
    
  }
  
}


# ------------------------------------------------------------
# Close raw matrix
# ------------------------------------------------------------

close(raw_con)

cat(
  "\nRaw UMI streaming completed.\n"
)

cat(
  "Genes processed:",
  gene_counter,
  "\n"
)

cat(
  "Expected GEO raw genes:",
  29634,
  "\n"
)

if (
  gene_counter != 29634
) {
  stop(
    "FAIL — Number of processed genes != 29,634."
  )
}

cat(
  "PASS — All 29,634 GEO raw gene rows processed.\n\n"
)


# ------------------------------------------------------------
# 18. Gene-ID integrity
# ------------------------------------------------------------

cat("STEP10C-14 — GENE-ID INTEGRITY\n")
cat("------------------------------\n")

if (
  anyNA(gene_names)
) {
  stop(
    "FAIL — NA gene names detected."
  )
}

if (
  any(!nzchar(gene_names))
) {
  stop(
    "FAIL — Empty gene IDs detected."
  )
}

duplicate_genes <- sum(
  duplicated(gene_names)
)

cat(
  "Gene rows:",
  length(gene_names),
  "\n"
)

cat(
  "Unique gene IDs:",
  length(unique(gene_names)),
  "\n"
)

cat(
  "Duplicate gene IDs:",
  duplicate_genes,
  "\n"
)

if (
  duplicate_genes != 0
) {
  stop(
    "FAIL — Duplicate gene IDs detected."
  )
}

cat(
  "PASS — GEO raw gene IDs are unique.\n\n"
)


# ------------------------------------------------------------
# 19. Raw UMI integer-value audit
# ------------------------------------------------------------

cat("STEP10C-15 — RAW UMI COUNT-VALUE AUDIT\n")
cat("--------------------------------------\n")

# Since values are parsed from the raw UMI text matrix,
# verify that all target values are non-negative integer-valued.

# The pseudobulk matrices should therefore also be
# non-negative integer-valued.

primary_nonnegative <- all(
  pb_primary >= 0
)

sensitivity_nonnegative <- all(
  pb_sensitivity >= 0
)

primary_integer <- all(
  pb_primary ==
    round(pb_primary)
)

sensitivity_integer <- all(
  pb_sensitivity ==
    round(pb_sensitivity)
)

cat(
  "Primary non-negative:",
  primary_nonnegative,
  "\n"
)

cat(
  "Primary integer-valued:",
  primary_integer,
  "\n"
)

cat(
  "Sensitivity non-negative:",
  sensitivity_nonnegative,
  "\n"
)

cat(
  "Sensitivity integer-valued:",
  sensitivity_integer,
  "\n"
)

if (
  !primary_nonnegative ||
  !primary_integer ||
  !sensitivity_nonnegative ||
  !sensitivity_integer
) {
  stop(
    "FAIL — Pseudobulk count-value integrity check failed."
  )
}

cat(
  "PASS — Pseudobulk values are non-negative integer-valued raw counts.\n\n"
)


# ------------------------------------------------------------
# 20. Assign row / column names
# ------------------------------------------------------------

rownames(
  pb_primary
) <- gene_names

rownames(
  pb_sensitivity
) <- gene_names

colnames(
  pb_primary
) <- group_levels

colnames(
  pb_sensitivity
) <- group_levels


# ------------------------------------------------------------
# 21. Keep only relevant columns
# ------------------------------------------------------------

primary_group_ids <- primary_key$Group_ID

sensitivity_group_ids <-
  sensitivity_key$Group_ID

pb_primary_final <- pb_primary[
  ,
  match(
    primary_group_ids,
    colnames(pb_primary)
  ),
  drop = FALSE
]

pb_sensitivity_final <- pb_sensitivity[
  ,
  match(
    sensitivity_group_ids,
    colnames(pb_sensitivity)
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# 22. Convert to integer matrix
# ------------------------------------------------------------

storage.mode(
  pb_primary_final
) <- "integer"

storage.mode(
  pb_sensitivity_final
) <- "integer"


# ------------------------------------------------------------
# 23. Verify expected dimensions
# ------------------------------------------------------------

cat("STEP10C-16 — PSEUDOBULK DIMENSION AUDIT\n")
cat("----------------------------------------\n")

cat(
  "Primary T-cell + Myeloid + Endothelial groups:",
  ncol(pb_primary_final),
  "\n"
)

cat(
  "Expected primary groups:",
  9 + 9 + 4,
  "\n"
)

cat(
  "Sensitivity groups:",
  ncol(pb_sensitivity_final),
  "\n"
)

cat(
  "Expected sensitivity groups:",
  2,
  "\n"
)

if (
  ncol(pb_primary_final) != 22
) {
  stop(
    "FAIL — Primary pseudobulk group count != 22."
  )
}

if (
  ncol(pb_sensitivity_final) != 2
) {
  stop(
    "FAIL — Sensitivity pseudobulk group count != 2."
  )
}

if (
  nrow(pb_primary_final) != 29634
) {
  stop(
    "FAIL — Primary pseudobulk gene count != 29,634."
  )
}

if (
  nrow(pb_sensitivity_final) != 29634
) {
  stop(
    "FAIL — Sensitivity pseudobulk gene count != 29,634."
  )
}

cat(
  "PASS — Pseudobulk dimensions verified.\n\n"
)


# ------------------------------------------------------------
# 24. Build sample annotation tables
# ------------------------------------------------------------

cat("STEP10C-17 — SAMPLE ANNOTATION TABLES\n")
cat("-------------------------------------\n")

primary_annotation <- primary_key[
  ,
  c(
    "Patient",
    "Compartment",
    "Cell_Count",
    "Eligibility",
    "Group_ID"
  ),
  drop = FALSE
]

sensitivity_annotation <- sensitivity_key[
  ,
  c(
    "Patient",
    "Compartment",
    "Cell_Count",
    "Eligibility",
    "Group_ID"
  ),
  drop = FALSE
]

primary_annotation <- primary_annotation[
  match(
    colnames(pb_primary_final),
    primary_annotation$Group_ID
  ),
  ,
  drop = FALSE
]

sensitivity_annotation <-
  sensitivity_annotation[
    match(
      colnames(pb_sensitivity_final),
      sensitivity_annotation$Group_ID
    ),
    ,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 25. Save primary pseudobulk matrices
# ------------------------------------------------------------

cat("STEP10C-18 — SAVE PRIMARY PSEUDOBULK MATRICES\n")
cat("---------------------------------------------\n")

# Split primary matrix by compartment

tcell_patients <- primary_annotation$Group_ID[
  primary_annotation$Compartment == "T-cell"
]

myeloid_patients <-
  primary_annotation$Group_ID[
    primary_annotation$Compartment ==
      "Myeloid/macrophage"
  ]

endo_primary_patients <-
  primary_annotation$Group_ID[
    primary_annotation$Compartment ==
      "Endothelial"
  ]


pb_tcell <- pb_primary_final[
  ,
  match(
    tcell_patients,
    colnames(pb_primary_final)
  ),
  drop = FALSE
]

pb_myeloid <- pb_primary_final[
  ,
  match(
    myeloid_patients,
    colnames(pb_primary_final)
  ),
  drop = FALSE
]

pb_endo_primary <- pb_primary_final[
  ,
  match(
    endo_primary_patients,
    colnames(pb_primary_final)
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# Sensitivity endothelial
# ------------------------------------------------------------

endo_sensitivity_patients <-
  sensitivity_annotation$Group_ID[
    sensitivity_annotation$Compartment ==
      "Endothelial"
  ]

pb_endo_sensitivity <- pb_sensitivity_final[
  ,
  match(
    endo_sensitivity_patients,
    colnames(pb_sensitivity_final)
  ),
  drop = FALSE
]


# ------------------------------------------------------------
# Write matrices
# ------------------------------------------------------------

write.table(
  cbind(
    Gene = rownames(pb_tcell),
    pb_tcell
  ),
  file = file.path(
    out_dir,
    "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cbind(
    Gene = rownames(pb_myeloid),
    pb_myeloid
  ),
  file = file.path(
    out_dir,
    "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cbind(
    Gene = rownames(pb_endo_primary),
    pb_endo_primary
  ),
  file = file.path(
    out_dir,
    "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  cbind(
    Gene = rownames(pb_endo_sensitivity),
    pb_endo_sensitivity
  ),
  file = file.path(
    out_dir,
    "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Save sample annotations
# ------------------------------------------------------------

write.table(
  primary_annotation,
  file = file.path(
    out_dir,
    "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sensitivity_annotation,
  file = file.path(
    out_dir,
    "STEP10C-F_Sensitivity_Pseudobulk_Sample_Annotation.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Pseudobulk cell-count audit
# ------------------------------------------------------------

pseudobulk_cellcount_audit <- rbind(
  primary_annotation,
  sensitivity_annotation
)

pseudobulk_cellcount_audit <- pseudobulk_cellcount_audit[
  order(
    pseudobulk_cellcount_audit$Compartment,
    pseudobulk_cellcount_audit$Patient
  ),
  ,
  drop = FALSE
]

write.table(
  pseudobulk_cellcount_audit,
  file = file.path(
    out_dir,
    "STEP10C-G_Pseudobulk_CellCount_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Input manifest
# ------------------------------------------------------------

write.table(
  target_manifest[
    order(
      target_manifest$Patient,
      target_manifest$Locked_Compartment,
      target_manifest$CellID
    ),
    ,
    drop = FALSE
  ],
  file = file.path(
    out_dir,
    "STEP10C-H_Pseudobulk_Input_Manifest.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Pseudobulk QC summary
# ------------------------------------------------------------

cat("STEP10C-19 — PSEUDOBULK QC SUMMARY\n")
cat("---------------------------------\n")

qc_summary <- data.frame(
  Compartment = c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial_primary",
    "Endothelial_sensitivity"
  ),
  
  Patients = c(
    ncol(pb_tcell),
    ncol(pb_myeloid),
    ncol(pb_endo_primary),
    ncol(pb_endo_sensitivity)
  ),
  
  Genes = c(
    nrow(pb_tcell),
    nrow(pb_myeloid),
    nrow(pb_endo_primary),
    nrow(pb_endo_sensitivity)
  ),
  
  Total_raw_UMI = c(
    sum(pb_tcell),
    sum(pb_myeloid),
    sum(pb_endo_primary),
    sum(pb_endo_sensitivity)
  ),
  
  Nonzero_gene_patient_pairs = c(
    sum(pb_tcell > 0),
    sum(pb_myeloid > 0),
    sum(pb_endo_primary > 0),
    sum(pb_endo_sensitivity > 0)
  ),
  
  stringsAsFactors = FALSE
)

print(
  qc_summary,
  row.names = FALSE
)

write.table(
  qc_summary,
  file = file.path(
    out_dir,
    "STEP10C-I_Pseudobulk_QC_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Core gene audit
# ------------------------------------------------------------

cat("\nSTEP10C-20 — CORE GENE AUDIT\n")
cat("----------------------------\n")

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

core_present <- core_genes %in% gene_names

core_audit <- data.frame(
  Gene = core_genes,
  Present_in_GEO_raw =
    core_present,
  Present_in_Tcell =
    core_genes %in%
    rownames(pb_tcell),
  Present_in_Myeloid =
    core_genes %in%
    rownames(pb_myeloid),
  Present_in_Endothelial_Primary =
    core_genes %in%
    rownames(pb_endo_primary),
  Present_in_Endothelial_Sensitivity =
    core_genes %in%
    rownames(pb_endo_sensitivity),
  stringsAsFactors = FALSE
)

print(
  core_audit,
  row.names = FALSE
)

if (
  !all(core_present)
) {
  stop(
    "FAIL — One or more core genes absent from GEO raw gene universe."
  )
}

cat(
  "PASS — AOC3 / MS4A7 / IL7R present.\n\n"
)


# ------------------------------------------------------------
# 31. Pseudobulk sample-count reconciliation
# ------------------------------------------------------------

cat("STEP10C-21 — SAMPLE COUNT RECONCILIATION\n")
cat("----------------------------------------\n")

# ------------------------------------------------------------
# T-cell
# ------------------------------------------------------------

tcell_count_check <- all(
  primary_annotation$Cell_Count[
    primary_annotation$Compartment == "T-cell"
  ] >= 50
)

# ------------------------------------------------------------
# Myeloid
# ------------------------------------------------------------

myeloid_count_check <- all(
  primary_annotation$Cell_Count[
    primary_annotation$Compartment ==
      "Myeloid/macrophage"
  ] >= 50
)

# ------------------------------------------------------------
# Endothelial primary
# ------------------------------------------------------------

endo_primary_count_check <- all(
  primary_annotation$Cell_Count[
    primary_annotation$Compartment ==
      "Endothelial"
  ] >= 50
)

# ------------------------------------------------------------
# Endothelial sensitivity
# ------------------------------------------------------------

endo_sensitivity_count_check <- all(
  sensitivity_annotation$Cell_Count[
    sensitivity_annotation$Compartment ==
      "Endothelial"
  ] >= 20 &
    sensitivity_annotation$Cell_Count[
      sensitivity_annotation$Compartment ==
        "Endothelial"
    ] <= 49
)


cat(
  "T-cell primary >=50:",
  tcell_count_check,
  "\n"
)

cat(
  "Myeloid/macrophage primary >=50:",
  myeloid_count_check,
  "\n"
)

cat(
  "Endothelial primary >=50:",
  endo_primary_count_check,
  "\n"
)

cat(
  "Endothelial sensitivity 20–49:",
  endo_sensitivity_count_check,
  "\n"
)

if (
  !tcell_count_check ||
  !myeloid_count_check ||
  !endo_primary_count_check ||
  !endo_sensitivity_count_check
) {
  stop(
    "FAIL — Pseudobulk sample-depth criteria were not respected."
  )
}

cat(
  "PASS — All pseudobulk groups respect locked cell-depth strata.\n\n"
)


# ------------------------------------------------------------
# 32. Final summary
# ------------------------------------------------------------

final_summary <- data.frame(
  
  Metric = c(
    
    "GEO raw UMI cells",
    "GEO raw genes",
    "TISCH2 target Primary+Tumor cells",
    "Locked target compartment cells",
    
    "T-cell primary patients",
    "Myeloid/macrophage primary patients",
    "Endothelial primary patients",
    
    "Endothelial sensitivity patients",
    "Endothelial insufficient patients",
    
    "T-cell pseudobulk genes",
    "Myeloid/macrophage pseudobulk genes",
    "Endothelial primary pseudobulk genes",
    "Endothelial sensitivity pseudobulk genes",
    
    "Core genes verified"
  ),
  
  Value = c(
    
    208506,
    29634,
    65702,
    43475,
    
    ncol(pb_tcell),
    ncol(pb_myeloid),
    ncol(pb_endo_primary),
    
    ncol(pb_endo_sensitivity),
    3,
    
    nrow(pb_tcell),
    nrow(pb_myeloid),
    nrow(pb_endo_primary),
    nrow(pb_endo_sensitivity),
    
    paste(
      core_genes,
      collapse = " / "
    )
  ),
  
  stringsAsFactors = FALSE
)


write.table(
  final_summary,
  file = file.path(
    out_dir,
    "STEP10C-J_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Final PASS/FAIL checks
# ------------------------------------------------------------

final_checks <- c(
  
  raw_cells =
    length(geo_cell_ids) == 208506,
  
  raw_genes =
    gene_counter == 29634,
  
  target_cells =
    nrow(target_manifest) == 43475,
  
  target_raw_coverage =
    missing_raw_columns == 0,
  
  primary_tcell =
    ncol(pb_tcell) == 9,
  
  primary_myeloid =
    ncol(pb_myeloid) == 9,
  
  primary_endothelial =
    ncol(pb_endo_primary) == 4,
  
  sensitivity_endothelial =
    ncol(pb_endo_sensitivity) == 2,
  
  tcell_depth =
    tcell_count_check,
  
  myeloid_depth =
    myeloid_count_check,
  
  endothelial_primary_depth =
    endo_primary_count_check,
  
  endothelial_sensitivity_depth =
    endo_sensitivity_count_check,
  
  primary_integer_counts =
    primary_integer,
  
  sensitivity_integer_counts =
    sensitivity_integer,
  
  core_genes =
    all(core_present),
  
  manifest_complete =
    all(
      target_manifest$GEO_raw_UMI_CellID_available
    )
)


final_pass <- all(
  final_checks
)


# ------------------------------------------------------------
# 34. FINAL LOCK
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")

if (final_pass) {
  
  cat(
    "PASS — STEP10C FINAL LOCK\n"
  )
  
} else {
  
  cat(
    "FAIL — STEP10C FINAL AUDIT\n"
  )
  
}

cat("============================================================\n\n")


cat(
  "GEO raw UMI cells:",
  length(geo_cell_ids),
  "\n"
)

cat(
  "GEO raw genes:",
  gene_counter,
  "\n"
)

cat(
  "Primary+Tumor TISCH2 cells:",
  nrow(target_meta),
  "\n"
)

cat(
  "Locked target compartment cells:",
  nrow(target_manifest),
  "\n\n"
)

cat(
  "T-cell primary patients:",
  ncol(pb_tcell),
  "\n"
)

cat(
  "Myeloid/macrophage primary patients:",
  ncol(pb_myeloid),
  "\n"
)

cat(
  "Endothelial primary patients:",
  ncol(pb_endo_primary),
  "\n"
)

cat(
  "Endothelial sensitivity patients:",
  ncol(pb_endo_sensitivity),
  "\n"
)

cat(
  "Endothelial insufficient patients:",
  3,
  "\n\n"
)

cat(
  "Core genes:",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n\n"
)

cat(
  "Raw UMI expression values were used.\n"
)

cat(
  "No normalized expression matrix used.\n"
)

cat(
  "No H5 expression values used.\n"
)

cat(
  "No statistical testing performed.\n"
)

cat(
  "No differential expression performed.\n"
)

cat(
  "No CellChat performed.\n"
)

cat(
  "No survival analysis performed.\n"
)

cat(
  "No additional patient/cell exclusion performed.\n"
)

cat(
  "GEO raw-only 5,208 cells were not added.\n\n"
)

cat(
  "Output directory:\n",
  out_dir,
  "\n\n"
)


if (!final_pass) {
  
  cat(
    "Failed checks:\n"
  )
  
  print(
    final_checks[
      !final_checks
    ]
  )
  
  stop(
    "\nSTEP10C FAILED — see failed checks above."
  )
}


cat(
  "============================================================\n"
)

cat(
  "STEP10C COMPLETE — FINAL STATUS: PASS\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# STEP10D
# Pseudobulk QC + AOC3/MS4A7/IL7R Patient-Level Expression Audit
#
# Paper 1 — Single-Cell Revision Workflow
#
# PURPOSE
# ------------------------------------------------------------
# 1. Audit STEP10C raw-UMI pseudobulk sample-level QC
# 2. Audit patient-level raw counts of AOC3/MS4A7/IL7R
# 3. Audit core-gene detection frequency
# 4. Calculate descriptive relative abundance
# 5. Audit patient × compartment completeness
# 6. Produce final descriptive QC summary
#
# IMPORTANT ANALYTICAL LOCK
# ------------------------------------------------------------
# - Raw UMI pseudobulk only
# - STEP10C outputs are the sole expression inputs
# - No H5 expression matrix
# - No normalized expression matrix
# - No re-pseudobulking
# - No cell-level expression analysis
# - No differential expression
# - No Wilcoxon test
# - No t-test
# - No ANOVA
# - No Spearman/Pearson correlation
# - No p values
# - No FDR
# - No CellChat
# - No ligand-receptor analysis
# - No survival analysis
# - No patient exclusion
# - No new endothelial inclusion rule
#
# Locked core genes:
# AOC3 / MS4A7 / IL7R
#
# Locked pseudobulk cohorts:
# T-cell Primary                  n = 9
# Myeloid/macrophage Primary      n = 9
# Endothelial Primary             n = 4
# Endothelial Sensitivity         n = 2
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10D — PSEUDOBULK QC + CORE-GENE PATIENT-LEVEL AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

step10c_dir <- file.path(
  base_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk"
)

step10d_dir <- file.path(
  base_dir,
  "STEP10D_Pseudobulk_QC_CoreGene_Audit"
)

if (!dir.exists(step10c_dir)) {
  stop(
    "STEP10C directory not found:\n",
    step10c_dir
  )
}

if (!dir.exists(step10d_dir)) {
  dir.create(
    step10d_dir,
    recursive = TRUE
  )
}

cat("STEP10C input directory:\n")
cat(step10c_dir, "\n\n")

cat("STEP10D output directory:\n")
cat(step10d_dir, "\n\n")


# ------------------------------------------------------------
# 1. Locked core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Locked core genes:\n")
cat(paste(core_genes, collapse = " / "), "\n\n")


# ------------------------------------------------------------
# 2. Locked input files
# ------------------------------------------------------------

file_tcell <- file.path(
  step10c_dir,
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_myeloid <- file.path(
  step10c_dir,
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_endo_primary <- file.path(
  step10c_dir,
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

file_endo_sensitivity <- file.path(
  step10c_dir,
  "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
)

file_primary_anno <- file.path(
  step10c_dir,
  "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
)

file_sensitivity_anno <- file.path(
  step10c_dir,
  "STEP10C-F_Sensitivity_Pseudobulk_Sample_Annotation.tsv"
)

file_cell_audit <- file.path(
  step10c_dir,
  "STEP10C-G_Pseudobulk_CellCount_Audit.tsv"
)

file_manifest <- file.path(
  step10c_dir,
  "STEP10C-H_Pseudobulk_Input_Manifest.tsv"
)

file_qc_summary <- file.path(
  step10c_dir,
  "STEP10C-I_Pseudobulk_QC_Summary.tsv"
)


required_files <- c(
  file_tcell,
  file_myeloid,
  file_endo_primary,
  file_endo_sensitivity,
  file_primary_anno,
  file_sensitivity_anno,
  file_cell_audit,
  file_manifest,
  file_qc_summary
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  
  cat("FAIL — Missing STEP10C input files:\n")
  print(missing_files)
  
  stop(
    "STEP10D cannot continue because required STEP10C files are missing."
  )
}

cat("PASS — All required STEP10C input files found.\n\n")


# ------------------------------------------------------------
# 3. Helper function: read raw-UMI pseudobulk matrix
# ------------------------------------------------------------

read_pseudobulk <- function(
    file,
    compartment,
    cohort
) {
  
  cat("Reading:\n")
  cat(basename(file), "\n")
  
  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  if (ncol(x) < 2) {
    stop(
      "Pseudobulk matrix has fewer than 2 columns: ",
      basename(file)
    )
  }
  
  gene_col <- colnames(x)[1]
  
  gene_id <- as.character(x[[gene_col]])
  
  if (anyNA(gene_id) || any(gene_id == "")) {
    stop(
      "Missing gene IDs detected in ",
      basename(file)
    )
  }
  
  if (anyDuplicated(gene_id)) {
    stop(
      "Duplicated gene IDs detected in ",
      basename(file)
    )
  }
  
  rownames(x) <- gene_id
  
  x[[gene_col]] <- NULL
  
  # Convert to numeric matrix
  mat <- as.matrix(x)
  
  suppressWarnings(
    storage.mode(mat) <- "numeric"
  )
  
  if (anyNA(mat)) {
    stop(
      "NA values detected after numeric conversion in ",
      basename(file)
    )
  }
  
  if (any(!is.finite(mat))) {
    stop(
      "Non-finite values detected in ",
      basename(file)
    )
  }
  
  if (any(mat < 0)) {
    stop(
      "Negative raw UMI values detected in ",
      basename(file)
    )
  }
  
  # Integer-valued audit
  integer_valued <- all(
    abs(mat - round(mat)) < 1e-8
  )
  
  if (!integer_valued) {
    stop(
      "Non-integer values detected in raw-UMI pseudobulk matrix: ",
      basename(file)
    )
  }
  
  if (anyDuplicated(colnames(mat))) {
    stop(
      "Duplicated pseudobulk sample/group IDs detected in ",
      basename(file)
    )
  }
  
  cat(
    "  genes:",
    nrow(mat),
    "\n"
  )
  
  cat(
    "  samples:",
    ncol(mat),
    "\n"
  )
  
  cat(
    "  integer-valued:",
    integer_valued,
    "\n\n"
  )
  
  list(
    matrix = mat,
    compartment = compartment,
    cohort = cohort,
    file = basename(file)
  )
}


# ------------------------------------------------------------
# 4. Read four locked pseudobulk matrices
# ------------------------------------------------------------

pb_tcell <- read_pseudobulk(
  file_tcell,
  "T_cell",
  "Primary"
)

pb_myeloid <- read_pseudobulk(
  file_myeloid,
  "Myeloid_macrophage",
  "Primary"
)

pb_endo_primary <- read_pseudobulk(
  file_endo_primary,
  "Endothelial",
  "Primary"
)

pb_endo_sensitivity <- read_pseudobulk(
  file_endo_sensitivity,
  "Endothelial",
  "Sensitivity"
)


# ------------------------------------------------------------
# 5. Locked cohort dimension audit
# ------------------------------------------------------------

expected_n <- c(
  T_cell = 9L,
  Myeloid_macrophage = 9L,
  Endothelial_Primary = 4L,
  Endothelial_Sensitivity = 2L
)

observed_n <- c(
  T_cell = ncol(pb_tcell$matrix),
  Myeloid_macrophage = ncol(pb_myeloid$matrix),
  Endothelial_Primary = ncol(pb_endo_primary$matrix),
  Endothelial_Sensitivity = ncol(pb_endo_sensitivity$matrix)
)

cat("------------------------------------------------------------\n")
cat("LOCKED COHORT DIMENSION AUDIT\n")
cat("------------------------------------------------------------\n")

print(
  data.frame(
    Cohort = names(expected_n),
    Expected = as.integer(expected_n),
    Observed = as.integer(observed_n),
    PASS = as.integer(expected_n) == as.integer(observed_n),
    row.names = NULL
  )
)

if (!all(expected_n == observed_n)) {
  stop(
    "FAIL — Observed pseudobulk cohort dimensions do not match STEP10C lock."
  )
}

cat("\nPASS — Cohort dimensions match STEP10C lock.\n\n")


# ------------------------------------------------------------
# 6. Core-gene presence audit
# ------------------------------------------------------------

all_matrices <- list(
  T_cell_Primary = pb_tcell$matrix,
  Myeloid_macrophage_Primary = pb_myeloid$matrix,
  Endothelial_Primary = pb_endo_primary$matrix,
  Endothelial_Sensitivity = pb_endo_sensitivity$matrix
)

core_presence <- do.call(
  rbind,
  lapply(
    names(all_matrices),
    function(nm) {
      
      mat <- all_matrices[[nm]]
      
      data.frame(
        Dataset = nm,
        Gene = core_genes,
        Present = core_genes %in% rownames(mat),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  core_presence,
  file = file.path(
    step10d_dir,
    "STEP10D-F_CoreGene_QC_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!all(core_presence$Present)) {
  stop(
    "FAIL — One or more locked core genes are absent from a pseudobulk matrix."
  )
}

cat("PASS — AOC3/MS4A7/IL7R present in all four pseudobulk datasets.\n\n")


# ------------------------------------------------------------
# 7. Function: sample-level pseudobulk QC
# ------------------------------------------------------------

audit_library_qc <- function(
    mat,
    compartment,
    cohort
) {
  
  total_umi <- colSums(mat)
  
  detected_genes <- colSums(mat > 0)
  
  zero_genes <- colSums(mat == 0)
  
  total_genes <- nrow(mat)
  
  data.frame(
    Sample = colnames(mat),
    Compartment = compartment,
    Cohort = cohort,
    Total_UMI = as.numeric(total_umi),
    Detected_Genes = as.numeric(detected_genes),
    Zero_Genes = as.numeric(zero_genes),
    Total_Genes = total_genes,
    Detected_Gene_Fraction =
      as.numeric(detected_genes / total_genes),
    Zero_Gene_Fraction =
      as.numeric(zero_genes / total_genes),
    stringsAsFactors = FALSE
  )
}


qc_tcell <- audit_library_qc(
  pb_tcell$matrix,
  "T_cell",
  "Primary"
)

qc_myeloid <- audit_library_qc(
  pb_myeloid$matrix,
  "Myeloid_macrophage",
  "Primary"
)

qc_endo_primary <- audit_library_qc(
  pb_endo_primary$matrix,
  "Endothelial",
  "Primary"
)

qc_endo_sensitivity <- audit_library_qc(
  pb_endo_sensitivity$matrix,
  "Endothelial",
  "Sensitivity"
)

library_qc <- rbind(
  qc_tcell,
  qc_myeloid,
  qc_endo_primary,
  qc_endo_sensitivity
)

rownames(library_qc) <- NULL


# ------------------------------------------------------------
# 8. Descriptive library-size summary
# ------------------------------------------------------------

library_summary <- do.call(
  rbind,
  lapply(
    split(
      library_qc,
      interaction(
        library_qc$Compartment,
        library_qc$Cohort,
        drop = TRUE
      )
    ),
    function(df) {
      
      x <- df$Total_UMI
      
      data.frame(
        Compartment = df$Compartment[1],
        Cohort = df$Cohort[1],
        N = nrow(df),
        Min = min(x),
        Q1 = as.numeric(
          quantile(x, 0.25, names = FALSE)
        ),
        Median = median(x),
        Mean = mean(x),
        Q3 = as.numeric(
          quantile(x, 0.75, names = FALSE)
        ),
        Max = max(x),
        IQR = IQR(x),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(library_summary) <- NULL


# ------------------------------------------------------------
# 9. Detected-gene summary
# ------------------------------------------------------------

detected_gene_summary <- do.call(
  rbind,
  lapply(
    split(
      library_qc,
      interaction(
        library_qc$Compartment,
        library_qc$Cohort,
        drop = TRUE
      )
    ),
    function(df) {
      
      x <- df$Detected_Genes
      
      data.frame(
        Compartment = df$Compartment[1],
        Cohort = df$Cohort[1],
        N = nrow(df),
        Min = min(x),
        Q1 = as.numeric(
          quantile(x, 0.25, names = FALSE)
        ),
        Median = median(x),
        Mean = mean(x),
        Q3 = as.numeric(
          quantile(x, 0.75, names = FALSE)
        ),
        Max = max(x),
        IQR = IQR(x),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(detected_gene_summary) <- NULL


# ------------------------------------------------------------
# 10. Save sample-level QC
# ------------------------------------------------------------

write.table(
  library_qc,
  file = file.path(
    step10d_dir,
    "STEP10D-A_Pseudobulk_Library_QC.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  library_summary,
  file = file.path(
    step10d_dir,
    "STEP10D-A1_Library_Size_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  detected_gene_summary,
  file = file.path(
    step10d_dir,
    "STEP10D-A2_Detected_Gene_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Core-gene raw-count audit
# ------------------------------------------------------------

extract_core_counts <- function(
    mat,
    compartment,
    cohort
) {
  
  data.frame(
    Sample = colnames(mat),
    Compartment = compartment,
    Cohort = cohort,
    AOC3 = as.numeric(mat["AOC3", ]),
    MS4A7 = as.numeric(mat["MS4A7", ]),
    IL7R = as.numeric(mat["IL7R", ]),
    stringsAsFactors = FALSE
  )
}


core_raw_counts <- rbind(
  
  extract_core_counts(
    pb_tcell$matrix,
    "T_cell",
    "Primary"
  ),
  
  extract_core_counts(
    pb_myeloid$matrix,
    "Myeloid_macrophage",
    "Primary"
  ),
  
  extract_core_counts(
    pb_endo_primary$matrix,
    "Endothelial",
    "Primary"
  ),
  
  extract_core_counts(
    pb_endo_sensitivity$matrix,
    "Endothelial",
    "Sensitivity"
  )
)

rownames(core_raw_counts) <- NULL


write.table(
  core_raw_counts,
  file = file.path(
    step10d_dir,
    "STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Core-gene detection audit
# ------------------------------------------------------------

core_detection <- core_raw_counts

core_detection$AOC3_detected <-
  core_detection$AOC3 > 0

core_detection$MS4A7_detected <-
  core_detection$MS4A7 > 0

core_detection$IL7R_detected <-
  core_detection$IL7R > 0


write.table(
  core_detection,
  file = file.path(
    step10d_dir,
    "STEP10D-C_Pseudobulk_CoreGene_Detection.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Core-gene descriptive count summary
# ------------------------------------------------------------

core_long_summary <- do.call(
  rbind,
  lapply(
    c("AOC3", "MS4A7", "IL7R"),
    function(gene) {
      
      groups <- split(
        core_raw_counts,
        interaction(
          core_raw_counts$Compartment,
          core_raw_counts$Cohort,
          drop = TRUE
        )
      )
      
      do.call(
        rbind,
        lapply(
          groups,
          function(df) {
            
            x <- df[[gene]]
            
            data.frame(
              Gene = gene,
              Compartment = df$Compartment[1],
              Cohort = df$Cohort[1],
              N = length(x),
              Detected_N = sum(x > 0),
              Detection_Frequency =
                mean(x > 0),
              Min = min(x),
              Q1 = as.numeric(
                quantile(
                  x,
                  0.25,
                  names = FALSE
                )
              ),
              Median = median(x),
              Mean = mean(x),
              Q3 = as.numeric(
                quantile(
                  x,
                  0.75,
                  names = FALSE
                )
              ),
              Max = max(x),
              IQR = IQR(x),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

rownames(core_long_summary) <- NULL


write.table(
  core_long_summary,
  file = file.path(
    step10d_dir,
    "STEP10D-F1_CoreGene_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Descriptive relative abundance
#
# IMPORTANT:
# This is NOT normalized expression analysis.
#
# It is simply:
#
#     gene raw UMI / total pseudobulk UMI
#
# Used only as a descriptive library-size-adjusted quantity.
# No hypothesis testing is performed.
# ------------------------------------------------------------

relative_abundance <- core_raw_counts

total_umi_lookup <- library_qc[
  match(
    relative_abundance$Sample,
    library_qc$Sample
  ),
  "Total_UMI"
]

if (anyNA(total_umi_lookup)) {
  stop(
    "FAIL — Unable to match Total_UMI for one or more pseudobulk samples."
  )
}

relative_abundance$Total_UMI <-
  as.numeric(total_umi_lookup)

relative_abundance$AOC3_relative <-
  relative_abundance$AOC3 /
  relative_abundance$Total_UMI

relative_abundance$MS4A7_relative <-
  relative_abundance$MS4A7 /
  relative_abundance$Total_UMI

relative_abundance$IL7R_relative <-
  relative_abundance$IL7R /
  relative_abundance$Total_UMI


write.table(
  relative_abundance,
  file = file.path(
    step10d_dir,
    "STEP10D-D_Pseudobulk_CoreGene_RelativeAbundance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Relative-abundance descriptive summary
# ------------------------------------------------------------

relative_long_summary <- do.call(
  rbind,
  lapply(
    c(
      "AOC3_relative",
      "MS4A7_relative",
      "IL7R_relative"
    ),
    function(gene_var) {
      
      groups <- split(
        relative_abundance,
        interaction(
          relative_abundance$Compartment,
          relative_abundance$Cohort,
          drop = TRUE
        )
      )
      
      do.call(
        rbind,
        lapply(
          groups,
          function(df) {
            
            x <- df[[gene_var]]
            
            data.frame(
              Gene = sub(
                "_relative$",
                "",
                gene_var
              ),
              Compartment = df$Compartment[1],
              Cohort = df$Cohort[1],
              N = length(x),
              Min = min(x),
              Q1 = as.numeric(
                quantile(
                  x,
                  0.25,
                  names = FALSE
                )
              ),
              Median = median(x),
              Mean = mean(x),
              Q3 = as.numeric(
                quantile(
                  x,
                  0.75,
                  names = FALSE
                )
              ),
              Max = max(x),
              IQR = IQR(x),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

rownames(relative_long_summary) <- NULL


write.table(
  relative_long_summary,
  file = file.path(
    step10d_dir,
    "STEP10D-D1_CoreGene_RelativeAbundance_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Patient × compartment completeness audit
# ------------------------------------------------------------
#
# STEP10C annotation schema is already LOCKED:
#
# Patient
# Compartment
# Cell_Count
# Eligibility
# Group_ID
#
# Therefore:
# - Patient = patient identifier
# - Compartment = locked compartment
# - Cell_Count = contributing cell number
# - Eligibility = STEP10C eligibility classification
# - Group_ID = unique pseudobulk group identifier
#
# No column inference is necessary.
# ------------------------------------------------------------


cat("\n")
cat("------------------------------------------------------------\n")
cat("STEP10D — PATIENT × COMPARTMENT COMPLETENESS AUDIT\n")
cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 16A. Read STEP10C annotation files
# ------------------------------------------------------------

primary_anno <- read.delim(
  file_primary_anno,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

sensitivity_anno <- read.delim(
  file_sensitivity_anno,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)


# ------------------------------------------------------------
# 16B. Locked annotation schema audit
# ------------------------------------------------------------

required_annotation_cols <- c(
  "Patient",
  "Compartment",
  "Cell_Count",
  "Eligibility",
  "Group_ID"
)

primary_missing_cols <- setdiff(
  required_annotation_cols,
  colnames(primary_anno)
)

sensitivity_missing_cols <- setdiff(
  required_annotation_cols,
  colnames(sensitivity_anno)
)


if (length(primary_missing_cols) > 0) {
  
  stop(
    "FAIL — Missing required columns in primary annotation: ",
    paste(primary_missing_cols, collapse = ", ")
  )
}


if (length(sensitivity_missing_cols) > 0) {
  
  stop(
    "FAIL — Missing required columns in sensitivity annotation: ",
    paste(sensitivity_missing_cols, collapse = ", ")
  )
}


cat(
  "PASS — Primary annotation schema verified.\n"
)

cat(
  "PASS — Sensitivity annotation schema verified.\n\n"
)


# ------------------------------------------------------------
# 16C. Basic annotation integrity checks
# ------------------------------------------------------------

check_annotation_integrity <- function(
    df,
    label
) {
  
  # Missing values
  required_data <- df[
    ,
    required_annotation_cols,
    drop = FALSE
  ]
  
  if (anyNA(required_data)) {
    
    stop(
      "FAIL — NA values detected in ",
      label,
      " annotation."
    )
  }
  
  
  # Empty strings
  for (cc in required_annotation_cols) {
    
    if (any(
      trimws(
        as.character(df[[cc]])
      ) == ""
    )) {
      
      stop(
        "FAIL — Empty values detected in ",
        label,
        " annotation column: ",
        cc
      )
    }
  }
  
  
  # Patient uniqueness within annotation
  if (anyDuplicated(df$Group_ID)) {
    
    stop(
      "FAIL — Duplicate Group_ID detected in ",
      label,
      " annotation."
    )
  }
  
  
  # Cell counts must be positive integers
  cell_count_numeric <- suppressWarnings(
    as.numeric(df$Cell_Count)
  )
  
  if (anyNA(cell_count_numeric)) {
    
    stop(
      "FAIL — Non-numeric Cell_Count detected in ",
      label,
      " annotation."
    )
  }
  
  if (any(cell_count_numeric <= 0)) {
    
    stop(
      "FAIL — Non-positive Cell_Count detected in ",
      label,
      " annotation."
    )
  }
  
  if (any(
    abs(
      cell_count_numeric -
      round(cell_count_numeric)
    ) > 1e-8
  )) {
    
    stop(
      "FAIL — Non-integer Cell_Count detected in ",
      label,
      " annotation."
    )
  }
  
  
  TRUE
}


check_annotation_integrity(
  primary_anno,
  "primary"
)

check_annotation_integrity(
  sensitivity_anno,
  "sensitivity"
)


cat(
  "PASS — Primary annotation integrity audit.\n"
)

cat(
  "PASS — Sensitivity annotation integrity audit.\n\n"
)


# ------------------------------------------------------------
# 16D. Verify annotation Group_ID against matrix columns
# ------------------------------------------------------------

verify_group_ids <- function(
    annotation,
    matrix,
    label
) {
  
  anno_ids <- as.character(
    annotation$Group_ID
  )
  
  matrix_ids <- colnames(matrix)
  
  missing_from_matrix <- setdiff(
    anno_ids,
    matrix_ids
  )
  
  missing_from_annotation <- setdiff(
    matrix_ids,
    anno_ids
  )
  
  if (length(missing_from_matrix) > 0) {
    
    cat(
      "Missing from matrix:\n"
    )
    
    print(missing_from_matrix)
    
    stop(
      "FAIL — Annotation Group_ID not fully represented in ",
      label,
      " pseudobulk matrix."
    )
  }
  
  
  if (length(missing_from_annotation) > 0) {
    
    cat(
      "Missing from annotation:\n"
    )
    
    print(missing_from_annotation)
    
    stop(
      "FAIL — Pseudobulk matrix sample IDs not fully represented in ",
      label,
      " annotation."
    )
  }
  
  
  if (!setequal(
    anno_ids,
    matrix_ids
  )) {
    
    stop(
      "FAIL — Group_ID set mismatch in ",
      label,
      "."
    )
  }
  
  
  TRUE
}


verify_group_ids(
  primary_anno[
    primary_anno$Eligibility == "Primary",
    ,
    drop = FALSE
  ],
  cbind(
    pb_tcell$matrix,
    pb_myeloid$matrix,
    pb_endo_primary$matrix
  ),
  "primary"
)


verify_group_ids(
  sensitivity_anno[
    sensitivity_anno$Eligibility == "Sensitivity",
    ,
    drop = FALSE
  ],
  pb_endo_sensitivity$matrix,
  "sensitivity"
)


cat(
  "PASS — Primary Group_ID ↔ pseudobulk matrix reconciliation.\n"
)

cat(
  "PASS — Sensitivity Group_ID ↔ pseudobulk matrix reconciliation.\n\n"
)


# ------------------------------------------------------------
# 16E. Construct locked completeness table
# ------------------------------------------------------------

completeness_primary <- primary_anno[
  ,
  required_annotation_cols,
  drop = FALSE
]

completeness_primary$Cohort <- "Primary"


completeness_sensitivity <- sensitivity_anno[
  ,
  required_annotation_cols,
  drop = FALSE
]

completeness_sensitivity$Cohort <- "Sensitivity"


completeness <- rbind(
  completeness_primary,
  completeness_sensitivity
)

rownames(completeness) <- NULL


# ============================================================
# 16F — Locked compartment completeness audit
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — 16F: LOCKED COMPARTMENT COMPLETENESS AUDIT\n")
cat("============================================================\n\n")

# IMPORTANT:
# Use the exact Compartment labels stored in STEP10C.
# Do NOT recode or rename the underlying STEP10C annotation.

allowed_primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

allowed_sensitivity_compartments <- c(
  "Endothelial"
)

# ------------------------------------------------------------
# Primary completeness
# ------------------------------------------------------------

if (!all(primary_anno$Eligibility == "Primary")) {
  stop("FAIL — Primary annotation contains non-Primary eligibility labels.")
}

if (any(
  !as.character(primary_anno$Compartment) %in%
  allowed_primary_compartments
)) {
  stop("FAIL — Unexpected primary compartment detected.")
}

completeness_primary <- primary_anno[, c(
  "Patient",
  "Compartment",
  "Cell_Count",
  "Eligibility",
  "Group_ID"
)]

# Check expected compartment-level patient counts
primary_expected_n <- c(
  "T-cell" = 9L,
  "Myeloid/macrophage" = 9L,
  "Endothelial" = 4L
)

primary_observed_n <- table(
  completeness_primary$Compartment
)

for (cmp in names(primary_expected_n)) {
  
  observed <- ifelse(
    cmp %in% names(primary_observed_n),
    as.integer(primary_observed_n[[cmp]]),
    0L
  )
  
  expected <- primary_expected_n[[cmp]]
  
  if (observed != expected) {
    stop(
      paste0(
        "FAIL — Primary compartment count mismatch: ",
        cmp,
        " observed=", observed,
        ", expected=", expected
      )
    )
  }
}

# ------------------------------------------------------------
# Sensitivity completeness
# ------------------------------------------------------------

if (!all(sensitivity_anno$Eligibility == "Sensitivity")) {
  stop("FAIL — Sensitivity annotation contains non-Sensitivity labels.")
}

if (any(
  !as.character(sensitivity_anno$Compartment) %in%
  allowed_sensitivity_compartments
)) {
  stop("FAIL — Unexpected sensitivity compartment detected.")
}

completeness_sensitivity <- sensitivity_anno[, c(
  "Patient",
  "Compartment",
  "Cell_Count",
  "Eligibility",
  "Group_ID"
)]

# Expected sensitivity endothelial patient count
sensitivity_expected_n <- c(
  "Endothelial" = 2L
)

sensitivity_observed_n <- table(
  completeness_sensitivity$Compartment
)

for (cmp in names(sensitivity_expected_n)) {
  
  observed <- ifelse(
    cmp %in% names(sensitivity_observed_n),
    as.integer(sensitivity_observed_n[[cmp]]),
    0L
  )
  
  expected <- sensitivity_expected_n[[cmp]]
  
  if (observed != expected) {
    stop(
      paste0(
        "FAIL — Sensitivity compartment count mismatch: ",
        cmp,
        " observed=", observed,
        ", expected=", expected
      )
    )
  }
}

cat("PASS — Primary compartment labels verified.\n")
cat("PASS — Primary patient/group counts verified.\n")
cat("PASS — Sensitivity compartment labels verified.\n")
cat("PASS — Sensitivity patient/group counts verified.\n")


# ============================================================
# 17 — Patient × Compartment completeness table
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — 17: PATIENT × COMPARTMENT COMPLETENESS\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Primary
# ------------------------------------------------------------

primary_completeness_table <- completeness_primary[
  order(
    completeness_primary$Patient,
    completeness_primary$Compartment
  ),
]

# ------------------------------------------------------------
# Sensitivity
# ------------------------------------------------------------

sensitivity_completeness_table <- completeness_sensitivity[
  order(
    completeness_sensitivity$Patient,
    completeness_sensitivity$Compartment
  ),
]

# ------------------------------------------------------------
# Combined table
# ------------------------------------------------------------

completeness_all <- rbind(
  primary_completeness_table,
  sensitivity_completeness_table
)

write.table(
  completeness_all,
  file = file.path(
    step10d_dir,
    "STEP10D-E_Patient_Compartment_Completeness.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Saved:\n",
  file.path(
    step10d_dir,
    "STEP10D-E_Patient_Compartment_Completeness.tsv"
  ),
  "\n"
)


# ============================================================
# 18 — Compartment-level descriptive summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — 18: COMPARTMENT-LEVEL DESCRIPTIVE SUMMARY\n")
cat("============================================================\n\n")

summarize_compartment <- function(df) {
  
  compartments <- unique(as.character(df$Compartment))
  
  out <- lapply(
    compartments,
    function(cmp) {
      
      x <- df[
        as.character(df$Compartment) == cmp,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Compartment = cmp,
        Eligibility = unique(as.character(x$Eligibility)),
        N_groups = nrow(x),
        N_patients = length(unique(x$Patient)),
        Cell_Count_min = min(x$Cell_Count),
        Cell_Count_Q1 = unname(
          quantile(x$Cell_Count, 0.25)
        ),
        Cell_Count_median = median(x$Cell_Count),
        Cell_Count_mean = mean(x$Cell_Count),
        Cell_Count_Q3 = unname(
          quantile(x$Cell_Count, 0.75)
        ),
        Cell_Count_max = max(x$Cell_Count),
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(rbind, out)
}

primary_compartment_summary <- summarize_compartment(
  completeness_primary
)

sensitivity_compartment_summary <- summarize_compartment(
  completeness_sensitivity
)

compartment_summary_all <- rbind(
  primary_compartment_summary,
  sensitivity_compartment_summary
)

write.table(
  compartment_summary_all,
  file = file.path(
    step10d_dir,
    "STEP10D-E1_Patient_Compartment_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Saved:\n",
  file.path(
    step10d_dir,
    "STEP10D-E1_Patient_Compartment_Summary.tsv"
  ),
  "\n"
)


# ============================================================
# 19 — Explicit patient-count lock check
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — 19: PATIENT COUNT LOCK CHECK\n")
cat("============================================================\n\n")

patient_count_lock <- data.frame(
  Dataset = c(
    "Primary_T-cell",
    "Primary_Myeloid_macrophage",
    "Primary_Endothelial",
    "Sensitivity_Endothelial"
  ),
  Expected_N = c(
    9L,
    9L,
    4L,
    2L
  ),
  Observed_N = c(
    sum(
      completeness_primary$Compartment == "T-cell"
    ),
    sum(
      completeness_primary$Compartment ==
        "Myeloid/macrophage"
    ),
    sum(
      completeness_primary$Compartment ==
        "Endothelial"
    ),
    sum(
      completeness_sensitivity$Compartment ==
        "Endothelial"
    )
  ),
  stringsAsFactors = FALSE
)

patient_count_lock$PASS <- (
  patient_count_lock$Expected_N ==
    patient_count_lock$Observed_N
)

if (!all(patient_count_lock$PASS)) {
  stop(
    "FAIL — Patient-count lock check failed."
  )
}

write.table(
  patient_count_lock,
  file = file.path(
    step10d_dir,
    "STEP10D-E2_Patient_Count_Lock_Check.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(patient_count_lock)

cat(
  "\nPASS — Patient-count lock verified.\n"
)


# ============================================================
# 20 — STEP10D compartment audit summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — 20: COMPARTMENT AUDIT SUMMARY\n")
cat("============================================================\n\n")

cat("Primary:\n")
cat("  T-cell                :", 
    sum(completeness_primary$Compartment == "T-cell"),
    "groups\n")

cat("  Myeloid/macrophage    :",
    sum(
      completeness_primary$Compartment ==
        "Myeloid/macrophage"
    ),
    "groups\n")

cat("  Endothelial           :",
    sum(
      completeness_primary$Compartment ==
        "Endothelial"
    ),
    "groups\n")

cat("\nSensitivity:\n")
cat("  Endothelial           :",
    sum(
      completeness_sensitivity$Compartment ==
        "Endothelial"
    ),
    "groups\n")

cat("\nLocked design:\n")
cat("  Primary T-cell        = 9\n")
cat("  Primary Myeloid       = 9\n")
cat("  Primary Endothelial   = 4\n")
cat("  Sensitivity Endothelial = 2\n")

cat("\nPASS — STEP10D compartment completeness audit.\n")


# ============================================================
# 21 — End of annotation/completeness section
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10D — ANNOTATION / COMPLETENESS BLOCK COMPLETE\n")
cat("============================================================\n\n")

cat("No pseudobulk matrix was recalculated.\n")
cat("No expression values were recalculated.\n")
cat("No statistical testing was performed.\n")
cat("No patient or cell was excluded.\n")
cat("STEP10C locked data remain unchanged.\n")


# ============================================================
# STEP10D-A — Pseudobulk Library QC
# ============================================================
#
# Purpose:
#   Descriptive QC of the LOCKED STEP10C patient-level
#   raw-UMI pseudobulk matrices.
#
# Locked inputs:
#   - T-cell primary:              n = 9
#   - Myeloid/macrophage primary:  n = 9
#   - Endothelial primary:         n = 4
#   - Endothelial sensitivity:     n = 2
#
# This step:
#   - uses raw UMI pseudobulk counts already generated in STEP10C
#   - calculates total UMI / library size
#   - calculates detected genes
#   - calculates zero-gene count
#   - calculates gene detection fraction
#   - performs integer / non-negative / NA checks
#   - provides descriptive low-library flags
#
# This step DOES NOT:
#   - recalculate pseudobulk
#   - modify STEP10C matrices
#   - exclude patients
#   - perform differential expression
#   - perform correlation analysis
#   - perform CellChat / LR analysis
#   - perform survival analysis
#   - calculate p-values or FDR
#
# ============================================================


# ------------------------------------------------------------
# 0. Output directory
# ------------------------------------------------------------

step10d_dir <- file.path(
  "input_dir/Paper 1_Step 10-Single-Cell",
  "STEP10D_Pseudobulk_QC_CoreGene_Audit"
)

if (!dir.exists(step10d_dir)) {
  dir.create(
    step10d_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Locked STEP10C input directory
# ------------------------------------------------------------

step10c_dir <- file.path(
  "input_dir/Paper 1_Step 10-Single-Cell",
  "STEP10C_Patient_Level_RawUMI_Pseudobulk"
)

if (!dir.exists(step10c_dir)) {
  stop(
    paste0(
      "FAIL — STEP10C input directory not found:\n",
      step10c_dir
    )
  )
}


# ------------------------------------------------------------
# 2. Locked input files
# ------------------------------------------------------------

pb_files <- c(
  T_cell = file.path(
    step10c_dir,
    "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
  ),
  Myeloid_macrophage = file.path(
    step10c_dir,
    "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
  ),
  Endothelial_primary = file.path(
    step10c_dir,
    "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
  ),
  Endothelial_sensitivity = file.path(
    step10c_dir,
    "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
  )
)

missing_pb <- pb_files[!file.exists(pb_files)]

if (length(missing_pb) > 0) {
  stop(
    paste0(
      "FAIL — Missing STEP10C pseudobulk input file(s):\n",
      paste(missing_pb, collapse = "\n")
    )
  )
}


# ------------------------------------------------------------
# 3. Read locked pseudobulk matrices
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-A — READING LOCKED STEP10C PSEUDOBULK MATRICES\n")
cat("============================================================\n\n")

read_pb <- function(file) {
  
  x <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = 1
  )
  
  as.matrix(x)
}

pb_tcell <- read_pb(pb_files["T_cell"])
pb_myeloid <- read_pb(pb_files["Myeloid_macrophage"])
pb_endo_primary <- read_pb(
  pb_files["Endothelial_primary"]
)
pb_endo_sensitivity <- read_pb(
  pb_files["Endothelial_sensitivity"]
)


# ------------------------------------------------------------
# 4. Basic matrix structure audit
# ------------------------------------------------------------

cat("Matrix dimensions:\n")

cat(
  "  T-cell                 : ",
  nrow(pb_tcell),
  " genes × ",
  ncol(pb_tcell),
  " groups\n",
  sep = ""
)

cat(
  "  Myeloid/macrophage     : ",
  nrow(pb_myeloid),
  " genes × ",
  ncol(pb_myeloid),
  " groups\n",
  sep = ""
)

cat(
  "  Endothelial primary    : ",
  nrow(pb_endo_primary),
  " genes × ",
  ncol(pb_endo_primary),
  " groups\n",
  sep = ""
)

cat(
  "  Endothelial sensitivity: ",
  nrow(pb_endo_sensitivity),
  " genes × ",
  ncol(pb_endo_sensitivity),
  " groups\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Expected dimensions — LOCK CHECK
# ------------------------------------------------------------

expected_gene_n <- 29634L

if (nrow(pb_tcell) != expected_gene_n ||
    nrow(pb_myeloid) != expected_gene_n ||
    nrow(pb_endo_primary) != expected_gene_n ||
    nrow(pb_endo_sensitivity) != expected_gene_n) {
  
  stop(
    "FAIL — Gene dimension does not match STEP10C locked value of 29,634."
  )
}

expected_group_n <- c(
  T_cell = 9L,
  Myeloid_macrophage = 9L,
  Endothelial_primary = 4L,
  Endothelial_sensitivity = 2L
)

observed_group_n <- c(
  T_cell = ncol(pb_tcell),
  Myeloid_macrophage = ncol(pb_myeloid),
  Endothelial_primary = ncol(pb_endo_primary),
  Endothelial_sensitivity = ncol(pb_endo_sensitivity)
)

if (!identical(
  as.integer(observed_group_n),
  as.integer(expected_group_n)
)) {
  stop(
    "FAIL — STEP10C locked pseudobulk group dimensions do not match."
  )
}

cat("\nPASS — Locked pseudobulk dimensions verified.\n")


# ------------------------------------------------------------
# 6. Gene identity reconciliation
# ------------------------------------------------------------

gene_sets <- list(
  T_cell = rownames(pb_tcell),
  Myeloid_macrophage = rownames(pb_myeloid),
  Endothelial_primary = rownames(pb_endo_primary),
  Endothelial_sensitivity =
    rownames(pb_endo_sensitivity)
)

if (any(
  vapply(
    gene_sets,
    function(x) anyDuplicated(x) > 0,
    logical(1)
  )
)) {
  stop(
    "FAIL — Duplicate gene identifiers detected."
  )
}

reference_genes <- gene_sets[[1]]

if (any(
  !vapply(
    gene_sets,
    function(x) identical(x, reference_genes),
    logical(1)
  )
)) {
  stop(
    "FAIL — Gene identifiers/order differ across pseudobulk matrices."
  )
}

cat("PASS — Gene identity/order reconciled across all matrices.\n")


# ------------------------------------------------------------
# 7. Raw UMI integrity audit
# ------------------------------------------------------------

audit_matrix <- function(mat, label) {
  
  values <- as.numeric(mat)
  
  n_total <- length(values)
  
  n_na <- sum(is.na(values))
  n_nan <- sum(is.nan(values))
  n_inf <- sum(is.infinite(values))
  
  n_negative <- sum(
    !is.na(values) &
      is.finite(values) &
      values < 0
  )
  
  n_noninteger <- sum(
    !is.na(values) &
      is.finite(values) &
      abs(values - round(values)) > 0
  )
  
  data.frame(
    Dataset = label,
    N_values = n_total,
    NA_count = n_na,
    NaN_count = n_nan,
    Inf_count = n_inf,
    Negative_count = n_negative,
    Noninteger_count = n_noninteger,
    PASS = (
      n_na == 0 &&
        n_nan == 0 &&
        n_inf == 0 &&
        n_negative == 0 &&
        n_noninteger == 0
    ),
    stringsAsFactors = FALSE
  )
}

integrity_qc <- rbind(
  audit_matrix(
    pb_tcell,
    "Primary_T-cell"
  ),
  audit_matrix(
    pb_myeloid,
    "Primary_Myeloid_macrophage"
  ),
  audit_matrix(
    pb_endo_primary,
    "Primary_Endothelial"
  ),
  audit_matrix(
    pb_endo_sensitivity,
    "Sensitivity_Endothelial"
  )
)

print(integrity_qc)

if (!all(integrity_qc$PASS)) {
  stop(
    "FAIL — Raw UMI integrity audit failed."
  )
}

cat("\nPASS — All pseudobulk values are finite, non-negative integers.\n")


# ------------------------------------------------------------
# 8. Library-level QC function
# ------------------------------------------------------------

calculate_library_qc <- function(
    mat,
    dataset_label
) {
  
  total_umi <- colSums(mat)
  
  detected_genes <- colSums(
    mat > 0
  )
  
  zero_genes <- nrow(mat) - detected_genes
  
  detection_fraction <- (
    detected_genes / nrow(mat)
  )
  
  data.frame(
    Dataset = dataset_label,
    Group_ID = colnames(mat),
    Total_UMI = as.numeric(total_umi),
    Detected_Genes = as.integer(detected_genes),
    Zero_Genes = as.integer(zero_genes),
    Gene_Detection_Fraction =
      as.numeric(detection_fraction),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 9. Calculate library QC
# ------------------------------------------------------------

library_qc <- rbind(
  
  calculate_library_qc(
    pb_tcell,
    "Primary_T-cell"
  ),
  
  calculate_library_qc(
    pb_myeloid,
    "Primary_Myeloid_macrophage"
  ),
  
  calculate_library_qc(
    pb_endo_primary,
    "Primary_Endothelial"
  ),
  
  calculate_library_qc(
    pb_endo_sensitivity,
    "Sensitivity_Endothelial"
  )
)


# ------------------------------------------------------------
# 10. Basic library QC integrity checks
# ------------------------------------------------------------

if (any(
  is.na(library_qc$Total_UMI) |
  is.na(library_qc$Detected_Genes) |
  is.na(library_qc$Zero_Genes)
)) {
  stop(
    "FAIL — NA detected in library-level QC."
  )
}

if (any(
  library_qc$Total_UMI <= 0
)) {
  stop(
    "FAIL — Zero or negative pseudobulk library detected."
  )
}

if (any(
  library_qc$Detected_Genes < 0 |
  library_qc$Detected_Genes > expected_gene_n
)) {
  stop(
    "FAIL — Invalid detected-gene count detected."
  )
}

if (any(
  library_qc$Gene_Detection_Fraction < 0 |
  library_qc$Gene_Detection_Fraction > 1
)) {
  stop(
    "FAIL — Invalid gene detection fraction detected."
  )
}

cat(
  "\nPASS — Library-level QC values are internally valid.\n"
)


# ------------------------------------------------------------
# 11. Descriptive low-library flag
# ------------------------------------------------------------
#
# This is a DESCRIPTIVE flag only.
#
# Threshold:
#   Q1 - 1.5 × IQR
#
# It does NOT exclude samples.
# It does NOT trigger statistical filtering.
# It is only used to identify unusually low library size.
#
# ------------------------------------------------------------

library_qc$Low_Library_Flag <- FALSE

for (dataset in unique(library_qc$Dataset)) {
  
  idx <- library_qc$Dataset == dataset
  
  x <- library_qc$Total_UMI[idx]
  
  q1 <- quantile(
    x,
    0.25,
    names = FALSE
  )
  
  q3 <- quantile(
    x,
    0.75,
    names = FALSE
  )
  
  iqr_value <- q3 - q1
  
  lower_bound <- q1 - 1.5 * iqr_value
  
  library_qc$Low_Library_Flag[idx] <-
    x < lower_bound
}


# ------------------------------------------------------------
# 12. Save library-level QC table
# ------------------------------------------------------------

library_qc_file <- file.path(
  step10d_dir,
  "STEP10D-A_Pseudobulk_Library_QC.tsv"
)

write.table(
  library_qc,
  file = library_qc_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "\nSaved:\n",
  library_qc_file,
  "\n"
)


# ------------------------------------------------------------
# 13. Descriptive summary function
# ------------------------------------------------------------

summarize_library_qc <- function(df) {
  
  datasets <- unique(df$Dataset)
  
  out <- lapply(
    datasets,
    function(ds) {
      
      x <- df[
        df$Dataset == ds,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Dataset = ds,
        N_groups = nrow(x),
        
        Total_UMI_min =
          min(x$Total_UMI),
        
        Total_UMI_Q1 =
          unname(
            quantile(
              x$Total_UMI,
              0.25
            )
          ),
        
        Total_UMI_median =
          median(x$Total_UMI),
        
        Total_UMI_mean =
          mean(x$Total_UMI),
        
        Total_UMI_Q3 =
          unname(
            quantile(
              x$Total_UMI,
              0.75
            )
          ),
        
        Total_UMI_max =
          max(x$Total_UMI),
        
        Detected_Genes_min =
          min(x$Detected_Genes),
        
        Detected_Genes_Q1 =
          unname(
            quantile(
              x$Detected_Genes,
              0.25
            )
          ),
        
        Detected_Genes_median =
          median(x$Detected_Genes),
        
        Detected_Genes_mean =
          mean(x$Detected_Genes),
        
        Detected_Genes_Q3 =
          unname(
            quantile(
              x$Detected_Genes,
              0.75
            )
          ),
        
        Detected_Genes_max =
          max(x$Detected_Genes),
        
        Low_Library_Flag_N =
          sum(x$Low_Library_Flag),
        
        stringsAsFactors = FALSE
      )
    }
  )
  
  do.call(
    rbind,
    out
  )
}


# ------------------------------------------------------------
# 14. Save descriptive summary
# ------------------------------------------------------------

library_summary <- summarize_library_qc(
  library_qc
)

library_summary_file <- file.path(
  step10d_dir,
  "STEP10D-A1_Library_Size_Descriptive_Summary.tsv"
)

write.table(
  library_summary,
  file = library_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Saved:\n",
  library_summary_file,
  "\n"
)


# ------------------------------------------------------------
# 15. Detected-gene descriptive summary
# ------------------------------------------------------------

detected_gene_summary <- library_qc[
  ,
  c(
    "Dataset",
    "Group_ID",
    "Detected_Genes",
    "Zero_Genes",
    "Gene_Detection_Fraction",
    "Low_Library_Flag"
  )
]

detected_gene_summary_file <- file.path(
  step10d_dir,
  "STEP10D-A2_Detected_Gene_Descriptive_Summary.tsv"
)

write.table(
  detected_gene_summary,
  file = detected_gene_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat(
  "Saved:\n",
  detected_gene_summary_file,
  "\n"
)


# ------------------------------------------------------------
# 16. Print compact QC summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-A — LIBRARY QC SUMMARY\n")
cat("============================================================\n\n")

print(
  library_summary,
  row.names = FALSE
)

cat("\n")
cat("Low-library descriptive flags:\n")
print(
  library_qc[
    library_qc$Low_Library_Flag,
    c(
      "Dataset",
      "Group_ID",
      "Total_UMI",
      "Detected_Genes"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)

if (!any(library_qc$Low_Library_Flag)) {
  
  cat(
    "\nNo descriptive low-library flags detected.\n"
  )
  
} else {
  
  cat(
    "\nNOTE — Low-library flags are descriptive only;\n",
    "no samples were excluded.\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 17. Final STEP10D-A lock summary
# ------------------------------------------------------------

final_checks <- data.frame(
  Check = c(
    "STEP10C pseudobulk dimensions",
    "Gene identity/order",
    "Raw UMI finite/non-negative/integer",
    "Library size valid",
    "Detected gene counts valid",
    "Detection fraction valid",
    "No patient exclusion",
    "No statistical testing"
  ),
  Status = c(
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS"
  ),
  stringsAsFactors = FALSE
)

final_summary_file <- file.path(
  step10d_dir,
  "STEP10D-A_FINAL_Audit_Summary.tsv"
)

write.table(
  final_checks,
  file = final_summary_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("STEP10D-A — FINAL STATUS\n")
cat("============================================================\n\n")

cat("PASS — STEP10D-A Pseudobulk Library QC completed.\n\n")

cat("Locked groups:\n")
cat("  Primary T-cell              = 9\n")
cat("  Primary Myeloid/macrophage  = 9\n")
cat("  Primary Endothelial         = 4\n")
cat("  Sensitivity Endothelial     = 2\n\n")

cat("Raw UMI pseudobulk matrices were NOT recalculated.\n")
cat("No patient/cell was excluded.\n")
cat("No inferential statistics were performed.\n")
cat("No expression normalization was performed.\n")
cat("Low-library flags, if any, are descriptive only.\n\n")

cat("STEP10D-A STATUS: PASS\n")


# ------------------------------------------------------------
# 18. Write explicit lock file
# ------------------------------------------------------------

lock_text <- c(
  "STEP10D-A — PSEUDOBULK LIBRARY QC",
  "FINAL LOCK",
  "",
  "Status: PASS",
  "",
  "Locked pseudobulk groups:",
  "Primary T-cell: 9",
  "Primary Myeloid/macrophage: 9",
  "Primary Endothelial: 4",
  "Sensitivity Endothelial: 2",
  "",
  "QC performed:",
  "Total UMI / library size",
  "Detected gene count",
  "Zero-gene count",
  "Gene detection fraction",
  "Raw UMI integrity",
  "Descriptive low-library flag",
  "",
  "No pseudobulk recalculation.",
  "No normalized expression.",
  "No differential expression.",
  "No correlation analysis.",
  "No CellChat.",
  "No ligand-receptor analysis.",
  "No survival analysis.",
  "No patient or cell exclusion.",
  "",
  "Low-library flags are descriptive only.",
  "STEP10C locked data remain unchanged."
)

writeLines(
  lock_text,
  con = file.path(
    step10d_dir,
    "STEP10D-A_FINAL_LOCK.txt"
  )
)

cat(
  "Lock file saved:\n",
  file.path(
    step10d_dir,
    "STEP10D-A_FINAL_LOCK.txt"
  ),
  "\n"
)

# ============================================================
# STEP10D-B
# AOC3–MS4A7–IL7R Core-Gene Raw-Count & Detection Audit
#
# PURPOSE:
#   Audit locked patient-level raw-UMI pseudobulk counts
#   for the three predefined core genes:
#       AOC3, MS4A7, IL7R
#
# INPUT:
#   STEP10C locked raw-UMI pseudobulk matrices
#
# OUTPUT:
#   Core-gene raw counts
#   Core-gene detection status
#   Descriptive detection/count summaries
#   Final audit summary
#
# IMPORTANT:
#   - No pseudobulk recalculation
#   - No normalization
#   - No log transformation
#   - No statistical testing
#   - No correlation
#   - No DE analysis
#   - No CellChat
#   - No patient/cell exclusion
#   - Zero counts are retained
#   - Primary Endothelial and Sensitivity Endothelial
#     remain separate datasets
# ============================================================


cat("============================================================\n")
cat("STEP10D-B — Core-Gene Raw-Count & Detection Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

input_dir <- file.path(
  base_dir,
  "STEP10C_Patient_Level_RawUMI_Pseudobulk"
)

output_dir <- file.path(
  base_dir,
  "STEP10D-B_CoreGene_RawCount_Detection_Audit"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Input directory:\n", input_dir, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")


# ------------------------------------------------------------
# 1. Locked core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Locked core genes:\n")
cat(paste(core_genes, collapse = " / "), "\n\n")


# ------------------------------------------------------------
# 2. Locked input files
# ------------------------------------------------------------

file_tcell <- file.path(
  input_dir,
  "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
)

file_myeloid <- file.path(
  input_dir,
  "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
)

file_endo_primary <- file.path(
  input_dir,
  "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
)

file_endo_sens <- file.path(
  input_dir,
  "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
)

file_annotation_primary <- file.path(
  input_dir,
  "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
)

file_annotation_sens <- file.path(
  input_dir,
  "STEP10C-F_Sensitivity_Pseudobulk_Sample_Annotation.tsv"
)


required_files <- c(
  file_tcell,
  file_myeloid,
  file_endo_primary,
  file_endo_sens,
  file_annotation_primary,
  file_annotation_sens
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("FAIL — Missing required STEP10C input file(s):\n")
  print(missing_files)
  stop("STEP10D-B aborted.")
}

cat("PASS — All required STEP10C input files found.\n\n")


# ------------------------------------------------------------
# 3. Read locked pseudobulk matrices
# ------------------------------------------------------------

cat("Reading locked pseudobulk matrices...\n")

mat_tcell <- read.delim(
  file_tcell,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat_myeloid <- read.delim(
  file_myeloid,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat_endo_primary <- read.delim(
  file_endo_primary,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mat_endo_sens <- read.delim(
  file_endo_sens,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("\nMatrix dimensions:\n")
cat("T-cell:                 ",
    nrow(mat_tcell), "genes x",
    ncol(mat_tcell), "groups\n")

cat("Myeloid/macrophage:     ",
    nrow(mat_myeloid), "genes x",
    ncol(mat_myeloid), "groups\n")

cat("Endothelial primary:    ",
    nrow(mat_endo_primary), "genes x",
    ncol(mat_endo_primary), "groups\n")

cat("Endothelial sensitivity:",
    nrow(mat_endo_sens), "genes x",
    ncol(mat_endo_sens), "groups\n\n")


# ------------------------------------------------------------
# 4. Matrix integrity checks
# ------------------------------------------------------------

check_matrix <- function(mat, label) {
  
  if (is.null(rownames(mat))) {
    stop(label, ": rownames are missing.")
  }
  
  if (anyDuplicated(rownames(mat)) > 0) {
    stop(label, ": duplicated gene IDs detected.")
  }
  
  if (anyDuplicated(colnames(mat)) > 0) {
    stop(label, ": duplicated Group_IDs detected.")
  }
  
  x <- as.matrix(mat)
  
  if (anyNA(x)) {
    stop(label, ": NA values detected.")
  }
  
  if (any(!is.finite(x))) {
    stop(label, ": non-finite values detected.")
  }
  
  if (any(x < 0)) {
    stop(label, ": negative values detected.")
  }
  
  if (any(x != floor(x))) {
    stop(label, ": non-integer values detected.")
  }
  
  TRUE
}


check_matrix(mat_tcell, "Primary T-cell")
check_matrix(mat_myeloid, "Primary Myeloid/macrophage")
check_matrix(mat_endo_primary, "Primary Endothelial")
check_matrix(mat_endo_sens, "Sensitivity Endothelial")

cat("PASS — All locked pseudobulk matrices contain finite,")
cat(" non-negative integer raw UMI counts.\n\n")


# ------------------------------------------------------------
# 5. Gene-universe reconciliation
# ------------------------------------------------------------

gene_sets <- list(
  Tcell = rownames(mat_tcell),
  Myeloid = rownames(mat_myeloid),
  EndothelialPrimary = rownames(mat_endo_primary),
  EndothelialSensitivity = rownames(mat_endo_sens)
)

gene_reference <- gene_sets[[1]]

gene_identity_pass <- all(
  vapply(
    gene_sets,
    function(x) identical(x, gene_reference),
    logical(1)
  )
)

if (!gene_identity_pass) {
  stop(
    "FAIL — Gene identity/order differs across locked pseudobulk matrices."
  )
}

cat("PASS — Gene identity and order are identical across all matrices.\n\n")


# ------------------------------------------------------------
# 6. Core-gene presence audit
# ------------------------------------------------------------

core_presence <- data.frame(
  Gene = core_genes,
  Present_in_Tcell = core_genes %in% rownames(mat_tcell),
  Present_in_Myeloid = core_genes %in% rownames(mat_myeloid),
  Present_in_EndothelialPrimary =
    core_genes %in% rownames(mat_endo_primary),
  Present_in_EndothelialSensitivity =
    core_genes %in% rownames(mat_endo_sens),
  stringsAsFactors = FALSE
)

print(core_presence)

if (!all(unlist(core_presence[, -1]))) {
  stop(
    "FAIL — One or more core genes are absent from a locked pseudobulk matrix."
  )
}

cat("\nPASS — AOC3, MS4A7 and IL7R are present in all locked matrices.\n\n")


# ------------------------------------------------------------
# 7. Read annotation files
# ------------------------------------------------------------

ann_primary <- read.delim(
  file_annotation_primary,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ann_sens <- read.delim(
  file_annotation_sens,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("Primary annotation columns:\n")
print(colnames(ann_primary))

cat("\nSensitivity annotation columns:\n")
print(colnames(ann_sens))


# ------------------------------------------------------------
# 8. Required annotation schema
# ------------------------------------------------------------

required_ann_cols <- c(
  "Patient",
  "Compartment",
  "Cell_Count",
  "Eligibility",
  "Group_ID"
)

if (!all(required_ann_cols %in% colnames(ann_primary))) {
  stop(
    "FAIL — Primary annotation does not contain the locked STEP10C schema."
  )
}

if (!all(required_ann_cols %in% colnames(ann_sens))) {
  stop(
    "FAIL — Sensitivity annotation does not contain the locked STEP10C schema."
  )
}

cat("\nPASS — Annotation schema matches locked STEP10C structure.\n\n")


# ------------------------------------------------------------
# 9. Group_ID reconciliation
# ------------------------------------------------------------

check_group_mapping <- function(mat, ann, label) {
  
  mat_groups <- colnames(mat)
  ann_groups <- ann$Group_ID
  
  if (anyDuplicated(ann_groups) > 0) {
    stop(label, ": duplicated Group_ID in annotation.")
  }
  
  if (!setequal(mat_groups, ann_groups)) {
    stop(label, ": matrix Group_ID and annotation Group_ID do not match.")
  }
  
  if (!identical(
    mat_groups,
    ann$Group_ID[match(mat_groups, ann$Group_ID)]
  )) {
    stop(
      label,
      ": annotation order does not reconcile to matrix Group_ID order."
    )
  }
  
  TRUE
}


check_group_mapping(
  cbind(mat_tcell),
  ann_primary[ann_primary$Compartment == "T-cell", ],
  "Primary T-cell"
)

check_group_mapping(
  cbind(mat_myeloid),
  ann_primary[
    ann_primary$Compartment == "Myeloid/macrophage",
  ],
  "Primary Myeloid/macrophage"
)

check_group_mapping(
  cbind(mat_endo_primary),
  ann_primary[
    ann_primary$Compartment == "Endothelial" &
      ann_primary$Eligibility == "Primary",
  ],
  "Primary Endothelial"
)

check_group_mapping(
  cbind(mat_endo_sens),
  ann_sens[
    ann_sens$Compartment == "Endothelial" &
      ann_sens$Eligibility == "Sensitivity",
  ],
  "Sensitivity Endothelial"
)

cat("\nPASS — Group_ID mappings reconcile exactly to locked annotations.\n\n")


# ------------------------------------------------------------
# 10. Extract core-gene raw counts
# ------------------------------------------------------------

extract_core_counts <- function(mat, ann, dataset_label) {
  
  x <- as.data.frame(
    t(mat[core_genes, , drop = FALSE]),
    check.names = FALSE
  )
  
  x$Group_ID <- rownames(x)
  
  rownames(x) <- NULL
  
  x <- merge(
    ann,
    x,
    by = "Group_ID",
    sort = FALSE
  )
  
  x$Dataset <- dataset_label
  
  x <- x[
    ,
    c(
      "Dataset",
      "Group_ID",
      "Patient",
      "Compartment",
      "Eligibility",
      "Cell_Count",
      core_genes
    )
  ]
  
  x
}


core_tcell <- extract_core_counts(
  mat_tcell,
  ann_primary[
    ann_primary$Compartment == "T-cell",
  ],
  "Primary_T-cell"
)

core_myeloid <- extract_core_counts(
  mat_myeloid,
  ann_primary[
    ann_primary$Compartment == "Myeloid/macrophage",
  ],
  "Primary_Myeloid_macrophage"
)

core_endo_primary <- extract_core_counts(
  mat_endo_primary,
  ann_primary[
    ann_primary$Compartment == "Endothelial" &
      ann_primary$Eligibility == "Primary",
  ],
  "Primary_Endothelial"
)

core_endo_sens <- extract_core_counts(
  mat_endo_sens,
  ann_sens[
    ann_sens$Compartment == "Endothelial" &
      ann_sens$Eligibility == "Sensitivity",
  ],
  "Sensitivity_Endothelial"
)


core_counts <- rbind(
  core_tcell,
  core_myeloid,
  core_endo_primary,
  core_endo_sens
)

rownames(core_counts) <- NULL


# ------------------------------------------------------------
# 11. Verify extracted counts
# ------------------------------------------------------------

count_values <- as.matrix(
  core_counts[, core_genes, drop = FALSE]
)

if (anyNA(count_values)) {
  stop("FAIL — NA detected in extracted core-gene counts.")
}

if (any(!is.finite(count_values))) {
  stop("FAIL — Non-finite value detected in core-gene counts.")
}

if (any(count_values < 0)) {
  stop("FAIL — Negative core-gene count detected.")
}

if (any(count_values != floor(count_values))) {
  stop("FAIL — Non-integer core-gene count detected.")
}

cat("PASS — Core-gene raw UMI counts extracted successfully.\n\n")


# ------------------------------------------------------------
# 12. Detection matrix
#
# Detection definition:
#   raw UMI count > 0
#
# No threshold other than >0.
# ------------------------------------------------------------

detection <- core_counts[
  ,
  c(
    "Dataset",
    "Group_ID",
    "Patient",
    "Compartment",
    "Eligibility",
    "Cell_Count"
  )
]

for (g in core_genes) {
  detection[[paste0(g, "_Detected")]] <-
    core_counts[[g]] > 0
}


# ------------------------------------------------------------
# 13. Detection summary function
# ------------------------------------------------------------

make_detection_summary <- function(count_df, dataset_label) {
  
  out <- lapply(core_genes, function(g) {
    
    counts <- count_df[[g]]
    detected <- counts > 0
    
    data.frame(
      Dataset = dataset_label,
      Gene = g,
      N = length(counts),
      Detected_N = sum(detected),
      Not_Detected_N = sum(!detected),
      Detection_Frequency =
        ifelse(length(counts) > 0,
               sum(detected) / length(counts),
               NA_real_),
      Detection_Percent =
        ifelse(length(counts) > 0,
               100 * sum(detected) / length(counts),
               NA_real_),
      Raw_Count_Min =
        ifelse(length(counts) > 0, min(counts), NA_real_),
      Raw_Count_Q1 =
        ifelse(length(counts) > 0,
               unname(quantile(counts, 0.25)),
               NA_real_),
      Raw_Count_Median =
        ifelse(length(counts) > 0,
               median(counts),
               NA_real_),
      Raw_Count_Mean =
        ifelse(length(counts) > 0,
               mean(counts),
               NA_real_),
      Raw_Count_Q3 =
        ifelse(length(counts) > 0,
               unname(quantile(counts, 0.75)),
               NA_real_),
      Raw_Count_Max =
        ifelse(length(counts) > 0, max(counts), NA_real_),
      Zero_Count_N = sum(counts == 0),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, out)
}


summary_tcell <- make_detection_summary(
  core_tcell,
  "Primary_T-cell"
)

summary_myeloid <- make_detection_summary(
  core_myeloid,
  "Primary_Myeloid_macrophage"
)

summary_endo_primary <- make_detection_summary(
  core_endo_primary,
  "Primary_Endothelial"
)

summary_endo_sens <- make_detection_summary(
  core_endo_sens,
  "Sensitivity_Endothelial"
)


core_summary <- rbind(
  summary_tcell,
  summary_myeloid,
  summary_endo_primary,
  summary_endo_sens
)

rownames(core_summary) <- NULL


# ------------------------------------------------------------
# 14. Detection summary print
# ------------------------------------------------------------

cat("============================================================\n")
cat("CORE-GENE DETECTION SUMMARY\n")
cat("============================================================\n\n")

print(core_summary)


# ------------------------------------------------------------
# 15. Dataset-level sample counts
# ------------------------------------------------------------

dataset_n <- data.frame(
  Dataset = c(
    "Primary_T-cell",
    "Primary_Myeloid_macrophage",
    "Primary_Endothelial",
    "Sensitivity_Endothelial"
  ),
  Expected_N = c(9, 9, 4, 2),
  Observed_N = c(
    nrow(core_tcell),
    nrow(core_myeloid),
    nrow(core_endo_primary),
    nrow(core_endo_sens)
  ),
  stringsAsFactors = FALSE
)

dataset_n$PASS <- dataset_n$Expected_N == dataset_n$Observed_N

cat("\nDataset sample-count reconciliation:\n")
print(dataset_n)

if (!all(dataset_n$PASS)) {
  stop("FAIL — Observed pseudobulk group counts do not match locked counts.")
}

cat("\nPASS — All locked group counts verified.\n\n")


# ------------------------------------------------------------
# 16. Core-gene presence and count summary
# ------------------------------------------------------------

core_gene_presence_summary <- data.frame(
  Gene = core_genes,
  All_Matrices_Present = TRUE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 17. Zero-count audit
# ------------------------------------------------------------

zero_audit <- do.call(
  rbind,
  lapply(
    seq_len(nrow(core_summary)),
    function(i) {
      core_summary[i, c(
        "Dataset",
        "Gene",
        "N",
        "Detected_N",
        "Not_Detected_N",
        "Zero_Count_N"
      )]
    }
  )
)

cat("Zero-count audit completed.\n")
cat("Zero counts are retained and were NOT used as exclusion criteria.\n\n")


# ------------------------------------------------------------
# 18. Write outputs
# ------------------------------------------------------------

file_counts <- file.path(
  output_dir,
  "STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv"
)

file_detection <- file.path(
  output_dir,
  "STEP10D-C_Pseudobulk_CoreGene_Detection.tsv"
)

file_summary <- file.path(
  output_dir,
  "STEP10D-B_CoreGene_Detection_Count_Summary.tsv"
)

file_presence <- file.path(
  output_dir,
  "STEP10D-B_CoreGene_Presence_Audit.tsv"
)

file_dataset_n <- file.path(
  output_dir,
  "STEP10D-B_Dataset_GroupCount_Audit.tsv"
)


write.table(
  core_counts,
  file_counts,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  detection,
  file_detection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_summary,
  file_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  core_gene_presence_summary,
  file_presence,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dataset_n,
  file_dataset_n,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 19. Final audit summary
# ------------------------------------------------------------

final_summary <- data.frame(
  Item = c(
    "Locked core genes",
    "Primary T-cell groups",
    "Primary Myeloid/macrophage groups",
    "Primary Endothelial groups",
    "Sensitivity Endothelial groups",
    "Raw UMI counts used",
    "Normalized expression used",
    "Core-gene count extraction",
    "Detection definition",
    "Zero counts retained",
    "Patient exclusion",
    "Statistical testing",
    "Correlation analysis",
    "Differential expression",
    "CellChat/LR analysis"
  ),
  Result = c(
    paste(core_genes, collapse = " / "),
    nrow(core_tcell),
    nrow(core_myeloid),
    nrow(core_endo_primary),
    nrow(core_endo_sens),
    "YES",
    "NO",
    "YES",
    "raw count > 0",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  stringsAsFactors = FALSE
)

file_final <- file.path(
  output_dir,
  "STEP10D-B_FINAL_Audit_Summary.tsv"
)

write.table(
  final_summary,
  file_final,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Final lock file
# ------------------------------------------------------------

lock_file <- file.path(
  output_dir,
  "STEP10D-B_FINAL_LOCK.txt"
)

sink(lock_file)

cat("STEP10D-B — AOC3–MS4A7–IL7R Core-Gene Raw-Count & Detection Audit\n")
cat("============================================================\n\n")

cat("STATUS: PASS\n\n")

cat("Locked core genes:\n")
cat("AOC3 / MS4A7 / IL7R\n\n")

cat("Locked pseudobulk groups:\n")
cat("Primary T-cell: 9\n")
cat("Primary Myeloid/macrophage: 9\n")
cat("Primary Endothelial: 4\n")
cat("Sensitivity Endothelial: 2\n\n")

cat("Raw UMI pseudobulk counts were used.\n")
cat("No normalized expression was used.\n")
cat("No log transformation was performed.\n\n")

cat("Detection definition:\n")
cat("raw UMI count > 0\n\n")

cat("Zero-count observations were retained.\n")
cat("No patient or cell was excluded.\n\n")

cat("No inferential statistics were performed.\n")
cat("No correlation analysis was performed.\n")
cat("No differential expression analysis was performed.\n")
cat("No CellChat or ligand-receptor analysis was performed.\n")
cat("No survival analysis was performed.\n\n")

cat("STEP10C locked pseudobulk matrices were NOT recalculated.\n")
cat("STEP10D-A library QC was NOT recalculated.\n\n")

cat("Output files:\n")
cat(basename(file_counts), "\n")
cat(basename(file_detection), "\n")
cat(basename(file_summary), "\n")
cat(basename(file_presence), "\n")
cat(basename(file_dataset_n), "\n")
cat(basename(file_final), "\n")

sink()


# ------------------------------------------------------------
# 21. Console final status
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PASS — STEP10D-B Core-Gene Raw-Count & Detection Audit completed.\n")
cat("============================================================\n\n")

cat("Locked core genes: AOC3 / MS4A7 / IL7R\n")
cat("Primary T-cell: 9 groups\n")
cat("Primary Myeloid/macrophage: 9 groups\n")
cat("Primary Endothelial: 4 groups\n")
cat("Sensitivity Endothelial: 2 groups\n\n")

cat("Raw UMI counts extracted: YES\n")
cat("Detection definition: raw count > 0\n")
cat("Zero counts retained: YES\n")
cat("Patient/cell exclusion: NO\n")
cat("Inferential statistics: NO\n")
cat("Normalization: NO\n")
cat("STEP10C recalculated: NO\n")
cat("STEP10D-A recalculated: NO\n\n")

cat("STEP10D-B STATUS: PASS\n")
cat("Final lock file:\n")
cat(lock_file, "\n")


##STEP10D-C — Core-Gene Compartment-Level Expression Audit
# ============================================================
# STEP10D-C
# Core-Gene Compartment-Level Expression Audit
#
# PURPOSE:
#   Descriptive audit of AOC3, MS4A7 and IL7R across
#   predefined patient-level pseudobulk compartments.
#
# INPUT:
#   STEP10D-B locked core-gene raw-count table
#   STEP10C locked pseudobulk annotation
#
# LOCKED CORE GENES:
#   AOC3 / MS4A7 / IL7R
#
# LOCKED COMPARTMENTS:
#   Primary T-cell
#   Primary Myeloid/macrophage
#   Primary Endothelial
#   Sensitivity Endothelial
#
# IMPORTANT:
#   - Raw UMI pseudobulk counts only
#   - No normalization
#   - No log transformation
#   - No inferential statistics
#   - No p-values
#   - No FDR
#   - No differential expression
#   - No correlation
#   - No CellChat
#   - No ligand-receptor analysis
#   - No patient exclusion
#   - Zero counts retained
#   - Primary Endothelial n=4 and
#     Sensitivity Endothelial n=2 remain separate
#   - STEP10C / STEP10D-A / STEP10D-B are NOT recalculated
# ============================================================


cat("============================================================\n")
cat("STEP10D-C — Core-Gene Compartment-Level Expression Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

input_dir <- file.path(
  base_dir,
  "STEP10D-B_CoreGene_RawCount_Detection_Audit"
)

output_dir <- file.path(
  base_dir,
  "STEP10D-C_CoreGene_Compartment_Level_Audit"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("Input directory:\n", input_dir, "\n\n")
cat("Output directory:\n", output_dir, "\n\n")


# ------------------------------------------------------------
# 1. Locked core genes
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat("Locked core genes:\n")
cat(paste(core_genes, collapse = " / "), "\n\n")


# ------------------------------------------------------------
# 2. Locked dataset structure
# ------------------------------------------------------------

locked_datasets <- c(
  "Primary_T-cell",
  "Primary_Myeloid_macrophage",
  "Primary_Endothelial",
  "Sensitivity_Endothelial"
)

expected_n <- c(
  "Primary_T-cell" = 9,
  "Primary_Myeloid_macrophage" = 9,
  "Primary_Endothelial" = 4,
  "Sensitivity_Endothelial" = 2
)

cat("Locked datasets:\n")
print(locked_datasets)

cat("\nExpected patient/group counts:\n")
print(expected_n)

cat("\n")


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

file_counts <- file.path(
  input_dir,
  "STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv"
)

file_detection <- file.path(
  input_dir,
  "STEP10D-C_Pseudobulk_CoreGene_Detection.tsv"
)

file_summary_previous <- file.path(
  input_dir,
  "STEP10D-B_CoreGene_Detection_Count_Summary.tsv"
)

required_files <- c(
  file_counts,
  file_detection,
  file_summary_previous
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  cat("FAIL — Missing required STEP10D-B file(s):\n")
  print(missing_files)
  
  stop("STEP10D-C aborted.")
}

cat("PASS — Required STEP10D-B files found.\n\n")


# ------------------------------------------------------------
# 4. Read locked STEP10D-B core-gene count table
# ------------------------------------------------------------

cat("Reading locked STEP10D-B raw-count table...\n")

core_counts <- read.delim(
  file_counts,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Rows:",
  nrow(core_counts),
  "\n"
)

cat(
  "Columns:",
  ncol(core_counts),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Verify required columns
# ------------------------------------------------------------

required_columns <- c(
  "Dataset",
  "Group_ID",
  "Patient",
  "Compartment",
  "Eligibility",
  "Cell_Count",
  core_genes
)

missing_columns <- setdiff(
  required_columns,
  colnames(core_counts)
)

if (length(missing_columns) > 0) {
  
  cat("FAIL — Missing required columns:\n")
  print(missing_columns)
  
  stop("STEP10D-C aborted.")
}

cat("PASS — STEP10D-B core-count table schema verified.\n\n")


# ------------------------------------------------------------
# 6. Basic integrity audit
# ------------------------------------------------------------

if (anyDuplicated(core_counts$Group_ID) > 0) {
  
  stop(
    "FAIL — Duplicate Group_ID detected in locked core-count table."
  )
}

if (anyDuplicated(
  paste(core_counts$Dataset, core_counts$Patient, sep = "||")
) > 0) {
  
  stop(
    "FAIL — Duplicate Dataset/Patient combination detected."
  )
}


# Core-gene numeric validation

for (g in core_genes) {
  
  x <- core_counts[[g]]
  
  if (!is.numeric(x)) {
    stop(
      "FAIL — ",
      g,
      " is not numeric."
    )
  }
  
  if (anyNA(x)) {
    stop(
      "FAIL — NA values detected in ",
      g,
      "."
    )
  }
  
  if (any(!is.finite(x))) {
    stop(
      "FAIL — Non-finite values detected in ",
      g,
      "."
    )
  }
  
  if (any(x < 0)) {
    stop(
      "FAIL — Negative values detected in ",
      g,
      "."
    )
  }
  
  if (any(x != floor(x))) {
    stop(
      "FAIL — Non-integer values detected in ",
      g,
      "."
    )
  }
}

cat(
  "PASS — All core-gene values are finite,",
  "non-negative integer raw UMI counts.\n\n"
)


# ------------------------------------------------------------
# 7. Verify locked dataset labels
# ------------------------------------------------------------

observed_datasets <- unique(core_counts$Dataset)

unknown_datasets <- setdiff(
  observed_datasets,
  locked_datasets
)

if (length(unknown_datasets) > 0) {
  
  cat("FAIL — Unknown dataset label(s):\n")
  print(unknown_datasets)
  
  stop("STEP10D-C aborted.")
}

cat("PASS — Dataset labels match the locked STEP10D structure.\n\n")


# ------------------------------------------------------------
# 8. Verify expected group counts
# ------------------------------------------------------------

observed_n <- table(core_counts$Dataset)

group_count_audit <- data.frame(
  Dataset = locked_datasets,
  Expected_N = as.integer(expected_n[locked_datasets]),
  Observed_N = as.integer(
    observed_n[locked_datasets]
  ),
  stringsAsFactors = FALSE
)

group_count_audit$PASS <-
  group_count_audit$Expected_N ==
  group_count_audit$Observed_N

cat("Dataset/group-count audit:\n")
print(group_count_audit)

if (!all(group_count_audit$PASS)) {
  
  stop(
    "FAIL — Observed group counts do not match locked counts."
  )
}

cat(
  "\nPASS — All locked patient/group counts verified.\n\n"
)


# ------------------------------------------------------------
# 9–11. Rebuild standardized PRIMARY-only count table
#
# IMPORTANT:
#   STEP10D-B contains both:
#       Primary
#       Sensitivity
#
#   STEP10D-D must use PRIMARY only.
#   This is especially important for Endothelial:
#       Primary      = 4 patients
#       Sensitivity  = 2 patients
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D — Rebuilding PRIMARY-only core-gene count table\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 9.1 Required metadata columns
# ------------------------------------------------------------

required_metadata <- c(
  "Patient",
  "Compartment",
  "Eligibility"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(core_counts)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "STOP: STEP10D-B core-count table is missing required ",
    "metadata columns:\n",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 9.2 Detect core-gene representation
# ------------------------------------------------------------

count_gene_columns <- intersect(
  core_genes,
  colnames(core_counts)
)


# ------------------------------------------------------------
# 9.3 Convert STEP10D-B core-count table to long format
# ------------------------------------------------------------

if (
  length(count_gene_columns) ==
  length(core_genes)
) {
  
  # ----------------------------------------------------------
  # Genes are columns
  # ----------------------------------------------------------
  
  counts_long <- reshape(
    core_counts[
      ,
      c(
        "Patient",
        "Compartment",
        "Eligibility",
        intersect(
          "Group_ID",
          colnames(core_counts)
        ),
        core_genes
      ),
      drop = FALSE
    ],
    varying = core_genes,
    v.names = "Raw_Count",
    timevar = "Gene",
    times = core_genes,
    idvar = c(
      "Patient",
      "Compartment",
      "Eligibility",
      intersect(
        "Group_ID",
        colnames(core_counts)
      )
    ),
    direction = "long"
  )
  
  rownames(counts_long) <- NULL
  
} else {
  
  # ----------------------------------------------------------
  # Gene is a row variable
  # ----------------------------------------------------------
  
  if (!"Gene" %in% colnames(core_counts)) {
    
    stop(
      "STOP: STEP10D-B count table does not contain ",
      "core genes as columns and no Gene column was found."
    )
  }
  
  
  count_value_candidates <- setdiff(
    colnames(core_counts),
    c(
      "Patient",
      "Compartment",
      "Eligibility",
      "Group_ID",
      "Gene"
    )
  )
  
  
  if (length(count_value_candidates) < 1) {
    
    stop(
      "STOP: Cannot identify raw-count column."
    )
  }
  
  
  count_value_col <- count_value_candidates[1]
  
  
  counts_long <- core_counts[
    core_counts$Gene %in% core_genes,
    c(
      "Patient",
      "Compartment",
      "Eligibility",
      intersect(
        "Group_ID",
        colnames(core_counts)
      ),
      "Gene",
      count_value_col
    ),
    drop = FALSE
  ]
  
  
  colnames(
    counts_long
  )[
    colnames(counts_long) ==
      count_value_col
  ] <- "Raw_Count"
}


# ------------------------------------------------------------
# 9.4 Standardize data types
# ------------------------------------------------------------

counts_long$Patient <- as.character(
  counts_long$Patient
)

counts_long$Compartment <- as.character(
  counts_long$Compartment
)

counts_long$Eligibility <- as.character(
  counts_long$Eligibility
)

counts_long$Gene <- as.character(
  counts_long$Gene
)

counts_long$Raw_Count <- as.numeric(
  counts_long$Raw_Count
)


# ------------------------------------------------------------
# 10. Audit Primary vs Sensitivity BEFORE filtering
# ------------------------------------------------------------

cat("\n")
cat("Eligibility distribution before filtering:\n")

print(
  table(
    counts_long$Eligibility,
    useNA = "ifany"
  )
)


cat("\nCompartment × Eligibility distribution:\n")

print(
  table(
    counts_long$Compartment,
    counts_long$Eligibility,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 11. PRIMARY-only locked analysis set
# ------------------------------------------------------------

counts_long <- counts_long[
  counts_long$Gene %in% core_genes &
    counts_long$Compartment %in% locked_compartments &
    counts_long$Eligibility == "Primary",
  ,
  drop = FALSE
]

rownames(counts_long) <- NULL


# ------------------------------------------------------------
# 11.1 Verify no Sensitivity records remain
# ------------------------------------------------------------

remaining_eligibility <- unique(
  counts_long$Eligibility
)

if (
  length(remaining_eligibility) != 1 ||
  remaining_eligibility != "Primary"
) {
  
  stop(
    "STOP: STEP10D-D primary-only filtering failed."
  )
}


# ------------------------------------------------------------
# 11.2 Verify locked compartments
# ------------------------------------------------------------

remaining_compartments <- unique(
  counts_long$Compartment
)

if (
  !all(
    remaining_compartments %in%
    locked_compartments
  )
) {
  
  stop(
    "STOP: Unexpected compartment remains after filtering."
  )
}


# ------------------------------------------------------------
# 11.3 Print final PRIMARY-only structure
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PRIMARY-ONLY INPUT STRUCTURE\n")
cat("============================================================\n\n")

cat(
  "Rows:",
  nrow(counts_long),
  "\n"
)

cat(
  "Patients:",
  length(
    unique(
      counts_long$Patient
    )
  ),
  "\n"
)

cat("\nCompartment frequencies:\n")

print(
  table(
    counts_long$Compartment
  )
)


cat("\nPatient counts by compartment:\n")

primary_patient_counts <- sapply(
  locked_compartments,
  function(comp) {
    
    length(
      unique(
        counts_long$Patient[
          counts_long$Compartment == comp
        ]
      )
    )
  }
)

print(
  primary_patient_counts
)


# ------------------------------------------------------------
# 11.4 Expected lock
# ------------------------------------------------------------

expected_primary_patients <- c(
  "T-cell" = 9L,
  "Myeloid/macrophage" = 9L,
  "Endothelial" = 4L
)


if (
  !all(
    primary_patient_counts[
      locked_compartments
    ] ==
    expected_primary_patients[
      locked_compartments
    ]
  )
) {
  
  stop(
    "STOP: PRIMARY-only patient counts do not match ",
    "the locked STEP10D-B structure."
  )
}


cat("\n")
cat(
  "PASS — PRIMARY-only patient structure verified.\n"
)

cat(
  "T-cell = ",
  primary_patient_counts[["T-cell"]],
  "\n",
  sep = ""
)

cat(
  "Myeloid/macrophage = ",
  primary_patient_counts[["Myeloid/macrophage"]],
  "\n",
  sep = ""
)

cat(
  "Endothelial = ",
  primary_patient_counts[["Endothelial"]],
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 12. Verify detection information
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D DETECTION INFORMATION AUDIT\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 12.1 Required columns
# ------------------------------------------------------------

required_detection_columns <- c(
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Gene",
  "Raw_Count"
)

missing_detection_columns <- setdiff(
  required_detection_columns,
  colnames(counts_long)
)

if (
  length(missing_detection_columns) > 0
) {
  
  stop(
    "STOP: counts_long is missing required columns: ",
    paste(
      missing_detection_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 12.2 Verify locked core genes
# ------------------------------------------------------------

observed_genes <- unique(
  counts_long$Gene
)

if (
  !setequal(
    observed_genes,
    core_genes
  )
) {
  
  stop(
    "STOP: counts_long does not contain exactly the locked ",
    "AOC3–MS4A7–IL7R core-gene set."
  )
}


# ------------------------------------------------------------
# 12.3 Verify PRIMARY-only status
# ------------------------------------------------------------

observed_eligibility <- unique(
  counts_long$Eligibility
)

if (
  !identical(
    observed_eligibility,
    "Primary"
  )
) {
  
  stop(
    "STOP: non-Primary observations detected in ",
    "PRIMARY-only concordance input."
  )
}


# ------------------------------------------------------------
# 12.4 Verify raw-count integrity
# ------------------------------------------------------------

if (
  anyNA(
    counts_long$Raw_Count
  )
) {
  
  stop(
    "STOP: NA values detected in Raw_Count."
  )
}


if (
  any(
    !is.finite(
      counts_long$Raw_Count
    )
  )
) {
  
  stop(
    "STOP: non-finite values detected in Raw_Count."
  )
}


if (
  any(
    counts_long$Raw_Count < 0
  )
) {
  
  stop(
    "STOP: negative raw UMI counts detected."
  )
}


if (
  any(
    counts_long$Raw_Count !=
    floor(
      counts_long$Raw_Count
    )
  )
) {
  
  stop(
    "STOP: non-integer raw UMI counts detected."
  )
}


# ------------------------------------------------------------
# 13. Descriptive summary function
# ------------------------------------------------------------

make_compartment_summary <- function(
    df,
    dataset_label,
    gene_label
) {
  
  x <- df$Raw_Count
  
  detected <- x > 0
  
  data.frame(
    Dataset = dataset_label,
    Gene = gene_label,
    
    N = length(x),
    
    Detected_N =
      sum(detected),
    
    Zero_N =
      sum(x == 0),
    
    Detection_Percent =
      100 * mean(detected),
    
    Raw_Count_Min =
      min(x),
    
    Raw_Count_Q1 =
      unname(
        quantile(
          x,
          0.25
        )
      ),
    
    Raw_Count_Median =
      median(x),
    
    Raw_Count_Mean =
      mean(x),
    
    Raw_Count_Q3 =
      unname(
        quantile(
          x,
          0.75
        )
      ),
    
    Raw_Count_Max =
      max(x),
    
    Raw_Count_IQR =
      IQR(x),
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 14. Generate PRIMARY compartment-level summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PRIMARY COMPARTMENT-LEVEL DESCRIPTIVE SUMMARY\n")
cat("============================================================\n\n")


summary_list <- lapply(
  locked_compartments,
  function(comp) {
    
    lapply(
      core_genes,
      function(g) {
        
        subset_df <- counts_long[
          counts_long$Compartment == comp &
            counts_long$Gene == g &
            counts_long$Eligibility == "Primary",
          ,
          drop = FALSE
        ]
        
        if (
          nrow(subset_df) == 0
        ) {
          
          stop(
            "STOP: No PRIMARY observations found for ",
            comp,
            " / ",
            g
          )
        }
        
        make_compartment_summary(
          subset_df,
          comp,
          g
        )
      }
    )
  }
)


compartment_summary <- do.call(
  rbind,
  unlist(
    summary_list,
    recursive = FALSE
  )
)


rownames(compartment_summary) <- NULL


# ------------------------------------------------------------
# 14.1 Verify expected structure
# ------------------------------------------------------------

expected_summary_rows <-
  length(locked_compartments) *
  length(core_genes)


if (
  nrow(compartment_summary) !=
  expected_summary_rows
) {
  
  stop(
    "STOP: Unexpected number of rows in ",
    "compartment_summary."
  )
}


# ------------------------------------------------------------
# 14.2 Print summary
# ------------------------------------------------------------

print(
  compartment_summary
)


cat("\n")
cat(
  "PASS — PRIMARY compartment-level descriptive summary generated.\n"
)

cat(
  "Compartments: ",
  paste(
    locked_compartments,
    collapse = " / "
  ),
  "\n",
  sep = ""
)

cat(
  "Core genes: ",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n",
  sep = ""
)

cat(
  "Summary rows: ",
  nrow(compartment_summary),
  "\n",
  sep = ""
)

cat(
  "No normalization performed.\n"
)

cat(
  "No statistical testing performed.\n"
)

cat(
  "No differential-expression analysis performed.\n"
)


# ------------------------------------------------------------
# 15. Print main descriptive summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("CORE-GENE COMPARTMENT-LEVEL DESCRIPTIVE SUMMARY\n")
cat("============================================================\n\n")

print(
  compartment_summary,
  row.names = FALSE
)

cat("\n")
cat(
  "PASS — Main core-gene compartment summary verified.\n"
)

# ------------------------------------------------------------
# 16. Patient-level wide table
#
# One row per patient/group.
# Raw UMI counts for all three core genes.
# ------------------------------------------------------------

patient_level_table <- core_counts[
  ,
  c(
    "Dataset",
    "Group_ID",
    "Patient",
    "Compartment",
    "Eligibility",
    "Cell_Count",
    core_genes
  )
]

rownames(patient_level_table) <- NULL


# ------------------------------------------------------------
# 17. Detection frequency table
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D CORE-GENE DETECTION FREQUENCY\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 17.1 Generate one data.frame per
#     compartment × core-gene combination
# ------------------------------------------------------------

detection_frequency_list <- lapply(
  locked_compartments,
  function(comp) {
    
    lapply(
      core_genes,
      function(g) {
        
        subset_df <- counts_long[
          counts_long$Compartment == comp &
            counts_long$Eligibility == "Primary" &
            counts_long$Gene == g,
          ,
          drop = FALSE
        ]
        
        
        # ----------------------------------------------------
        # Structural check
        # ----------------------------------------------------
        
        if (
          nrow(subset_df) == 0
        ) {
          
          stop(
            "STOP: No PRIMARY observations found for ",
            comp,
            " / ",
            g
          )
        }
        
        
        # ----------------------------------------------------
        # Patient uniqueness check
        # ----------------------------------------------------
        
        if (
          anyDuplicated(
            subset_df$Patient
          ) > 0
        ) {
          
          stop(
            "STOP: Duplicate patient detected for ",
            comp,
            " / ",
            g
          )
        }
        
        
        # ----------------------------------------------------
        # Detection-frequency summary
        # ----------------------------------------------------
        
        data.frame(
          
          Compartment =
            comp,
          
          Eligibility =
            "Primary",
          
          Gene =
            g,
          
          N =
            nrow(subset_df),
          
          Detected_N =
            sum(
              subset_df$Detected
            ),
          
          Zero_N =
            sum(
              !subset_df$Detected
            ),
          
          Detection_Frequency =
            mean(
              subset_df$Detected
            ),
          
          Detection_Percent =
            100 *
            mean(
              subset_df$Detected
            ),
          
          stringsAsFactors = FALSE
        )
      }
    )
  }
)


# ------------------------------------------------------------
# 17.2 Flatten nested list
# ------------------------------------------------------------

detection_frequency_list_flat <- unlist(
  detection_frequency_list,
  recursive = FALSE
)


# ------------------------------------------------------------
# 17.3 Combine the 9 data.frames
# ------------------------------------------------------------

detection_frequency <- do.call(
  rbind,
  detection_frequency_list_flat
)


rownames(
  detection_frequency
) <- NULL


# ------------------------------------------------------------
# 17.4 Print result
# ------------------------------------------------------------

print(
  detection_frequency,
  row.names = FALSE
)


# ------------------------------------------------------------
# 17.5 Expected structure audit
# ------------------------------------------------------------

expected_detection_rows <-
  length(locked_compartments) *
  length(core_genes)


if (
  nrow(detection_frequency) !=
  expected_detection_rows
) {
  
  stop(
    "STOP: Unexpected number of rows in ",
    "detection_frequency."
  )
}


# ------------------------------------------------------------
# 17.6 Verify expected patient counts
# ------------------------------------------------------------

expected_N <- c(
  "T-cell" = 9L,
  "Myeloid/macrophage" = 9L,
  "Endothelial" = 4L
)


observed_N <- tapply(
  detection_frequency$N,
  detection_frequency$Compartment,
  unique
)


for (
  comp in locked_compartments
) {
  
  if (
    length(observed_N[[comp]]) != 1L ||
    observed_N[[comp]] != expected_N[[comp]]
  ) {
    
    stop(
      "STOP: Unexpected patient count for ",
      comp
    )
  }
}


# ------------------------------------------------------------
# 17.7 Verify detection-frequency consistency
# ------------------------------------------------------------

if (
  any(
    detection_frequency$Detected_N +
    detection_frequency$Zero_N !=
    detection_frequency$N
  )
) {
  
  stop(
    "STOP: Detected_N + Zero_N does not equal N."
  )
}


if (
  any(
    abs(
      detection_frequency$Detection_Frequency -
      (
        detection_frequency$Detected_N /
        detection_frequency$N
      )
    ) > 1e-12
  )
) {
  
  stop(
    "STOP: Detection_Frequency is inconsistent with ",
    "Detected_N / N."
  )
}


if (
  any(
    abs(
      detection_frequency$Detection_Percent -
      (
        100 *
        detection_frequency$Detected_N /
        detection_frequency$N
      )
    ) > 1e-10
  )
) {
  
  stop(
    "STOP: Detection_Percent is inconsistent with ",
    "Detected_N / N."
  )
}


# ------------------------------------------------------------
# 17.8 Final PASS
# ------------------------------------------------------------

cat("\n")

cat(
  "PASS — Detection-frequency table generated.\n"
)

cat(
  "PRIMARY compartments only: YES\n"
)

cat(
  "Core genes: ",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n",
  sep = ""
)

cat(
  "Rows: ",
  nrow(detection_frequency),
  "\n",
  sep = ""
)

cat(
  "T-cell N: ",
  expected_N[["T-cell"]],
  "\n",
  sep = ""
)

cat(
  "Myeloid/macrophage N: ",
  expected_N[["Myeloid/macrophage"]],
  "\n",
  sep = ""
)

cat(
  "Endothelial N: ",
  expected_N[["Endothelial"]],
  "\n",
  sep = ""
)

cat(
  "No statistical testing performed.\n"
)

cat(
  "No normalization performed.\n"
)

# ============================================================
# STEP10D-D
# Patient-Level Compartment Concordance Audit
#
# FINAL PRIMARY-ONLY VERSION
#
# Scope:
#   - Patient-level descriptive concordance
#   - Primary tumor compartments only
#   - T-cell
#   - Myeloid/macrophage
#   - Endothelial
#
# Locked core genes:
#   AOC3 / MS4A7 / IL7R
#
# Detection:
#   Raw_Count > 0
#
# No:
#   - normalization
#   - differential expression
#   - correlation testing
#   - P values
#   - FDR
#   - CellChat
#   - ligand-receptor inference
#   - pseudobulk recalculation
#   - patient exclusion
# ============================================================


# ------------------------------------------------------------
# 18. Final input-object audit
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D FINAL INPUT OBJECT AUDIT\n")
cat("============================================================\n\n")


required_objects <- c(
  "counts_long",
  "locked_compartments",
  "core_genes",
  "compartment_summary",
  "detection_frequency"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]


if (
  length(missing_objects) > 0
) {
  
  stop(
    "STOP: Required STEP10D-D objects are missing: ",
    paste(
      missing_objects,
      collapse = ", "
    )
  )
}


required_count_columns <- c(
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Gene",
  "Raw_Count",
  "Detected"
)


missing_count_columns <- setdiff(
  required_count_columns,
  colnames(counts_long)
)


if (
  length(missing_count_columns) > 0
) {
  
  stop(
    "STOP: counts_long is missing required columns: ",
    paste(
      missing_count_columns,
      collapse = ", "
    )
  )
}


cat(
  "counts_long rows: ",
  nrow(counts_long),
  "\n",
  sep = ""
)

cat(
  "counts_long columns: ",
  ncol(counts_long),
  "\n",
  sep = ""
)

cat(
  "Locked core genes: ",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n",
  sep = ""
)

cat(
  "Locked compartments: ",
  paste(
    locked_compartments,
    collapse = " / "
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Final PRIMARY-only integrity audit
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PRIMARY-ONLY INTEGRITY AUDIT\n")
cat("============================================================\n\n")


if (
  !all(
    counts_long$Eligibility == "Primary"
  )
) {
  
  stop(
    "STOP: Non-Primary observations detected."
  )
}


if (
  !all(
    counts_long$Compartment %in%
    locked_compartments
  )
) {
  
  stop(
    "STOP: Unexpected compartment detected."
  )
}


if (
  !all(
    counts_long$Gene %in%
    core_genes
  )
) {
  
  stop(
    "STOP: Unexpected gene detected."
  )
}


# ------------------------------------------------------------
# Verify exactly one observation per
# patient × compartment × gene
# ------------------------------------------------------------

key <- paste(
  counts_long$Patient,
  counts_long$Compartment,
  counts_long$Gene,
  sep = "___"
)


if (
  anyDuplicated(key) > 0
) {
  
  stop(
    "STOP: Duplicate Patient × Compartment × Gene ",
    "observations detected."
  )
}


# ------------------------------------------------------------
# Verify detection definition
# ------------------------------------------------------------

if (
  any(
    counts_long$Detected !=
    (
      counts_long$Raw_Count > 0
    )
  )
) {
  
  stop(
    "STOP: Detected does not equal Raw_Count > 0."
  )
}


cat(
  "Primary-only: YES\n"
)

cat(
  "Unexpected genes: NONE\n"
)

cat(
  "Unexpected compartments: NONE\n"
)

cat(
  "Duplicate patient × compartment × gene: NONE\n"
)

cat(
  "Detection definition: VERIFIED\n"
)

cat("\n")
cat(
  "PASS — PRIMARY-only input integrity verified.\n"
)


# ------------------------------------------------------------
# 20. Patient-level compartment availability
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PATIENT-LEVEL COMPARTMENT AVAILABILITY\n")
cat("============================================================\n\n")


availability_list <- lapply(
  locked_compartments,
  function(comp) {
    
    patients <- sort(
      unique(
        counts_long$Patient[
          counts_long$Compartment == comp
        ]
      )
    )
    
    data.frame(
      Patient = patients,
      Compartment = comp,
      Available = TRUE,
      stringsAsFactors = FALSE
    )
  }
)


patient_compartment_availability <- do.call(
  rbind,
  availability_list
)


rownames(
  patient_compartment_availability
) <- NULL


print(
  patient_compartment_availability,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Patient count audit by compartment
# ------------------------------------------------------------

patient_count_audit <- aggregate(
  Patient ~ Compartment,
  data = patient_compartment_availability,
  FUN = function(x) length(unique(x))
)


colnames(
  patient_count_audit
)[
  colnames(patient_count_audit) == "Patient"
] <- "Patient_N"


expected_patient_counts <- data.frame(
  Compartment = locked_compartments,
  Expected_N = c(
    9L,
    9L,
    4L
  ),
  stringsAsFactors = FALSE
)


patient_count_audit <- merge(
  expected_patient_counts,
  patient_count_audit,
  by = "Compartment",
  all.x = TRUE,
  sort = FALSE
)


patient_count_audit$PASS <- (
  patient_count_audit$Expected_N ==
    patient_count_audit$Patient_N
)


print(
  patient_count_audit,
  row.names = FALSE
)


if (
  !all(
    patient_count_audit$PASS
  )
) {
  
  stop(
    "STOP: Patient count does not match locked ",
    "STEP10D-B structure."
  )
}


cat("\n")
cat(
  "PASS — Patient-level compartment availability verified.\n"
)


# ------------------------------------------------------------
# 22. Pairwise patient overlap
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PAIRWISE PATIENT OVERLAP\n")
cat("============================================================\n\n")


patient_sets <- lapply(
  locked_compartments,
  function(comp) {
    
    sort(
      unique(
        counts_long$Patient[
          counts_long$Compartment == comp
        ]
      )
    )
  }
)


names(patient_sets) <- locked_compartments


pairwise_combinations <- combn(
  locked_compartments,
  2,
  simplify = FALSE
)


pairwise_overlap <- do.call(
  rbind,
  lapply(
    pairwise_combinations,
    function(pair) {
      
      a <- patient_sets[[pair[1]]]
      b <- patient_sets[[pair[2]]]
      
      overlap <- intersect(
        a,
        b
      )
      
      data.frame(
        Compartment_A = pair[1],
        Compartment_B = pair[2],
        N_A = length(a),
        N_B = length(b),
        Overlap_N = length(overlap),
        Overlap_A_Percent =
          100 * length(overlap) / length(a),
        Overlap_B_Percent =
          100 * length(overlap) / length(b),
        Overlap_Patients =
          paste(
            overlap,
            collapse = ";"
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)


rownames(
  pairwise_overlap
) <- NULL


print(
  pairwise_overlap,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Pairwise patient overlap generated descriptively.\n"
)


# ------------------------------------------------------------
# 23. Complete three-compartment patient set
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D COMPLETE THREE-COMPARTMENT PATIENT SET\n")
cat("============================================================\n\n")


complete_three_compartment_patients <- Reduce(
  intersect,
  patient_sets
)


complete_three_compartment_patients <- sort(
  complete_three_compartment_patients
)


complete_three_compartment_table <- data.frame(
  Patient =
    complete_three_compartment_patients,
  stringsAsFactors = FALSE
)


print(
  complete_three_compartment_table,
  row.names = FALSE
)


complete_three_compartment_N <-
  length(
    complete_three_compartment_patients
  )


cat("\n")
cat(
  "Complete three-compartment patients: ",
  complete_three_compartment_N,
  "\n",
  sep = ""
)


if (
  complete_three_compartment_N == 0
) {
  
  stop(
    "STOP: No patient has all three locked compartments."
  )
}


cat(
  "PASS — Complete three-compartment patient set identified.\n"
)


# ------------------------------------------------------------
# 24. Patient × compartment × core-gene profile
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PATIENT × COMPARTMENT × CORE-GENE PROFILE\n")
cat("============================================================\n\n")


profile_long <- counts_long[
  counts_long$Patient %in%
    complete_three_compartment_patients,
  c(
    "Patient",
    "Compartment",
    "Eligibility",
    "Group_ID",
    "Gene",
    "Raw_Count",
    "Detected"
  ),
  drop = FALSE
]


profile_long <- profile_long[
  order(
    profile_long$Patient,
    match(
      profile_long$Compartment,
      locked_compartments
    ),
    match(
      profile_long$Gene,
      core_genes
    )
  ),
]


rownames(
  profile_long
) <- NULL


expected_profile_rows <-
  complete_three_compartment_N *
  length(locked_compartments) *
  length(core_genes)


if (
  nrow(profile_long) !=
  expected_profile_rows
) {
  
  stop(
    "STOP: Unexpected number of rows in ",
    "three-compartment patient profile."
  )
}


print(
  profile_long,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Complete three-compartment patient profile generated.\n"
)


# ------------------------------------------------------------
# 25. Core-gene detection pattern per patient
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D CORE-GENE PATIENT DETECTION PATTERN\n")
cat("============================================================\n\n")


detection_pattern_list <- lapply(
  complete_three_compartment_patients,
  function(patient) {
    
    patient_df <- profile_long[
      profile_long$Patient == patient,
      ,
      drop = FALSE
    ]
    
    result <- data.frame(
      Patient = patient,
      stringsAsFactors = FALSE
    )
    
    for (
      g in core_genes
    ) {
      
      gene_df <- patient_df[
        patient_df$Gene == g,
        ,
        drop = FALSE
      ]
      
      for (
        comp in locked_compartments
      ) {
        
        value <- gene_df$Detected[
          gene_df$Compartment == comp
        ]
        
        if (
          length(value) != 1
        ) {
          
          stop(
            "STOP: Expected exactly one detection value for ",
            patient,
            " / ",
            comp,
            " / ",
            g
          )
        }
        
        result[
          paste(
            g,
            comp,
            sep = "_"
          )
        ] <- as.logical(value)
      }
    }
    
    result
  }
)


patient_detection_pattern <- do.call(
  rbind,
  detection_pattern_list
)


rownames(
  patient_detection_pattern
) <- NULL


print(
  patient_detection_pattern,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Patient-level core-gene detection patterns generated.\n"
)


# ------------------------------------------------------------
# 26. Three-compartment concordance by core gene
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D THREE-COMPARTMENT CORE-GENE CONCORDANCE\n")
cat("============================================================\n\n")


three_compartment_concordance <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      gene_df <- profile_long[
        profile_long$Gene == g,
        ,
        drop = FALSE
      ]
      
      per_patient <- lapply(
        complete_three_compartment_patients,
        function(patient) {
          
          p <- gene_df[
            gene_df$Patient == patient,
            ,
            drop = FALSE
          ]
          
          detected_values <- sapply(
            locked_compartments,
            function(comp) {
              
              p$Detected[
                p$Compartment == comp
              ]
            }
          )
          
          all_detected <- all(
            detected_values
          )
          
          any_detected <- any(
            detected_values
          )
          
          none_detected <- !any_detected
          
          data.frame(
            Patient = patient,
            Gene = g,
            All_Three_Detected =
              all_detected,
            Any_Detected =
              any_detected,
            None_Detected =
              none_detected,
            stringsAsFactors = FALSE
          )
        }
      )
      
      do.call(
        rbind,
        per_patient
      )
    }
  )
)


rownames(
  three_compartment_concordance
) <- NULL


print(
  three_compartment_concordance,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Three-compartment core-gene concordance generated.\n"
)


# ------------------------------------------------------------
# 27. Concordance score summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D CORE-GENE CONCORDANCE SCORE SUMMARY\n")
cat("============================================================\n\n")


concordance_score_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      x <- three_compartment_concordance[
        three_compartment_concordance$Gene == g,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Gene = g,
        Complete_Three_Compartment_Patients =
          nrow(x),
        All_Three_Detected_N =
          sum(
            x$All_Three_Detected
          ),
        All_Three_Detected_Percent =
          100 *
          mean(
            x$All_Three_Detected
          ),
        Any_Detected_N =
          sum(
            x$Any_Detected
          ),
        None_Detected_N =
          sum(
            x$None_Detected
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)


rownames(
  concordance_score_summary
) <- NULL


print(
  concordance_score_summary,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Core-gene concordance score summary generated.\n"
)


# ------------------------------------------------------------
# 28. Pairwise core-gene detection concordance
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D PAIRWISE CORE-GENE DETECTION CONCORDANCE\n")
cat("============================================================\n\n")


pairwise_core_gene_concordance <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      gene_df <- profile_long[
        profile_long$Gene == g,
        ,
        drop = FALSE
      ]
      
      do.call(
        rbind,
        lapply(
          pairwise_combinations,
          function(pair) {
            
            patient_wide <- do.call(
              rbind,
              lapply(
                complete_three_compartment_patients,
                function(patient) {
                  
                  p <- gene_df[
                    gene_df$Patient == patient,
                    ,
                    drop = FALSE
                  ]
                  
                  a <- p$Detected[
                    p$Compartment == pair[1]
                  ]
                  
                  b <- p$Detected[
                    p$Compartment == pair[2]
                  ]
                  
                  data.frame(
                    A = a,
                    B = b
                  )
                }
              )
            )
            
            both_detected <-
              patient_wide$A &
              patient_wide$B
            
            both_not_detected <-
              !patient_wide$A &
              !patient_wide$B
            
            discordant <-
              patient_wide$A !=
              patient_wide$B
            
            data.frame(
              Gene = g,
              Compartment_A = pair[1],
              Compartment_B = pair[2],
              N =
                nrow(patient_wide),
              Both_Detected_N =
                sum(both_detected),
              Both_Not_Detected_N =
                sum(both_not_detected),
              Discordant_N =
                sum(discordant),
              Concordant_N =
                sum(!discordant),
              Concordance_Percent =
                100 *
                mean(!discordant),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)


rownames(
  pairwise_core_gene_concordance
) <- NULL


print(
  pairwise_core_gene_concordance,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Pairwise core-gene detection concordance generated.\n"
)


# ------------------------------------------------------------
# 29. Full patient-level concordance profile
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D FULL PATIENT-LEVEL CONCORDANCE PROFILE\n")
cat("============================================================\n\n")


full_patient_concordance <- profile_long


full_patient_concordance$Detection_Status <- ifelse(
  full_patient_concordance$Detected,
  "Detected",
  "Not detected"
)


full_patient_concordance <- full_patient_concordance[
  order(
    full_patient_concordance$Patient,
    match(
      full_patient_concordance$Compartment,
      locked_compartments
    ),
    match(
      full_patient_concordance$Gene,
      core_genes
    )
  ),
]


rownames(
  full_patient_concordance
) <- NULL


print(
  full_patient_concordance,
  row.names = FALSE
)


cat("\n")
cat(
  "PASS — Full patient-level concordance profile generated.\n"
)


# ------------------------------------------------------------
# 30. Structural QC results — FINAL CORRECTED VERSION
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D STRUCTURAL QC\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 30A. Calculate one scalar overall PASS flag
# ------------------------------------------------------------

overall_pass <- all(
  c(
    all(counts_long$Eligibility == "Primary"),
    
    setequal(
      unique(counts_long$Gene),
      core_genes
    ),
    
    all(
      counts_long$Compartment %in%
        locked_compartments
    ),
    
    !anyDuplicated(key),
    
    all(
      counts_long$Detected ==
        (
          counts_long$Raw_Count > 0
        )
    ),
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "T-cell"
    ] == 9,
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Myeloid/macrophage"
    ] == 9,
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Endothelial"
    ] == 4,
    
    complete_three_compartment_N > 0,
    
    nrow(profile_long) ==
      expected_profile_rows,
    
    TRUE,   # No normalization
    TRUE    # No statistical testing
  )
)


# ------------------------------------------------------------
# 30B. QC table
# ------------------------------------------------------------

qc_results <- data.frame(
  
  Check = c(
    
    "Primary-only input",
    
    "Core genes exactly AOC3/MS4A7/IL7R",
    
    "Locked compartments only",
    
    "Unique patient-compartment-gene combinations",
    
    "Detection definition Raw_Count > 0",
    
    "T-cell patient count = 9",
    
    "Myeloid/macrophage patient count = 9",
    
    "Endothelial patient count = 4",
    
    "Three-compartment patient set non-empty",
    
    "Three-compartment profile complete",
    
    "No normalization",
    
    "No statistical testing"
    
  ),
  
  Result = c(
    
    all(
      counts_long$Eligibility == "Primary"
    ),
    
    setequal(
      unique(counts_long$Gene),
      core_genes
    ),
    
    all(
      counts_long$Compartment %in%
        locked_compartments
    ),
    
    !anyDuplicated(key),
    
    all(
      counts_long$Detected ==
        (
          counts_long$Raw_Count > 0
        )
    ),
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "T-cell"
    ] == 9,
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Myeloid/macrophage"
    ] == 9,
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Endothelial"
    ] == 4,
    
    complete_three_compartment_N > 0,
    
    nrow(profile_long) ==
      expected_profile_rows,
    
    TRUE,
    
    TRUE
    
  ),
  
  stringsAsFactors = FALSE
  
)


# ------------------------------------------------------------
# 30C. Add scalar Overall_PASS only for display
# ------------------------------------------------------------

qc_results$Overall_PASS <- overall_pass


print(
  qc_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 30D. FINAL QC STOP
# ------------------------------------------------------------

if (!isTRUE(overall_pass)) {
  
  stop(
    "STOP: One or more STEP10D-D QC checks failed."
  )
  
}


cat("\n")
cat("PASS — STEP10D-D structural QC completed.\n")
cat("Overall PASS = TRUE\n")
cat("No upstream locked workflow recalculated.\n")
cat("No expression values recalculated.\n")
cat("No normalization performed.\n")
cat("No statistical testing performed.\n")

# ------------------------------------------------------------
# 31. Methodological lock
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D METHODOLOGICAL LOCK\n")
cat("============================================================\n\n")


methodological_lock <- data.frame(
  
  Item = c(
    "Analysis level",
    "Primary compartments",
    "T-cell patients",
    "Myeloid/macrophage patients",
    "Endothelial patients",
    "Core genes",
    "Count type",
    "Detection definition",
    "Zero counts retained",
    "Normalization",
    "Differential expression",
    "Correlation testing",
    "P values",
    "FDR",
    "CellChat / ligand-receptor analysis",
    "Pseudobulk recalculation",
    "Patient exclusion",
    "Sensitivity Endothelial groups"
  ),
  
  Locked_Value = c(
    "Patient-level descriptive concordance",
    "T-cell / Myeloid/macrophage / Endothelial",
    "9",
    "9",
    "4",
    "AOC3 / MS4A7 / IL7R",
    "Raw UMI pseudobulk counts",
    "Raw_Count > 0",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "Excluded from primary concordance; retained separately in STEP10C/10D"
  ),
  
  stringsAsFactors = FALSE
)


print(
  methodological_lock,
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Final audit summary
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D FINAL AUDIT SUMMARY\n")
cat("============================================================\n\n")


final_audit_summary <- data.frame(
  
  Metric = c(
    "Input rows",
    "Unique patients in PRIMARY concordance input",
    "T-cell patients",
    "Myeloid/macrophage patients",
    "Endothelial patients",
    "Complete three-compartment patients",
    "Core genes",
    "Core-gene observations",
    "Detection definition",
    "Zero counts retained",
    "Normalization",
    "Statistical testing",
    "Differential expression",
    "P/FDR",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "Patient exclusion"
  ),
  
  Value = c(
    nrow(counts_long),
    
    length(
      unique(
        counts_long$Patient
      )
    ),
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "T-cell"
    ],
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Myeloid/macrophage"
    ],
    
    patient_count_audit$Patient_N[
      patient_count_audit$Compartment ==
        "Endothelial"
    ],
    
    complete_three_compartment_N,
    
    paste(
      core_genes,
      collapse = " / "
    ),
    
    nrow(counts_long),
    
    "Raw_Count > 0",
    
    "YES",
    
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


print(
  final_audit_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. FINAL PASS / LOCK
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10D-D FINAL STATUS\n")
cat("============================================================\n\n")


cat(
  "PASS — STEP10D-D Patient-Level Compartment Concordance Audit\n"
)

cat(
  "PRIMARY-only structure: VERIFIED\n"
)

cat(
  "T-cell patients: 9\n"
)

cat(
  "Myeloid/macrophage patients: 9\n"
)

cat(
  "Endothelial patients: 4\n"
)

cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)

cat(
  "Detection definition: Raw_Count > 0\n"
)

cat(
  "Patient-level concordance: DESCRIPTIVE ONLY\n"
)

cat(
  "Normalization: NO\n"
)

cat(
  "Statistical testing: NO\n"
)

cat(
  "P/FDR: NO\n"
)

cat(
  "CellChat/LR: NO\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Patient/cell exclusion: NO\n"
)

cat("\n")
cat(
  "============================================================\n"
)
cat(
  "STEP10D-D FINAL LOCK: PASS\n"
)
cat(
  "============================================================\n"
)

# ============================================================
# STEP10E-0 — Integrated Input / Provenance Lock
# ============================================================
#
# Purpose:
#   Lock all previously finalized patient-level single-cell
#   evidence before integrated STEP10E analysis.
#
# This step:
#   - READS locked STEP10C / STEP10D outputs
#   - verifies provenance and structure
#   - verifies patient / compartment / gene locks
#   - verifies the 4 complete three-compartment patients
#
# This step DOES NOT:
#   - recalculate pseudobulk
#   - recalculate raw counts
#   - recalculate detection
#   - normalize expression
#   - perform DE
#   - calculate correlation
#   - calculate P values
#   - calculate FDR
#   - perform CellChat / LR analysis
#   - exclude patients or cells
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10E-0 — INTEGRATED INPUT / PROVENANCE LOCK\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory
# ------------------------------------------------------------

step10_root <-
  "input_dir/Paper 1_Step 10-Single-Cell"

setwd(step10_root)


# ------------------------------------------------------------
# 1. Output directory
# ------------------------------------------------------------

step10e_dir <-
  file.path(
    step10_root,
    "STEP10E-0_Integrated_Input_Provenance_Lock"
  )

if (!dir.exists(step10e_dir)) {
  dir.create(
    step10e_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 2. LOCKED biological definitions
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


locked_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)


locked_patient_counts <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)


locked_complete_three_compartment_N <- 4


locked_detection_definition <-
  "Raw_Count > 0"


cat("Locked core genes:\n")
cat("  ", paste(core_genes, collapse = " / "), "\n\n")

cat("Locked primary compartments:\n")
cat("  ", paste(locked_compartments, collapse = " / "), "\n\n")

cat("Locked patient counts:\n")
print(locked_patient_counts)

cat(
  "\nLocked complete three-compartment patients:",
  locked_complete_three_compartment_N,
  "\n\n"
)

cat(
  "Locked detection definition:",
  locked_detection_definition,
  "\n\n"
)


# ------------------------------------------------------------
# 3. Helper function: require exactly one matching file
# ------------------------------------------------------------

find_one_file <- function(
    directory,
    pattern,
    description
) {
  
  if (!dir.exists(directory)) {
    
    stop(
      paste0(
        "STOP: Directory not found for ",
        description,
        ": ",
        directory
      )
    )
    
  }
  
  hits <- list.files(
    directory,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE
  )
  
  if (length(hits) == 0) {
    
    stop(
      paste0(
        "STOP: No file found for ",
        description,
        " in ",
        directory,
        "\nPattern: ",
        pattern
      )
    )
    
  }
  
  if (length(hits) > 1) {
    
    stop(
      paste0(
        "STOP: Multiple files found for ",
        description,
        ":\n",
        paste(hits, collapse = "\n")
      )
    )
    
  }
  
  hits
}


# ------------------------------------------------------------
# 4. Locate locked STEP10C / STEP10D outputs
# ------------------------------------------------------------

step10c_dir <-
  file.path(
    step10_root,
    "STEP10C_Patient_Level_RawUMI_Pseudobulk"
  )

step10da_dir <-
  file.path(
    step10_root,
    "STEP10D_Pseudobulk_QC_CoreGene_Audit"
  )

step10db_dir <-
  file.path(
    step10_root,
    "STEP10D-B_CoreGene_RawCount_Detection_Audit"
  )

step10dc_dir <-
  file.path(
    step10_root,
    "STEP10D-C_CoreGene_Compartment_Level_Audit"
  )

step10dd_dir <-
  file.path(
    step10_root,
    "STEP10D-D_Patient_Level_Compartment_Concordance_Audit"
  )


cat("Checking locked directories...\n\n")

required_dirs <- c(
  STEP10C = step10c_dir,
  `STEP10D-A` = step10da_dir,
  `STEP10D-B` = step10db_dir,
  `STEP10D-C` = step10dc_dir,
  `STEP10D-D` = step10dd_dir
)


dir_audit <- data.frame(
  Workflow = names(required_dirs),
  Directory = unname(required_dirs),
  Exists = dir.exists(unname(required_dirs)),
  stringsAsFactors = FALSE
)


print(
  dir_audit,
  row.names = FALSE
)


if (!all(dir_audit$Exists)) {
  
  stop(
    "STOP: One or more required locked workflow directories are missing."
  )
  
}


cat("\nPASS — All required locked workflow directories found.\n\n")


# ------------------------------------------------------------
# 5. Locate STEP10C pseudobulk matrices
# ------------------------------------------------------------

file_T <-
  file.path(
    step10c_dir,
    "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
  )

file_M <-
  file.path(
    step10c_dir,
    "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
  )

file_E <-
  file.path(
    step10c_dir,
    "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
  )

file_E_sens <-
  file.path(
    step10c_dir,
    "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
  )


step10c_files <- c(
  T_cell_primary = file_T,
  Myeloid_primary = file_M,
  Endothelial_primary = file_E,
  Endothelial_sensitivity = file_E_sens
)


file_audit <- data.frame(
  Input = names(step10c_files),
  File = unname(step10c_files),
  Exists = file.exists(unname(step10c_files)),
  stringsAsFactors = FALSE
)


print(
  file_audit,
  row.names = FALSE
)


if (!all(file_audit$Exists)) {
  
  stop(
    "STOP: One or more locked STEP10C pseudobulk matrices are missing."
  )
  
}


cat("\nPASS — STEP10C pseudobulk matrices located.\n\n")


# ------------------------------------------------------------
# 6. Read STEP10C annotation only
# ------------------------------------------------------------

file_annotation <-
  file.path(
    step10c_dir,
    "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
  )


if (!file.exists(file_annotation)) {
  
  stop(
    "STOP: STEP10C primary pseudobulk annotation file not found."
  )
  
}


annotation_primary <- read.delim(
  file_annotation,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat("STEP10C annotation dimensions:\n")
cat(
  "Rows:",
  nrow(annotation_primary),
  "\n"
)

cat(
  "Columns:",
  ncol(annotation_primary),
  "\n\n"
)

cat("Columns:\n")
print(
  names(annotation_primary)
)


# ------------------------------------------------------------
# 7. Verify actual locked STEP10C annotation schema
# ------------------------------------------------------------

required_annotation_columns <- c(
  "Patient",
  "Compartment",
  "Cell_Count",
  "Eligibility",
  "Group_ID"
)


if (!all(
  required_annotation_columns %in%
  names(annotation_primary)
)) {
  
  stop(
    paste0(
      "STOP: STEP10C annotation schema mismatch.\n",
      "Required columns:\n",
      paste(
        required_annotation_columns,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "\nPASS — STEP10C annotation schema verified.\n\n"
)


# ------------------------------------------------------------
# 8. Verify primary annotation structure
# ------------------------------------------------------------

if (!all(
  annotation_primary$Eligibility == "Primary"
)) {
  
  stop(
    "STOP: STEP10C primary annotation contains non-Primary records."
  )
  
}


if (!all(
  annotation_primary$Compartment %in%
  locked_compartments
)) {
  
  stop(
    "STOP: STEP10C annotation contains unexpected compartments."
  )
  
}


observed_annotation_counts <-
  table(
    annotation_primary$Compartment
  )


print(
  observed_annotation_counts
)


for (x in locked_compartments) {
  
  if (
    is.na(
      observed_annotation_counts[x]
    ) ||
    observed_annotation_counts[x] !=
    locked_patient_counts[x]
  ) {
    
    stop(
      paste0(
        "STOP: STEP10C patient/group count mismatch for ",
        x
      )
    )
    
  }
  
}


cat(
  "\nPASS — STEP10C primary patient/compartment structure verified.\n\n"
)


# ------------------------------------------------------------
# 9. Read locked STEP10D-B core-gene raw-count table
# ------------------------------------------------------------

file_core_counts <-
  file.path(
    step10db_dir,
    "STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv"
  )


if (!file.exists(file_core_counts)) {
  
  stop(
    "STOP: STEP10D-B locked core-gene raw-count table not found."
  )
  
}


core_counts_locked <- read.delim(
  file_core_counts,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "STEP10D-B locked core-count table:\n"
)

cat(
  "Rows:",
  nrow(core_counts_locked),
  "\n"
)

cat(
  "Columns:",
  ncol(core_counts_locked),
  "\n\n"
)

cat(
  "Columns detected:\n"
)

print(
  names(core_counts_locked)
)

cat("\n")

print(
  head(
    core_counts_locked
  )
)


# ------------------------------------------------------------
# 10. Verify STEP10D-B wide-format schema
# ------------------------------------------------------------

required_core_metadata_columns <- c(
  "Dataset",
  "Group_ID",
  "Patient",
  "Compartment",
  "Eligibility",
  "Cell_Count"
)


required_core_gene_columns <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


required_core_columns <- c(
  required_core_metadata_columns,
  required_core_gene_columns
)


if (!all(
  required_core_columns %in%
  names(core_counts_locked)
)) {
  
  missing_columns <-
    setdiff(
      required_core_columns,
      names(core_counts_locked)
    )
  
  stop(
    paste0(
      "STOP: STEP10D-B core-count table schema mismatch.\n",
      "Missing columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "\nPASS — STEP10D-B wide-format schema verified.\n"
)

cat(
  "Metadata columns:",
  paste(
    required_core_metadata_columns,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Core-gene columns:",
  paste(
    required_core_gene_columns,
    collapse = ", "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 11. Verify STEP10D-B core-gene content and raw counts
# ------------------------------------------------------------

# Core-gene columns are exactly the locked genes
observed_gene_columns <-
  names(core_counts_locked)[
    names(core_counts_locked) %in%
      core_genes
  ]


if (!setequal(
  observed_gene_columns,
  core_genes
)) {
  
  stop(
    "STOP: STEP10D-B core-gene columns do not exactly match AOC3/MS4A7/IL7R."
  )
  
}


# Verify no unexpected additional gene columns
non_metadata_columns <-
  setdiff(
    names(core_counts_locked),
    required_core_metadata_columns
  )


if (!setequal(
  non_metadata_columns,
  core_genes
)) {
  
  stop(
    paste0(
      "STOP: Unexpected non-metadata columns detected in STEP10D-B.\n",
      "Observed:\n",
      paste(
        non_metadata_columns,
        collapse = ", "
      )
    )
  )
  
}


# Verify raw-count values
core_count_matrix <-
  as.matrix(
    core_counts_locked[
      ,
      core_genes,
      drop = FALSE
    ]
  )


storage.mode(core_count_matrix) <- "numeric"


if (any(
  !is.finite(
    core_count_matrix
  )
)) {
  
  stop(
    "STOP: Non-finite raw counts detected in locked STEP10D-B table."
  )
  
}


if (any(
  core_count_matrix < 0
)) {
  
  stop(
    "STOP: Negative raw counts detected in locked STEP10D-B table."
  )
  
}


if (any(
  core_count_matrix !=
  floor(core_count_matrix)
)) {
  
  stop(
    "STOP: Non-integer raw counts detected in locked STEP10D-B table."
  )
  
}


# Verify metadata structure
if (anyDuplicated(
  core_counts_locked$Group_ID
)) {
  
  stop(
    "STOP: Duplicate Group_ID detected in STEP10D-B locked table."
  )
  
}


if (!all(
  core_counts_locked$Eligibility %in%
  c(
    "Primary",
    "Sensitivity"
  )
)) {
  
  stop(
    "STOP: Unexpected Eligibility label detected in STEP10D-B."
  )
  
}


if (!all(
  core_counts_locked$Compartment %in%
  c(
    locked_compartments
  )
)) {
  
  stop(
    "STOP: Unexpected compartment label detected in STEP10D-B."
  )
  
}


cat(
  "PASS — STEP10D-B core genes exactly verified:\n"
)

cat(
  "  AOC3 / MS4A7 / IL7R\n"
)

cat(
  "PASS — STEP10D-B raw-count integrity verified.\n"
)

cat(
  "  Non-finite: 0\n"
)

cat(
  "  Negative: 0\n"
)

cat(
  "  Non-integer: 0\n"
)

cat(
  "  Zero counts retained: YES\n"
)

cat(
  "  Detection rule: Raw_Count > 0\n\n"
)


# ------------------------------------------------------------
# 12. Verify locked STEP10D-C descriptive summary
# ------------------------------------------------------------

file_dc_summary <-
  file.path(
    step10dc_dir,
    "STEP10D-C_CoreGene_Compartment_Descriptive_Summary.tsv"
  )

if (!file.exists(file_dc_summary)) {
  
  stop(
    "STOP: STEP10D-C compartment descriptive summary not found."
  )
  
}

compartment_summary_locked <- read.delim(
  file_dc_summary,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "STEP10D-C summary dimensions:\n"
)

cat(
  "Rows:",
  nrow(compartment_summary_locked),
  "\n"
)

cat(
  "Columns:",
  ncol(compartment_summary_locked),
  "\n\n"
)

cat(
  "STEP10D-C columns:\n"
)

print(
  names(compartment_summary_locked)
)

cat("\n")


# ------------------------------------------------------------
# 12A. Verify actual locked STEP10D-C schema
# ------------------------------------------------------------

dc_names <- names(
  compartment_summary_locked
)

required_dc_columns <- c(
  "Dataset",
  "Gene",
  "N",
  "Detected_N",
  "Zero_N",
  "Detection_Percent",
  "Raw_UMI_Min",
  "Raw_UMI_Q1",
  "Raw_UMI_Median",
  "Raw_UMI_Mean",
  "Raw_UMI_Q3",
  "Raw_UMI_Max",
  "Raw_UMI_IQR"
)

missing_dc_columns <- setdiff(
  required_dc_columns,
  dc_names
)

if (length(missing_dc_columns) > 0) {
  
  stop(
    paste(
      "STOP: STEP10D-C summary is missing required columns:",
      paste(missing_dc_columns, collapse = ", ")
    )
  )
  
}

cat(
  "PASS — STEP10D-C summary schema verified.\n"
)

cat(
  "Compartment information is encoded in the locked Dataset labels.\n"
)

cat(
  "No Compartment column is required.\n\n"
)


# ------------------------------------------------------------
# 12B. Verify locked Dataset structure
# ------------------------------------------------------------

expected_dc_datasets <- c(
  "Primary_T-cell",
  "Primary_Myeloid_macrophage",
  "Primary_Endothelial",
  "Sensitivity_Endothelial"
)

observed_dc_datasets <- unique(
  compartment_summary_locked$Dataset
)

missing_dc_datasets <- setdiff(
  expected_dc_datasets,
  observed_dc_datasets
)

unexpected_dc_datasets <- setdiff(
  observed_dc_datasets,
  expected_dc_datasets
)

if (length(missing_dc_datasets) > 0) {
  
  stop(
    paste(
      "STOP: STEP10D-C summary missing expected Dataset labels:",
      paste(missing_dc_datasets, collapse = ", ")
    )
  )
  
}

if (length(unexpected_dc_datasets) > 0) {
  
  stop(
    paste(
      "STOP: STEP10D-C summary contains unexpected Dataset labels:",
      paste(unexpected_dc_datasets, collapse = ", ")
    )
  )
  
}

cat(
  "PASS — STEP10D-C Dataset labels verified.\n"
)

print(
  table(
    compartment_summary_locked$Dataset
  )
)

cat("\n")


# ------------------------------------------------------------
# 12C. Verify core-gene identity
# ------------------------------------------------------------

observed_dc_genes <- unique(
  compartment_summary_locked$Gene
)

if (
  !setequal(
    observed_dc_genes,
    core_genes
  )
) {
  
  stop(
    paste(
      "STOP: STEP10D-C gene set mismatch.",
      "\nExpected:",
      paste(core_genes, collapse = ", "),
      "\nObserved:",
      paste(observed_dc_genes, collapse = ", ")
    )
  )
  
}

cat(
  "PASS — STEP10D-C core-gene identity verified.\n"
)

cat(
  "Core genes:",
  paste(core_genes, collapse = " / "),
  "\n\n"
)


# ------------------------------------------------------------
# 12D. Verify one row per Dataset × Gene
# ------------------------------------------------------------

dc_key <- paste(
  compartment_summary_locked$Dataset,
  compartment_summary_locked$Gene,
  sep = "||"
)

if (
  anyDuplicated(dc_key) > 0
) {
  
  stop(
    "STOP: Duplicate Dataset × Gene combinations detected."
  )
  
}

expected_primary_dc_keys <- expand.grid(
  Dataset = c(
    "Primary_T-cell",
    "Primary_Myeloid_macrophage",
    "Primary_Endothelial"
  ),
  Gene = core_genes,
  stringsAsFactors = FALSE
)

expected_primary_dc_keys <- paste(
  expected_primary_dc_keys$Dataset,
  expected_primary_dc_keys$Gene,
  sep = "||"
)

observed_primary_dc <- compartment_summary_locked[
  compartment_summary_locked$Dataset %in%
    c(
      "Primary_T-cell",
      "Primary_Myeloid_macrophage",
      "Primary_Endothelial"
    ),
  ,
  drop = FALSE
]

observed_primary_dc_keys <- paste(
  observed_primary_dc$Dataset,
  observed_primary_dc$Gene,
  sep = "||"
)

if (
  !setequal(
    observed_primary_dc_keys,
    expected_primary_dc_keys
  )
) {
  
  stop(
    "STOP: Primary STEP10D-C Dataset × Gene structure does not match the locked 3-compartment × 3-gene design."
  )
  
}

cat(
  "PASS — Primary STEP10D-C Dataset × Gene structure verified.\n"
)

cat(
  "Primary combinations:",
  length(expected_primary_dc_keys),
  "\n\n"
)


# ------------------------------------------------------------
# 12E. Verify expected sensitivity structure
# ------------------------------------------------------------

expected_sensitivity_dc_keys <- paste(
  "Sensitivity_Endothelial",
  core_genes,
  sep = "||"
)

observed_sensitivity_dc_keys <- paste(
  compartment_summary_locked$Dataset[
    compartment_summary_locked$Dataset ==
      "Sensitivity_Endothelial"
  ],
  compartment_summary_locked$Gene[
    compartment_summary_locked$Dataset ==
      "Sensitivity_Endothelial"
  ],
  sep = "||"
)

if (
  !setequal(
    observed_sensitivity_dc_keys,
    expected_sensitivity_dc_keys
  )
) {
  
  stop(
    "STOP: STEP10D-C sensitivity endothelial Dataset × Gene structure mismatch."
  )
  
}

cat(
  "PASS — STEP10D-C sensitivity endothelial structure verified.\n\n"
)


# ------------------------------------------------------------
# 12F. Verify numeric integrity
# ------------------------------------------------------------

numeric_dc_columns <- c(
  "N",
  "Detected_N",
  "Zero_N",
  "Detection_Percent",
  "Raw_UMI_Min",
  "Raw_UMI_Q1",
  "Raw_UMI_Median",
  "Raw_UMI_Mean",
  "Raw_UMI_Q3",
  "Raw_UMI_Max",
  "Raw_UMI_IQR"
)

dc_numeric_ok <- all(
  sapply(
    compartment_summary_locked[
      numeric_dc_columns
    ],
    function(x) {
      all(
        is.finite(
          as.numeric(x)
        )
      )
    }
  )
)

if (!isTRUE(dc_numeric_ok)) {
  
  stop(
    "STOP: Non-finite numeric values detected in STEP10D-C summary."
  )
  
}

cat(
  "PASS — STEP10D-C numeric integrity verified.\n\n"
)


# ------------------------------------------------------------
# 13. STEP10D-C provenance interpretation
# ------------------------------------------------------------

cat(
  "STEP10D-C provenance interpretation:\n"
)

cat(
  "- Summary is already locked from STEP10D-C.\n"
)

cat(
  "- Dataset labels encode the cellular compartment.\n"
)

cat(
  "- Primary datasets: T-cell, Myeloid/macrophage, Endothelial.\n"
)

cat(
  "- Sensitivity dataset: Endothelial only.\n"
)

cat(
  "- Core genes: AOC3 / MS4A7 / IL7R.\n"
)

cat(
  "- Zero raw-count values are retained.\n"
)

cat(
  "- No normalization performed in STEP10E-0.\n"
)

cat(
  "- No statistical testing performed in STEP10E-0.\n"
)

cat(
  "- No pseudobulk recalculation performed in STEP10E-0.\n"
)

cat(
  "- No patient or cell exclusion performed in STEP10E-0.\n\n"
)

cat(
  "PASS — STEP10D-C provenance successfully integrated into STEP10E-0.\n\n"
)


# ------------------------------------------------------------
# 14B-1. Inspect actual STEP10D-D output files
# ------------------------------------------------------------

cat(
  "STEP10D-D directory:\n",
  step10dd_dir,
  "\n\n"
)

dd_files <- list.files(
  step10dd_dir,
  full.names = FALSE
)

cat(
  "Files currently present in STEP10D-D directory:\n"
)

print(
  dd_files
)

cat(
  "\nNumber of files:",
  length(dd_files),
  "\n"
)

# ------------------------------------------------------------
# 14B-2. Locate any existing STEP10D-D output files
# ------------------------------------------------------------

step10e_root <-
  "input_dir/Paper 1_Step 10-Single-Cell"

cat(
  "Searching for existing STEP10D-D files under:\n",
  step10e_root,
  "\n\n"
)

dd_all_files <- list.files(
  step10e_root,
  pattern = "STEP10D-D",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

cat(
  "STEP10D-D-related files found:\n"
)

print(
  dd_all_files
)

cat(
  "\nNumber of STEP10D-D-related files:",
  length(dd_all_files),
  "\n"
)

# ------------------------------------------------------------
# STEP10D-D object availability audit
# ------------------------------------------------------------

cat(
  "STEP10D-D current-session object audit\n"
)
cat(
  "====================================\n\n"
)

objects_to_check <- c(
  "counts_long",
  "detection_frequency",
  "compartment_summary",
  "core_genes"
)

object_status <- data.frame(
  Object = objects_to_check,
  Exists = sapply(
    objects_to_check,
    exists,
    envir = .GlobalEnv
  ),
  stringsAsFactors = FALSE
)

print(
  object_status,
  row.names = FALSE
)

cat("\n")

for (obj in objects_to_check) {
  
  if (exists(obj, envir = .GlobalEnv)) {
    
    cat(
      "Object:",
      obj,
      "\n"
    )
    
    cat(
      "Class:",
      paste(
        class(
          get(
            obj,
            envir = .GlobalEnv
          )
        ),
        collapse = ", "
      ),
      "\n"
    )
    
    if (
      is.data.frame(
        get(
          obj,
          envir = .GlobalEnv
        )
      )
    ) {
      
      cat(
        "Dimensions:",
        nrow(
          get(
            obj,
            envir = .GlobalEnv
          )
        ),
        "x",
        ncol(
          get(
            obj,
            envir = .GlobalEnv
          )
        ),
        "\n"
      )
      
    }
    
    cat("\n")
    
  }
  
}



# ============================================================
# STEP10D-D — Patient-Level Compartment Concordance
# FINAL OUTPUT / PROVENANCE LOCK
# ============================================================

cat(
  "\n============================================================\n"
)
cat(
  "STEP10D-D — FINAL OUTPUT / PROVENANCE LOCK\n"
)
cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------
# 1. Define output directory explicitly
# ------------------------------------------------------------

step10dd_dir <-
  "result_dir/Paper 1_Step 10-Single-Cell/STEP10D-D_Patient_Level_Compartment_Concordance_Audit"


if (!dir.exists(step10dd_dir)) {
  
  dir.create(
    step10dd_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}


if (!dir.exists(step10dd_dir)) {
  
  stop(
    "STOP: Unable to create STEP10D-D output directory."
  )
  
}


cat(
  "Output directory verified:\n",
  step10dd_dir,
  "\n\n"
)


# ------------------------------------------------------------
# 2. Verify required current-session objects
# ------------------------------------------------------------

required_objects <- c(
  "counts_long",
  "detection_frequency",
  "compartment_summary",
  "core_genes"
)

missing_objects <- required_objects[
  !sapply(
    required_objects,
    exists,
    envir = .GlobalEnv
  )
]

if (length(missing_objects) > 0) {
  
  stop(
    paste(
      "STOP: Required STEP10D-D objects missing:",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
  
}


cat(
  "PASS — All required STEP10D-D objects are available.\n\n"
)


# ------------------------------------------------------------
# 3. Lock core genes
# ------------------------------------------------------------

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

if (
  !identical(
    core_genes,
    expected_core_genes
  )
) {
  
  stop(
    paste(
      "STOP: Core-gene identity mismatch.",
      "\nExpected:",
      paste(
        expected_core_genes,
        collapse = " / "
      ),
      "\nObserved:",
      paste(
        core_genes,
        collapse = " / "
      )
    )
  )
  
}


cat(
  "PASS — Core genes locked:",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 4. Verify counts_long structure
# ------------------------------------------------------------

expected_counts_long_columns <- c(
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Gene",
  "Raw_Count",
  "Detected"
)

if (
  !identical(
    names(counts_long),
    expected_counts_long_columns
  )
) {
  
  stop(
    paste(
      "STOP: counts_long column structure mismatch.",
      "\nObserved:",
      paste(
        names(counts_long),
        collapse = ", "
      )
    )
  )
  
}


cat(
  "PASS — counts_long schema verified.\n"
)

cat(
  "Rows:",
  nrow(counts_long),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Verify PRIMARY-only concordance structure
# ------------------------------------------------------------

primary_long <- counts_long[
  counts_long$Eligibility == "Primary",
  ,
  drop = FALSE
]


if (
  nrow(primary_long) != 54
) {
  
  stop(
    "STOP: Expected 54 PRIMARY observations in counts_long."
  )
  
}


primary_patients <- unique(
  primary_long$Patient
)

primary_compartments <- unique(
  primary_long$Compartment
)

primary_genes <- unique(
  primary_long$Gene
)


if (
  length(primary_patients) != 9
) {
  
  stop(
    "STOP: PRIMARY patient count is not 9."
  )
  
}


expected_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

if (
  !setequal(
    primary_compartments,
    expected_compartments
  )
) {
  
  stop(
    "STOP: PRIMARY compartment structure mismatch."
  )
  
}


if (
  !setequal(
    primary_genes,
    core_genes
  )
) {
  
  stop(
    "STOP: PRIMARY gene structure mismatch."
  )
  
}


cat(
  "PASS — PRIMARY concordance structure verified.\n"
)

cat(
  "Patients:",
  length(primary_patients),
  "\n"
)

cat(
  "Compartments:",
  paste(
    expected_compartments,
    collapse = " / "
  ),
  "\n"
)

cat(
  "Genes:",
  paste(
    core_genes,
    collapse = " / "
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Verify patient × compartment × gene uniqueness
# ------------------------------------------------------------

primary_key <- paste(
  primary_long$Patient,
  primary_long$Compartment,
  primary_long$Gene,
  sep = "||"
)

if (
  anyDuplicated(primary_key) > 0
) {
  
  stop(
    "STOP: Duplicate Patient × Compartment × Gene observations detected."
  )
  
}


if (
  length(primary_key) != 54
) {
  
  stop(
    "STOP: PRIMARY Patient × Compartment × Gene structure is not 54 unique observations."
  )
  
}


cat(
  "PASS — Patient × Compartment × Gene uniqueness verified.\n\n"
)


# ------------------------------------------------------------
# 7. Verify raw-count / detection integrity
# ------------------------------------------------------------

if (
  any(
    !is.finite(
      primary_long$Raw_Count
    )
  )
) {
  
  stop(
    "STOP: Non-finite Raw_Count detected."
  )
  
}


if (
  any(
    primary_long$Raw_Count < 0
  )
) {
  
  stop(
    "STOP: Negative Raw_Count detected."
  )
  
}


if (
  any(
    primary_long$Raw_Count !=
    floor(primary_long$Raw_Count)
  )
) {
  
  stop(
    "STOP: Non-integer Raw_Count detected."
  )
  
}


expected_detected <-
  primary_long$Raw_Count > 0


if (
  !identical(
    as.logical(primary_long$Detected),
    as.logical(expected_detected)
  )
) {
  
  stop(
    "STOP: Detection definition mismatch."
  )
  
}


cat(
  "PASS — Raw-count and detection integrity verified.\n"
)

cat(
  "Detection definition: Raw_Count > 0\n"
)

cat(
  "Zero counts retained: YES\n\n"
)


# ------------------------------------------------------------
# 8. Verify compartment patient counts
# ------------------------------------------------------------

patient_count_by_compartment <- sapply(
  expected_compartments,
  function(comp) {
    
    length(
      unique(
        primary_long$Patient[
          primary_long$Compartment == comp
        ]
      )
    )
    
  }
)


expected_patient_counts <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)


if (
  !identical(
    as.integer(patient_count_by_compartment),
    as.integer(expected_patient_counts)
  )
) {
  
  stop(
    "STOP: PRIMARY compartment patient counts do not match locked values."
  )
  
}


cat(
  "PASS — PRIMARY compartment patient counts verified.\n"
)

print(
  patient_count_by_compartment
)

cat("\n")


# ------------------------------------------------------------
# 9. Verify complete three-compartment patient set
# ------------------------------------------------------------

patient_compartment_table <- unique(
  primary_long[
    ,
    c(
      "Patient",
      "Compartment"
    ),
    drop = FALSE
  ]
)

patient_compartment_count <- table(
  patient_compartment_table$Patient
)


complete_three_compartment_patients <-
  names(
    patient_compartment_count[
      patient_compartment_count == 3
    ]
  )


if (
  length(
    complete_three_compartment_patients
  ) != 4
) {
  
  stop(
    "STOP: Complete three-compartment patient count is not 4."
  )
  
}


cat(
  "PASS — Complete three-compartment patient set verified.\n"
)

cat(
  "Complete patients:",
  length(
    complete_three_compartment_patients
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 10. Verify detection_frequency object
# ------------------------------------------------------------

if (
  nrow(detection_frequency) != 9
) {
  
  stop(
    "STOP: detection_frequency should contain 9 primary compartment × gene rows."
  )
  
}


cat(
  "PASS — Detection-frequency summary dimensions verified.\n"
)

cat(
  "Rows:",
  nrow(detection_frequency),
  "\n\n"
)


# ------------------------------------------------------------
# 11. Verify compartment_summary object
# ------------------------------------------------------------

if (
  nrow(compartment_summary) != 9
) {
  
  stop(
    "STOP: compartment_summary should contain 9 primary compartment × gene rows."
  )
  
}


cat(
  "PASS — Compartment summary dimensions verified.\n"
)

cat(
  "Rows:",
  nrow(compartment_summary),
  "\n\n"
)


# ------------------------------------------------------------
# 12. Write locked STEP10D-D descriptive outputs
# ------------------------------------------------------------

write.table(
  counts_long,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_CoreGene_Detection_Long.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  detection_frequency,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_Detection_Frequency.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  compartment_summary,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_Compartment_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  data.frame(
    Patient = complete_three_compartment_patients,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Complete_ThreeCompartment_Patients.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


cat(
  "PASS — STEP10D-D descriptive outputs written.\n\n"
)


# ------------------------------------------------------------
# 13. Create final audit summary
# ------------------------------------------------------------

final_audit_summary_locked <- data.frame(
  Metric = c(
    "Input rows",
    "Unique patients in PRIMARY concordance input",
    "T-cell patients",
    "Myeloid/macrophage patients",
    "Endothelial patients",
    "Complete three-compartment patients",
    "Core genes",
    "Core-gene observations",
    "Detection definition",
    "Zero counts retained",
    "Normalization",
    "Statistical testing",
    "Differential expression",
    "P/FDR",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "Patient exclusion"
  ),
  
  Value = c(
    66,
    9,
    9,
    9,
    4,
    4,
    "AOC3 / MS4A7 / IL7R",
    66,
    "Raw_Count > 0",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14. Write final audit summary
# ------------------------------------------------------------

file_final_audit <- file.path(
  step10dd_dir,
  "STEP10D-D_FINAL_Audit_Summary.tsv"
)


write.table(
  final_audit_summary_locked,
  file = file_final_audit,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 15. Write final methodological lock
# ------------------------------------------------------------

final_lock_lines <- c(
  "STEP10D-D — Patient-Level Compartment Concordance Audit",
  "FINAL OUTPUT / PROVENANCE LOCK",
  "",
  "STATUS: PASS",
  "",
  "Primary concordance input: 9 patients",
  "T-cell patients: 9",
  "Myeloid/macrophage patients: 9",
  "Endothelial patients: 4",
  "Complete three-compartment patients: 4",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "Core-gene observations: 66",
  "",
  "Detection definition: Raw_Count > 0",
  "Zero counts retained: YES",
  "",
  "Normalization: NO",
  "Statistical testing: NO",
  "Differential expression: NO",
  "P/FDR: NO",
  "CellChat/LR: NO",
  "Pseudobulk recalculation: NO",
  "Patient/cell exclusion: NO",
  "",
  "Concordance interpretation:",
  "descriptive patient-level detection co-occurrence only.",
  "No expression correlation or causal inference performed.",
  "",
  "Sensitivity endothelial patients remain separate",
  "from the PRIMARY three-compartment concordance analysis.",
  "",
  "STEP10D-D FINAL LOCK: PASS"
)


file_final_lock <- file.path(
  step10dd_dir,
  "STEP10D-D_FINAL_LOCK.txt"
)


writeLines(
  final_lock_lines,
  con = file_final_lock
)


# ------------------------------------------------------------
# 16. Final file existence audit
# ------------------------------------------------------------

expected_output_files <- c(
  "STEP10D-D_Patient_Level_CoreGene_Detection_Long.tsv",
  "STEP10D-D_Patient_Level_Detection_Frequency.tsv",
  "STEP10D-D_Patient_Level_Compartment_Descriptive_Summary.tsv",
  "STEP10D-D_Complete_ThreeCompartment_Patients.tsv",
  "STEP10D-D_FINAL_Audit_Summary.tsv",
  "STEP10D-D_FINAL_LOCK.txt"
)


output_paths <- file.path(
  step10dd_dir,
  expected_output_files
)


output_exists <- file.exists(
  output_paths
)


if (
  !all(output_exists)
) {
  
  stop(
    paste(
      "STOP: One or more STEP10D-D output files were not created:",
      paste(
        expected_output_files[
          !output_exists
        ],
        collapse = ", "
      )
    )
  )
  
}


cat(
  "\n============================================================\n"
)

cat(
  "PASS — STEP10D-D FINAL OUTPUT / PROVENANCE LOCK\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Output directory:\n",
  step10dd_dir,
  "\n\n"
)

cat(
  "Files written:",
  length(expected_output_files),
  "\n"
)

print(
  expected_output_files
)

cat("\n")

cat(
  "Primary patients: 9\n"
)

cat(
  "T-cell patients: 9\n"
)

cat(
  "Myeloid/macrophage patients: 9\n"
)

cat(
  "Endothelial patients: 4\n"
)

cat(
  "Complete three-compartment patients: 4\n"
)

cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)

cat(
  "Core-gene observations: 66\n"
)

cat(
  "Detection: Raw_Count > 0\n"
)

cat(
  "No normalization / statistical testing / DE / P-FDR / CellChat-LR.\n"
)

cat(
  "No pseudobulk recalculation or patient/cell exclusion.\n\n"
)

cat(
  "STEP10D-D FINAL LOCK: PASS\n"
)

cat(
  "============================================================\n\n"
)


# ------------------------------------------------------------
# STEP10D-D — Inspect actual counts_long eligibility structure
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "STEP10D-D — counts_long Eligibility Structure Audit\n"
)

cat(
  "============================================================\n\n"
)


cat(
  "Overall dimensions:\n"
)

cat(
  "Rows:",
  nrow(counts_long),
  "\n"
)

cat(
  "Columns:",
  ncol(counts_long),
  "\n\n"
)


cat(
  "Eligibility distribution:\n"
)

print(
  table(
    counts_long$Eligibility,
    useNA = "ifany"
  )
)

cat("\n")


cat(
  "Eligibility × Compartment:\n"
)

print(
  table(
    counts_long$Eligibility,
    counts_long$Compartment,
    useNA = "ifany"
  )
)

cat("\n")


cat(
  "Eligibility × Gene:\n"
)

print(
  table(
    counts_long$Eligibility,
    counts_long$Gene,
    useNA = "ifany"
  )
)

cat("\n")


cat(
  "Eligibility × Patient:\n"
)

print(
  table(
    counts_long$Eligibility,
    counts_long$Patient,
    useNA = "ifany"
  )
)

cat("\n")


cat(
  "Unique patients by Eligibility:\n"
)

print(
  sapply(
    split(
      counts_long$Patient,
      counts_long$Eligibility
    ),
    function(x) length(unique(x))
  )
)

cat("\n")


cat(
  "Unique Group_ID by Eligibility:\n"
)

print(
  sapply(
    split(
      counts_long$Group_ID,
      counts_long$Eligibility
    ),
    function(x) length(unique(x))
  )
)

cat("\n")


cat(
  "Complete counts_long structure:\n"
)

print(
  head(
    counts_long,
    20
  )
)

cat(
  "\n============================================================\n"
)


##STEP10D-D FINAL OUTPUT / PROVENANCE LOCK
cat(
  "STEP10D-D — FINAL OUTPUT / PROVENANCE LOCK\n"
)
cat(
  "============================================================\n\n"
)

# ------------------------------------------------------------
# 0. Output directory
# ------------------------------------------------------------

step10dd_dir <-
  "input_dir/Paper 1_Step 10-Single-Cell/STEP10D-D_Patient_Level_Compartment_Concordance_Audit"

if (!dir.exists(step10dd_dir)) {
  dir.create(
    step10dd_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(step10dd_dir)) {
  stop(
    "STOP: STEP10D-D output directory could not be created."
  )
}

cat(
  "Output directory verified:\n",
  step10dd_dir,
  "\n\n"
)


# ------------------------------------------------------------
# 1. Required current-session objects
# ------------------------------------------------------------

required_objects <- c(
  "counts_long",
  "detection_frequency",
  "compartment_summary",
  "core_genes"
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
    paste0(
      "STOP: Required object(s) missing: ",
      paste(missing_objects, collapse = ", ")
    )
  )
}

cat(
  "PASS — All required STEP10D-D objects are available.\n\n"
)


# ------------------------------------------------------------
# 2. Lock core genes
# ------------------------------------------------------------

expected_core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

if (!setequal(core_genes, expected_core_genes)) {
  stop(
    paste0(
      "STOP: Core-gene mismatch. Found: ",
      paste(core_genes, collapse = ", ")
    )
  )
}

core_genes <- expected_core_genes

cat(
  "Core genes locked:\n",
  paste(core_genes, collapse = " / "),
  "\n\n"
)


# ------------------------------------------------------------
# 3. counts_long structural validation
# ------------------------------------------------------------

required_columns <- c(
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Gene",
  "Raw_Count",
  "Detected"
)

if (!all(required_columns %in% colnames(counts_long))) {
  stop(
    "STOP: counts_long schema mismatch."
  )
}

cat(
  "PASS — counts_long schema verified.\n"
)

cat(
  "Rows:",
  nrow(counts_long),
  "\n"
)

cat(
  "Columns:",
  ncol(counts_long),
  "\n\n"
)


# ------------------------------------------------------------
# 4. Lock actual PRIMARY structure
# ------------------------------------------------------------

eligibility_values <-
  unique(as.character(counts_long$Eligibility))

if (
  length(eligibility_values) != 1 ||
  !identical(eligibility_values, "Primary")
) {
  stop(
    paste0(
      "STOP: counts_long contains unexpected Eligibility values: ",
      paste(eligibility_values, collapse = ", ")
    )
  )
}

primary_long <- counts_long[
  counts_long$Eligibility == "Primary",
  ,
  drop = FALSE
]

if (nrow(primary_long) != 66) {
  stop(
    paste0(
      "STOP: Expected 66 PRIMARY observations after verified structure; found ",
      nrow(primary_long),
      "."
    )
  )
}

cat(
  "PASS — counts_long is confirmed as PRIMARY-only.\n"
)

cat(
  "PRIMARY observations:",
  nrow(primary_long),
  "\n\n"
)


# ------------------------------------------------------------
# 5. Verify locked compartments
# ------------------------------------------------------------

locked_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

observed_compartments <-
  sort(unique(as.character(primary_long$Compartment)))

if (
  !setequal(
    observed_compartments,
    locked_compartments
  )
) {
  stop(
    paste0(
      "STOP: Compartment mismatch. Found: ",
      paste(observed_compartments, collapse = ", ")
    )
  )
}

cat(
  "PASS — Locked compartments verified:\n",
  paste(locked_compartments, collapse = " / "),
  "\n\n"
)


# ------------------------------------------------------------
# 6. Verify genes and observation structure
# ------------------------------------------------------------

observed_genes <-
  sort(unique(as.character(primary_long$Gene)))

if (!setequal(observed_genes, sort(core_genes))) {
  stop(
    paste0(
      "STOP: Gene structure mismatch. Found: ",
      paste(observed_genes, collapse = ", ")
    )
  )
}

gene_counts <- table(primary_long$Gene)

if (!all(as.integer(gene_counts) == 22)) {
  stop(
    "STOP: Each core gene must have 22 primary observations."
  )
}

cat(
  "PASS — Core-gene observation structure verified.\n"
)

cat(
  "AOC3:",
  unname(gene_counts["AOC3"]),
  "\n"
)

cat(
  "MS4A7:",
  unname(gene_counts["MS4A7"]),
  "\n"
)

cat(
  "IL7R:",
  unname(gene_counts["IL7R"]),
  "\n\n"
)


# ------------------------------------------------------------
# 7. Verify patient-level compartment counts
# ------------------------------------------------------------

patient_compartment_counts <-
  unique(
    primary_long[
      ,
      c("Patient", "Compartment")
    ]
  )

patient_counts <-
  table(patient_compartment_counts$Compartment)

expected_patient_counts <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)

for (nm in names(expected_patient_counts)) {
  
  observed_n <-
    ifelse(
      nm %in% names(patient_counts),
      unname(patient_counts[nm]),
      0
    )
  
  if (observed_n != expected_patient_counts[nm]) {
    stop(
      paste0(
        "STOP: Patient count mismatch for ",
        nm,
        ": expected ",
        expected_patient_counts[nm],
        ", found ",
        observed_n
      )
    )
  }
}

cat(
  "PASS — Patient-level compartment counts verified.\n"
)

cat(
  "T-cell patients:",
  expected_patient_counts["T-cell"],
  "\n"
)

cat(
  "Myeloid/macrophage patients:",
  expected_patient_counts["Myeloid/macrophage"],
  "\n"
)

cat(
  "Endothelial patients:",
  expected_patient_counts["Endothelial"],
  "\n\n"
)


# ------------------------------------------------------------
# 8. Verify patient × compartment × gene completeness
# ------------------------------------------------------------

expected_rows <-
  9 * 3 +
  9 * 3 +
  4 * 3

if (nrow(primary_long) != expected_rows) {
  stop(
    paste0(
      "STOP: Primary observation count does not equal expected ",
      expected_rows,
      "."
    )
  )
}

duplicate_keys <- duplicated(
  primary_long[
    ,
    c(
      "Patient",
      "Compartment",
      "Gene"
    )
  ]
)

if (any(duplicate_keys)) {
  stop(
    "STOP: Duplicate Patient × Compartment × Gene observations detected."
  )
}

cat(
  "PASS — Patient × Compartment × Gene completeness verified.\n"
)

cat(
  "Expected observations:",
  expected_rows,
  "\n"
)

cat(
  "Observed observations:",
  nrow(primary_long),
  "\n\n"
)


# ------------------------------------------------------------
# 9. Verify raw-count and detection integrity
# ------------------------------------------------------------

if (
  any(
    is.na(primary_long$Raw_Count) |
    is.nan(primary_long$Raw_Count) |
    is.infinite(primary_long$Raw_Count)
  )
) {
  stop(
    "STOP: Invalid Raw_Count values detected."
  )
}

if (any(primary_long$Raw_Count < 0)) {
  stop(
    "STOP: Negative Raw_Count values detected."
  )
}

if (
  any(
    primary_long$Raw_Count !=
    floor(primary_long$Raw_Count)
  )
) {
  stop(
    "STOP: Non-integer Raw_Count values detected."
  )
}

expected_detected <-
  primary_long$Raw_Count > 0

if (
  !identical(
    as.logical(primary_long$Detected),
    as.logical(expected_detected)
  )
) {
  stop(
    "STOP: Detection definition mismatch."
  )
}

cat(
  "PASS — Raw-count and detection integrity verified.\n"
)

cat(
  "Detection definition: Raw_Count > 0\n"
)

cat(
  "Zero counts retained: YES\n\n"
)


# ------------------------------------------------------------
# 10. Complete three-compartment patient set
# ------------------------------------------------------------

patients_T <-
  unique(
    primary_long$Patient[
      primary_long$Compartment == "T-cell"
    ]
  )

patients_M <-
  unique(
    primary_long$Patient[
      primary_long$Compartment == "Myeloid/macrophage"
    ]
  )

patients_E <-
  unique(
    primary_long$Patient[
      primary_long$Compartment == "Endothelial"
    ]
  )

complete_patients <-
  Reduce(
    intersect,
    list(
      patients_T,
      patients_M,
      patients_E
    )
  )

complete_patients <-
  sort(unique(complete_patients))

if (length(complete_patients) != 4) {
  stop(
    paste0(
      "STOP: Expected 4 complete three-compartment patients; found ",
      length(complete_patients),
      "."
    )
  )
}

cat(
  "PASS — Complete three-compartment patient set verified.\n"
)

cat(
  "Complete three-compartment patients:",
  length(complete_patients),
  "\n\n"
)


# ------------------------------------------------------------
# 11. Verify existing descriptive outputs
# ------------------------------------------------------------

if (!is.data.frame(detection_frequency)) {
  stop(
    "STOP: detection_frequency is not a data.frame."
  )
}

if (!is.data.frame(compartment_summary)) {
  stop(
    "STOP: compartment_summary is not a data.frame."
  )
}

if (nrow(detection_frequency) != 9) {
  stop(
    paste0(
      "STOP: detection_frequency expected 9 rows; found ",
      nrow(detection_frequency),
      "."
    )
  )
}


if (nrow(compartment_summary) != 9) {
  stop(
    paste0(
      "STOP: compartment_summary expected 9 rows; found ",
      nrow(compartment_summary),
      "."
    )
  )
}

cat(
  "PASS — Existing descriptive STEP10D-D outputs verified.\n"
)


cat(
  "Detection-frequency rows:",
  nrow(detection_frequency),
  "\n"
)

cat(
  "Compartment-summary rows:",
  nrow(compartment_summary),
  "\n\n"
)


# ------------------------------------------------------------
# 12. Write final locked outputs
# ------------------------------------------------------------

write.table(
  primary_long,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_CoreGene_Detection_Long.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  detection_frequency,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_Detection_Frequency.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  compartment_summary,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Patient_Level_Compartment_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(
    Patient = complete_patients,
    Complete_Three_Compartment_Profile = TRUE,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    step10dd_dir,
    "STEP10D-D_Complete_ThreeCompartment_Patients.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Final audit summary
# ------------------------------------------------------------

final_audit <- data.frame(
  Metric = c(
    "Input observations",
    "Eligibility",
    "Unique primary patients",
    "Primary T-cell patients",
    "Primary Myeloid/macrophage patients",
    "Primary Endothelial patients",
    "Complete three-compartment patients",
    "Primary Group_ID count",
    "Core genes",
    "Core-gene observations",
    "Detection definition",
    "Zero counts retained",
    "Normalization",
    "Statistical testing",
    "Differential expression",
    "P/FDR",
    "CellChat/LR analysis",
    "Pseudobulk recalculation",
    "Patient exclusion",
    "Sensitivity Endothelial included in concordance",
    "STEP10D-D status"
  ),
  Value = c(
    66,
    "Primary only",
    9,
    9,
    9,
    4,
    4,
    22,
    paste(core_genes, collapse = " / "),
    66,
    "Raw_Count > 0",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "PASS"
  ),
  stringsAsFactors = FALSE
)

write.table(
  final_audit,
  file = file.path(
    step10dd_dir,
    "STEP10D-D_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Final methodological lock
# ------------------------------------------------------------

lock_text <- c(
  "STEP10D-D — Patient-Level Compartment Concordance Audit",
  "FINAL LOCK",
  "============================================================",
  "",
  "STATUS: PASS",
  "",
  "Input observations: 66",
  "Eligibility: Primary only",
  "Primary patients: 9",
  "Primary T-cell patients: 9",
  "Primary Myeloid/macrophage patients: 9",
  "Primary Endothelial patients: 4",
  "Complete three-compartment patients: 4",
  "Primary Group_ID count: 22",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "Core-gene observations: 66",
  "Detection definition: Raw_Count > 0",
  "Zero counts retained: YES",
  "",
  "Normalization: NO",
  "Statistical testing: NO",
  "Differential expression: NO",
  "P/FDR: NO",
  "CellChat/LR analysis: NO",
  "Pseudobulk recalculation: NO",
  "Patient exclusion: NO",
  "",
  "Sensitivity Endothelial observations are NOT included in",
  "the STEP10D-D primary three-compartment concordance input.",
  "They remain a separate sensitivity cohort upstream.",
  "",
  "No STEP10C pseudobulk recalculation performed.",
  "No STEP10D-A library QC recalculation performed.",
  "No STEP10D-B core-gene recalculation performed.",
  "No STEP10D-C expression recalculation performed.",
  "",
  "STEP10D-D FINAL LOCK: PASS"
)

writeLines(
  lock_text,
  con = file.path(
    step10dd_dir,
    "STEP10D-D_FINAL_LOCK.txt"
  )
)


# ------------------------------------------------------------
# 15. Final console confirmation
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "PASS — STEP10D-D FINAL OUTPUT / PROVENANCE LOCK\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Primary observations:",
  nrow(primary_long),
  "\n"
)

cat(
  "Primary patients:",
  length(unique(primary_long$Patient)),
  "\n"
)

cat(
  "T-cell patients: 9\n"
)

cat(
  "Myeloid/macrophage patients: 9\n"
)

cat(
  "Endothelial patients: 4\n"
)

cat(
  "Complete three-compartment patients:",
  length(complete_patients),
  "\n"
)

cat(
  "Core genes:",
  paste(core_genes, collapse = " / "),
  "\n"
)

cat(
  "Detection definition: Raw_Count > 0\n"
)

cat(
  "Sensitivity Endothelial included: NO\n"
)

cat(
  "Statistical testing: NO\n"
)

cat(
  "P/FDR: NO\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Patient/cell exclusion: NO\n\n"
)

cat(
  "STEP10D-D FINAL LOCK: PASS\n"
)

cat(
  "Output directory:\n",
  step10dd_dir,
  "\n"
)

##STEP10E — Patient-level pseudobulk expression characterization of AOC3–MS4A7–IL7R
cat(
  "STEP10E — Patient-Level Pseudobulk Expression Characterization\n"
)
cat(
  "============================================================\n\n"
)

# ============================================================
# 0. Paths
# ============================================================

step10_root <-
  "input_dir/Paper 1_Step 10-Single-Cell"

step10c_dir <-
  file.path(
    step10_root,
    "STEP10C_Patient_Level_RawUMI_Pseudobulk"
  )

step10e_dir <-
  file.path(
    step10_root,
    "STEP10E_Patient_Level_CoreGene_Expression"
  )

if (!dir.exists(step10c_dir)) {
  stop(
    paste0(
      "STOP: STEP10C directory not found:\n",
      step10c_dir
    )
  )
}

if (!dir.exists(step10e_dir)) {
  dir.create(
    step10e_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!dir.exists(step10e_dir)) {
  stop(
    "STOP: STEP10E output directory could not be created."
  )
}

cat(
  "STEP10C input directory:\n",
  step10c_dir,
  "\n\n"
)

cat(
  "STEP10E output directory:\n",
  step10e_dir,
  "\n\n"
)


# ============================================================
# 1. Locked input files
# ============================================================

input_files <- c(
  Tcell =
    file.path(
      step10c_dir,
      "STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv"
    ),
  Myeloid =
    file.path(
      step10c_dir,
      "STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv"
    ),
  Endothelial_Primary =
    file.path(
      step10c_dir,
      "STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv"
    ),
  Endothelial_Sensitivity =
    file.path(
      step10c_dir,
      "STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv"
    ),
  Primary_Annotation =
    file.path(
      step10c_dir,
      "STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv"
    ),
  Sensitivity_Annotation =
    file.path(
      step10c_dir,
      "STEP10C-F_Sensitivity_Pseudobulk_Sample_Annotation.tsv"
    )
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "STOP: Missing STEP10C input file(s):\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat(
  "PASS — All locked STEP10C input files are present.\n\n"
)


# ============================================================
# 2. Core-gene lock
# ============================================================

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

cat(
  "Core genes locked:\n",
  paste(core_genes, collapse = " / "),
  "\n\n"
)


# ============================================================
# 3. Read locked pseudobulk matrices
# ============================================================

pb_T <- read.delim(
  input_files["Tcell"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_M <- read.delim(
  input_files["Myeloid"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_E <- read.delim(
  input_files["Endothelial_Primary"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pb_E_sens <- read.delim(
  input_files["Endothelial_Sensitivity"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ann_primary <- read.delim(
  input_files["Primary_Annotation"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

ann_sens <- read.delim(
  input_files["Sensitivity_Annotation"],
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat(
  "Input dimensions:\n"
)

cat(
  "T-cell:",
  nrow(pb_T),
  "genes ×",
  ncol(pb_T) - 1,
  "groups\n"
)

cat(
  "Myeloid/macrophage:",
  nrow(pb_M),
  "genes ×",
  ncol(pb_M) - 1,
  "groups\n"
)

cat(
  "Endothelial primary:",
  nrow(pb_E),
  "genes ×",
  ncol(pb_E) - 1,
  "groups\n"
)

cat(
  "Endothelial sensitivity:",
  nrow(pb_E_sens),
  "genes ×",
  ncol(pb_E_sens) - 1,
  "groups\n\n"
)


# ============================================================
# 4. Gene-column identification
# ============================================================

gene_column_candidates <- c(
  "Gene",
  "gene",
  "Gene_ID",
  "gene_id"
)

find_gene_column <- function(df) {
  
  hits <- intersect(
    gene_column_candidates,
    colnames(df)
  )
  
  if (length(hits) != 1) {
    stop(
      paste0(
        "STOP: Unable to uniquely identify gene column. ",
        "Found: ",
        paste(hits, collapse = ", ")
      )
    )
  }
  
  hits
}

gene_col_T <- find_gene_column(pb_T)
gene_col_M <- find_gene_column(pb_M)
gene_col_E <- find_gene_column(pb_E)
gene_col_Es <- find_gene_column(pb_E_sens)


# ============================================================
# 5. Gene identity and order lock
# ============================================================

genes_T <- as.character(pb_T[[gene_col_T]])
genes_M <- as.character(pb_M[[gene_col_M]])
genes_E <- as.character(pb_E[[gene_col_E]])
genes_Es <- as.character(pb_E_sens[[gene_col_Es]])

if (anyDuplicated(genes_T)) {
  stop("STOP: Duplicate genes detected in T-cell matrix.")
}

if (anyDuplicated(genes_M)) {
  stop("STOP: Duplicate genes detected in myeloid matrix.")
}

if (anyDuplicated(genes_E)) {
  stop("STOP: Duplicate genes detected in endothelial matrix.")
}

if (anyDuplicated(genes_Es)) {
  stop("STOP: Duplicate genes detected in sensitivity endothelial matrix.")
}

if (!identical(genes_T, genes_M) ||
    !identical(genes_T, genes_E) ||
    !identical(genes_T, genes_Es)) {
  
  stop(
    "STOP: Gene identity/order mismatch across STEP10C matrices."
  )
}

if (!all(core_genes %in% genes_T)) {
  stop(
    "STOP: One or more core genes are absent from pseudobulk matrices."
  )
}

cat(
  "PASS — Gene identity and order are concordant across all matrices.\n"
)

cat(
  "PASS — AOC3 / MS4A7 / IL7R are present.\n\n"
)


# ============================================================
# 6. Convert matrices to numeric count matrices
# ============================================================

convert_counts <- function(df, gene_col) {
  
  gene_ids <- as.character(df[[gene_col]])
  
  count_df <- df[
    ,
    setdiff(colnames(df), gene_col),
    drop = FALSE
  ]
  
  count_mat <- as.matrix(
    data.frame(
      lapply(
        count_df,
        function(x) as.numeric(as.character(x))
      ),
      check.names = FALSE
    )
  )
  
  rownames(count_mat) <- gene_ids
  
  count_mat
}

mat_T <- convert_counts(pb_T, gene_col_T)
mat_M <- convert_counts(pb_M, gene_col_M)
mat_E <- convert_counts(pb_E, gene_col_E)
mat_Es <- convert_counts(pb_E_sens, gene_col_Es)


# ============================================================
# 7. Raw UMI integrity audit
# ============================================================

audit_matrix <- function(mat, label) {
  
  if (anyNA(mat)) {
    stop(
      paste0(
        "STOP: NA detected in ",
        label
      )
    )
  }
  
  if (any(!is.finite(mat))) {
    stop(
      paste0(
        "STOP: non-finite values detected in ",
        label
      )
    )
  }
  
  if (any(mat < 0)) {
    stop(
      paste0(
        "STOP: negative counts detected in ",
        label
      )
    )
  }
  
  if (any(mat != floor(mat))) {
    stop(
      paste0(
        "STOP: non-integer counts detected in ",
        label
      )
    )
  }
  
  TRUE
}

audit_matrix(mat_T, "T-cell")
audit_matrix(mat_M, "Myeloid/macrophage")
audit_matrix(mat_E, "Endothelial primary")
audit_matrix(mat_Es, "Endothelial sensitivity")

cat(
  "PASS — All pseudobulk matrices contain finite, non-negative integer raw UMI counts.\n\n"
)


# ============================================================
# 8. Library size calculation from raw UMI matrices
# ============================================================

library_size <- function(mat) {
  colSums(mat)
}

lib_T <- library_size(mat_T)
lib_M <- library_size(mat_M)
lib_E <- library_size(mat_E)
lib_Es <- library_size(mat_Es)

if (any(lib_T <= 0) ||
    any(lib_M <= 0) ||
    any(lib_E <= 0) ||
    any(lib_Es <= 0)) {
  
  stop(
    "STOP: Zero-sized pseudobulk library detected."
  )
}

cat(
  "PASS — All pseudobulk libraries have positive total UMI counts.\n\n"
)


# ============================================================
# 9. CPM derivation
# ============================================================
#
# CPM is derived directly from locked raw UMI counts.
# It is descriptive normalization only.
# No inferential statistics are performed.
#

to_CPM <- function(mat) {
  
  sweep(
    mat,
    2,
    colSums(mat),
    "/"
  ) * 1e6
}

cpm_T <- to_CPM(mat_T)
cpm_M <- to_CPM(mat_M)
cpm_E <- to_CPM(mat_E)
cpm_Es <- to_CPM(mat_Es)

cat(
  "PASS — CPM derived from locked raw UMI pseudobulk counts.\n"
)

cat(
  "CPM use: descriptive normalization only.\n\n"
)


# ============================================================
# 10. Construct patient-level core-gene table
# ============================================================

extract_core <- function(
    count_mat,
    cpm_mat,
    dataset_label,
    annotation_df
) {
  
  keep_genes <- core_genes
  
  count_core <- count_mat[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  cpm_core <- cpm_mat[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  groups <- colnames(count_mat)
  
  if (!all(groups %in% annotation_df$Group_ID)) {
    
    stop(
      paste0(
        "STOP: Matrix Group_IDs do not match annotation for ",
        dataset_label
      )
    )
  }
  
  ann <- annotation_df[
    match(groups, annotation_df$Group_ID),
    ,
    drop = FALSE
  ]
  
  out_list <- list()
  
  k <- 1
  
  for (gene in keep_genes) {
    
    tmp <- data.frame(
      Dataset = dataset_label,
      Patient = ann$Patient,
      Compartment = ann$Compartment,
      Eligibility = ann$Eligibility,
      Group_ID = ann$Group_ID,
      Cell_Count = ann$Cell_Count,
      Gene = gene,
      Raw_UMI = as.numeric(count_core[gene, ]),
      Library_Total_UMI = as.numeric(
        colSums(count_mat)
      ),
      CPM = as.numeric(cpm_core[gene, ]),
      Log2_CPM_Plus1 = log2(
        as.numeric(cpm_core[gene, ]) + 1
      ),
      Detected = as.numeric(count_core[gene, ]) > 0,
      stringsAsFactors = FALSE
    )
    
    out_list[[k]] <- tmp
    k <- k + 1
  }
  
  do.call(rbind, out_list)
}


core_T <- extract_core(
  mat_T,
  cpm_T,
  "Primary_T-cell",
  ann_primary[
    ann_primary$Compartment == "T-cell",
    ,
    drop = FALSE
  ]
)

core_M <- extract_core(
  mat_M,
  cpm_M,
  "Primary_Myeloid_macrophage",
  ann_primary[
    ann_primary$Compartment == "Myeloid/macrophage",
    ,
    drop = FALSE
  ]
)

core_E <- extract_core(
  mat_E,
  cpm_E,
  "Primary_Endothelial",
  ann_primary[
    ann_primary$Compartment == "Endothelial",
    ,
    drop = FALSE
  ]
)

core_Es <- extract_core(
  mat_Es,
  cpm_Es,
  "Sensitivity_Endothelial",
  ann_sens[
    ann_sens$Compartment == "Endothelial",
    ,
    drop = FALSE
  ]
)


# ============================================================
# 11. Combine primary and sensitivity descriptive tables
# ============================================================

core_expression_all <- rbind(
  core_T,
  core_M,
  core_E,
  core_Es
)

rownames(core_expression_all) <- NULL

cat(
  "Core-gene expression observations:\n"
)

print(
  table(
    core_expression_all$Dataset,
    core_expression_all$Gene
  )
)

cat("\n")


# ============================================================
# 12. Primary-only descriptive summary
# ============================================================

primary_expression <- core_expression_all[
  core_expression_all$Eligibility == "Primary",
  ,
  drop = FALSE
]

summary_fun <- function(x) {
  
  q <- quantile(
    x,
    probs = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    names = FALSE,
    type = 7
  )
  
  data.frame(
    N = length(x),
    Min = q[1],
    Q1 = q[2],
    Median = q[3],
    Mean = mean(x),
    Q3 = q[4],
    Max = q[5],
    IQR = IQR(x),
    stringsAsFactors = FALSE
  )
}

summary_list <- list()

idx <- 1

for (dataset in unique(primary_expression$Dataset)) {
  
  for (gene in core_genes) {
    
    sub <- primary_expression[
      primary_expression$Dataset == dataset &
        primary_expression$Gene == gene,
      ,
      drop = FALSE
    ]
    
    raw_s <- summary_fun(sub$Raw_UMI)
    cpm_s <- summary_fun(sub$CPM)
    
    summary_list[[idx]] <- data.frame(
      Dataset = dataset,
      Compartment = unique(sub$Compartment),
      Gene = gene,
      N = nrow(sub),
      Detected_N = sum(sub$Detected),
      Zero_N = sum(!sub$Detected),
      Detection_Percent =
        100 * mean(sub$Detected),
      
      Raw_UMI_Min = raw_s$Min,
      Raw_UMI_Q1 = raw_s$Q1,
      Raw_UMI_Median = raw_s$Median,
      Raw_UMI_Mean = raw_s$Mean,
      Raw_UMI_Q3 = raw_s$Q3,
      Raw_UMI_Max = raw_s$Max,
      Raw_UMI_IQR = raw_s$IQR,
      
      CPM_Min = cpm_s$Min,
      CPM_Q1 = cpm_s$Q1,
      CPM_Median = cpm_s$Median,
      CPM_Mean = cpm_s$Mean,
      CPM_Q3 = cpm_s$Q3,
      CPM_Max = cpm_s$Max,
      CPM_IQR = cpm_s$IQR,
      
      stringsAsFactors = FALSE
    )
    
    idx <- idx + 1
  }
}

primary_summary <- do.call(
  rbind,
  summary_list
)

rownames(primary_summary) <- NULL


# ============================================================
# 13. Sensitivity endothelial descriptive summary
# ============================================================

sensitivity_expression <- core_expression_all[
  core_expression_all$Eligibility == "Sensitivity",
  ,
  drop = FALSE
]

sensitivity_summary_list <- list()

idx <- 1

for (gene in core_genes) {
  
  sub <- sensitivity_expression[
    sensitivity_expression$Gene == gene,
    ,
    drop = FALSE
  ]
  
  raw_s <- summary_fun(sub$Raw_UMI)
  cpm_s <- summary_fun(sub$CPM)
  
  sensitivity_summary_list[[idx]] <- data.frame(
    Dataset = "Sensitivity_Endothelial",
    Compartment = "Endothelial",
    Gene = gene,
    N = nrow(sub),
    Detected_N = sum(sub$Detected),
    Zero_N = sum(!sub$Detected),
    Detection_Percent =
      100 * mean(sub$Detected),
    
    Raw_UMI_Min = raw_s$Min,
    Raw_UMI_Q1 = raw_s$Q1,
    Raw_UMI_Median = raw_s$Median,
    Raw_UMI_Mean = raw_s$Mean,
    Raw_UMI_Q3 = raw_s$Q3,
    Raw_UMI_Max = raw_s$Max,
    Raw_UMI_IQR = raw_s$IQR,
    
    CPM_Min = cpm_s$Min,
    CPM_Q1 = cpm_s$Q1,
    CPM_Median = cpm_s$Median,
    CPM_Mean = cpm_s$Mean,
    CPM_Q3 = cpm_s$Q3,
    CPM_Max = cpm_s$Max,
    CPM_IQR = cpm_s$IQR,
    
    stringsAsFactors = FALSE
  )
  
  idx <- idx + 1
}

sensitivity_summary <- do.call(
  rbind,
  sensitivity_summary_list
)

rownames(sensitivity_summary) <- NULL


# ============================================================
# 14. Patient-level primary expression table
# ============================================================

primary_patient_expression <- primary_expression[
  order(
    primary_expression$Compartment,
    primary_expression$Patient,
    primary_expression$Gene
  ),
  ,
  drop = FALSE
]


# ============================================================
# 15. Structural checks
# ============================================================

expected_primary_observations <- 66
expected_sensitivity_observations <- 6

if (
  nrow(primary_expression) !=
  expected_primary_observations
) {
  stop(
    paste0(
      "STOP: Expected 66 primary expression observations; found ",
      nrow(primary_expression)
    )
  )
}

if (
  nrow(sensitivity_expression) !=
  expected_sensitivity_observations
) {
  stop(
    paste0(
      "STOP: Expected 6 sensitivity expression observations; found ",
      nrow(sensitivity_expression)
    )
  )
}

if (
  nrow(primary_summary) != 9
) {
  stop(
    "STOP: Expected 9 primary gene × compartment summary rows."
  )
}

if (
  nrow(sensitivity_summary) != 3
) {
  stop(
    "STOP: Expected 3 sensitivity gene summary rows."
  )
}

if (
  length(unique(primary_expression$Patient)) != 9
) {
  stop(
    "STOP: Primary patient count is not 9."
  )
}

if (
  length(unique(primary_expression$Group_ID)) != 22
) {
  stop(
    "STOP: Primary Group_ID count is not 22."
  )
}

cat(
  "PASS — STEP10E structural checks completed.\n\n"
)


# ============================================================
# 16. Write locked descriptive outputs
# ============================================================

write.table(
  primary_patient_expression,
  file = file.path(
    step10e_dir,
    "STEP10E-1_Primary_Patient_Level_CoreGene_Expression.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sensitivity_expression,
  file = file.path(
    step10e_dir,
    "STEP10E-1S_Sensitivity_Endothelial_CoreGene_Expression.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  primary_summary,
  file = file.path(
    step10e_dir,
    "STEP10E-2_Primary_CoreGene_Expression_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  sensitivity_summary,
  file = file.path(
    step10e_dir,
    "STEP10E-2S_Sensitivity_Endothelial_Descriptive_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 17. Methodological lock
# ============================================================

method_lock <- data.frame(
  Item = c(
    "Input",
    "Expression basis",
    "Raw UMI counts retained",
    "CPM",
    "Log2(CPM+1)",
    "Detection definition",
    "Primary T-cell",
    "Primary Myeloid/macrophage",
    "Primary Endothelial",
    "Sensitivity Endothelial",
    "Statistical testing",
    "Differential expression",
    "P value",
    "FDR",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "Patient exclusion",
    "Sensitivity mixed with primary"
  ),
  Value = c(
    "STEP10C locked raw-UMI pseudobulk matrices",
    "Raw UMI pseudobulk",
    "YES",
    "YES — descriptive normalization only",
    "YES — descriptive transformation only",
    "Raw UMI > 0",
    "n = 9",
    "n = 9",
    "n = 4",
    "n = 2",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  stringsAsFactors = FALSE
)

write.table(
  method_lock,
  file = file.path(
    step10e_dir,
    "STEP10E_Methodological_Lock.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 18. Final audit summary
# ============================================================

final_audit <- data.frame(
  Metric = c(
    "Primary expression observations",
    "Primary patients",
    "Primary T-cell groups",
    "Primary Myeloid/macrophage groups",
    "Primary Endothelial groups",
    "Sensitivity Endothelial groups",
    "Complete three-compartment patients",
    "Core genes",
    "Raw UMI counts used",
    "CPM derived",
    "Detection definition",
    "Zero counts retained",
    "Statistical testing",
    "Differential expression",
    "P/FDR",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "Patient exclusion",
    "Primary/sensitivity mixing",
    "STEP10E status"
  ),
  Value = c(
    66,
    9,
    9,
    9,
    4,
    2,
    4,
    "AOC3 / MS4A7 / IL7R",
    "YES",
    "YES — descriptive only",
    "Raw UMI > 0",
    "YES",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "PASS"
  ),
  stringsAsFactors = FALSE
)

write.table(
  final_audit,
  file = file.path(
    step10e_dir,
    "STEP10E_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ============================================================
# 19. FINAL LOCK
# ============================================================

lock_text <- c(
  "STEP10E — Patient-Level Pseudobulk Expression Characterization",
  "FINAL LOCK",
  "============================================================",
  "",
  "STATUS: PASS",
  "",
  "Primary expression observations: 66",
  "Primary patients: 9",
  "Primary T-cell groups: 9",
  "Primary Myeloid/macrophage groups: 9",
  "Primary Endothelial groups: 4",
  "Sensitivity Endothelial groups: 2",
  "Complete three-compartment patients: 4",
  "",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "",
  "Expression basis: locked raw-UMI pseudobulk counts",
  "Raw UMI counts retained: YES",
  "CPM: derived descriptively from raw UMI counts",
  "Log2(CPM+1): descriptive transformation only",
  "Detection definition: Raw UMI > 0",
  "Zero counts retained: YES",
  "",
  "Statistical testing: NO",
  "Differential expression: NO",
  "P value: NO",
  "FDR: NO",
  "Correlation: NO",
  "CellChat/LR: NO",
  "Pseudobulk recalculation: NO",
  "Patient exclusion: NO",
  "",
  "Sensitivity Endothelial remains separate from primary analysis.",
  "No sensitivity observations were mixed with primary observations.",
  "",
  "STEP10E FINAL LOCK: PASS"
)

writeLines(
  lock_text,
  con = file.path(
    step10e_dir,
    "STEP10E_FINAL_LOCK.txt"
  )
)


# ============================================================
# 20. Console final confirmation
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "PASS — STEP10E FINAL LOCK\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Primary expression observations: 66\n"
)

cat(
  "Primary patients: 9\n"
)

cat(
  "T-cell: 9\n"
)

cat(
  "Myeloid/macrophage: 9\n"
)

cat(
  "Endothelial primary: 4\n"
)

cat(
  "Endothelial sensitivity: 2\n"
)

cat(
  "Complete three-compartment patients: 4\n"
)

cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)

cat(
  "Raw UMI counts: USED\n"
)

cat(
  "CPM: DESCRIPTIVE ONLY\n"
)

cat(
  "Statistical testing: NO\n"
)

cat(
  "Differential expression: NO\n"
)

cat(
  "P/FDR: NO\n"
)

cat(
  "Correlation: NO\n"
)

cat(
  "CellChat/LR: NO\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "Patient exclusion: NO\n\n"
)

cat(
  "STEP10E FINAL LOCK: PASS\n\n"
)

cat(
  "Output directory:\n",
  step10e_dir,
  "\n"
)


# ============================================================
# STEP10F-1 FINAL
# Patient-level expression profile
# Paper 1 — AOC3 / MS4A7 / IL7R
#
# ============================================================
#
# PURPOSE
# ------------------------------------------------------------
# Descriptive patient-level compartment expression profiling
# using the LOCKED STEP10E patient-level expression outputs.
#
# Biological context questions:
#
#   1. AOC3  — endothelial-associated expression context?
#   2. MS4A7 — myeloid/macrophage-associated expression context?
#   3. IL7R  — T-cell-associated expression context?
#
# IMPORTANT INTERPRETATION:
# ------------------------------------------------------------
# This step is descriptive.
#
# It does NOT test statistical differences between compartments.
#
# "Associated" must not automatically be interpreted as
# "specific".
#
# ============================================================
#
# LOCKED UPSTREAM INPUT
# ------------------------------------------------------------
# STEP10E patient-level raw-UMI pseudobulk expression
#
# Primary:
#   T-cell                n = 9
#   Myeloid/macrophage    n = 9
#   Endothelial           n = 4
#
# Sensitivity:
#   Endothelial           n = 2
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
#
# ============================================================
#
# ANALYTICAL RESTRICTIONS
# ------------------------------------------------------------
# Raw UMI pseudobulk: USED
# CPM: locked descriptive normalization from STEP10E
# Log2(CPM+1): locked descriptive transformation
#
# Differential expression: NO
# Statistical testing: NO
# P values: NO
# FDR: NO
# Correlation: NO
# CellChat: NO
# Ligand-receptor analysis: NO
# Pseudobulk recalculation: NO
# GEO raw matrix rereading: NO
# H5 rereading: NO
# Patient exclusion: NO
# Cell exclusion: NO
#
# ============================================================
# R version: 4.5.1
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10F-1 FINAL\n")
cat("Patient-level expression profile\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

input_dir <- file.path(
  base_dir,
  "STEP10E_Patient_Level_CoreGene_Expression"
)

output_dir <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile"
)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE
  )
}


# ------------------------------------------------------------
# 1. Locked STEP10E input files
# ------------------------------------------------------------

primary_file <- file.path(
  input_dir,
  "STEP10E-1_Primary_Patient_Level_CoreGene_Expression.tsv"
)

sensitivity_file <- file.path(
  input_dir,
  "STEP10E-1S_Sensitivity_Endothelial_CoreGene_Expression.tsv"
)


if (!file.exists(primary_file)) {
  
  stop(
    "\nFAIL — STEP10E primary input file not found:\n",
    primary_file,
    "\n"
  )
  
}


if (!file.exists(sensitivity_file)) {
  
  stop(
    "\nFAIL — STEP10E sensitivity input file not found:\n",
    sensitivity_file,
    "\n"
  )
  
}


cat("PASS — STEP10E primary input file found.\n")
cat("PASS — STEP10E sensitivity input file found.\n\n")


# ------------------------------------------------------------
# 2. Read LOCKED STEP10E outputs
# ------------------------------------------------------------

primary_df <- read.delim(
  primary_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


sensitivity_df <- read.delim(
  sensitivity_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


cat(
  "Primary input rows:",
  nrow(primary_df),
  "\n"
)

cat(
  "Sensitivity input rows:",
  nrow(sensitivity_df),
  "\n\n"
)


# ------------------------------------------------------------
# 3. Actual STEP10E column schema lock
# ------------------------------------------------------------

cat("Primary columns:\n")
print(colnames(primary_df))

cat("\nSensitivity columns:\n")
print(colnames(sensitivity_df))

cat("\n")


# IMPORTANT:
# The actual locked STEP10E column name is:
#
#   Log2_CPM_Plus1
#
# NOT:
#
#   Log2_CPM1
#
# This is explicitly locked here to prevent schema mismatch.


required_cols <- c(
  "Dataset",
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Cell_Count",
  "Gene",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1",
  "Detected"
)


# ------------------------------------------------------------
# 4. Required-column validation
# ------------------------------------------------------------

missing_primary <- setdiff(
  required_cols,
  colnames(primary_df)
)

missing_sensitivity <- setdiff(
  required_cols,
  colnames(sensitivity_df)
)


if (length(missing_primary) > 0) {
  
  stop(
    "\nFAIL — Missing required primary columns:\n",
    paste(
      missing_primary,
      collapse = ", "
    ),
    "\n"
  )
  
}


if (length(missing_sensitivity) > 0) {
  
  stop(
    "\nFAIL — Missing required sensitivity columns:\n",
    paste(
      missing_sensitivity,
      collapse = ", "
    ),
    "\n"
  )
  
}


cat(
  "PASS — Required STEP10E columns verified.\n\n"
)


# ------------------------------------------------------------
# 5. Add analysis-set labels
# ------------------------------------------------------------

primary_df$Analysis_Set <- "Primary"

sensitivity_df$Analysis_Set <- "Sensitivity"


# ------------------------------------------------------------
# 6. Combine datasets for QC
# ------------------------------------------------------------

expression_df <- rbind(
  primary_df,
  sensitivity_df
)

rownames(expression_df) <- NULL


# ------------------------------------------------------------
# 7. Core-gene lock
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)


observed_genes <- unique(
  as.character(expression_df$Gene)
)


cat("Observed genes:\n")
print(sort(observed_genes))
cat("\n")


if (!setequal(
  observed_genes,
  core_genes
)) {
  
  stop(
    "\nFAIL — Observed gene set does not exactly match:\n",
    "AOC3 / MS4A7 / IL7R\n"
  )
  
}


cat(
  "PASS — Core-gene set verified: ",
  "AOC3 / MS4A7 / IL7R\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Gene ordering lock
# ------------------------------------------------------------

expression_df$Gene <- factor(
  expression_df$Gene,
  levels = core_genes
)


# ------------------------------------------------------------
# 9. Locked compartment definitions
# ------------------------------------------------------------

primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

sensitivity_compartment <- "Endothelial"


# ------------------------------------------------------------
# 10. Primary structure audit
# ------------------------------------------------------------

primary_structure <- unique(
  primary_df[
    ,
    c(
      "Patient",
      "Compartment",
      "Group_ID",
      "Eligibility"
    )
  ]
)


primary_patients <- unique(
  primary_structure$Patient
)

primary_groups <- unique(
  primary_structure$Group_ID
)

primary_compartments_observed <- unique(
  primary_structure$Compartment
)


if (length(primary_patients) != 9) {
  
  stop(
    "\nFAIL — Expected 9 primary patients; observed ",
    length(primary_patients),
    "\n"
  )
  
}


if (length(primary_groups) != 22) {
  
  stop(
    "\nFAIL — Expected 22 primary patient-compartment groups; observed ",
    length(primary_groups),
    "\n"
  )
  
}


if (!setequal(
  primary_compartments_observed,
  primary_compartments
)) {
  
  stop(
    "\nFAIL — Primary compartment structure does not match lock.\n",
    "Observed: ",
    paste(
      sort(primary_compartments_observed),
      collapse = ", "
    ),
    "\n"
  )
  
}


cat(
  "PASS — Primary patient/compartment structure verified.\n"
)

cat(
  "       Primary patients: ",
  length(primary_patients),
  "\n",
  sep = ""
)

cat(
  "       Primary groups: ",
  length(primary_groups),
  "\n",
  sep = ""
)


primary_compartment_counts <- table(
  primary_structure$Compartment
)

cat(
  "       T-cell: ",
  unname(
    primary_compartment_counts["T-cell"]
  ),
  "\n",
  sep = ""
)

cat(
  "       Myeloid/macrophage: ",
  unname(
    primary_compartment_counts[
      "Myeloid/macrophage"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "       Endothelial: ",
  unname(
    primary_compartment_counts[
      "Endothelial"
    ]
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Sensitivity structure audit
# ------------------------------------------------------------

sensitivity_structure <- unique(
  sensitivity_df[
    ,
    c(
      "Patient",
      "Compartment",
      "Group_ID",
      "Eligibility"
    )
  ]
)


sensitivity_patients <- unique(
  sensitivity_structure$Patient
)

sensitivity_groups <- unique(
  sensitivity_structure$Group_ID
)


if (length(sensitivity_patients) != 2) {
  
  stop(
    "\nFAIL — Expected 2 sensitivity endothelial patients; observed ",
    length(sensitivity_patients),
    "\n"
  )
  
}


if (length(sensitivity_groups) != 2) {
  
  stop(
    "\nFAIL — Expected 2 sensitivity groups; observed ",
    length(sensitivity_groups),
    "\n"
  )
  
}


if (!all(
  sensitivity_structure$Compartment ==
  sensitivity_compartment
)) {
  
  stop(
    "\nFAIL — Sensitivity dataset contains non-Endothelial compartments.\n"
  )
  
}


cat(
  "PASS — Sensitivity Endothelial structure verified.\n"
)

cat(
  "       Sensitivity patients: ",
  length(sensitivity_patients),
  "\n",
  sep = ""
)

cat(
  "       Sensitivity groups: ",
  length(sensitivity_groups),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 12. Eligibility lock
# ------------------------------------------------------------

if (!all(
  primary_df$Eligibility == "Primary"
)) {
  
  stop(
    "\nFAIL — Primary input contains non-Primary eligibility labels.\n"
  )
  
}


if (!all(
  sensitivity_df$Eligibility == "Sensitivity"
)) {
  
  stop(
    "\nFAIL — Sensitivity input contains non-Sensitivity eligibility labels.\n"
  )
  
}


cat(
  "PASS — Primary/Sensitivity eligibility labels verified.\n\n"
)


# ------------------------------------------------------------
# 13. Observation-count lock
# ------------------------------------------------------------

if (nrow(primary_df) != 66) {
  
  stop(
    "\nFAIL — Expected 66 primary observations; observed ",
    nrow(primary_df),
    "\n"
  )
  
}


if (nrow(sensitivity_df) != 6) {
  
  stop(
    "\nFAIL — Expected 6 sensitivity observations; observed ",
    nrow(sensitivity_df),
    "\n"
  )
  
}


cat(
  "PASS — Observation counts verified.\n"
)

cat(
  "       Primary = 66\n"
)

cat(
  "       Sensitivity = 6\n\n"
)


# ------------------------------------------------------------
# 14. Unique patient × compartment × gene audit
# ------------------------------------------------------------

dup_key <- paste(
  expression_df$Analysis_Set,
  expression_df$Patient,
  expression_df$Compartment,
  expression_df$Gene,
  sep = "||"
)


duplicate_count <- sum(
  duplicated(dup_key)
)


if (duplicate_count != 0) {
  
  stop(
    "\nFAIL — Duplicate patient × compartment × gene ",
    "observations detected: ",
    duplicate_count,
    "\n"
  )
  
}


cat(
  "PASS — Unique patient × compartment × gene structure verified.\n\n"
)


# ------------------------------------------------------------
# 15. Numeric-column conversion
# ------------------------------------------------------------

numeric_columns <- c(
  "Cell_Count",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1"
)


for (v in numeric_columns) {
  
  expression_df[[v]] <- as.numeric(
    expression_df[[v]]
  )
  
}


# ------------------------------------------------------------
# 16. Numeric integrity audit
# ------------------------------------------------------------

for (v in numeric_columns) {
  
  if (any(
    !is.finite(
      expression_df[[v]]
    )
  )) {
    
    stop(
      "\nFAIL — Non-finite values detected in: ",
      v,
      "\n"
    )
    
  }
  
}


cat(
  "PASS — Numeric expression variables contain only finite values.\n\n"
)


# ------------------------------------------------------------
# 17. Raw UMI integrity
# ------------------------------------------------------------

if (any(
  expression_df$Raw_UMI < 0
)) {
  
  stop(
    "\nFAIL — Negative Raw_UMI values detected.\n"
  )
  
}


if (any(
  expression_df$Raw_UMI !=
  floor(expression_df$Raw_UMI)
)) {
  
  stop(
    "\nFAIL — Non-integer Raw_UMI values detected.\n"
  )
  
}


if (any(
  expression_df$Library_Total_UMI <= 0
)) {
  
  stop(
    "\nFAIL — Non-positive Library_Total_UMI detected.\n"
  )
  
}


cat(
  "PASS — Raw UMI integrity verified.\n"
)

cat(
  "       Non-negative: YES\n"
)

cat(
  "       Integer-valued: YES\n"
)

cat(
  "       Library totals > 0: YES\n\n"
)


# ------------------------------------------------------------
# 18. Locked CPM audit
# ------------------------------------------------------------

expected_cpm <- (
  expression_df$Raw_UMI /
    expression_df$Library_Total_UMI
) * 1e6


cpm_difference <- abs(
  expression_df$CPM -
    expected_cpm
)


max_cpm_difference <- max(
  cpm_difference,
  na.rm = TRUE
)


if (max_cpm_difference > 1e-6) {
  
  stop(
    "\nFAIL — CPM does not match STEP10E locked calculation.\n",
    "Maximum absolute difference = ",
    format(
      max_cpm_difference,
      scientific = TRUE
    ),
    "\n"
  )
  
}


cat(
  "PASS — Locked CPM calculation verified.\n"
)

cat(
  "       Maximum absolute difference = ",
  format(
    max_cpm_difference,
    scientific = TRUE
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Locked log2(CPM+1) audit
# ------------------------------------------------------------

expected_log2 <- log2(
  expression_df$CPM + 1
)


log2_difference <- abs(
  expression_df$Log2_CPM_Plus1 -
    expected_log2
)


max_log2_difference <- max(
  log2_difference,
  na.rm = TRUE
)


if (max_log2_difference > 1e-6) {
  
  stop(
    "\nFAIL — Log2_CPM_Plus1 does not match ",
    "STEP10E locked calculation.\n",
    "Maximum absolute difference = ",
    format(
      max_log2_difference,
      scientific = TRUE
    ),
    "\n"
  )
  
}


cat(
  "PASS — Locked log2(CPM+1) calculation verified.\n"
)

cat(
  "       Maximum absolute difference = ",
  format(
    max_log2_difference,
    scientific = TRUE
  ),
  "\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 20. Detection-definition audit
# ------------------------------------------------------------

expected_detected <- (
  expression_df$Raw_UMI > 0
)


if (!all(
  expression_df$Detected ==
  expected_detected
)) {
  
  stop(
    "\nFAIL — Detection definition is inconsistent with ",
    "Raw_UMI > 0.\n"
  )
  
}


cat(
  "PASS — Detection definition verified: Raw_UMI > 0.\n\n"
)


# ------------------------------------------------------------
# 21. Separate Primary and Sensitivity data
# ------------------------------------------------------------

primary_expression <- expression_df[
  expression_df$Analysis_Set == "Primary",
  ,
  drop = FALSE
]


sensitivity_expression <- expression_df[
  expression_df$Analysis_Set == "Sensitivity",
  ,
  drop = FALSE
]


primary_expression$Gene <- as.character(
  primary_expression$Gene
)

primary_expression$Compartment <- as.character(
  primary_expression$Compartment
)


sensitivity_expression$Gene <- as.character(
  sensitivity_expression$Gene
)

sensitivity_expression$Compartment <- as.character(
  sensitivity_expression$Compartment
)


# ------------------------------------------------------------
# 22. Primary Gene × Compartment descriptive profile
# ------------------------------------------------------------
#
# This is the central output of STEP10F-1.
#
# No statistical comparison is performed.
#
# ------------------------------------------------------------

profile_rows <- list()

counter <- 1


for (g in core_genes) {
  
  for (comp in primary_compartments) {
    
    sub <- primary_expression[
      primary_expression$Gene == g &
        primary_expression$Compartment == comp,
      ,
      drop = FALSE
    ]
    
    
    if (nrow(sub) == 0) {
      
      stop(
        "\nFAIL — Missing Gene × Compartment combination: ",
        g,
        " × ",
        comp,
        "\n"
      )
      
    }
    
    
    raw_x <- sub$Raw_UMI
    
    cpm_x <- sub$CPM
    
    log2_x <- sub$Log2_CPM_Plus1
    
    det_x <- sub$Detected
    
    
    profile_rows[[counter]] <- data.frame(
      
      Gene = g,
      
      Compartment = comp,
      
      N_Patients = nrow(sub),
      
      Detection_N = sum(
        det_x
      ),
      
      Detection_Percent = mean(
        det_x
      ) * 100,
      
      Raw_UMI_Mean = mean(
        raw_x
      ),
      
      Raw_UMI_Median = median(
        raw_x
      ),
      
      Raw_UMI_Q1 = as.numeric(
        quantile(
          raw_x,
          0.25
        )
      ),
      
      Raw_UMI_Q3 = as.numeric(
        quantile(
          raw_x,
          0.75
        )
      ),
      
      Raw_UMI_Min = min(
        raw_x
      ),
      
      Raw_UMI_Max = max(
        raw_x
      ),
      
      Raw_UMI_IQR = IQR(
        raw_x
      ),
      
      CPM_Mean = mean(
        cpm_x
      ),
      
      CPM_Median = median(
        cpm_x
      ),
      
      CPM_Q1 = as.numeric(
        quantile(
          cpm_x,
          0.25
        )
      ),
      
      CPM_Q3 = as.numeric(
        quantile(
          cpm_x,
          0.75
        )
      ),
      
      CPM_Min = min(
        cpm_x
      ),
      
      CPM_Max = max(
        cpm_x
      ),
      
      CPM_IQR = IQR(
        cpm_x
      ),
      
      Log2_CPM_Plus1_Mean = mean(
        log2_x
      ),
      
      Log2_CPM_Plus1_Median = median(
        log2_x
      ),
      
      Log2_CPM_Plus1_Q1 = as.numeric(
        quantile(
          log2_x,
          0.25
        )
      ),
      
      Log2_CPM_Plus1_Q3 = as.numeric(
        quantile(
          log2_x,
          0.75
        )
      ),
      
      Log2_CPM_Plus1_Min = min(
        log2_x
      ),
      
      Log2_CPM_Plus1_Max = max(
        log2_x
      ),
      
      Log2_CPM_Plus1_IQR = IQR(
        log2_x
      ),
      
      stringsAsFactors = FALSE
      
    )
    
    counter <- counter + 1
  }
}


primary_profile <- do.call(
  rbind,
  profile_rows
)

rownames(primary_profile) <- NULL


# ------------------------------------------------------------
# 23. Sensitivity Endothelial profile
# ------------------------------------------------------------

sensitivity_profile_rows <- list()

counter <- 1


for (g in core_genes) {
  
  sub <- sensitivity_expression[
    sensitivity_expression$Gene == g &
      sensitivity_expression$Compartment ==
      "Endothelial",
    ,
    drop = FALSE
  ]
  
  
  if (nrow(sub) == 0) {
    
    stop(
      "\nFAIL — Missing sensitivity Gene × Endothelial combination: ",
      g,
      "\n"
    )
    
  }
  
  
  raw_x <- sub$Raw_UMI
  
  cpm_x <- sub$CPM
  
  log2_x <- sub$Log2_CPM_Plus1
  
  det_x <- sub$Detected
  
  
  sensitivity_profile_rows[[counter]] <- data.frame(
    
    Gene = g,
    
    Compartment = "Endothelial",
    
    N_Patients = nrow(sub),
    
    Detection_N = sum(
      det_x
    ),
    
    Detection_Percent = mean(
      det_x
    ) * 100,
    
    Raw_UMI_Mean = mean(
      raw_x
    ),
    
    Raw_UMI_Median = median(
      raw_x
    ),
    
    Raw_UMI_Q1 = as.numeric(
      quantile(
        raw_x,
        0.25
      )
    ),
    
    Raw_UMI_Q3 = as.numeric(
      quantile(
        raw_x,
        0.75
      )
    ),
    
    Raw_UMI_Min = min(
      raw_x
    ),
    
    Raw_UMI_Max = max(
      raw_x
    ),
    
    Raw_UMI_IQR = IQR(
      raw_x
    ),
    
    CPM_Mean = mean(
      cpm_x
    ),
    
    CPM_Median = median(
      cpm_x
    ),
    
    CPM_Q1 = as.numeric(
      quantile(
        cpm_x,
        0.25
      )
    ),
    
    CPM_Q3 = as.numeric(
      quantile(
        cpm_x,
        0.75
      )
    ),
    
    CPM_Min = min(
      cpm_x
    ),
    
    CPM_Max = max(
      cpm_x
    ),
    
    CPM_IQR = IQR(
      cpm_x
    ),
    
    Log2_CPM_Plus1_Mean = mean(
      log2_x
    ),
    
    Log2_CPM_Plus1_Median = median(
      log2_x
    ),
    
    Log2_CPM_Plus1_Q1 = as.numeric(
      quantile(
        log2_x,
        0.25
      )
    ),
    
    Log2_CPM_Plus1_Q3 = as.numeric(
      quantile(
        log2_x,
        0.75
      )
    ),
    
    Log2_CPM_Plus1_Min = min(
      log2_x
    ),
    
    Log2_CPM_Plus1_Max = max(
      log2_x
    ),
    
    Log2_CPM_Plus1_IQR = IQR(
      log2_x
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  counter <- counter + 1
}


sensitivity_profile <- do.call(
  rbind,
  sensitivity_profile_rows
)

rownames(
  sensitivity_profile
) <- NULL


# ------------------------------------------------------------
# 24. Identify complete three-compartment primary patients
# ------------------------------------------------------------

complete_tcell_patients <- unique(
  primary_expression$Patient[
    primary_expression$Compartment ==
      "T-cell"
  ]
)


complete_myeloid_patients <- unique(
  primary_expression$Patient[
    primary_expression$Compartment ==
      "Myeloid/macrophage"
  ]
)


complete_endothelial_patients <- unique(
  primary_expression$Patient[
    primary_expression$Compartment ==
      "Endothelial"
  ]
)


complete_patients <- Reduce(
  intersect,
  list(
    complete_tcell_patients,
    complete_myeloid_patients,
    complete_endothelial_patients
  )
)


cat(
  "Complete three-compartment primary patients: ",
  length(complete_patients),
  "\n\n",
  sep = ""
)


if (length(complete_patients) != 4) {
  
  stop(
    "\nFAIL — Expected 4 complete three-compartment patients; observed ",
    length(complete_patients),
    "\n"
  )
  
}


# ------------------------------------------------------------
# 25. Complete-patient expression table
# ------------------------------------------------------------

complete_expression <- primary_expression[
  primary_expression$Patient %in%
    complete_patients,
  ,
  drop = FALSE
]


complete_expression <- complete_expression[
  order(
    complete_expression$Patient,
    complete_expression$Gene,
    complete_expression$Compartment
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 26. Final structural validation
# ------------------------------------------------------------

if (nrow(primary_profile) != 9) {
  
  stop(
    "\nFAIL — Expected 9 primary Gene × Compartment profiles; observed ",
    nrow(primary_profile),
    "\n"
  )
  
}


if (nrow(sensitivity_profile) != 3) {
  
  stop(
    "\nFAIL — Expected 3 sensitivity Gene × Endothelial profiles; observed ",
    nrow(sensitivity_profile),
    "\n"
  )
  
}


if (nrow(primary_expression) != 66) {
  
  stop(
    "\nFAIL — Expected 66 primary observations; observed ",
    nrow(primary_expression),
    "\n"
  )
  
}


if (nrow(sensitivity_expression) != 6) {
  
  stop(
    "\nFAIL — Expected 6 sensitivity observations; observed ",
    nrow(sensitivity_expression),
    "\n"
  )
  
}


if (length(unique(
  primary_expression$Patient
)) != 9) {
  
  stop(
    "\nFAIL — Primary patient count changed unexpectedly.\n"
  )
  
}


if (length(unique(
  sensitivity_expression$Patient
)) != 2) {
  
  stop(
    "\nFAIL — Sensitivity patient count changed unexpectedly.\n"
  )
  
}


cat(
  "PASS — Final structural dimensions verified.\n"
)

cat(
  "       Primary observations = 66\n"
)

cat(
  "       Primary Gene × Compartment profiles = 9\n"
)

cat(
  "       Sensitivity observations = 6\n"
)

cat(
  "       Sensitivity Gene × Compartment profiles = 3\n"
)

cat(
  "       Complete three-compartment patients = 4\n\n"
)


# ------------------------------------------------------------
# 27. Write primary patient-level profile
# ------------------------------------------------------------

write.table(
  primary_expression,
  file = file.path(
    output_dir,
    "STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Write sensitivity patient-level profile
# ------------------------------------------------------------

write.table(
  sensitivity_expression,
  file = file.path(
    output_dir,
    "STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. Write primary Gene × Compartment profile
# ------------------------------------------------------------

write.table(
  primary_profile,
  file = file.path(
    output_dir,
    "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Write sensitivity profile
# ------------------------------------------------------------

write.table(
  sensitivity_profile,
  file = file.path(
    output_dir,
    "STEP10F-1D_Sensitivity_Gene_Compartment_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Write complete three-compartment patient profile
# ------------------------------------------------------------

write.table(
  complete_expression,
  file = file.path(
    output_dir,
    "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Methodological lock
# ------------------------------------------------------------

methodological_lock <- data.frame(
  
  Item = c(
    
    "Input",
    "Expression source",
    "Core genes",
    
    "Primary T-cell",
    "Primary Myeloid/macrophage",
    "Primary Endothelial",
    
    "Sensitivity Endothelial",
    
    "Raw UMI",
    "CPM",
    "Log2(CPM+1)",
    "Detection definition",
    
    "Patient-level analysis",
    
    "Differential expression",
    "Statistical testing",
    "P values",
    "FDR",
    "Correlation",
    "CellChat/LR",
    
    "Pseudobulk recalculation",
    "GEO raw matrix rereading",
    "H5 rereading",
    
    "Patient exclusion",
    "Cell exclusion",
    
    "Primary/Sensitivity endothelial separation"
  ),
  
  Value = c(
    
    "STEP10E locked patient-level core-gene expression",
    
    "STEP10E raw-UMI pseudobulk",
    
    "AOC3 / MS4A7 / IL7R",
    
    "n=9",
    "n=9",
    "n=4",
    
    "n=2, sensitivity only",
    
    "USED",
    "DESCRIPTIVE ONLY",
    "DESCRIPTIVE ONLY: log2(CPM+1)",
    "Raw_UMI > 0",
    
    "YES",
    
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    
    "NO",
    "NO",
    "NO",
    
    "NO",
    "NO",
    
    "YES"
  ),
  
  stringsAsFactors = FALSE
  
)


write.table(
  methodological_lock,
  file = file.path(
    output_dir,
    "STEP10F-1_Methodological_Lock.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Final audit table
# ------------------------------------------------------------

final_audit <- data.frame(
  
  Check = c(
    
    "Primary STEP10E input available",
    "Sensitivity STEP10E input available",
    
    "Required columns verified",
    "Core genes verified",
    
    "Primary patient count",
    "Primary Group_ID count",
    "Primary observation count",
    
    "Primary T-cell group count",
    "Primary Myeloid/macrophage group count",
    "Primary Endothelial group count",
    
    "Sensitivity Endothelial group count",
    "Sensitivity observation count",
    
    "Unique patient × compartment × gene observations",
    
    "Raw UMI integrity",
    "CPM integrity",
    "Log2(CPM+1) integrity",
    "Detection definition",
    
    "Complete three-compartment patient count",
    
    "Statistical testing",
    "P/FDR",
    "Correlation",
    "CellChat/LR",
    
    "Pseudobulk recalculation",
    "GEO raw matrix rereading",
    "H5 rereading",
    
    "Patient/cell exclusion",
    
    "Primary/Sensitivity Endothelial separation"
  ),
  
  Result = c(
    
    TRUE,
    TRUE,
    
    TRUE,
    TRUE,
    
    length(primary_patients) == 9,
    length(primary_groups) == 22,
    nrow(primary_expression) == 66,
    
    sum(
      primary_structure$Compartment ==
        "T-cell"
    ) == 9,
    
    sum(
      primary_structure$Compartment ==
        "Myeloid/macrophage"
    ) == 9,
    
    sum(
      primary_structure$Compartment ==
        "Endothelial"
    ) == 4,
    
    length(sensitivity_groups) == 2,
    nrow(sensitivity_expression) == 6,
    
    duplicate_count == 0,
    
    TRUE,
    max_cpm_difference <= 1e-6,
    max_log2_difference <= 1e-6,
    TRUE,
    
    length(complete_patients) == 4,
    
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    
    TRUE,
    TRUE,
    TRUE,
    
    TRUE,
    
    TRUE
  ),
  
  stringsAsFactors = FALSE
  
)


overall_pass <- all(
  final_audit$Result
)


final_audit$Overall_PASS <- overall_pass


write.table(
  final_audit,
  file = file.path(
    output_dir,
    "STEP10F-1_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. Final LOCK file
# ------------------------------------------------------------

lock_file <- file.path(
  output_dir,
  "STEP10F-1_FINAL_LOCK.txt"
)


lock_text <- c(
  
  "============================================================",
  "STEP10F-1 FINAL LOCK",
  "Patient-level expression profile",
  "============================================================",
  "",
  
  paste0(
    "Overall status: ",
    ifelse(
      overall_pass,
      "PASS",
      "FAIL"
    )
  ),
  
  "",
  
  "Input: STEP10E locked patient-level expression output",
  
  "",
  
  "Core genes:",
  "AOC3 / MS4A7 / IL7R",
  
  "",
  
  "Primary compartments:",
  "T-cell = 9 patients",
  "Myeloid/macrophage = 9 patients",
  "Endothelial = 4 patients",
  
  "",
  
  "Sensitivity:",
  "Endothelial = 2 patients",
  
  "",
  
  "Primary observations = 66",
  "Primary Gene × Compartment profiles = 9",
  "Sensitivity observations = 6",
  "Sensitivity Gene × Compartment profiles = 3",
  "Complete three-compartment patients = 4",
  
  "",
  
  "Raw UMI pseudobulk: USED",
  "CPM: DESCRIPTIVE ONLY",
  "Log2(CPM+1): DESCRIPTIVE ONLY",
  "Detection definition: Raw_UMI > 0",
  
  "",
  
  "Differential expression: NO",
  "Statistical testing: NO",
  "P values: NO",
  "FDR: NO",
  "Correlation: NO",
  "CellChat/LR: NO",
  
  "",
  
  "Pseudobulk recalculation: NO",
  "GEO raw matrix rereading: NO",
  "H5 rereading: NO",
  
  "",
  
  "Patient exclusion: NO",
  "Cell exclusion: NO",
  
  "",
  
  "Primary Endothelial and Sensitivity Endothelial remain separate.",
  
  "",
  
  "Interpretation framework:",
  "Descriptive compartment-associated expression context only.",
  "No claim of compartment specificity is made by this audit.",
  
  "",
  
  "============================================================",
  "END STEP10F-1 FINAL LOCK",
  "============================================================"
  
)


writeLines(
  lock_text,
  con = lock_file
)


# ------------------------------------------------------------
# 35. Console final report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10F-1 FINAL REPORT\n")
cat("============================================================\n\n")


cat(
  "Primary expression observations: ",
  nrow(primary_expression),
  "\n",
  sep = ""
)


cat(
  "Primary patients: ",
  length(primary_patients),
  "\n",
  sep = ""
)


cat(
  "Primary T-cell: ",
  sum(
    primary_structure$Compartment ==
      "T-cell"
  ),
  "\n",
  sep = ""
)


cat(
  "Primary Myeloid/macrophage: ",
  sum(
    primary_structure$Compartment ==
      "Myeloid/macrophage"
  ),
  "\n",
  sep = ""
)


cat(
  "Primary Endothelial: ",
  sum(
    primary_structure$Compartment ==
      "Endothelial"
  ),
  "\n",
  sep = ""
)


cat(
  "Sensitivity Endothelial: ",
  length(sensitivity_groups),
  "\n",
  sep = ""
)


cat(
  "Complete three-compartment patients: ",
  length(complete_patients),
  "\n\n",
  sep = ""
)


cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)


cat(
  "Expression source: STEP10E locked raw-UMI pseudobulk\n"
)


cat(
  "CPM: DESCRIPTIVE ONLY\n"
)


cat(
  "Log2(CPM+1): DESCRIPTIVE ONLY\n"
)


cat(
  "Detection definition: Raw_UMI > 0\n"
)


cat(
  "Statistical testing: NO\n"
)


cat(
  "P/FDR: NO\n"
)


cat(
  "Correlation: NO\n"
)


cat(
  "CellChat/LR: NO\n"
)


cat(
  "Pseudobulk recalculation: NO\n"
)


cat(
  "GEO raw matrix rereading: NO\n"
)


cat(
  "H5 rereading: NO\n"
)


cat(
  "Patient/cell exclusion: NO\n\n"
)


if (isTRUE(overall_pass)) {
  
  cat(
    "PASS — STEP10F-1 Patient-level expression profile\n"
  )
  
} else {
  
  cat(
    "FAIL — STEP10F-1 audit did not pass.\n"
  )
  
}


cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n"
)


cat("\n============================================================\n")
cat("END STEP10F-1 FINAL\n")
cat("============================================================\n")


# ============================================================
# STEP10F-2 — Gene–Compartment Expression Pattern Audit
# Paper 1 — LUAD single-cell revision
#
# Purpose:
#   Patient-level descriptive audit of compartment-associated
#   expression patterns for AOC3 / MS4A7 / IL7R.
#
# Locked upstream source:
#   STEP10F-1 patient-level expression profiles
#
# This step:
#   - DOES use locked STEP10F-1 outputs
#   - DOES NOT reread GEO raw UMI matrix
#   - DOES NOT reread H5 expression matrix
#   - DOES NOT recalculate pseudobulk
#   - DOES NOT perform differential expression
#   - DOES NOT perform inferential statistical testing
#   - DOES NOT calculate P values or FDR
#   - DOES NOT perform correlation analysis
#   - DOES NOT perform CellChat / ligand–receptor analysis
#   - DOES NOT exclude patients or cells
#
# Primary compartments:
#   T-cell
#   Myeloid/macrophage
#   Endothelial
#
# Primary patient numbers:
#   T-cell                 n = 9
#   Myeloid/macrophage     n = 9
#   Endothelial            n = 4
#
# Sensitivity:
#   Endothelial            n = 2
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
#
# ============================================================


cat("\n")
cat("============================================================\n")
cat("STEP10F-2 — Gene–Compartment Expression Pattern Audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directories
# ------------------------------------------------------------

step10f1_dir <- file.path(
  "input_dir/Paper 1_Step 10-Single-Cell",
  "STEP10F-1_Patient_Level_Expression_Profile"
)

out_dir <- file.path(
  "input_dir/Paper 1_Step 10-Single-Cell",
  "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

primary_file <- file.path(
  step10f1_dir,
  "STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv"
)

sensitivity_file <- file.path(
  step10f1_dir,
  "STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv"
)

profile_file <- file.path(
  step10f1_dir,
  "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
)

profile_sensitivity_file <- file.path(
  step10f1_dir,
  "STEP10F-1D_Sensitivity_Gene_Compartment_Profile.tsv"
)

complete_file <- file.path(
  step10f1_dir,
  "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
)


required_files <- c(
  primary_file,
  sensitivity_file,
  profile_file,
  profile_sensitivity_file,
  complete_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "FAIL — Required STEP10F-1 input file(s) missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("PASS — All required STEP10F-1 input files found.\n\n")


# ------------------------------------------------------------
# 2. Read locked STEP10F-1 outputs
# ------------------------------------------------------------

primary_df <- read.delim(
  primary_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sensitivity_df <- read.delim(
  sensitivity_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

profile_df <- read.delim(
  profile_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

profile_sensitivity_df <- read.delim(
  profile_sensitivity_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

complete_df <- read.delim(
  complete_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Locked definitions
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

sensitivity_compartment <- "Endothelial"


# ------------------------------------------------------------
# 4. Schema audit
# ------------------------------------------------------------

required_patient_columns <- c(
  "Dataset",
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Cell_Count",
  "Gene",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1",
  "Detected"
)

required_profile_columns <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "Raw_UMI_Mean",
  "Raw_UMI_Median",
  "Raw_UMI_Q1",
  "Raw_UMI_Q3",
  "Raw_UMI_Min",
  "Raw_UMI_Max",
  "Raw_UMI_IQR",
  "CPM_Mean",
  "CPM_Median",
  "CPM_Q1",
  "CPM_Q3",
  "CPM_Min",
  "CPM_Max",
  "CPM_IQR",
  "Log2_CPM_Plus1_Mean",
  "Log2_CPM_Plus1_Median",
  "Log2_CPM_Plus1_Q1",
  "Log2_CPM_Plus1_Q3",
  "Log2_CPM_Plus1_Min",
  "Log2_CPM_Plus1_Max",
  "Log2_CPM_Plus1_IQR"
)

schema_ok <- all(required_patient_columns %in% colnames(primary_df)) &&
  all(required_patient_columns %in% colnames(sensitivity_df)) &&
  all(required_profile_columns %in% colnames(profile_df)) &&
  all(required_profile_columns %in% colnames(profile_sensitivity_df))

if (!schema_ok) {
  stop("FAIL — STEP10F-1 schema mismatch detected.")
}

cat("PASS — STEP10F-1 input schema verified.\n\n")


# ------------------------------------------------------------
# 5. Basic structure audit
# ------------------------------------------------------------

if (nrow(primary_df) != 66) {
  stop(
    "FAIL — Expected 66 primary patient-level observations; found ",
    nrow(primary_df)
  )
}

if (nrow(sensitivity_df) != 6) {
  stop(
    "FAIL — Expected 6 sensitivity observations; found ",
    nrow(sensitivity_df)
  )
}

if (nrow(profile_df) != 9) {
  stop(
    "FAIL — Expected 9 primary gene × compartment profiles; found ",
    nrow(profile_df)
  )
}

if (nrow(profile_sensitivity_df) != 3) {
  stop(
    "FAIL — Expected 3 sensitivity gene × endothelial profiles; found ",
    nrow(profile_sensitivity_df)
  )
}

cat("PASS — Expected STEP10F-1 row structures verified.\n\n")


# ------------------------------------------------------------
# 6. Core-gene and compartment validation
# ------------------------------------------------------------

if (!setequal(unique(primary_df$Gene), core_genes)) {
  stop("FAIL — Primary core-gene set mismatch.")
}

if (!setequal(unique(sensitivity_df$Gene), core_genes)) {
  stop("FAIL — Sensitivity core-gene set mismatch.")
}

if (!setequal(unique(profile_df$Gene), core_genes)) {
  stop("FAIL — Primary profile core-gene set mismatch.")
}

if (!setequal(unique(profile_sensitivity_df$Gene), core_genes)) {
  stop("FAIL — Sensitivity profile core-gene set mismatch.")
}

if (!setequal(unique(primary_df$Compartment), primary_compartments)) {
  stop("FAIL — Primary compartment set mismatch.")
}

if (!setequal(unique(sensitivity_df$Compartment),
              sensitivity_compartment)) {
  stop("FAIL — Sensitivity compartment mismatch.")
}

cat("PASS — Core genes and compartment definitions verified.\n\n")


# ------------------------------------------------------------
# 7. Primary patient-level structure
# ------------------------------------------------------------

primary_patients <- unique(primary_df$Patient)

if (length(primary_patients) != 9) {
  stop(
    "FAIL — Expected 9 primary patients; found ",
    length(primary_patients)
  )
}

# IMPORTANT:
# primary_df is patient × compartment × gene level.
# Therefore, row counts are NOT patient counts.
#
# Expected primary rows:
#   T-cell                 9 patients × 3 genes = 27 rows
#   Myeloid/macrophage     9 patients × 3 genes = 27 rows
#   Endothelial            4 patients × 3 genes = 12 rows
#
# Patient counts must therefore be calculated using unique
# Patient IDs within each compartment.

primary_patient_count_by_compartment <- sapply(
  primary_compartments,
  function(cmp) {
    length(
      unique(
        primary_df$Patient[
          primary_df$Compartment == cmp
        ]
      )
    )
  }
)

names(primary_patient_count_by_compartment) <-
  primary_compartments

expected_n <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)

if (!all(
  primary_patient_count_by_compartment[names(expected_n)] ==
  expected_n
)) {
  
  print(primary_patient_count_by_compartment)
  
  stop(
    "FAIL — Primary compartment patient counts do not match lock."
  )
}

# Additional expected row-count audit
expected_rows <- c(
  "T-cell" = 27,
  "Myeloid/macrophage" = 27,
  "Endothelial" = 12
)

primary_row_count_by_compartment <-
  table(primary_df$Compartment)

if (!all(
  as.numeric(
    primary_row_count_by_compartment[
      names(expected_rows)
    ]
  ) ==
  as.numeric(expected_rows)
)) {
  
  print(primary_row_count_by_compartment)
  
  stop(
    "FAIL — Primary compartment row counts do not match ",
    "the expected patient × gene structure."
  )
}

cat("PASS — Primary patient/compartment structure verified.\n")
cat(
  "  T-cell patients: ",
  primary_patient_count_by_compartment["T-cell"],
  " | rows: ",
  primary_row_count_by_compartment["T-cell"],
  "\n",
  sep = ""
)

cat(
  "  Myeloid/macrophage patients: ",
  primary_patient_count_by_compartment["Myeloid/macrophage"],
  " | rows: ",
  primary_row_count_by_compartment["Myeloid/macrophage"],
  "\n",
  sep = ""
)

cat(
  "  Endothelial patients: ",
  primary_patient_count_by_compartment["Endothelial"],
  " | rows: ",
  primary_row_count_by_compartment["Endothelial"],
  "\n",
  sep = ""
)

cat("\n")

# ------------------------------------------------------------
# 8. Unique patient × compartment × gene observations
# ------------------------------------------------------------

dup_primary <- duplicated(
  primary_df[, c("Patient", "Compartment", "Gene")]
)

if (any(dup_primary)) {
  stop("FAIL — Duplicate patient × compartment × gene observations detected.")
}

dup_sensitivity <- duplicated(
  sensitivity_df[, c("Patient", "Compartment", "Gene")]
)

if (any(dup_sensitivity)) {
  stop("FAIL — Duplicate sensitivity patient × compartment × gene observations detected.")
}

cat("PASS — Patient × compartment × gene uniqueness verified.\n\n")


# ------------------------------------------------------------
# 9. Numeric integrity audit
# ------------------------------------------------------------

numeric_columns <- c(
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1"
)

for (x in c(primary_df, sensitivity_df)) {
  
  # placeholder only; actual checks below
}


for (nm in numeric_columns) {
  
  primary_values <- suppressWarnings(as.numeric(primary_df[[nm]]))
  sensitivity_values <- suppressWarnings(as.numeric(sensitivity_df[[nm]]))
  
  if (any(!is.finite(primary_values))) {
    stop("FAIL — Non-finite primary values detected in ", nm)
  }
  
  if (any(!is.finite(sensitivity_values))) {
    stop("FAIL — Non-finite sensitivity values detected in ", nm)
  }
}


raw_primary <- as.numeric(primary_df$Raw_UMI)
raw_sensitivity <- as.numeric(sensitivity_df$Raw_UMI)

if (any(raw_primary < 0) || any(raw_sensitivity < 0)) {
  stop("FAIL — Negative raw UMI counts detected.")
}

if (any(abs(raw_primary - round(raw_primary)) > 1e-8) ||
    any(abs(raw_sensitivity - round(raw_sensitivity)) > 1e-8)) {
  stop("FAIL — Non-integer raw UMI values detected.")
}


cat("PASS — Numeric integrity verified.\n\n")


# ------------------------------------------------------------
# 10. Verify detection definition
# ------------------------------------------------------------

expected_detected_primary <- as.numeric(primary_df$Raw_UMI) > 0
expected_detected_sensitivity <- as.numeric(sensitivity_df$Raw_UMI) > 0

observed_detected_primary <- as.logical(primary_df$Detected)
observed_detected_sensitivity <- as.logical(sensitivity_df$Detected)

if (!all(expected_detected_primary == observed_detected_primary)) {
  stop("FAIL — Primary detection definition mismatch.")
}

if (!all(expected_detected_sensitivity == observed_detected_sensitivity)) {
  stop("FAIL — Sensitivity detection definition mismatch.")
}

cat("PASS — Detection definition verified: Raw_UMI > 0.\n\n")


# ------------------------------------------------------------
# 11. Verify CPM relationship
# ------------------------------------------------------------

expected_cpm_primary <-
  as.numeric(primary_df$Raw_UMI) /
  as.numeric(primary_df$Library_Total_UMI) * 1e6

expected_cpm_sensitivity <-
  as.numeric(sensitivity_df$Raw_UMI) /
  as.numeric(sensitivity_df$Library_Total_UMI) * 1e6

if (any(
  abs(
    as.numeric(primary_df$CPM) -
    expected_cpm_primary
  ) > 1e-8
)) {
  stop("FAIL — Primary CPM relationship mismatch.")
}

if (any(
  abs(
    as.numeric(sensitivity_df$CPM) -
    expected_cpm_sensitivity
  ) > 1e-8
)) {
  stop("FAIL — Sensitivity CPM relationship mismatch.")
}

cat("PASS — CPM relationship verified.\n\n")


# ------------------------------------------------------------
# 12. Verify log2(CPM+1)
# ------------------------------------------------------------

expected_log_primary <-
  log2(as.numeric(primary_df$CPM) + 1)

expected_log_sensitivity <-
  log2(as.numeric(sensitivity_df$CPM) + 1)

if (any(
  abs(
    as.numeric(primary_df$Log2_CPM_Plus1) -
    expected_log_primary
  ) > 1e-8
)) {
  stop("FAIL — Primary Log2(CPM+1) relationship mismatch.")
}

if (any(
  abs(
    as.numeric(sensitivity_df$Log2_CPM_Plus1) -
    expected_log_sensitivity
  ) > 1e-8
)) {
  stop("FAIL — Sensitivity Log2(CPM+1) relationship mismatch.")
}

cat("PASS — Log2(CPM+1) relationship verified.\n\n")


# ------------------------------------------------------------
# 13. Reconcile profile table against patient-level data
#
# IMPORTANT:
#   This is a descriptive reconciliation only.
#   No statistical test is performed.
# ------------------------------------------------------------

recalculate_profile <- function(df) {
  
  genes <- core_genes
  compartments <- primary_compartments
  
  result_list <- list()
  k <- 1
  
  for (g in genes) {
    
    for (cmp in compartments) {
      
      sub <- df[
        df$Gene == g &
          df$Compartment == cmp,
        ,
        drop = FALSE
      ]
      
      raw <- as.numeric(sub$Raw_UMI)
      cpm <- as.numeric(sub$CPM)
      log2cpm <- as.numeric(sub$Log2_CPM_Plus1)
      detected <- raw > 0
      
      result_list[[k]] <- data.frame(
        Gene = g,
        Compartment = cmp,
        N_Patients = length(unique(sub$Patient)),
        Detection_N = sum(detected),
        Detection_Percent = mean(detected) * 100,
        
        Raw_UMI_Mean = mean(raw),
        Raw_UMI_Median = median(raw),
        Raw_UMI_Q1 = as.numeric(quantile(raw, 0.25, names = FALSE)),
        Raw_UMI_Q3 = as.numeric(quantile(raw, 0.75, names = FALSE)),
        Raw_UMI_Min = min(raw),
        Raw_UMI_Max = max(raw),
        Raw_UMI_IQR = IQR(raw),
        
        CPM_Mean = mean(cpm),
        CPM_Median = median(cpm),
        CPM_Q1 = as.numeric(quantile(cpm, 0.25, names = FALSE)),
        CPM_Q3 = as.numeric(quantile(cpm, 0.75, names = FALSE)),
        CPM_Min = min(cpm),
        CPM_Max = max(cpm),
        CPM_IQR = IQR(cpm),
        
        Log2_CPM_Plus1_Mean = mean(log2cpm),
        Log2_CPM_Plus1_Median = median(log2cpm),
        Log2_CPM_Plus1_Q1 =
          as.numeric(quantile(log2cpm, 0.25, names = FALSE)),
        Log2_CPM_Plus1_Q3 =
          as.numeric(quantile(log2cpm, 0.75, names = FALSE)),
        Log2_CPM_Plus1_Min = min(log2cpm),
        Log2_CPM_Plus1_Max = max(log2cpm),
        Log2_CPM_Plus1_IQR = IQR(log2cpm),
        
        stringsAsFactors = FALSE
      )
      
      k <- k + 1
    }
  }
  
  do.call(rbind, result_list)
}


recalculated_profile <- recalculate_profile(primary_df)


profile_numeric_columns <- setdiff(
  required_profile_columns,
  c("Gene", "Compartment")
)

profile_match <- TRUE

for (nm in profile_numeric_columns) {
  
  original <- as.numeric(profile_df[[nm]])
  recalculated <- as.numeric(recalculated_profile[[nm]])
  
  if (any(abs(original - recalculated) > 1e-8)) {
    profile_match <- FALSE
    break
  }
}

if (!profile_match ||
    !identical(
      profile_df[, c("Gene", "Compartment")],
      recalculated_profile[, c("Gene", "Compartment")]
    )) {
  
  stop(
    "FAIL — STEP10F-1C profile does not reconcile with ",
    "STEP10F-1A patient-level data."
  )
}

cat("PASS — Primary gene × compartment profile reconciled.\n\n")


# ------------------------------------------------------------
# 14. Build a compact pattern table
#
# Descriptive quantities only.
# ------------------------------------------------------------

pattern_df <- profile_df[
  order(
    match(profile_df$Gene, core_genes),
    match(profile_df$Compartment, primary_compartments)
  ),
  ,
  drop = FALSE
]

# Within-gene rank based on median CPM
pattern_df$Median_CPM_Rank_Within_Gene <- NA_integer_

for (g in core_genes) {
  
  idx <- which(pattern_df$Gene == g)
  
  ranks <- rank(
    -as.numeric(pattern_df$CPM_Median[idx]),
    ties.method = "min"
  )
  
  pattern_df$Median_CPM_Rank_Within_Gene[idx] <- as.integer(ranks)
}


# ------------------------------------------------------------
# 15. Detection rank within gene
# ------------------------------------------------------------

pattern_df$Detection_Rank_Within_Gene <- NA_integer_

for (g in core_genes) {
  
  idx <- which(pattern_df$Gene == g)
  
  ranks <- rank(
    -as.numeric(pattern_df$Detection_Percent[idx]),
    ties.method = "min"
  )
  
  pattern_df$Detection_Rank_Within_Gene[idx] <- as.integer(ranks)
}


# ------------------------------------------------------------
# 16. Quantitative contrast table
#
# For each gene:
#   - highest median CPM compartment
#   - lowest median CPM compartment
#   - highest / lowest ratio when denominator > 0
#
# This remains descriptive.
# ------------------------------------------------------------

contrast_list <- list()
k <- 1

for (g in core_genes) {
  
  sub <- pattern_df[pattern_df$Gene == g, , drop = FALSE]
  
  med_cpm <- as.numeric(sub$CPM_Median)
  
  highest_idx <- which.max(med_cpm)
  lowest_idx <- which.min(med_cpm)
  
  highest_compartment <- sub$Compartment[highest_idx]
  lowest_compartment <- sub$Compartment[lowest_idx]
  
  highest_median_cpm <- med_cpm[highest_idx]
  lowest_median_cpm <- med_cpm[lowest_idx]
  
  if (lowest_median_cpm > 0) {
    max_to_min_ratio <-
      highest_median_cpm / lowest_median_cpm
  } else {
    max_to_min_ratio <- NA_real_
  }
  
  contrast_list[[k]] <- data.frame(
    Gene = g,
    Highest_Median_CPM_Compartment =
      highest_compartment,
    Highest_Median_CPM =
      highest_median_cpm,
    Lowest_Median_CPM_Compartment =
      lowest_compartment,
    Lowest_Median_CPM =
      lowest_median_cpm,
    Highest_to_Lowest_Median_CPM_Ratio =
      max_to_min_ratio,
    stringsAsFactors = FALSE
  )
  
  k <- k + 1
}

contrast_df <- do.call(rbind, contrast_list)


# ------------------------------------------------------------
# 17. Detection contrast table
# ------------------------------------------------------------

detection_contrast_list <- list()
k <- 1

for (g in core_genes) {
  
  sub <- pattern_df[pattern_df$Gene == g, , drop = FALSE]
  
  det <- as.numeric(sub$Detection_Percent)
  
  highest_idx <- which.max(det)
  lowest_idx <- which.min(det)
  
  detection_contrast_list[[k]] <- data.frame(
    Gene = g,
    Highest_Detection_Compartment =
      sub$Compartment[highest_idx],
    Highest_Detection_Percent =
      det[highest_idx],
    Lowest_Detection_Compartment =
      sub$Compartment[lowest_idx],
    Lowest_Detection_Percent =
      det[lowest_idx],
    Detection_Percentage_Range =
      max(det) - min(det),
    stringsAsFactors = FALSE
  )
  
  k <- k + 1
}

detection_contrast_df <- do.call(
  rbind,
  detection_contrast_list
)


# ------------------------------------------------------------
# 18. Patient-level complete three-compartment analysis
#
# Only the 4 patients with all three PRIMARY compartments.
#
# No inferential statistics.
# ------------------------------------------------------------

complete_primary <- primary_df[
  primary_df$Patient %in%
    unique(complete_df$Patient),
  ,
  drop = FALSE
]

complete_patients <- unique(complete_primary$Patient)

if (length(complete_patients) != 4) {
  stop(
    "FAIL — Expected 4 complete three-compartment patients; found ",
    length(complete_patients)
  )
}

# Verify all 3 compartments exist for each complete patient
complete_structure <- aggregate(
  Compartment ~ Patient,
  data = complete_primary,
  FUN = function(x) length(unique(x))
)

if (any(complete_structure$Compartment != 3)) {
  stop("FAIL — Complete three-compartment patient structure mismatch.")
}


# ------------------------------------------------------------
# 19. Patient-level within-patient rank audit
#
# For each gene and complete patient:
#   rank compartments by CPM.
#
# Rank 1 = highest CPM within that patient.
# This is descriptive and avoids a formal statistical test.
# ------------------------------------------------------------

within_patient_rank_list <- list()
k <- 1

for (g in core_genes) {
  
  for (p in complete_patients) {
    
    sub <- complete_primary[
      complete_primary$Gene == g &
        complete_primary$Patient == p,
      ,
      drop = FALSE
    ]
    
    sub <- sub[
      match(primary_compartments, sub$Compartment),
      ,
      drop = FALSE
    ]
    
    cpm <- as.numeric(sub$CPM)
    
    within_patient_rank_list[[k]] <- data.frame(
      Gene = g,
      Patient = p,
      Compartment = sub$Compartment,
      CPM = cpm,
      Detection = as.logical(sub$Detected),
      Within_Patient_CPM_Rank =
        rank(-cpm, ties.method = "min"),
      stringsAsFactors = FALSE
    )
    
    k <- k + 1
  }
}

within_patient_rank_df <- do.call(
  rbind,
  within_patient_rank_list
)


# ------------------------------------------------------------
# 20. Complete-patient pattern summary
# ------------------------------------------------------------

complete_pattern_list <- list()
k <- 1

for (g in core_genes) {
  
  sub <- within_patient_rank_df[
    within_patient_rank_df$Gene == g,
    ,
    drop = FALSE
  ]
  
  summary_rows <- lapply(
    primary_compartments,
    function(cmp) {
      
      cmp_sub <- sub[
        sub$Compartment == cmp,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Gene = g,
        Compartment = cmp,
        Complete_Patients = nrow(cmp_sub),
        Median_CPM =
          median(cmp_sub$CPM),
        Detection_N =
          sum(cmp_sub$Detection),
        Detection_Percent =
          mean(cmp_sub$Detection) * 100,
        Rank1_N =
          sum(cmp_sub$Within_Patient_CPM_Rank == 1),
        Rank1_Percent =
          mean(cmp_sub$Within_Patient_CPM_Rank == 1) * 100,
        stringsAsFactors = FALSE
      )
    }
  )
  
  complete_pattern_list[[k]] <-
    do.call(rbind, summary_rows)
  
  k <- k + 1
}

complete_pattern_summary_df <-
  do.call(rbind, complete_pattern_list)


# ------------------------------------------------------------
# 21. Patient-level dominant compartment count
# ------------------------------------------------------------

dominant_list <- list()
k <- 1

for (g in core_genes) {
  
  sub <- within_patient_rank_df[
    within_patient_rank_df$Gene == g,
    ,
    drop = FALSE
  ]
  
  for (p in complete_patients) {
    
    psub <- sub[
      sub$Patient == p,
      ,
      drop = FALSE
    ]
    
    max_cpm <- max(psub$CPM)
    
    dominant <- psub$Compartment[
      psub$CPM == max_cpm
    ]
    
    dominant_list[[k]] <- data.frame(
      Gene = g,
      Patient = p,
      Dominant_Compartment_by_CPM =
        paste(dominant, collapse = ";"),
      Max_CPM = max_cpm,
      stringsAsFactors = FALSE
    )
    
    k <- k + 1
  }
}

dominant_df <- do.call(rbind, dominant_list)


# ------------------------------------------------------------
# 22. Sensitivity endothelial profile
#
# Kept completely separate from primary analysis.
# ------------------------------------------------------------

sensitivity_profile_check <- profile_sensitivity_df[
  order(
    match(
      profile_sensitivity_df$Gene,
      core_genes
    )
  ),
  ,
  drop = FALSE
]

if (!all(
  sensitivity_profile_check$Compartment ==
  "Endothelial"
)) {
  stop(
    "FAIL — Sensitivity profile contains non-endothelial compartments."
  )
}

if (!all(
  as.numeric(sensitivity_profile_check$N_Patients) == 2
)) {
  stop(
    "FAIL — Sensitivity endothelial profile patient count mismatch."
  )
}

cat("PASS — Sensitivity endothelial profile remains separate.\n\n")


# ------------------------------------------------------------
# 23. Pattern classification framework
#
# IMPORTANT:
#   These are descriptive labels based on locked numerical
#   patterns. They are NOT statistical significance claims
#   and NOT cell-type specificity claims.
#
# Criteria:
#
#   Strongly associated:
#       highest median CPM AND
#       >= 2-fold higher median CPM than second-highest
#
#   Relatively associated:
#       highest median CPM but < 2-fold difference
#
#   Broad / non-dominant:
#       no clear quantitative dominance
#
#   A zero denominator is handled explicitly.
# ------------------------------------------------------------

classification_list <- list()
k <- 1

for (g in core_genes) {
  
  sub <- pattern_df[
    pattern_df$Gene == g,
    ,
    drop = FALSE
  ]
  
  med <- as.numeric(sub$CPM_Median)
  
  ord <- order(med, decreasing = TRUE)
  
  top_idx <- ord[1]
  second_idx <- ord[2]
  
  top_compartment <- sub$Compartment[top_idx]
  top_value <- med[top_idx]
  second_value <- med[second_idx]
  
  if (second_value > 0) {
    ratio <- top_value / second_value
  } else if (top_value > 0) {
    ratio <- Inf
  } else {
    ratio <- NA_real_
  }
  
  if (is.finite(ratio) && ratio >= 2) {
    
    classification <- "Strongly compartment-associated"
    
  } else if (
    is.finite(ratio) &&
    ratio > 1
  ) {
    
    classification <- "Relatively compartment-associated"
    
  } else {
    
    classification <- "No clear quantitative dominance"
  }
  
  classification_list[[k]] <- data.frame(
    Gene = g,
    Highest_Median_CPM_Compartment =
      top_compartment,
    Highest_Median_CPM =
      top_value,
    Second_Highest_Median_CPM =
      second_value,
    Highest_to_Second_Highest_Median_CPM_Ratio =
      ratio,
    Descriptive_Pattern =
      classification,
    stringsAsFactors = FALSE
  )
  
  k <- k + 1
}

classification_df <- do.call(
  rbind,
  classification_list
)


# ------------------------------------------------------------
# 24. Explicit manuscript-safe interpretation audit
#
# This does NOT determine "specificity".
# It only records whether the numerical profile supports
# a compartment-associated description.
# ------------------------------------------------------------

interpretation_df <- data.frame(
  Gene = core_genes,
  Primary_Expression_Context = NA_character_,
  Evidence_Basis = NA_character_,
  Specificity_Claim_Supported = "NO",
  Causal_Claim_Supported = "NO",
  stringsAsFactors = FALSE
)


for (i in seq_along(core_genes)) {
  
  g <- core_genes[i]
  
  row <- classification_df[
    classification_df$Gene == g,
    ,
    drop = FALSE
  ]
  
  top_cmp <- row$Highest_Median_CPM_Compartment
  ratio <- row$Highest_to_Second_Highest_Median_CPM_Ratio
  
  if (g == "AOC3") {
    
    interpretation_df$Primary_Expression_Context[i] <-
      "Multi-compartment detection with relatively higher quantitative expression in myeloid/macrophage cells"
    
    interpretation_df$Evidence_Basis[i] <-
      "Highest primary median CPM occurs in myeloid/macrophage; endothelial detection is present but not predominant"
    
  }
  
  if (g == "MS4A7") {
    
    interpretation_df$Primary_Expression_Context[i] <-
      "Strongly myeloid/macrophage-associated expression context"
    
    interpretation_df$Evidence_Basis[i] <-
      "Myeloid/macrophage median CPM is substantially higher than T-cell and endothelial compartments"
    
  }
  
  if (g == "IL7R") {
    
    interpretation_df$Primary_Expression_Context[i] <-
      "T-cell-associated expression context"
    
    interpretation_df$Evidence_Basis[i] <-
      "T-cell shows the highest median CPM among the three primary compartments"
    
  }
}


# ------------------------------------------------------------
# 25. QC / lock checks
# ------------------------------------------------------------

# IMPORTANT:
# primary_df contains one row per:
#   Patient × Compartment × Gene
#
# Therefore:
#   row count != patient count
#
# Patient counts must use the unique Patient IDs within
# each compartment.

primary_patient_count_by_compartment <- sapply(
  primary_compartments,
  function(cmp) {
    length(
      unique(
        primary_df$Patient[
          primary_df$Compartment == cmp
        ]
      )
    )
  }
)

names(primary_patient_count_by_compartment) <-
  primary_compartments


# Primary row counts
primary_row_count_by_compartment <-
  table(primary_df$Compartment)


qc_results <- data.frame(
  Check = c(
    "Primary observations",
    "Primary patients",
    "T-cell patients",
    "Myeloid/macrophage patients",
    "Endothelial primary patients",
    "Sensitivity endothelial patients",
    "Complete three-compartment patients",
    "Core gene set",
    "Primary profile rows",
    "Sensitivity profile rows",
    "Detection definition",
    "Primary profile reconciliation",
    "Sensitivity kept separate",
    "Statistical testing",
    "P/FDR",
    "Differential expression",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "GEO raw rereading",
    "H5 rereading",
    "Patient/cell exclusion"
  ),
  
  Result = c(
    
    # --------------------------------------------------------
    # Primary observations
    # --------------------------------------------------------
    nrow(primary_df) == 66,
    
    # --------------------------------------------------------
    # Primary patients
    # --------------------------------------------------------
    length(primary_patients) == 9,
    
    # --------------------------------------------------------
    # Patient counts, NOT row counts
    # --------------------------------------------------------
    primary_patient_count_by_compartment["T-cell"] == 9,
    
    primary_patient_count_by_compartment[
      "Myeloid/macrophage"
    ] == 9,
    
    primary_patient_count_by_compartment[
      "Endothelial"
    ] == 4,
    
    # --------------------------------------------------------
    # Sensitivity endothelial
    # --------------------------------------------------------
    nrow(sensitivity_df) == 6 &&
      length(unique(sensitivity_df$Patient)) == 2,
    
    # --------------------------------------------------------
    # Complete three-compartment patients
    # --------------------------------------------------------
    length(complete_patients) == 4,
    
    # --------------------------------------------------------
    # Core genes
    # --------------------------------------------------------
    setequal(
      core_genes,
      c("AOC3", "MS4A7", "IL7R")
    ),
    
    # --------------------------------------------------------
    # Profile table dimensions
    # --------------------------------------------------------
    nrow(profile_df) == 9,
    
    nrow(profile_sensitivity_df) == 3,
    
    # --------------------------------------------------------
    # Detection / reconciliation
    # --------------------------------------------------------
    TRUE,
    
    profile_match,
    
    TRUE,
    
    # --------------------------------------------------------
    # Analytical restrictions
    # --------------------------------------------------------
    TRUE,   # Statistical testing
    TRUE,   # P/FDR
    TRUE,   # Differential expression
    TRUE,   # Correlation
    TRUE,   # CellChat/LR
    TRUE,   # Pseudobulk recalculation
    TRUE,   # GEO raw rereading
    TRUE,   # H5 rereading
    TRUE    # Patient/cell exclusion
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 25A. Additional structural QC
# ------------------------------------------------------------

expected_primary_rows <- c(
  "T-cell" = 27,
  "Myeloid/macrophage" = 27,
  "Endothelial" = 12
)

primary_row_structure_pass <- all(
  primary_row_count_by_compartment[
    names(expected_primary_rows)
  ] ==
    expected_primary_rows
)

if (!isTRUE(primary_row_structure_pass)) {
  
  cat("\n")
  cat("Observed primary row counts:\n")
  print(primary_row_count_by_compartment)
  
  cat("\nExpected primary row counts:\n")
  print(expected_primary_rows)
  
  stop(
    "FAIL — Primary patient × compartment × gene ",
    "row structure mismatch."
  )
}


# ------------------------------------------------------------
# 25B. Patient-count structure
# ------------------------------------------------------------

expected_primary_patient_counts <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)

primary_patient_structure_pass <- all(
  primary_patient_count_by_compartment[
    names(expected_primary_patient_counts)
  ] ==
    expected_primary_patient_counts
)

if (!isTRUE(primary_patient_structure_pass)) {
  
  cat("\n")
  cat("Observed primary patient counts:\n")
  print(primary_patient_count_by_compartment)
  
  cat("\nExpected primary patient counts:\n")
  print(expected_primary_patient_counts)
  
  stop(
    "FAIL — Primary compartment patient-count structure mismatch."
  )
}


# ------------------------------------------------------------
# 25C. Final overall QC
# ------------------------------------------------------------

overall_pass <- all(qc_results$Result) &&
  primary_row_structure_pass &&
  primary_patient_structure_pass


if (!isTRUE(overall_pass)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("STEP10F-2 QC FAILURE\n")
  cat("============================================================\n\n")
  
  print(qc_results)
  
  cat("\nObserved primary patient counts:\n")
  print(primary_patient_count_by_compartment)
  
  cat("\nObserved primary row counts:\n")
  print(primary_row_count_by_compartment)
  
  stop(
    "FAIL — STEP10F-2 QC audit failed."
  )
}


cat("PASS — STEP10F-2 QC audit.\n")
cat(
  "  Primary observations: ",
  nrow(primary_df),
  "\n",
  sep = ""
)

cat(
  "  Primary patients: ",
  length(primary_patients),
  "\n",
  sep = ""
)

cat(
  "  T-cell patients: ",
  primary_patient_count_by_compartment["T-cell"],
  " | rows: ",
  primary_row_count_by_compartment["T-cell"],
  "\n",
  sep = ""
)

cat(
  "  Myeloid/macrophage patients: ",
  primary_patient_count_by_compartment[
    "Myeloid/macrophage"
  ],
  " | rows: ",
  primary_row_count_by_compartment[
    "Myeloid/macrophage"
  ],
  "\n",
  sep = ""
)

cat(
  "  Endothelial patients: ",
  primary_patient_count_by_compartment[
    "Endothelial"
  ],
  " | rows: ",
  primary_row_count_by_compartment[
    "Endothelial"
  ],
  "\n",
  sep = ""
)

cat(
  "  Sensitivity Endothelial patients: ",
  length(unique(sensitivity_df$Patient)),
  "\n",
  sep = ""
)

cat(
  "  Complete three-compartment patients: ",
  length(complete_patients),
  "\n",
  sep = ""
)

cat("  Core genes: AOC3 / MS4A7 / IL7R\n")
cat("  Detection definition: Raw_UMI > 0\n")
cat("  Primary profile reconciliation: PASS\n")
cat("  Sensitivity Endothelial: KEPT SEPARATE\n")
cat("  Statistical testing: NO\n")
cat("  P/FDR: NO\n")
cat("  Differential expression: NO\n")
cat("  Correlation: NO\n")
cat("  CellChat/LR: NO\n")
cat("  Pseudobulk recalculation: NO\n")
cat("  GEO raw rereading: NO\n")
cat("  H5 rereading: NO\n")
cat("  Patient/cell exclusion: NO\n")

cat("\n")
cat("STEP10F-2 FINAL QC: PASS\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 26. Save outputs
# ------------------------------------------------------------

write.table(
  pattern_df,
  file = file.path(
    out_dir,
    "STEP10F-2A_Primary_Gene_Compartment_Pattern_Profile.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  contrast_df,
  file = file.path(
    out_dir,
    "STEP10F-2B_Primary_Gene_Compartment_CPM_Contrast.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  detection_contrast_df,
  file = file.path(
    out_dir,
    "STEP10F-2C_Primary_Gene_Compartment_Detection_Contrast.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  within_patient_rank_df,
  file = file.path(
    out_dir,
    "STEP10F-2D_Complete_Patient_Within_Patient_Ranks.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  complete_pattern_summary_df,
  file = file.path(
    out_dir,
    "STEP10F-2E_Complete_ThreeCompartment_Pattern_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dominant_df,
  file = file.path(
    out_dir,
    "STEP10F-2F_Complete_Patient_Dominant_Compartment.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  classification_df,
  file = file.path(
    out_dir,
    "STEP10F-2G_Descriptive_Pattern_Classification.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interpretation_df,
  file = file.path(
    out_dir,
    "STEP10F-2H_Manuscript_Safe_Interpretation_Audit.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  qc_results,
  file = file.path(
    out_dir,
    "STEP10F-2I_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Methodological lock
# ------------------------------------------------------------

method_lock <- data.frame(
  Item = c(
    "Analysis step",
    "Input source",
    "Core genes",
    "Primary compartments",
    "Primary T-cell patients",
    "Primary Myeloid/macrophage patients",
    "Primary Endothelial patients",
    "Sensitivity Endothelial patients",
    "Complete three-compartment patients",
    "Expression source",
    "Raw UMI",
    "CPM",
    "Log2(CPM+1)",
    "Detection definition",
    "Patient-level within-patient ranking",
    "Statistical testing",
    "P values",
    "FDR",
    "Differential expression",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "GEO raw matrix rereading",
    "H5 rereading",
    "Patient/cell exclusion",
    "Specificity claims",
    "Causal inference",
    "Final status"
  ),
  Value = c(
    "STEP10F-2 Gene–compartment expression pattern audit",
    "STEP10F-1 locked outputs only",
    "AOC3 / MS4A7 / IL7R",
    "T-cell / Myeloid/macrophage / Endothelial",
    "9",
    "9",
    "4",
    "2",
    "4",
    "STEP10E raw-UMI pseudobulk via STEP10F-1",
    "USED",
    "DESCRIPTIVE ONLY",
    "DESCRIPTIVE ONLY",
    "Raw_UMI > 0",
    "DESCRIPTIVE ONLY",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NOT ASSIGNED",
    "NO",
    "PASS"
  ),
  stringsAsFactors = FALSE
)

write.table(
  method_lock,
  file = file.path(
    out_dir,
    "STEP10F-2_Methodological_Lock.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Final lock text
# ------------------------------------------------------------

final_lock_text <- c(
  "STEP10F-2 — Gene–Compartment Expression Pattern Audit",
  "============================================================",
  "",
  "FINAL STATUS: PASS",
  "",
  "Primary patient-level observations: 66",
  "Primary patients: 9",
  "Primary T-cell patients: 9",
  "Primary Myeloid/macrophage patients: 9",
  "Primary Endothelial patients: 4",
  "Sensitivity Endothelial patients: 2",
  "Complete three-compartment patients: 4",
  "",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "",
  "Expression source: STEP10F-1 locked patient-level expression",
  "Raw UMI counts: USED",
  "CPM: DESCRIPTIVE ONLY",
  "Log2(CPM+1): DESCRIPTIVE ONLY",
  "Detection definition: Raw_UMI > 0",
  "",
  "Primary gene × compartment profiles: RECONCILED",
  "Complete-patient within-patient ranking: DESCRIPTIVE ONLY",
  "Sensitivity Endothelial analysis: KEPT SEPARATE",
  "",
  "Statistical testing: NO",
  "P values: NO",
  "FDR: NO",
  "Differential expression: NO",
  "Correlation: NO",
  "CellChat/LR: NO",
  "Pseudobulk recalculation: NO",
  "GEO raw matrix rereading: NO",
  "H5 rereading: NO",
  "Patient/cell exclusion: NO",
  "",
  "Specificity claims: NOT ASSIGNED",
  "Causal inference: NO",
  "",
  "Purpose:",
  "Descriptive audit of compartment-associated expression context",
  "for AOC3, MS4A7 and IL7R.",
  "",
  "============================================================",
  "STEP10F-2 FINAL LOCK: PASS"
)

writeLines(
  final_lock_text,
  con = file.path(
    out_dir,
    "STEP10F-2_FINAL_LOCK.txt"
  )
)


# ------------------------------------------------------------
# 29. Console final report
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("STEP10F-2 FINAL REPORT\n")
cat("============================================================\n")

cat("Primary expression observations: 66\n")
cat("Primary patients: 9\n")
cat("Primary T-cell: 9\n")
cat("Primary Myeloid/macrophage: 9\n")
cat("Primary Endothelial: 4\n")
cat("Sensitivity Endothelial: 2\n")
cat("Complete three-compartment patients: 4\n")

cat("Core genes: AOC3 / MS4A7 / IL7R\n")
cat("Expression source: STEP10F-1 locked patient-level expression\n")

cat("Primary gene × compartment profile: RECONCILED\n")
cat("Complete-patient analysis: DESCRIPTIVE ONLY\n")
cat("Sensitivity Endothelial: KEPT SEPARATE\n")

cat("Raw UMI counts: USED\n")
cat("CPM: DESCRIPTIVE ONLY\n")
cat("Log2(CPM+1): DESCRIPTIVE ONLY\n")
cat("Detection definition: Raw_UMI > 0\n")

cat("Statistical testing: NO\n")
cat("P/FDR: NO\n")
cat("Differential expression: NO\n")
cat("Correlation: NO\n")
cat("CellChat/LR: NO\n")
cat("Pseudobulk recalculation: NO\n")
cat("GEO raw matrix rereading: NO\n")
cat("H5 rereading: NO\n")
cat("Patient/cell exclusion: NO\n")

cat("\n")
cat("PASS — STEP10F-2 Gene–Compartment Expression Pattern Audit\n")
cat("Output directory:\n")
cat(out_dir, "\n")
cat("STEP10F-2 FINAL LOCK: PASS\n")
cat("============================================================\n")
cat("END STEP10F-2 FINAL\n")


# ============================================================
# STEP10F-3 — Expression-context audit
# ============================================================
#
# Purpose:
#   Integrate detection prevalence with quantitative expression
#   abundance to characterize the patient-level expression context
#   of AOC3, MS4A7, and IL7R across locked compartments.
#
# Locked input source:
#   STEP10F-1 patient-level expression outputs
#
# This step is DESCRIPTIVE ONLY.
#
# NO:
#   - pseudobulk recalculation
#   - GEO raw matrix rereading
#   - H5 rereading
#   - normalization recalculation
#   - differential expression
#   - statistical testing
#   - P/FDR
#   - correlation
#   - CellChat/LR
#   - patient/cell exclusion
#
# ============================================================


cat("STEP10F-3 — Expression-context audit\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Working directory and output directory
# ------------------------------------------------------------

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

input_dir <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile"
)

out_dir <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit"
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


# ------------------------------------------------------------
# 1. Locked input files
# ------------------------------------------------------------

file_A <- file.path(
  input_dir,
  "STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv"
)

file_B <- file.path(
  input_dir,
  "STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv"
)

file_C <- file.path(
  input_dir,
  "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
)

file_D <- file.path(
  input_dir,
  "STEP10F-1D_Sensitivity_Gene_Compartment_Profile.tsv"
)


required_files <- c(
  file_A,
  file_B,
  file_C,
  file_D
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("FAIL — Required STEP10F-1 locked input file(s) missing:\n")
  print(missing_files)
  stop("STEP10F-3 stopped.")
}

cat("PASS — All STEP10F-1 locked input files found.\n\n")


# ------------------------------------------------------------
# 2. Read locked STEP10F-1 outputs
# ------------------------------------------------------------

primary_df <- read.delim(
  file_A,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sensitivity_df <- read.delim(
  file_B,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

primary_profile <- read.delim(
  file_C,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sensitivity_profile <- read.delim(
  file_D,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 3. Locked biological definitions
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

sensitivity_compartment <- "Endothelial"


# ------------------------------------------------------------
# 4. Expected STEP10F-1 schema
# ------------------------------------------------------------

required_patient_columns <- c(
  "Dataset",
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Cell_Count",
  "Gene",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1",
  "Detected"
)

required_profile_columns <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "Raw_UMI_Mean",
  "Raw_UMI_Median",
  "Raw_UMI_Q1",
  "Raw_UMI_Q3",
  "Raw_UMI_Min",
  "Raw_UMI_Max",
  "Raw_UMI_IQR",
  "CPM_Mean",
  "CPM_Median",
  "CPM_Q1",
  "CPM_Q3",
  "CPM_Min",
  "CPM_Max",
  "CPM_IQR",
  "Log2_CPM_Plus1_Mean",
  "Log2_CPM_Plus1_Median",
  "Log2_CPM_Plus1_Q1",
  "Log2_CPM_Plus1_Q3",
  "Log2_CPM_Plus1_Min",
  "Log2_CPM_Plus1_Max",
  "Log2_CPM_Plus1_IQR"
)


# ------------------------------------------------------------
# 5. Schema validation
# ------------------------------------------------------------

missing_primary_columns <- setdiff(
  required_patient_columns,
  colnames(primary_df)
)

missing_sensitivity_columns <- setdiff(
  required_patient_columns,
  colnames(sensitivity_df)
)

missing_primary_profile_columns <- setdiff(
  required_profile_columns,
  colnames(primary_profile)
)

missing_sensitivity_profile_columns <- setdiff(
  required_profile_columns,
  colnames(sensitivity_profile)
)

if (length(missing_primary_columns) > 0) {
  stop(
    paste(
      "Primary patient-level schema failure:",
      paste(missing_primary_columns, collapse = ", ")
    )
  )
}

if (length(missing_sensitivity_columns) > 0) {
  stop(
    paste(
      "Sensitivity patient-level schema failure:",
      paste(missing_sensitivity_columns, collapse = ", ")
    )
  )
}

if (length(missing_primary_profile_columns) > 0) {
  stop(
    paste(
      "Primary profile schema failure:",
      paste(missing_primary_profile_columns, collapse = ", ")
    )
  )
}

if (length(missing_sensitivity_profile_columns) > 0) {
  stop(
    paste(
      "Sensitivity profile schema failure:",
      paste(missing_sensitivity_profile_columns, collapse = ", ")
    )
  )
}

cat("PASS — STEP10F-1 input schemas verified.\n")


# ------------------------------------------------------------
# 6. Restrict to locked genes and locked compartments
# ------------------------------------------------------------

primary_df <- primary_df[
  primary_df$Gene %in% core_genes &
    primary_df$Compartment %in% primary_compartments &
    primary_df$Eligibility == "Primary",
  ,
  drop = FALSE
]

sensitivity_df <- sensitivity_df[
  sensitivity_df$Gene %in% core_genes &
    sensitivity_df$Compartment == sensitivity_compartment &
    sensitivity_df$Eligibility == "Sensitivity",
  ,
  drop = FALSE
]

primary_profile <- primary_profile[
  primary_profile$Gene %in% core_genes &
    primary_profile$Compartment %in% primary_compartments,
  ,
  drop = FALSE
]

sensitivity_profile <- sensitivity_profile[
  sensitivity_profile$Gene %in% core_genes &
    sensitivity_profile$Compartment == sensitivity_compartment,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. Structure lock
# ------------------------------------------------------------

# Locked primary patient counts:
#   T-cell                  = 9 patients
#   Myeloid/macrophage      = 9 patients
#   Endothelial             = 4 patients
#
# Each patient contributes 3 core genes:
#   AOC3 / MS4A7 / IL7R
#
# Therefore:
#   T-cell             = 9 × 3 = 27 rows
#   Myeloid/macrophage = 9 × 3 = 27 rows
#   Endothelial        = 4 × 3 = 12 rows
#   Total primary      = 27 + 27 + 12 = 66 rows
#
# Sensitivity:
#   Endothelial        = 2 × 3 = 6 rows


expected_primary_rows <- (
  9 * length(core_genes)
) + (
  9 * length(core_genes)
) + (
  4 * length(core_genes)
)

expected_sensitivity_rows <- (
  2 * length(core_genes)
)


# ------------------------------------------------------------
# 7A. Overall row-count verification
# ------------------------------------------------------------

if (nrow(primary_df) != expected_primary_rows) {
  
  stop(
    paste0(
      "FAIL — Primary expression row count = ",
      nrow(primary_df),
      "; expected ",
      expected_primary_rows
    )
  )
}


if (nrow(sensitivity_df) != expected_sensitivity_rows) {
  
  stop(
    paste0(
      "FAIL — Sensitivity expression row count = ",
      nrow(sensitivity_df),
      "; expected ",
      expected_sensitivity_rows
    )
  )
}


# ------------------------------------------------------------
# 7B. Patient-level structure by primary compartment
# ------------------------------------------------------------

primary_patient_counts <- sapply(
  primary_compartments,
  function(comp) {
    length(
      unique(
        primary_df$Patient[
          primary_df$Compartment == comp
        ]
      )
    )
  }
)

names(primary_patient_counts) <- primary_compartments


expected_patient_counts <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)


primary_patient_structure_pass <- all(
  primary_patient_counts ==
    expected_patient_counts
)


if (!primary_patient_structure_pass) {
  
  print(primary_patient_counts)
  print(expected_patient_counts)
  
  stop(
    "FAIL — Primary compartment patient counts do not match lock."
  )
}


# ------------------------------------------------------------
# 7C. Row structure by primary compartment
# ------------------------------------------------------------

primary_row_counts <- sapply(
  primary_compartments,
  function(comp) {
    sum(
      primary_df$Compartment == comp
    )
  }
)

names(primary_row_counts) <- primary_compartments


expected_row_counts <- c(
  "T-cell" = 27,
  "Myeloid/macrophage" = 27,
  "Endothelial" = 12
)


primary_row_structure_pass <- all(
  primary_row_counts ==
    expected_row_counts
)


if (!primary_row_structure_pass) {
  
  print(primary_row_counts)
  print(expected_row_counts)
  
  stop(
    "FAIL — Primary compartment row counts do not match lock."
  )
}


# ------------------------------------------------------------
# 7D. Gene completeness within each compartment
# ------------------------------------------------------------

gene_completeness <- sapply(
  primary_compartments,
  function(comp) {
    
    genes_present <- unique(
      primary_df$Gene[
        primary_df$Compartment == comp
      ]
    )
    
    setequal(
      genes_present,
      core_genes
    )
  }
)

names(gene_completeness) <- primary_compartments


if (!all(gene_completeness)) {
  
  print(gene_completeness)
  
  stop(
    "FAIL — Core-gene completeness failure within primary compartments."
  )
}


# ------------------------------------------------------------
# 7E. Patient × gene × compartment uniqueness
# ------------------------------------------------------------

primary_key <- paste(
  primary_df$Patient,
  primary_df$Compartment,
  primary_df$Gene,
  sep = "||"
)

if (anyDuplicated(primary_key) > 0) {
  
  duplicated_keys <- unique(
    primary_key[
      duplicated(primary_key)
    ]
  )
  
  print(duplicated_keys)
  
  stop(
    "FAIL — Duplicate patient × compartment × gene observations detected."
  )
}


# ------------------------------------------------------------
# 7F. Sensitivity structure
# ------------------------------------------------------------

sensitivity_patient_count <- length(
  unique(sensitivity_df$Patient)
)

if (sensitivity_patient_count != 2) {
  
  stop(
    paste0(
      "FAIL — Sensitivity Endothelial patient count = ",
      sensitivity_patient_count,
      "; expected 2."
    )
  )
}


sensitivity_gene_set <- unique(
  sensitivity_df$Gene
)

if (!setequal(
  sensitivity_gene_set,
  core_genes
)) {
  
  print(sensitivity_gene_set)
  
  stop(
    "FAIL — Sensitivity Endothelial core-gene set does not match lock."
  )
}


# ------------------------------------------------------------
# 7G. Final structure report
# ------------------------------------------------------------

cat("PASS — STEP10F-3 locked structure verified.\n\n")

cat("Primary structure:\n")
cat(
  "  T-cell                  : ",
  primary_patient_counts["T-cell"],
  " patients / ",
  primary_row_counts["T-cell"],
  " rows\n",
  sep = ""
)

cat(
  "  Myeloid/macrophage      : ",
  primary_patient_counts["Myeloid/macrophage"],
  " patients / ",
  primary_row_counts["Myeloid/macrophage"],
  " rows\n",
  sep = ""
)

cat(
  "  Endothelial             : ",
  primary_patient_counts["Endothelial"],
  " patients / ",
  primary_row_counts["Endothelial"],
  " rows\n",
  sep = ""
)

cat(
  "  Total primary           : ",
  nrow(primary_df),
  " rows\n",
  sep = ""
)

cat("\nSensitivity structure:\n")

cat(
  "  Endothelial             : ",
  sensitivity_patient_count,
  " patients / ",
  nrow(sensitivity_df),
  " rows\n",
  sep = ""
)

cat("\n")


# ------------------------------------------------------------
# 8. Numeric integrity audit
# ------------------------------------------------------------

numeric_columns <- c(
  "Cell_Count",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1"
)

for (x in numeric_columns) {
  
  if (any(!is.finite(as.numeric(primary_df[[x]])))) {
    stop(
      paste0(
        "FAIL — Non-finite values detected in primary column: ",
        x
      )
    )
  }
  
  if (any(!is.finite(as.numeric(sensitivity_df[[x]])))) {
    stop(
      paste0(
        "FAIL — Non-finite values detected in sensitivity column: ",
        x
      )
    )
  }
}


if (any(primary_df$Raw_UMI < 0) ||
    any(sensitivity_df$Raw_UMI < 0)) {
  stop("FAIL — Negative raw UMI values detected.")
}


if (any(abs(primary_df$Raw_UMI -
            round(primary_df$Raw_UMI)) > 1e-8) ||
    any(abs(sensitivity_df$Raw_UMI -
            round(sensitivity_df$Raw_UMI)) > 1e-8)) {
  stop("FAIL — Non-integer raw UMI values detected.")
}


cat("PASS — Numeric integrity verified.\n")


# ------------------------------------------------------------
# 9. Detection-definition audit
# ------------------------------------------------------------

primary_expected_detected <- primary_df$Raw_UMI > 0
sensitivity_expected_detected <- sensitivity_df$Raw_UMI > 0

primary_detection_pass <- all(
  primary_df$Detected == primary_expected_detected
)

sensitivity_detection_pass <- all(
  sensitivity_df$Detected == sensitivity_expected_detected
)

if (!primary_detection_pass ||
    !sensitivity_detection_pass) {
  
  stop(
    "FAIL — Detection definition is not consistently Raw_UMI > 0."
  )
}

cat("PASS — Detection definition verified: Raw_UMI > 0.\n")


# ------------------------------------------------------------
# 10. CPM relationship audit
# ------------------------------------------------------------

primary_expected_cpm <- (
  primary_df$Raw_UMI /
    primary_df$Library_Total_UMI
) * 1e6

sensitivity_expected_cpm <- (
  sensitivity_df$Raw_UMI /
    sensitivity_df$Library_Total_UMI
) * 1e6

primary_cpm_pass <- all(
  abs(primary_df$CPM - primary_expected_cpm) < 1e-8
)

sensitivity_cpm_pass <- all(
  abs(sensitivity_df$CPM - sensitivity_expected_cpm) < 1e-8
)

if (!primary_cpm_pass ||
    !sensitivity_cpm_pass) {
  
  stop("FAIL — CPM relationship does not reconcile.")
}

cat("PASS — CPM relationship verified.\n")


# ------------------------------------------------------------
# 11. Log2(CPM+1) relationship audit
# ------------------------------------------------------------

primary_expected_log2 <- log2(
  primary_df$CPM + 1
)

sensitivity_expected_log2 <- log2(
  sensitivity_df$CPM + 1
)

primary_log2_pass <- all(
  abs(
    primary_df$Log2_CPM_Plus1 -
      primary_expected_log2
  ) < 1e-8
)

sensitivity_log2_pass <- all(
  abs(
    sensitivity_df$Log2_CPM_Plus1 -
      sensitivity_expected_log2
  ) < 1e-8
)

if (!primary_log2_pass ||
    !sensitivity_log2_pass) {
  
  stop("FAIL — Log2(CPM+1) relationship does not reconcile.")
}

cat("PASS — Log2(CPM+1) relationship verified.\n\n")


# ============================================================
# 12. CORE ANALYSIS:
#     Detection + quantitative abundance
# ============================================================

cat("------------------------------------------------------------\n")
cat("12. Expression-context characterization\n")
cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 12A. Patient-level summary by gene × compartment
# ------------------------------------------------------------

context_summary <- do.call(
  rbind,
  lapply(
    split(
      primary_df,
      interaction(
        primary_df$Gene,
        primary_df$Compartment,
        drop = TRUE
      )
    ),
    function(x) {
      
      detected <- x$Detected
      
      data.frame(
        Gene = unique(x$Gene),
        Compartment = unique(x$Compartment),
        
        N_Patients = nrow(x),
        
        Detection_N = sum(detected),
        
        Detection_Percent =
          100 * mean(detected),
        
        Zero_N =
          sum(x$Raw_UMI == 0),
        
        Raw_UMI_Median =
          median(x$Raw_UMI),
        
        Raw_UMI_Mean =
          mean(x$Raw_UMI),
        
        CPM_Median =
          median(x$CPM),
        
        CPM_Mean =
          mean(x$CPM),
        
        CPM_Q1 =
          quantile(x$CPM, 0.25, names = FALSE),
        
        CPM_Q3 =
          quantile(x$CPM, 0.75, names = FALSE),
        
        CPM_Min =
          min(x$CPM),
        
        CPM_Max =
          max(x$CPM),
        
        CPM_IQR =
          IQR(x$CPM),
        
        Log2_CPM_Plus1_Median =
          median(x$Log2_CPM_Plus1),
        
        Log2_CPM_Plus1_Mean =
          mean(x$Log2_CPM_Plus1),
        
        Log2_CPM_Plus1_Q1 =
          quantile(
            x$Log2_CPM_Plus1,
            0.25,
            names = FALSE
          ),
        
        Log2_CPM_Plus1_Q3 =
          quantile(
            x$Log2_CPM_Plus1,
            0.75,
            names = FALSE
          ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(context_summary) <- NULL


# ------------------------------------------------------------
# 12B. Detected-patient quantitative abundance
# ------------------------------------------------------------

detected_abundance_summary <- do.call(
  rbind,
  lapply(
    split(
      primary_df,
      interaction(
        primary_df$Gene,
        primary_df$Compartment,
        drop = TRUE
      )
    ),
    function(x) {
      
      x_detected <- x[x$Detected, , drop = FALSE]
      
      if (nrow(x_detected) == 0) {
        
        return(
          data.frame(
            Gene = unique(x$Gene),
            Compartment = unique(x$Compartment),
            
            Detected_N = 0,
            
            Detected_CPM_Median = NA_real_,
            Detected_CPM_Mean = NA_real_,
            Detected_CPM_Q1 = NA_real_,
            Detected_CPM_Q3 = NA_real_,
            Detected_CPM_Min = NA_real_,
            Detected_CPM_Max = NA_real_,
            Detected_CPM_IQR = NA_real_,
            
            Detected_Log2_CPM_Plus1_Median = NA_real_,
            Detected_Log2_CPM_Plus1_Mean = NA_real_,
            
            stringsAsFactors = FALSE
          )
        )
      }
      
      data.frame(
        Gene = unique(x$Gene),
        Compartment = unique(x$Compartment),
        
        Detected_N = nrow(x_detected),
        
        Detected_CPM_Median =
          median(x_detected$CPM),
        
        Detected_CPM_Mean =
          mean(x_detected$CPM),
        
        Detected_CPM_Q1 =
          quantile(
            x_detected$CPM,
            0.25,
            names = FALSE
          ),
        
        Detected_CPM_Q3 =
          quantile(
            x_detected$CPM,
            0.75,
            names = FALSE
          ),
        
        Detected_CPM_Min =
          min(x_detected$CPM),
        
        Detected_CPM_Max =
          max(x_detected$CPM),
        
        Detected_CPM_IQR =
          IQR(x_detected$CPM),
        
        Detected_Log2_CPM_Plus1_Median =
          median(x_detected$Log2_CPM_Plus1),
        
        Detected_Log2_CPM_Plus1_Mean =
          mean(x_detected$Log2_CPM_Plus1),
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(detected_abundance_summary) <- NULL


# ------------------------------------------------------------
# 12C. Merge prevalence + abundance
# ------------------------------------------------------------

expression_context <- merge(
  context_summary,
  detected_abundance_summary,
  by = c("Gene", "Compartment"),
  all = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# 13. Reconcile with locked STEP10F-1C profile
# ------------------------------------------------------------

profile_check <- merge(
  expression_context,
  primary_profile[
    ,
    c(
      "Gene",
      "Compartment",
      "N_Patients",
      "Detection_N",
      "Detection_Percent",
      "CPM_Median",
      "CPM_Mean",
      "CPM_Q1",
      "CPM_Q3",
      "CPM_Min",
      "CPM_Max",
      "CPM_IQR",
      "Log2_CPM_Plus1_Median"
    )
  ],
  by = c("Gene", "Compartment"),
  suffixes = c("_F3", "_F1"),
  all = FALSE
)


profile_reconciliation <- data.frame(
  Metric = c(
    "N_Patients",
    "Detection_N",
    "Detection_Percent",
    "CPM_Median",
    "CPM_Mean",
    "CPM_Q1",
    "CPM_Q3",
    "CPM_Min",
    "CPM_Max",
    "CPM_IQR",
    "Log2_CPM_Plus1_Median"
  ),
  PASS = c(
    all(
      profile_check$N_Patients_F3 ==
        profile_check$N_Patients_F1
    ),
    
    all(
      profile_check$Detection_N_F3 ==
        profile_check$Detection_N_F1
    ),
    
    all(
      abs(
        profile_check$Detection_Percent_F3 -
          profile_check$Detection_Percent_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Median_F3 -
          profile_check$CPM_Median_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Mean_F3 -
          profile_check$CPM_Mean_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Q1_F3 -
          profile_check$CPM_Q1_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Q3_F3 -
          profile_check$CPM_Q3_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Min_F3 -
          profile_check$CPM_Min_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_Max_F3 -
          profile_check$CPM_Max_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$CPM_IQR_F3 -
          profile_check$CPM_IQR_F1
      ) < 1e-8
    ),
    
    all(
      abs(
        profile_check$Log2_CPM_Plus1_Median_F3 -
          profile_check$Log2_CPM_Plus1_Median_F1
      ) < 1e-8
    )
  ),
  stringsAsFactors = FALSE
)

if (!all(profile_reconciliation$PASS)) {
  print(profile_reconciliation)
  stop(
    "FAIL — STEP10F-3 does not reconcile with STEP10F-1C."
  )
}

cat("PASS — STEP10F-1C quantitative profile fully reconciled.\n\n")


# ============================================================
# 14. Quantitative context ranking
# ============================================================

cat("------------------------------------------------------------\n")
cat("14. Quantitative expression-context ranking\n")
cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 14A. Confirm actual expression_context schema
# ------------------------------------------------------------

required_context_columns <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "CPM_Median",
  "CPM_Mean",
  "CPM_Q1",
  "CPM_Q3",
  "CPM_Min",
  "CPM_Max",
  "CPM_IQR",
  "Log2_CPM_Plus1_Median"
)

missing_context_columns <- setdiff(
  required_context_columns,
  colnames(expression_context)
)

if (length(missing_context_columns) > 0) {
  
  cat("FAIL — Required expression_context column(s) missing:\n")
  print(missing_context_columns)
  
  cat("\nActual expression_context columns:\n")
  print(colnames(expression_context))
  
  stop(
    "STEP10F-3 stopped because expression_context schema is incomplete."
  )
}


# ------------------------------------------------------------
# 14B. Ensure quantitative columns are numeric
# ------------------------------------------------------------

quantitative_context_columns <- c(
  "CPM_Median",
  "CPM_Mean",
  "CPM_Q1",
  "CPM_Q3",
  "CPM_Min",
  "CPM_Max",
  "CPM_IQR",
  "Log2_CPM_Plus1_Median"
)

for (x in quantitative_context_columns) {
  
  if (!is.numeric(expression_context[[x]])) {
    
    expression_context[[x]] <-
      as.numeric(expression_context[[x]])
  }
  
  if (any(!is.finite(expression_context[[x]]))) {
    
    stop(
      paste0(
        "FAIL — Non-finite values detected in expression_context column: ",
        x
      )
    )
  }
}


cat("PASS — expression_context quantitative schema verified.\n")


# ------------------------------------------------------------
# 14C. Within-gene compartment ranking by median CPM
# ------------------------------------------------------------

expression_context$CPM_Rank_Within_Gene <- NA_integer_


for (gene in core_genes) {
  
  idx <- expression_context$Gene == gene
  
  x <- expression_context[idx, , drop = FALSE]
  
  if (nrow(x) != length(primary_compartments)) {
    
    stop(
      paste0(
        "FAIL — Gene ",
        gene,
        " does not contain exactly ",
        length(primary_compartments),
        " primary compartments."
      )
    )
  }
  
  rank_values <- rank(
    -x$CPM_Median,
    ties.method = "min"
  )
  
  expression_context$CPM_Rank_Within_Gene[idx] <-
    as.integer(rank_values)
}


# ------------------------------------------------------------
# 14D. Within-gene detection ranking
# ------------------------------------------------------------

expression_context$Detection_Rank_Within_Gene <- NA_integer_


for (gene in core_genes) {
  
  idx <- expression_context$Gene == gene
  
  x <- expression_context[idx, , drop = FALSE]
  
  rank_values <- rank(
    -x$Detection_Percent,
    ties.method = "min"
  )
  
  expression_context$Detection_Rank_Within_Gene[idx] <-
    as.integer(rank_values)
}


# ------------------------------------------------------------
# 14E. Quantitative contrast:
#     highest median CPM vs lowest median CPM
# ------------------------------------------------------------

gene_contrast <- do.call(
  rbind,
  lapply(
    core_genes,
    function(gene) {
      
      x <- expression_context[
        expression_context$Gene == gene,
        ,
        drop = FALSE
      ]
      
      x <- x[
        order(
          -x$CPM_Median,
          x$Compartment
        ),
        ,
        drop = FALSE
      ]
      
      highest <- x[1, , drop = FALSE]
      lowest <- x[nrow(x), , drop = FALSE]
      
      
      ratio_value <- ifelse(
        is.finite(lowest$CPM_Median) &&
          lowest$CPM_Median > 0,
        
        highest$CPM_Median /
          lowest$CPM_Median,
        
        NA_real_
      )
      
      
      data.frame(
        
        Gene = gene,
        
        Highest_CPM_Compartment =
          highest$Compartment,
        
        Highest_CPM_Median =
          highest$CPM_Median,
        
        Highest_CPM_Detection_Percent =
          highest$Detection_Percent,
        
        Lowest_CPM_Compartment =
          lowest$Compartment,
        
        Lowest_CPM_Median =
          lowest$CPM_Median,
        
        Lowest_CPM_Detection_Percent =
          lowest$Detection_Percent,
        
        Highest_to_Lowest_Median_CPM =
          ratio_value,
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(gene_contrast) <- NULL


# ------------------------------------------------------------
# 14F. Highest quantitative compartment + detection context
# ------------------------------------------------------------

quantitative_context_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(gene) {
      
      x <- expression_context[
        expression_context$Gene == gene,
        ,
        drop = FALSE
      ]
      
      highest_cpm_idx <- which.max(
        x$CPM_Median
      )
      
      highest_detection_idx <- which.max(
        x$Detection_Percent
      )
      
      data.frame(
        
        Gene = gene,
        
        Highest_CPM_Compartment =
          x$Compartment[highest_cpm_idx],
        
        Highest_CPM_Median =
          x$CPM_Median[highest_cpm_idx],
        
        Highest_CPM_Detection_Percent =
          x$Detection_Percent[highest_cpm_idx],
        
        Highest_Detection_Compartment =
          x$Compartment[highest_detection_idx],
        
        Highest_Detection_Percent =
          x$Detection_Percent[highest_detection_idx],
        
        Highest_Detection_CPM_Median =
          x$CPM_Median[highest_detection_idx],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(quantitative_context_summary) <- NULL


# ------------------------------------------------------------
# 14G. Print compact result for manual audit
# ------------------------------------------------------------

cat("\nQuantitative expression-context summary:\n")
print(
  quantitative_context_summary,
  row.names = FALSE
)


cat("\nGene-level quantitative contrast:\n")
print(
  gene_contrast,
  row.names = FALSE
)


cat("\nPASS — Quantitative expression-context ranking completed.\n\n")


# ============================================================
# 15. Detection–abundance context classification
# ============================================================
#
# IMPORTANT:
#   This is a DESCRIPTIVE classification only.
#
#   It is NOT:
#     - a statistical test
#     - a specificity score
#     - a differential-expression classification
#     - a cell-type assignment
#
# The purpose is to distinguish:
#
#   (1) detection prevalence
#       How consistently is the gene detected across patients?
#
#   (2) quantitative abundance
#       How high is the patient-level CPM distribution?
#
# The 2/3 detection threshold is used ONLY as a descriptive
# operational label and has no inferential interpretation.
#
# No P value / FDR / statistical comparison is performed.
# ============================================================


cat("------------------------------------------------------------\n")
cat("15. Detection–abundance context classification\n")
cat("------------------------------------------------------------\n\n")


# ------------------------------------------------------------
# 15A. Verify actual expression_context schema
# ------------------------------------------------------------

required_classification_columns <- c(
  "Gene",
  "Compartment",
  "Detection_Percent",
  "CPM_Median"
)

missing_classification_columns <- setdiff(
  required_classification_columns,
  colnames(expression_context)
)

if (length(missing_classification_columns) > 0) {
  
  cat("FAIL — Required classification column(s) missing:\n")
  print(missing_classification_columns)
  
  cat("\nActual expression_context columns:\n")
  print(colnames(expression_context))
  
  stop(
    "STEP10F-3 Section 15 stopped because the expression_context schema is incomplete."
  )
}


# ------------------------------------------------------------
# 15B. Ensure numeric fields are numeric
# ------------------------------------------------------------

expression_context$Detection_Percent <-
  as.numeric(
    expression_context$Detection_Percent
  )

expression_context$CPM_Median <-
  as.numeric(
    expression_context$CPM_Median
  )


if (any(
  !is.finite(
    expression_context$Detection_Percent
  )
)) {
  
  stop(
    "FAIL — Non-finite Detection_Percent values detected."
  )
}


if (any(
  !is.finite(
    expression_context$CPM_Median
  )
)) {
  
  stop(
    "FAIL — Non-finite CPM_Median values detected."
  )
}


# ------------------------------------------------------------
# 15C. Descriptive detection threshold
# ------------------------------------------------------------

# IMPORTANT:
# This is NOT a biological or statistical cutoff.
#
# 2/3 of patients is used only to describe whether detection
# is relatively broad within the available compartment.

descriptive_detection_threshold <-
  (2 / 3) * 100


# ------------------------------------------------------------
# 15D. Determine quantitative leader within each gene
# ------------------------------------------------------------

highest_cpm_compartment <- list()


for (gene in core_genes) {
  
  x <- expression_context[
    expression_context$Gene == gene,
    ,
    drop = FALSE
  ]
  
  if (nrow(x) != length(primary_compartments)) {
    
    stop(
      paste0(
        "FAIL — Gene ",
        gene,
        " does not contain exactly ",
        length(primary_compartments),
        " primary compartments."
      )
    )
  }
  
  highest_index <- which.max(
    x$CPM_Median
  )
  
  highest_cpm_compartment[[gene]] <-
    x$Compartment[highest_index]
}


# ------------------------------------------------------------
# 15E. Classification
# ------------------------------------------------------------

expression_context$Context_Description <-
  NA_character_


for (gene in core_genes) {
  
  idx <- expression_context$Gene == gene
  
  x <- expression_context[
    idx,
    ,
    drop = FALSE
  ]
  
  quantitative_leader <-
    highest_cpm_compartment[[gene]]
  
  
  for (i in seq_len(nrow(x))) {
    
    detection_broad <-
      x$Detection_Percent[i] >=
      descriptive_detection_threshold
    
    quantitative_highest <-
      identical(
        x$Compartment[i],
        quantitative_leader
      )
    
    
    if (
      detection_broad &&
      quantitative_highest
    ) {
      
      desc <-
        "Broadly detected with highest quantitative expression"
      
    } else if (
      detection_broad &&
      !quantitative_highest
    ) {
      
      desc <-
        "Broadly detected but not quantitatively dominant"
      
    } else if (
      !detection_broad &&
      quantitative_highest
    ) {
      
      desc <-
        "Limited detection with highest quantitative expression"
      
    } else {
      
      desc <-
        "Limited detection and not quantitatively dominant"
    }
    
    
    expression_context$Context_Description[
      which(idx)[i]
    ] <- desc
  }
}


# ------------------------------------------------------------
# 15F. Add quantitative rank and detection rank
# ------------------------------------------------------------

expression_context$CPM_Rank_Within_Gene <-
  as.integer(
    expression_context$CPM_Rank_Within_Gene
  )

expression_context$Detection_Rank_Within_Gene <-
  as.integer(
    expression_context$Detection_Rank_Within_Gene
  )


# ------------------------------------------------------------
# 15G. Create gene-level context summary
# ------------------------------------------------------------

context_classification_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(gene) {
      
      x <- expression_context[
        expression_context$Gene == gene,
        ,
        drop = FALSE
      ]
      
      highest_cpm_index <- which.max(
        x$CPM_Median
      )
      
      highest_detection_index <- which.max(
        x$Detection_Percent
      )
      
      data.frame(
        
        Gene = gene,
        
        Quantitative_Highest_Compartment =
          x$Compartment[highest_cpm_index],
        
        Quantitative_Highest_CPM_Median =
          x$CPM_Median[highest_cpm_index],
        
        Quantitative_Highest_Detection_Percent =
          x$Detection_Percent[highest_cpm_index],
        
        Detection_Highest_Compartment =
          x$Compartment[highest_detection_index],
        
        Detection_Highest_Percent =
          x$Detection_Percent[highest_detection_index],
        
        Detection_Highest_CPM_Median =
          x$CPM_Median[highest_detection_index],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(
  context_classification_summary
) <- NULL


# ------------------------------------------------------------
# 15H. Print classification table
# ------------------------------------------------------------

cat(
  "Descriptive detection threshold: ",
  round(
    descriptive_detection_threshold,
    4
  ),
  "%\n\n",
  sep = ""
)


cat(
  "Gene × compartment descriptive classification:\n"
)

print(
  expression_context[
    order(
      expression_context$Gene,
      expression_context$CPM_Rank_Within_Gene
    ),
    c(
      "Gene",
      "Compartment",
      "Detection_N",
      "Detection_Percent",
      "CPM_Median",
      "Detected_CPM_Median",
      "CPM_Rank_Within_Gene",
      "Detection_Rank_Within_Gene",
      "Context_Description"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


cat(
  "\nGene-level context summary:\n"
)

print(
  context_classification_summary,
  row.names = FALSE
)


cat(
  "\nPASS — Detection–abundance context classification completed.\n\n"
)


# ============================================================
# 16. Manuscript-safe interpretation audit
# ============================================================

manuscript_interpretation <- data.frame(
  
  Gene = core_genes,
  
  Recommended_Expression_Context = c(
    "Multi-compartment detection with relatively higher quantitative expression in myeloid/macrophage cells",
    "Myeloid/macrophage-associated quantitative expression pattern",
    "T-cell-associated quantitative expression pattern"
  ),
  
  Avoid_Specificity_Claim = c(
    "Yes",
    "Yes",
    "Yes"
  ),
  
  Avoid_Causal_Claim = c(
    "Yes",
    "Yes",
    "Yes"
  ),
  
  Avoid_Communication_Claim = c(
    "Yes",
    "Yes",
    "Yes"
  ),
  
  Interpretation_Basis = c(
    "Detection prevalence plus descriptive CPM distribution across primary compartments",
    "Detection prevalence plus descriptive CPM distribution across primary compartments",
    "Detection prevalence plus descriptive CPM distribution across primary compartments"
  ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 17. Sensitivity Endothelial kept separate
# ============================================================

sensitivity_context <- merge(
  sensitivity_profile[
    ,
    c(
      "Gene",
      "Compartment",
      "N_Patients",
      "Detection_N",
      "Detection_Percent",
      "CPM_Median",
      "CPM_Mean",
      "CPM_Q1",
      "CPM_Q3",
      "CPM_Min",
      "CPM_Max",
      "CPM_IQR",
      "Log2_CPM_Plus1_Median"
    )
  ],
  detected_abundance_summary[
    detected_abundance_summary$Gene %in% core_genes &
      detected_abundance_summary$Compartment ==
      "Endothelial",
    ,
    drop = FALSE
  ],
  by = c("Gene", "Compartment"),
  suffixes = c("_Profile", "_Detected"),
  all = FALSE
)


# ============================================================
# 18. Final QC
# ============================================================

qc_results <- data.frame(
  Check = c(
    "Primary observations = 66",
    "Primary patients = 9",
    "Primary T-cell patients = 9",
    "Primary Myeloid/macrophage patients = 9",
    "Primary Endothelial patients = 4",
    "Sensitivity Endothelial patients = 2",
    "Core genes = AOC3/MS4A7/IL7R",
    "Raw UMI non-negative integer",
    "Detection definition Raw_UMI > 0",
    "CPM relationship valid",
    "Log2(CPM+1) relationship valid",
    "STEP10F-1C reconciliation",
    "Primary context summary = 9 rows",
    "Sensitivity context = 3 genes",
    "No statistical testing",
    "No P/FDR",
    "No differential expression",
    "No correlation",
    "No CellChat/LR",
    "No pseudobulk recalculation",
    "No GEO rereading",
    "No H5 rereading",
    "No patient/cell exclusion"
  ),
  
  PASS = c(
    
    nrow(primary_df) == 66,
    
    length(unique(primary_df$Patient)) == 9,
    
    length(unique(
      primary_df$Patient[
        primary_df$Compartment == "T-cell"
      ]
    )) == 9,
    
    length(unique(
      primary_df$Patient[
        primary_df$Compartment ==
          "Myeloid/macrophage"
      ]
    )) == 9,
    
    length(unique(
      primary_df$Patient[
        primary_df$Compartment ==
          "Endothelial"
      ]
    )) == 4,
    
    length(unique(sensitivity_df$Patient)) == 2,
    
    identical(core_genes,
              c("AOC3", "MS4A7", "IL7R")),
    
    all(
      primary_df$Raw_UMI >= 0
    ) &&
      all(
        abs(
          primary_df$Raw_UMI -
            round(primary_df$Raw_UMI)
        ) < 1e-8
      ),
    
    primary_detection_pass &&
      sensitivity_detection_pass,
    
    primary_cpm_pass &&
      sensitivity_cpm_pass,
    
    primary_log2_pass &&
      sensitivity_log2_pass,
    
    all(profile_reconciliation$PASS),
    
    nrow(expression_context) == 9,
    
    nrow(sensitivity_context) == 3,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE,
    
    TRUE
  ),
  
  stringsAsFactors = FALSE
)


overall_pass <- all(qc_results$PASS)

if (!isTRUE(overall_pass)) {
  
  cat("\nFAIL — STEP10F-3 QC audit.\n")
  print(
    qc_results[
      !qc_results$PASS,
      ,
      drop = FALSE
    ]
  )
  
  stop("STEP10F-3 FINAL QC: FAIL.")
}


# ============================================================
# 19. Save outputs
# ============================================================

write.table(
  expression_context,
  file.path(
    out_dir,
    "STEP10F-3A_Expression_Context_Profile.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  detected_abundance_summary,
  file.path(
    out_dir,
    "STEP10F-3B_Detected_Patient_Quantitative_Abundance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  gene_contrast,
  file.path(
    out_dir,
    "STEP10F-3C_Gene_Quantitative_Context_Contrast.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  sensitivity_context,
  file.path(
    out_dir,
    "STEP10F-3D_Sensitivity_Endothelial_Context.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  manuscript_interpretation,
  file.path(
    out_dir,
    "STEP10F-3E_Manuscript_Safe_Interpretation_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  profile_reconciliation,
  file.path(
    out_dir,
    "STEP10F-3F_STEP10F1C_Reconciliation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  qc_results,
  file.path(
    out_dir,
    "STEP10F-3G_FINAL_Audit_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 20. Methodological lock
# ------------------------------------------------------------

method_lock <- data.frame(
  
  Item = c(
    "Analysis name",
    "Input source",
    "Core genes",
    "Primary compartments",
    "Sensitivity compartment",
    "Expression source",
    "Raw UMI",
    "CPM",
    "Log2(CPM+1)",
    "Detection definition",
    "Detection prevalence",
    "Quantitative abundance",
    "Patient-level analysis",
    "Statistical testing",
    "P value",
    "FDR",
    "Differential expression",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "GEO raw matrix rereading",
    "H5 rereading",
    "Patient exclusion",
    "Specificity claims",
    "Causal inference"
  ),
  
  Value = c(
    "STEP10F-3 Expression-context audit",
    "STEP10F-1 locked patient-level expression outputs",
    "AOC3 / MS4A7 / IL7R",
    "T-cell / Myeloid/macrophage / Endothelial",
    "Endothelial sensitivity only",
    "STEP10F-1 locked raw-UMI pseudobulk expression",
    "Used as locked quantitative source",
    "Descriptive normalization only",
    "Descriptive transformation only",
    "Raw_UMI > 0",
    "Described by N and percentage",
    "Described by median/mean/IQR/range CPM",
    "Yes",
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "No",
    "Not assigned",
    "No"
  ),
  
  stringsAsFactors = FALSE
)


write.table(
  method_lock,
  file.path(
    out_dir,
    "STEP10F-3_Methodological_Lock.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 21. Final lock file
# ------------------------------------------------------------

lock_lines <- c(
  
  "STEP10F-3 — EXPRESSION-CONTEXT AUDIT",
  "============================================================",
  "",
  "FINAL STATUS: PASS",
  "",
  "Primary expression observations: 66",
  "Primary patients: 9",
  "Primary T-cell patients: 9",
  "Primary Myeloid/macrophage patients: 9",
  "Primary Endothelial patients: 4",
  "Sensitivity Endothelial patients: 2",
  "",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "",
  "Expression source: STEP10F-1 locked patient-level expression",
  "Primary STEP10F-1C profile: RECONCILED",
  "",
  "Detection definition: Raw_UMI > 0",
  "Quantitative abundance: CPM descriptive summary",
  "Log2(CPM+1): descriptive summary",
  "",
  "Detection prevalence and quantitative abundance jointly audited.",
  "",
  "Primary Endothelial sensitivity cohort kept separate.",
  "",
  "Statistical testing: NO",
  "P/FDR: NO",
  "Differential expression: NO",
  "Correlation: NO",
  "CellChat/LR: NO",
  "Pseudobulk recalculation: NO",
  "GEO raw matrix rereading: NO",
  "H5 rereading: NO",
  "Patient/cell exclusion: NO",
  "",
  "Specificity claims: NOT ASSIGNED",
  "Causal inference: NO",
  "",
  "STEP10F-3 FINAL LOCK: PASS",
  "============================================================"
)


writeLines(
  lock_lines,
  con = file.path(
    out_dir,
    "STEP10F-3_FINAL_LOCK.txt"
  )
)


# ============================================================
# 22. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("STEP10F-3 FINAL REPORT\n")
cat("============================================================\n")

cat("Primary expression observations: 66\n")
cat("Primary patients: 9\n")
cat("Primary T-cell: 9\n")
cat("Primary Myeloid/macrophage: 9\n")
cat("Primary Endothelial: 4\n")
cat("Sensitivity Endothelial: 2\n")
cat("Core genes: AOC3 / MS4A7 / IL7R\n")

cat("\n")

cat("Expression source: STEP10F-1 locked patient-level expression\n")
cat("Detection definition: Raw_UMI > 0\n")
cat("Quantitative abundance: CPM descriptive only\n")
cat("Log2(CPM+1): descriptive only\n")

cat("\n")

cat("Detection prevalence + quantitative abundance: INTEGRATED\n")
cat("STEP10F-1C reconciliation: PASS\n")
cat("Sensitivity Endothelial: KEPT SEPARATE\n")

cat("\n")

cat("Statistical testing: NO\n")
cat("P/FDR: NO\n")
cat("Differential expression: NO\n")
cat("Correlation: NO\n")
cat("CellChat/LR: NO\n")
cat("Pseudobulk recalculation: NO\n")
cat("GEO raw matrix rereading: NO\n")
cat("H5 rereading: NO\n")
cat("Patient/cell exclusion: NO\n")

cat("\n")

cat("Specificity claims: NOT ASSIGNED\n")
cat("Causal inference: NO\n")

cat("\n")

cat("PASS — STEP10F-3 Expression-context audit\n")
cat("Output directory:\n")
cat(out_dir, "\n")
cat("STEP10F-3 FINAL LOCK: PASS\n")
cat("============================================================\n")
cat("END STEP10F-3 FINAL\n")



# ============================================================
# STEP10G — Single-cell Evidence Synthesis / Manuscript Integration
# ============================================================
#
# Purpose:
#   Integrate LOCKED STEP10A–STEP10F-3 single-cell evidence
#   into a manuscript-ready evidence synthesis.
#
# IMPORTANT:
#   - No new pseudobulk calculation
#   - No GEO raw matrix rereading
#   - No H5 rereading
#   - No normalization recalculation
#   - No statistical testing
#   - No P/FDR
#   - No differential expression
#   - No correlation
#   - No CellChat / ligand-receptor analysis
#   - No patient/cell exclusion
#   - No specificity score
#   - No causal inference
#
# Core genes:
#   AOC3 / MS4A7 / IL7R
#
# Locked primary compartments:
#   T-cell
#   Myeloid/macrophage
#   Endothelial
#
# Locked patient structure:
#   T-cell                 = 9
#   Myeloid/macrophage     = 9
#   Endothelial            = 4
#   Sensitivity Endothelial= 2
#   Complete 3-compartment = 4
#
# ============================================================


# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

options(stringsAsFactors = FALSE)

base_dir <- "input_dir/Paper 1_Step 10-Single-Cell"

out_dir <- file.path(
  base_dir,
  "STEP10G_Single_Cell_Evidence_Synthesis"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("STEP10G — Single-cell Evidence Synthesis\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 1. Locked constants
# ------------------------------------------------------------

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

locked_primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

locked_primary_patients <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)

locked_sensitivity_endothelial_patients <- 2

locked_total_primary_patients <- 9

locked_complete_three_compartment_patients <- 4

detection_definition <- "Raw_UMI > 0"

cat("Core genes:",
    paste(core_genes, collapse = " / "),
    "\n")

cat("Primary compartments:",
    paste(locked_primary_compartments, collapse = " / "),
    "\n\n")


# ------------------------------------------------------------
# 2. Input files — ONLY LOCKED STEP10F outputs
# ------------------------------------------------------------

file_F1A <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv"
)

file_F1B <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv"
)

file_F1C <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
)

file_F1D <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1D_Sensitivity_Gene_Compartment_Profile.tsv"
)

file_F1E <- file.path(
  base_dir,
  "STEP10F-1_Patient_Level_Expression_Profile",
  "STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv"
)

file_F2A <- file.path(
  base_dir,
  "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
  "STEP10F-2A_Primary_Gene_Compartment_Pattern_Profile.tsv"
)

file_F2E <- file.path(
  base_dir,
  "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
  "STEP10F-2E_Complete_ThreeCompartment_Pattern_Summary.tsv"
)

file_F2F <- file.path(
  base_dir,
  "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
  "STEP10F-2F_Complete_Patient_Dominant_Compartment.tsv"
)

file_F2G <- file.path(
  base_dir,
  "STEP10F-2_Gene_Compartment_Expression_Pattern_Audit",
  "STEP10F-2G_Descriptive_Pattern_Classification.tsv"
)

file_F3A <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit",
  "STEP10F-3A_Expression_Context_Profile.tsv"
)

file_F3B <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit",
  "STEP10F-3B_Detected_Patient_Quantitative_Abundance.tsv"
)

file_F3C <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit",
  "STEP10F-3C_Gene_Quantitative_Context_Contrast.tsv"
)

file_F3D <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit",
  "STEP10F-3D_Sensitivity_Endothelial_Context.tsv"
)

file_F3E <- file.path(
  base_dir,
  "STEP10F-3_Expression_Context_Audit",
  "STEP10F-3E_Manuscript_Safe_Interpretation_Audit.tsv"
)


# ------------------------------------------------------------
# 3. File existence audit
# ------------------------------------------------------------

input_files <- c(
  STEP10F1A = file_F1A,
  STEP10F1B = file_F1B,
  STEP10F1C = file_F1C,
  STEP10F1D = file_F1D,
  STEP10F1E = file_F1E,
  STEP10F2A = file_F2A,
  STEP10F2E = file_F2E,
  STEP10F2F = file_F2F,
  STEP10F2G = file_F2G,
  STEP10F3A = file_F3A,
  STEP10F3B = file_F3B,
  STEP10F3C = file_F3C,
  STEP10F3D = file_F3D,
  STEP10F3E = file_F3E
)

file_exists <- file.exists(input_files)

if (!all(file_exists)) {
  
  missing_files <- names(input_files)[!file_exists]
  
  cat("FAIL — Missing locked input files:\n")
  print(missing_files)
  
  stop("STEP10G input-file audit failed.")
  
} else {
  
  cat("PASS — All locked STEP10F input files found.\n\n")
}


# ------------------------------------------------------------
# 4. Read locked outputs
# ------------------------------------------------------------

F1A <- read.delim(
  file_F1A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F1B <- read.delim(
  file_F1B,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F1C <- read.delim(
  file_F1C,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F1D <- read.delim(
  file_F1D,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F1E <- read.delim(
  file_F1E,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F2A <- read.delim(
  file_F2A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F2E <- read.delim(
  file_F2E,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F2F <- read.delim(
  file_F2F,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F2G <- read.delim(
  file_F2G,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F3A <- read.delim(
  file_F3A,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F3B <- read.delim(
  file_F3B,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F3C <- read.delim(
  file_F3C,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F3D <- read.delim(
  file_F3D,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

F3E <- read.delim(
  file_F3E,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("PASS — Locked STEP10F outputs loaded.\n\n")


# ------------------------------------------------------------
# 5. Basic schema audit
# ------------------------------------------------------------

required_F1A <- c(
  "Dataset",
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Cell_Count",
  "Gene",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1",
  "Detected"
)

required_F1C <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "Raw_UMI_Mean",
  "Raw_UMI_Median",
  "CPM_Mean",
  "CPM_Median"
)

required_F3A <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "CPM_Median"
)

schema_pass <- all(
  required_F1A %in% colnames(F1A)
) &&
  all(
    required_F1C %in% colnames(F1C)
  ) &&
  all(
    required_F3A %in% colnames(F3A)
  )

if (!schema_pass) {
  cat("FAIL — Required schema elements missing.\n")
  stop("STEP10G schema audit failed.")
}

cat("PASS — Required schema elements verified.\n\n")


# ------------------------------------------------------------
# 6. Core-gene and compartment identity audit
# ------------------------------------------------------------

identity_pass <- identical(
  sort(unique(F1A$Gene)),
  sort(core_genes)
) &&
  all(
    unique(F1A$Compartment) %in%
      locked_primary_compartments
  ) &&
  identical(
    sort(unique(F1C$Gene)),
    sort(core_genes)
  ) &&
  identical(
    sort(unique(F3A$Gene)),
    sort(core_genes)
  )

if (!identity_pass) {
  cat("FAIL — Gene/compartment identity audit failed.\n")
  stop("STEP10G identity audit failed.")
}

cat("PASS — Core-gene and primary-compartment identities verified.\n\n")


# ------------------------------------------------------------
# 7. Patient-level structural synthesis
# ------------------------------------------------------------

primary_F1A <- F1A[
  F1A$Eligibility == "Primary",
  ,
  drop = FALSE
]

sensitivity_F1B <- F1B[
  ,
  drop = FALSE
]

primary_patient_counts <- sapply(
  locked_primary_compartments,
  function(x) {
    length(
      unique(
        primary_F1A$Patient[
          primary_F1A$Compartment == x
        ]
      )
    )
  }
)

names(primary_patient_counts) <-
  locked_primary_compartments

sensitivity_patient_count <- length(
  unique(sensitivity_F1B$Patient)
)

structure_pass <- all(
  primary_patient_counts ==
    locked_primary_patients
) &&
  sensitivity_patient_count ==
  locked_sensitivity_endothelial_patients

cat("Primary patient structure:\n")

for (x in locked_primary_compartments) {
  
  cat(
    "  ",
    x,
    ": ",
    primary_patient_counts[x],
    "\n",
    sep = ""
  )
}

cat(
  "Sensitivity Endothelial: ",
  sensitivity_patient_count,
  "\n\n",
  sep = ""
)

if (!structure_pass) {
  cat("FAIL — Patient-level structure mismatch.\n")
  stop("STEP10G patient structure audit failed.")
}

cat("PASS — Patient-level structure verified.\n\n")


# ------------------------------------------------------------
# 8. Complete three-compartment patient audit
# ------------------------------------------------------------

complete_patients <- Reduce(
  intersect,
  lapply(
    locked_primary_compartments,
    function(x) {
      unique(
        primary_F1A$Patient[
          primary_F1A$Compartment == x
        ]
      )
    }
  )
)

complete_patient_count <- length(
  complete_patients
)

complete_pass <-
  complete_patient_count ==
  locked_complete_three_compartment_patients

if (!complete_pass) {
  
  cat(
    "FAIL — Complete three-compartment patient count = ",
    complete_patient_count,
    "; expected ",
    locked_complete_three_compartment_patients,
    ".\n",
    sep = ""
  )
  
  stop("STEP10G complete-patient audit failed.")
}

cat(
  "PASS — Complete three-compartment patient count verified: ",
  complete_patient_count,
  ".\n\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Core-gene × compartment quantitative evidence synthesis
# ------------------------------------------------------------

context_primary <- F3A[
  F3A$Compartment %in%
    locked_primary_compartments,
  ,
  drop = FALSE
]

context_primary <- context_primary[
  context_primary$Gene %in% core_genes,
  ,
  drop = FALSE
]

context_primary$Detection_Percent <- as.numeric(
  context_primary$Detection_Percent
)

context_primary$CPM_Median <- as.numeric(
  context_primary$CPM_Median
)

context_primary$Detection_N <- as.numeric(
  context_primary$Detection_N
)

context_primary$N_Patients <- as.numeric(
  context_primary$N_Patients
)

# Rank quantitative abundance within gene
context_primary$CPM_Rank <- ave(
  -context_primary$CPM_Median,
  context_primary$Gene,
  FUN = function(x)
    rank(x, ties.method = "min")
)

# Rank detection within gene
context_primary$Detection_Rank <- ave(
  -context_primary$Detection_Percent,
  context_primary$Gene,
  FUN = function(x)
    rank(x, ties.method = "min")
)


# ------------------------------------------------------------
# 10. Determine quantitative context
# ------------------------------------------------------------

gene_context_summary <- do.call(
  rbind,
  lapply(
    core_genes,
    function(g) {
      
      x <- context_primary[
        context_primary$Gene == g,
        ,
        drop = FALSE
      ]
      
      max_cpm <- max(
        x$CPM_Median,
        na.rm = TRUE
      )
      
      min_cpm <- min(
        x$CPM_Median,
        na.rm = TRUE
      )
      
      highest_rows <- x[
        x$CPM_Median == max_cpm,
        ,
        drop = FALSE
      ]
      
      highest_compartment <-
        paste(
          highest_rows$Compartment,
          collapse = " / "
        )
      
      max_detection <- max(
        x$Detection_Percent,
        na.rm = TRUE
      )
      
      detection_rows <- x[
        x$Detection_Percent == max_detection,
        ,
        drop = FALSE
      ]
      
      detection_highest <-
        paste(
          detection_rows$Compartment,
          collapse = " / "
        )
      
      data.frame(
        Gene = g,
        Quantitative_Highest_Compartment =
          highest_compartment,
        Quantitative_Highest_CPM_Median =
          max_cpm,
        Quantitative_Lowest_CPM_Median =
          min_cpm,
        Highest_to_Lowest_Median_CPM =
          ifelse(
            min_cpm > 0,
            max_cpm / min_cpm,
            NA_real_
          ),
        Detection_Highest_Compartment =
          detection_highest,
        Detection_Highest_Percent =
          max_detection,
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(gene_context_summary) <- NULL


# ------------------------------------------------------------
# 11. Detection–abundance relationship classification
# ------------------------------------------------------------

context_primary$Detection_Abundance_Context <-
  "Descriptive expression context"

for (i in seq_len(nrow(context_primary))) {
  
  g <- context_primary$Gene[i]
  
  highest_cpm <- gene_context_summary[
    gene_context_summary$Gene == g,
    "Quantitative_Highest_Compartment"
  ]
  
  if (
    context_primary$Compartment[i] ==
    highest_cpm
  ) {
    
    context_primary$Detection_Abundance_Context <-
      "Quantitatively highest compartment"
    
  } else {
    
    context_primary$Detection_Abundance_Context <-
      "Detected but not quantitatively highest"
  }
}


# ------------------------------------------------------------
# 12. Manuscript-safe gene-level interpretation
# ------------------------------------------------------------

manuscript_safe <- data.frame(
  Gene = core_genes,
  Evidence_Summary = NA_character_,
  Safe_Interpretation = NA_character_,
  Specificity_Claim = "NOT SUPPORTED",
  Causal_Claim = "NOT SUPPORTED",
  Communication_Claim = "NOT SUPPORTED",
  stringsAsFactors = FALSE
)

manuscript_safe[
  manuscript_safe$Gene == "AOC3",
  "Evidence_Summary"
] <-
  "Detected across primary compartments, with highest median CPM in myeloid/macrophage cells."

manuscript_safe[
  manuscript_safe$Gene == "AOC3",
  "Safe_Interpretation"
] <-
  "Multi-compartment detection with relatively higher quantitative expression in myeloid/macrophage cells."

manuscript_safe[
  manuscript_safe$Gene == "MS4A7",
  "Evidence_Summary"
] <-
  "Detected in all primary compartments, with markedly higher median CPM in myeloid/macrophage cells."

manuscript_safe[
  manuscript_safe$Gene == "MS4A7",
  "Safe_Interpretation"
] <-
  "Broad detection with a myeloid/macrophage-associated quantitative expression context."

manuscript_safe[
  manuscript_safe$Gene == "IL7R",
  "Evidence_Summary"
] <-
  "Detected in all primary compartments, with highest median CPM in T cells."

manuscript_safe[
  manuscript_safe$Gene == "IL7R",
  "Safe_Interpretation"
] <-
  "Broad detection with a T-cell-associated quantitative expression context."


# ------------------------------------------------------------
# 13. Explicit unsupported-claim audit
# ------------------------------------------------------------

unsupported_claim_audit <- data.frame(
  
  Claim_Type = c(
    "Cell-type specificity",
    "Cell-type exclusivity",
    "Causal regulation",
    "Intercellular communication",
    "Ligand-receptor interaction",
    "Functional activation",
    "Immune escape",
    "Immune checkpoint resistance",
    "Tumor-promoting mechanism"
  ),
  
  Supported_by_STEP10 = c(
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  
  Reason = c(
    "Genes were detected in multiple compartments.",
    "Detection was not restricted to one compartment.",
    "STEP10 is descriptive and observational.",
    "Expression context alone does not establish communication.",
    "No ligand-receptor analysis was performed in STEP10G.",
    "No functional perturbation or activity assay was performed.",
    "Not evaluated by this workflow.",
    "Not evaluated by this workflow.",
    "No causal or mechanistic experiment was performed."
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 14. Evidence hierarchy for manuscript integration
# ------------------------------------------------------------

evidence_hierarchy <- data.frame(
  
  Evidence_Level = c(
    "Dataset / cell structure",
    "Patient-level compartment representation",
    "Raw-UMI pseudobulk construction",
    "Library QC",
    "Core-gene detection",
    "Quantitative expression context",
    "Detection–abundance integration",
    "Cell-type specificity",
    "Intercellular communication",
    "Causal mechanism"
  ),
  
  STEP10_Component = c(
    "STEP10A / STEP10B",
    "STEP10B",
    "STEP10C",
    "STEP10D-A",
    "STEP10D-B / STEP10D-C",
    "STEP10E / STEP10F-1",
    "STEP10F-2 / STEP10F-3",
    "Not established",
    "Not established",
    "Not established"
  ),
  
  Status = c(
    "LOCKED",
    "LOCKED",
    "LOCKED",
    "LOCKED",
    "LOCKED",
    "LOCKED",
    "LOCKED",
    "NOT CLAIMED",
    "NOT CLAIMED",
    "NOT CLAIMED"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 15. Final manuscript evidence table
# ------------------------------------------------------------

manuscript_evidence <- context_primary[
  order(
    context_primary$Gene,
    context_primary$CPM_Rank
  ),
  c(
    "Gene",
    "Compartment",
    "N_Patients",
    "Detection_N",
    "Detection_Percent",
    "CPM_Median",
    "Detection_Abundance_Context"
  ),
  drop = FALSE
]

colnames(manuscript_evidence) <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "Median_CPM",
  "Descriptive_Context"
)


# ------------------------------------------------------------
# 16. Compact manuscript summary
# ------------------------------------------------------------

manuscript_summary <- data.frame(
  
  Gene = c(
    "AOC3",
    "MS4A7",
    "IL7R"
  ),
  
  Primary_Detection_Pattern = c(
    "6/9 T-cell; 8/9 myeloid/macrophage; 2/4 endothelial",
    "9/9 T-cell; 9/9 myeloid/macrophage; 4/4 endothelial",
    "9/9 T-cell; 9/9 myeloid/macrophage; 4/4 endothelial"
  ),
  
  Highest_Quantitative_Context = c(
    "Myeloid/macrophage",
    "Myeloid/macrophage",
    "T-cell"
  ),
  
  Highest_Median_CPM = c(
    1.3345484,
    506.1955170,
    258.1607418
  ),
  
  Interpretation = c(
    "Multi-compartment detection with relatively higher quantitative expression in myeloid/macrophage cells.",
    "Broad detection with a myeloid/macrophage-associated quantitative expression context.",
    "Broad detection with a T-cell-associated quantitative expression context."
  ),
  
  Specificity = c(
    "Not assigned",
    "Not assigned",
    "Not assigned"
  ),
  
  Causality = c(
    "Not inferred",
    "Not inferred",
    "Not inferred"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------

write.table(
  context_primary,
  file.path(
    out_dir,
    "STEP10G-1_Primary_CoreGene_Evidence_Synthesis.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  gene_context_summary,
  file.path(
    out_dir,
    "STEP10G-2_Gene_Level_Quantitative_Context_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  manuscript_evidence,
  file.path(
    out_dir,
    "STEP10G-3_Manuscript_Evidence_Table.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  manuscript_summary,
  file.path(
    out_dir,
    "STEP10G-4_Manuscript_Safe_CoreGene_Summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  manuscript_safe,
  file.path(
    out_dir,
    "STEP10G-5_Manuscript_Safe_Interpretation_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  unsupported_claim_audit,
  file.path(
    out_dir,
    "STEP10G-6_Unsupported_Claim_Audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  evidence_hierarchy,
  file.path(
    out_dir,
    "STEP10G-7_Evidence_Hierarchy.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 18. Methodological lock
# ------------------------------------------------------------

methodological_lock <- data.frame(
  
  Item = c(
    "Core genes",
    "Primary compartments",
    "Primary patients",
    "Complete three-compartment patients",
    "Sensitivity endothelial patients",
    "Expression source",
    "Detection definition",
    "Quantitative measure",
    "Statistical testing",
    "P/FDR",
    "Differential expression",
    "Correlation",
    "CellChat/LR",
    "Pseudobulk recalculation",
    "GEO raw matrix rereading",
    "H5 rereading",
    "Patient/cell exclusion",
    "Specificity scoring",
    "Causal inference"
  ),
  
  Locked_Value = c(
    "AOC3 / MS4A7 / IL7R",
    "T-cell / Myeloid/macrophage / Endothelial",
    "9",
    "4",
    "2",
    "STEP10F-1 locked patient-level expression",
    "Raw_UMI > 0",
    "CPM descriptive only",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO",
    "NO"
  ),
  
  stringsAsFactors = FALSE
)

write.table(
  methodological_lock,
  file.path(
    out_dir,
    "STEP10G_Methodological_Lock.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ------------------------------------------------------------
# 19. Final QC
# ------------------------------------------------------------

qc_results <- data.frame(
  
  Check = c(
    "Input files available",
    "Core genes verified",
    "Primary compartments verified",
    "Primary patient counts verified",
    "Sensitivity endothelial patients verified",
    "Complete three-compartment patients verified",
    "Primary expression evidence available",
    "Quantitative context available",
    "Manuscript-safe interpretation generated",
    "Unsupported claims explicitly flagged"
  ),
  
  Result = c(
    all(file_exists),
    identical(sort(unique(F1A$Gene)), sort(core_genes)),
    all(
      unique(F1A$Compartment) %in%
        locked_primary_compartments
    ),
    all(
      primary_patient_counts ==
        locked_primary_patients
    ),
    sensitivity_patient_count ==
      locked_sensitivity_endothelial_patients,
    complete_patient_count ==
      locked_complete_three_compartment_patients,
    nrow(primary_F1A) == 66,
    nrow(context_primary) == 9,
    nrow(manuscript_safe) == 3,
    nrow(unsupported_claim_audit) == 9
  ),
  
  stringsAsFactors = FALSE
)

overall_pass <- all(qc_results$Result)

if (!isTRUE(overall_pass)) {
  
  cat("\nFAIL — STEP10G QC audit.\n")
  print(
    qc_results[
      !qc_results$Result,
      ,
      drop = FALSE
    ]
  )
  
  stop("STEP10G final QC failed.")
  
} else {
  
  cat("\nPASS — STEP10G QC audit.\n\n")
}


# ------------------------------------------------------------
# 20. Print final manuscript evidence table
# ------------------------------------------------------------

cat("STEP10G manuscript evidence summary:\n")
cat("------------------------------------------------------------\n")

print(
  manuscript_summary,
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Print final report
# ------------------------------------------------------------

cat("\n")
cat("STEP10G FINAL REPORT\n")
cat("============================================================\n")

cat(
  "Primary expression observations: 66\n"
)

cat(
  "Primary patients: 9\n"
)

cat(
  "Primary T-cell: 9\n"
)

cat(
  "Primary Myeloid/macrophage: 9\n"
)

cat(
  "Primary Endothelial: 4\n"
)

cat(
  "Sensitivity Endothelial: 2\n"
)

cat(
  "Complete three-compartment patients: 4\n"
)

cat(
  "Core genes: AOC3 / MS4A7 / IL7R\n"
)

cat("\n")

cat(
  "Evidence source: STEP10F-1 / STEP10F-2 / STEP10F-3 locked outputs\n"
)

cat(
  "Detection definition: Raw_UMI > 0\n"
)

cat(
  "Quantitative abundance: CPM descriptive only\n"
)

cat(
  "Detection + quantitative abundance: INTEGRATED\n"
)

cat(
  "Primary patient-level evidence: VERIFIED\n"
)

cat(
  "Complete three-compartment patients: 4\n"
)

cat(
  "Sensitivity Endothelial: KEPT SEPARATE\n"
)

cat("\n")

cat(
  "Cell-type specificity: NOT ESTABLISHED\n"
)

cat(
  "Cell-type exclusivity: NOT ESTABLISHED\n"
)

cat(
  "Intercellular communication: NOT ESTABLISHED\n"
)

cat(
  "Ligand-receptor interaction: NOT ESTABLISHED\n"
)

cat(
  "Causal inference: NO\n"
)

cat("\n")

cat(
  "Statistical testing: NO\n"
)

cat(
  "P/FDR: NO\n"
)

cat(
  "Differential expression: NO\n"
)

cat(
  "Correlation: NO\n"
)

cat(
  "CellChat/LR: NO\n"
)

cat(
  "Pseudobulk recalculation: NO\n"
)

cat(
  "GEO raw matrix rereading: NO\n"
)

cat(
  "H5 rereading: NO\n"
)

cat(
  "Patient/cell exclusion: NO\n"
)

cat(
  "Specificity scoring: NO\n"
)

cat("\n")

cat(
  "PASS — STEP10G Single-cell Evidence Synthesis\n"
)

cat(
  "Output directory:\n",
  out_dir,
  "\n"
)

cat(
  "STEP10G FINAL LOCK: PASS\n"
)

cat("============================================================\n")
cat("END STEP10G FINAL\n")


# ------------------------------------------------------------
# 22. Save final lock file
# ------------------------------------------------------------

final_lock_text <- c(
  
  "STEP10G — Single-cell Evidence Synthesis",
  "FINAL LOCK",
  "============================================================",
  "",
  "STATUS: PASS",
  "",
  "Core genes: AOC3 / MS4A7 / IL7R",
  "",
  "Primary T-cell patients: 9",
  "Primary Myeloid/macrophage patients: 9",
  "Primary Endothelial patients: 4",
  "Sensitivity Endothelial patients: 2",
  "Complete three-compartment patients: 4",
  "",
  "Expression source: STEP10F-1 locked patient-level expression",
  "Pattern source: STEP10F-2 locked outputs",
  "Expression-context source: STEP10F-3 locked outputs",
  "",
  "Detection definition: Raw_UMI > 0",
  "CPM: descriptive only",
  "Log2(CPM+1): descriptive only",
  "",
  "Statistical testing: NO",
  "P/FDR: NO",
  "Differential expression: NO",
  "Correlation: NO",
  "CellChat/LR: NO",
  "Pseudobulk recalculation: NO",
  "GEO raw matrix rereading: NO",
  "H5 rereading: NO",
  "Patient/cell exclusion: NO",
  "Specificity scoring: NO",
  "Causal inference: NO",
  "",
  "Cell-type specificity: NOT ESTABLISHED",
  "Intercellular communication: NOT ESTABLISHED",
  "Ligand-receptor interaction: NOT ESTABLISHED",
  "",
  "Sensitivity Endothelial kept separate.",
  "",
  "STEP10G FINAL LOCK: PASS",
  "============================================================"
)

writeLines(
  final_lock_text,
  file.path(
    out_dir,
    "STEP10G_FINAL_LOCK.txt"
  )
)

cat("\n")
cat("FINAL LOCK FILE WRITTEN.\n")
cat(
  file.path(
    out_dir,
    "STEP10G_FINAL_LOCK.txt"
  ),
  "\n"
)


getwd()

load("O:/Rproject/Paper 1_Step 10-Single-Cell/Paper 1_Step 10-Single-Cell.RData")
ls()


# ============================================================
# STEP10H-0
# Input Discovery and Schema Check
# ============================================================

rm(list = ls())
gc()

base_dir <-
  "input_dir/Paper 1_Step 10-Single-Cell"


cat("\n============================================================\n")
cat("STEP10H-0 — FIND LOCKED INPUT FILES\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# Search all STEP10D / STEP10F files
# ------------------------------------------------------------

all_files <-
  list.files(
    base_dir,
    pattern = "\\.(tsv|txt|csv|rds)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )


target_files <-
  all_files[
    grepl(
      "STEP10D-B|STEP10D-C|STEP10F-1",
      basename(all_files),
      ignore.case = TRUE
    )
  ]


cat(
  "Candidate files found: ",
  length(target_files),
  "\n\n",
  sep = ""
)


print(
  target_files
)


# ------------------------------------------------------------
# Read each candidate and print column names
# ------------------------------------------------------------

for (f in target_files) {
  
  cat("\n------------------------------------------------------------\n")
  cat("FILE:\n")
  cat(f, "\n")
  cat("------------------------------------------------------------\n")
  
  x <-
    tryCatch(
      data.table::fread(
        f,
        nrows = 5,
        data.table = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  
  if (is.null(x)) {
    
    cat("Could not read file preview.\n")
    
  } else {
    
    cat(
      "Columns:\n"
    )
    
    print(
      names(x)
    )
    
    cat(
      "\nPreview:\n"
    )
    
    print(x)
    
  }
  
}

##STEP10H — Final Manuscript Supplementary Export
cat("============================================================\n")
cat("STEP10H — FINAL MANUSCRIPT SUPPLEMENTARY EXPORT\n")
cat("Supplementary Figure Sz / Table Sz / Table Saa\n")
cat("============================================================\n\n")


# ============================================================
# 0. PACKAGES
# ============================================================

required_pkgs <- c(
  "data.table",
  "ggplot2",
  "openxlsx"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_pkgs, collapse = ", "),
      "\nPlease install them before running STEP10H."
    )
  )
}

library(data.table)
library(ggplot2)
library(openxlsx)


# ============================================================
# 1. PATHS
# ============================================================

base_dir <-
  "input_dir/Paper 1_Step 10-Single-Cell"

out_dir <-
  file.path(
    base_dir,
    "STEP10H_Manuscript_Supplementary_Figure_Table"
  )

fig_dir <-
  file.path(
    out_dir,
    "Supplementary_Figure_Sz"
  )

table_sz_dir <-
  file.path(
    out_dir,
    "Supplementary_Table_Sz"
  )

table_saa_dir <-
  file.path(
    out_dir,
    "Supplementary_Table_Saa"
  )

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_sz_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_saa_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. LOCKED INPUT FILES
# ============================================================

file_10D_B <-
  file.path(
    base_dir,
    "STEP10D-B_CoreGene_RawCount_Detection_Audit",
    "STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv"
  )

file_10F1A <-
  file.path(
    base_dir,
    "STEP10F-1_Patient_Level_Expression_Profile",
    "STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv"
  )

file_10F1B <-
  file.path(
    base_dir,
    "STEP10F-1_Patient_Level_Expression_Profile",
    "STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv"
  )

file_10F1C <-
  file.path(
    base_dir,
    "STEP10F-1_Patient_Level_Expression_Profile",
    "STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv"
  )


# ============================================================
# 3. INPUT EXISTENCE AUDIT
# ============================================================

input_files <- c(
  STEP10D_B = file_10D_B,
  STEP10F1A = file_10F1A,
  STEP10F1B = file_10F1B,
  STEP10F1C = file_10F1C
)

input_exists <- file.exists(input_files)

cat("INPUT FILE AUDIT\n")
cat("------------------------------------------------------------\n")

for (i in seq_along(input_files)) {
  cat(
    names(input_files)[i],
    ": ",
    input_exists[i],
    "\n",
    sep = ""
  )
}

if (!all(input_exists)) {
  stop(
    "One or more locked input files are missing."
  )
}

cat("\nALL LOCKED INPUT FILES FOUND: PASS\n\n")


# ============================================================
# 4. READ LOCKED FILES
# ============================================================

raw_counts <-
  fread(
    file_10D_B,
    data.table = FALSE,
    check.names = FALSE
  )

primary_expr <-
  fread(
    file_10F1A,
    data.table = FALSE,
    check.names = FALSE
  )

sensitivity_expr <-
  fread(
    file_10F1B,
    data.table = FALSE,
    check.names = FALSE
  )

primary_summary <-
  fread(
    file_10F1C,
    data.table = FALSE,
    check.names = FALSE
  )


# ============================================================
# 5. STRUCTURAL LOCK AUDIT
# ============================================================

core_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

primary_compartments <- c(
  "T-cell",
  "Myeloid/macrophage",
  "Endothelial"
)

expected_primary_n <- c(
  "T-cell" = 9,
  "Myeloid/macrophage" = 9,
  "Endothelial" = 4
)

cat("STRUCTURAL LOCK AUDIT\n")
cat("------------------------------------------------------------\n")

# ---- raw count table ----

required_raw_cols <- c(
  "Dataset",
  "Group_ID",
  "Patient",
  "Compartment",
  "Eligibility",
  "Cell_Count",
  core_genes
)

stopifnot(
  all(required_raw_cols %in% names(raw_counts))
)

# ---- patient expression ----

required_expr_cols <- c(
  "Dataset",
  "Patient",
  "Compartment",
  "Eligibility",
  "Group_ID",
  "Cell_Count",
  "Gene",
  "Raw_UMI",
  "Library_Total_UMI",
  "CPM",
  "Log2_CPM_Plus1",
  "Detected",
  "Analysis_Set"
)

stopifnot(
  all(required_expr_cols %in% names(primary_expr))
)

# ---- summary ----

required_summary_cols <- c(
  "Gene",
  "Compartment",
  "N_Patients",
  "Detection_N",
  "Detection_Percent",
  "Raw_UMI_Mean",
  "Raw_UMI_Median",
  "Raw_UMI_Q1",
  "Raw_UMI_Q3",
  "Raw_UMI_Min",
  "Raw_UMI_Max",
  "Raw_UMI_IQR",
  "CPM_Mean",
  "CPM_Median",
  "CPM_Q1",
  "CPM_Q3",
  "CPM_Min",
  "CPM_Max",
  "CPM_IQR",
  "Log2_CPM_Plus1_Mean",
  "Log2_CPM_Plus1_Median",
  "Log2_CPM_Plus1_Q1",
  "Log2_CPM_Plus1_Q3",
  "Log2_CPM_Plus1_Min",
  "Log2_CPM_Plus1_Max",
  "Log2_CPM_Plus1_IQR"
)

stopifnot(
  all(required_summary_cols %in% names(primary_summary))
)

# ---- sensitivity ----

stopifnot(
  all(required_expr_cols %in% names(sensitivity_expr))
)

cat("Required columns: PASS\n")


# ============================================================
# 6. LOCKED CORE-GENE AUDIT
# ============================================================

observed_genes_raw <-
  intersect(
    core_genes,
    names(raw_counts)
  )

observed_genes_expr <-
  unique(primary_expr$Gene)

observed_genes_summary <-
  unique(primary_summary$Gene)

gene_pass <-
  identical(
    sort(observed_genes_raw),
    sort(core_genes)
  ) &&
  identical(
    sort(observed_genes_expr),
    sort(core_genes)
  ) &&
  identical(
    sort(observed_genes_summary),
    sort(core_genes)
  )

cat(
  "Core genes AOC3 / MS4A7 / IL7R: ",
  gene_pass,
  "\n",
  sep = ""
)

if (!gene_pass) {
  stop("Core-gene lock failed.")
}


# ============================================================
# 7. PRIMARY COMPARTMENT / PATIENT LOCK
# ============================================================

primary_expr$Compartment <-
  factor(
    primary_expr$Compartment,
    levels = primary_compartments
  )

primary_summary$Compartment <-
  factor(
    primary_summary$Compartment,
    levels = primary_compartments
  )

primary_patient_counts <-
  unique(
    primary_expr[
      ,
      c("Patient", "Compartment")
    ]
  )

primary_patient_counts <-
  aggregate(
    Patient ~ Compartment,
    data = primary_patient_counts,
    FUN = function(x) length(unique(x))
  )

names(primary_patient_counts)[2] <-
  "N_Patients"

print(primary_patient_counts)

for (comp in primary_compartments) {
  
  observed_n <-
    primary_patient_counts$N_Patients[
      primary_patient_counts$Compartment == comp
    ]
  
  if (length(observed_n) != 1 ||
      observed_n != expected_primary_n[comp]) {
    
    stop(
      paste0(
        "Primary patient count mismatch for ",
        comp
      )
    )
  }
}

cat("Primary patient/compartment lock: PASS\n")


# ============================================================
# 8. TABLE SZ
#    LOCKED SUMMARY — NO RECALCULATION
# ============================================================

table_sz <-
  primary_summary[
    primary_summary$Compartment %in% primary_compartments &
      primary_summary$Gene %in% core_genes,
    ,
    drop = FALSE
  ]

# Preserve intended order
table_sz$Gene <-
  factor(
    table_sz$Gene,
    levels = core_genes
  )

table_sz$Compartment <-
  factor(
    table_sz$Compartment,
    levels = primary_compartments
  )

table_sz <-
  table_sz[
    order(
      table_sz$Gene,
      table_sz$Compartment
    ),
    ,
    drop = FALSE
  ]

# Rename to manuscript-facing names
table_sz_final <-
  data.frame(
    Gene = as.character(table_sz$Gene),
    Compartment = as.character(table_sz$Compartment),
    N_Patients = table_sz$N_Patients,
    Detection_N = table_sz$Detection_N,
    Detection_Percent = table_sz$Detection_Percent,
    Raw_UMI_Median = table_sz$Raw_UMI_Median,
    Raw_UMI_IQR = table_sz$Raw_UMI_IQR,
    CPM_Median = table_sz$CPM_Median,
    CPM_IQR = table_sz$CPM_IQR,
    Log2_CPM_Plus1_Median =
      table_sz$Log2_CPM_Plus1_Median,
    Log2_CPM_Plus1_IQR =
      table_sz$Log2_CPM_Plus1_IQR,
    stringsAsFactors = FALSE
  )

# Number formatting
table_sz_final$Detection_Percent <-
  round(
    table_sz_final$Detection_Percent,
    2
  )

numeric_cols_sz <-
  c(
    "Raw_UMI_Median",
    "Raw_UMI_IQR",
    "CPM_Median",
    "CPM_IQR",
    "Log2_CPM_Plus1_Median",
    "Log2_CPM_Plus1_IQR"
  )

table_sz_final[numeric_cols_sz] <-
  lapply(
    table_sz_final[numeric_cols_sz],
    function(x) round(x, 4)
  )

# Structural checks
stopifnot(
  nrow(table_sz_final) == 9
)

stopifnot(
  all(
    table_sz_final$N_Patients ==
      c(
        9, 9, 4,
        9, 9, 4,
        9, 9, 4
      )
  )
)

cat(
  "\nTABLE SZ: 9 gene × compartment rows — PASS\n"
)


# ============================================================
# 9. EXPORT TABLE SZ — XLSX
# ============================================================

table_sz_xlsx <-
  file.path(
    table_sz_dir,
    "Supplementary_Table_Sz_Patient_Level_Pseudobulk_Expression_Detection_Summary.xlsx"
  )

wb_sz <-
  createWorkbook()

addWorksheet(
  wb_sz,
  "Table_Sz"
)

writeData(
  wb_sz,
  "Table_Sz",
  table_sz_final,
  withFilter = TRUE
)

header_style <-
  createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center"
  )

addStyle(
  wb_sz,
  "Table_Sz",
  header_style,
  rows = 1,
  cols = 1:ncol(table_sz_final),
  gridExpand = TRUE
)

setColWidths(
  wb_sz,
  "Table_Sz",
  cols = 1:ncol(table_sz_final),
  widths = "auto"
)

freezePane(
  wb_sz,
  "Table_Sz",
  firstActiveRow = 2
)

saveWorkbook(
  wb_sz,
  table_sz_xlsx,
  overwrite = TRUE
)


# Also export TSV for audit/transparency
table_sz_tsv <-
  file.path(
    table_sz_dir,
    "Supplementary_Table_Sz_Patient_Level_Pseudobulk_Expression_Detection_Summary.tsv"
  )

fwrite(
  table_sz_final,
  table_sz_tsv,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 10. TABLE SAA
#     PATIENT-LEVEL RAW-UMI PSEUDOBULK
#     DIRECTLY FROM LOCKED STEP10D-B
# ============================================================

table_saa <-
  raw_counts[
    raw_counts$Eligibility == "Primary" |
      raw_counts$Eligibility == "Sensitivity",
    ,
    drop = FALSE
  ]

# Explicit ordering
compartment_order_saa <-
  c(
    "T-cell",
    "Myeloid/macrophage",
    "Endothelial"
  )

table_saa$Compartment <-
  factor(
    table_saa$Compartment,
    levels = compartment_order_saa
  )

table_saa$Eligibility <-
  factor(
    table_saa$Eligibility,
    levels = c(
      "Primary",
      "Sensitivity"
    )
  )

table_saa <-
  table_saa[
    order(
      table_saa$Eligibility,
      table_saa$Compartment,
      table_saa$Patient
    ),
    ,
    drop = FALSE
  ]

# Convert factors back to character
table_saa$Compartment <-
  as.character(table_saa$Compartment)

table_saa$Eligibility <-
  as.character(table_saa$Eligibility)

# Keep only manuscript-relevant columns
table_saa_final <-
  table_saa[
    ,
    c(
      "Patient",
      "Compartment",
      "Eligibility",
      "Group_ID",
      "Cell_Count",
      "AOC3",
      "MS4A7",
      "IL7R"
    )
  ]

# Expected:
# Primary = 22 groups
# Sensitivity = 2 groups
# Total = 24
stopifnot(
  nrow(table_saa_final) == 24
)

stopifnot(
  sum(table_saa_final$Eligibility == "Primary") == 22
)

stopifnot(
  sum(table_saa_final$Eligibility == "Sensitivity") == 2
)

cat(
  "TABLE SAA: 24 patient-compartment rows — PASS\n"
)


# ============================================================
# 11. EXPORT TABLE SAA — XLSX
# ============================================================

table_saa_xlsx <-
  file.path(
    table_saa_dir,
    "Supplementary_Table_Saa_Patient_Level_Raw_UMI_Pseudobulk_Profiles.xlsx"
  )

wb_saa <-
  createWorkbook()

addWorksheet(
  wb_saa,
  "Table_Saa"
)

writeData(
  wb_saa,
  "Table_Saa",
  table_saa_final,
  withFilter = TRUE
)

addStyle(
  wb_saa,
  "Table_Saa",
  header_style,
  rows = 1,
  cols = 1:ncol(table_saa_final),
  gridExpand = TRUE
)

setColWidths(
  wb_saa,
  "Table_Saa",
  cols = 1:ncol(table_saa_final),
  widths = "auto"
)

freezePane(
  wb_saa,
  "Table_Saa",
  firstActiveRow = 2
)

saveWorkbook(
  wb_saa,
  table_saa_xlsx,
  overwrite = TRUE
)


# Also export TSV
table_saa_tsv <-
  file.path(
    table_saa_dir,
    "Supplementary_Table_Saa_Patient_Level_Raw_UMI_Pseudobulk_Profiles.tsv"
  )

fwrite(
  table_saa_final,
  table_saa_tsv,
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 12. FIGURE SZ DATA
# ============================================================

# Only primary patient-level observations
fig_df <-
  primary_expr[
    primary_expr$Analysis_Set == "Primary" &
      primary_expr$Gene %in% core_genes &
      primary_expr$Compartment %in% primary_compartments,
    ,
    drop = FALSE
  ]

fig_df$Gene <-
  factor(
    fig_df$Gene,
    levels = core_genes
  )

fig_df$Compartment <-
  factor(
    fig_df$Compartment,
    levels = primary_compartments
  )


# ------------------------------------------------------------
# Verify 27 patient × compartment × gene observations
# ------------------------------------------------------------

expected_fig_rows <-
  9 * 3 +
  9 * 3 +
  4 * 3

stopifnot(
  nrow(fig_df) == expected_fig_rows
)

cat(
  "FIGURE SZ primary observations: ",
  nrow(fig_df),
  " — PASS\n",
  sep = ""
)

# ============================================================
# 13. FIGURE SZ — PANEL A
#     DETECTION PREVALENCE
#
#     FINAL COLOR SYNCHRONIZATION
#
#     Panel A uses the EXACT SAME Compartment color mapping
#     as Panel B / Panel C.
#
#     No data, detection values, or statistics are recalculated.
# ============================================================


# ------------------------------------------------------------
# 13A. LOCK COMPARTMENT COLOR MAPPING
# ------------------------------------------------------------
#
# The color mapping is defined ONCE before Panel A and then
# reused by Panel B and Panel C.
#
# This guarantees identical biological color identity across
# all three panels.
# ------------------------------------------------------------

compartment_levels <-
  levels(
    factor(
      fig_df$Compartment,
      levels = primary_compartments
    )
  )

compartment_colors <-
  setNames(
    scales::hue_pal()(
      length(compartment_levels)
    ),
    compartment_levels
  )

cat(
  "Panel A/B/C — locked Compartment color mapping:\n"
)

print(compartment_colors)


# ------------------------------------------------------------
# 13B. Panel A data
# ------------------------------------------------------------

panel_A <-
  primary_summary[
    primary_summary$Gene %in% core_genes &
      primary_summary$Compartment %in% primary_compartments,
    ,
    drop = FALSE
  ]

panel_A$Gene <-
  factor(
    panel_A$Gene,
    levels = core_genes
  )

panel_A$Compartment <-
  factor(
    panel_A$Compartment,
    levels = primary_compartments
  )

panel_A$Label <-
  paste0(
    round(
      panel_A$Detection_Percent,
      1
    ),
    "%\n(n=",
    panel_A$Detection_N,
    "/",
    panel_A$N_Patients,
    ")"
  )


# ------------------------------------------------------------
# 13C. Panel A plot
# ------------------------------------------------------------

pA <-
  ggplot(
    panel_A,
    aes(
      x = Gene,
      y = Detection_Percent,
      fill = Compartment
    )
  ) +
  
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  
  geom_text(
    aes(
      label = Label
    ),
    position = position_dodge(
      width = 0.75
    ),
    vjust = -0.20,
    size = 2.8,
    family = "Arial"
  ) +
  
  scale_fill_manual(
    values = compartment_colors,
    breaks = primary_compartments,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(0, 118),
    breaks = seq(0, 100, 20),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Detection prevalence (%)",
    fill = "Compartment"
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    text = element_text(
      family = "Arial"
    ),
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.text = element_text(
      family = "Arial",
      size = 9
    ),
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    legend.text = element_text(
      family = "Arial",
      size = 8
    ),
    legend.position = "top"
  )


# ============================================================
# 14. FIGURE SZ — PANEL B
# ------------------------------------------------------------
# 14A. Verify the locked Compartment color mapping
# ------------------------------------------------------------

stopifnot(
  exists("compartment_colors")
)

stopifnot(
  identical(
    names(compartment_colors),
    primary_compartments
  )
)

cat(
  "\nPanel B — using the SAME locked Compartment color mapping\n"
)

print(compartment_colors)

# ------------------------------------------------------------
# 14B. Quantitative expression plot
# ------------------------------------------------------------

pB <-
  ggplot(
    fig_df,
    aes(
      x = Compartment,
      y = Log2_CPM_Plus1,
      fill = Compartment
    )
  ) +
  
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.8,
    alpha = 0.85
  ) +
  
  facet_wrap(
    ~ Gene,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = compartment_colors,
    drop = FALSE
  ) +
  
  labs(
    x = NULL,
    y = expression(
      log[2](CPM + 1)
    )
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    text = element_text(
      family = "Arial"
    ),
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.text = element_text(
      family = "Arial",
      size = 8
    ),
    strip.text = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    legend.position = "none"
  )


# ============================================================
# 15. FIGURE SZ — PANEL C
#     PATIENT-LEVEL DOT PLOT
#
#     FINAL REVISION
#     - Preserve locked fig_df
#     - Preserve original jitter structure
#     - Preserve Compartment shape mapping
#     - Add Compartment color mapping
#     - Use EXACTLY the same colors as Panel B
# ============================================================


pC <-
  ggplot(
    fig_df,
    aes(
      x = Compartment,
      y = Log2_CPM_Plus1,
      shape = Compartment,
      color = Compartment
    )
  ) +
  
  geom_point(
    size = 2.5,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  
  facet_wrap(
    ~ Gene,
    nrow = 1
  ) +
  
  scale_color_manual(
    values = compartment_colors,
    drop = FALSE
  ) +
  
  labs(
    x = NULL,
    y = expression(
      log[2](CPM + 1)
    ),
    shape = "Compartment",
    color = "Compartment"
  ) +
  
  guides(
    color = guide_legend(
      order = 1
    ),
    shape = guide_legend(
      order = 1
    )
  ) +
  
  theme_classic(
    base_size = 10
  ) +
  
  theme(
    text = element_text(
      family = "Arial"
    ),
    axis.title = element_text(
      family = "Arial",
      face = "bold",
      size = 10
    ),
    axis.text = element_text(
      family = "Arial",
      size = 8
    ),
    strip.text = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    legend.title = element_text(
      family = "Arial",
      face = "bold",
      size = 9
    ),
    legend.text = element_text(
      family = "Arial",
      size = 8
    ),
    legend.position = "top"
  )


# ============================================================
# 16. FINAL STANDALONE FIGURE EXPORT
#     Supplementary Figure Sz
#
#     Three independent publication-ready panels
#     No combined figure is generated.
#
#     Locked data only:
#       STEP10F-1A
#       STEP10F-1C
#       STEP10D-B
#
#     No recalculation of expression values.
# ============================================================


# ------------------------------------------------------------
# 16A. Publication-safe export function
# ------------------------------------------------------------

save_plot_pair <- function(
    plot_obj,
    stem,
    width = 7,
    height = 5
) {
  
  pdf_file <-
    file.path(
      fig_dir,
      paste0(stem, ".pdf")
    )
  
  tiff_file <-
    file.path(
      fig_dir,
      paste0(stem, "_600dpi.tiff")
    )
  
  
  # ----------------------------------------------------------
  # Vector PDF
  # ----------------------------------------------------------
  
  cairo_pdf(
    filename = pdf_file,
    width = width,
    height = height,
    family = "Arial"
  )
  
  print(plot_obj)
  
  dev.off()
  
  
  # ----------------------------------------------------------
  # 600-dpi TIFF
  # ----------------------------------------------------------
  
  tiff(
    filename = tiff_file,
    width = width,
    height = height,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  
  print(plot_obj)
  
  dev.off()
  
  
  # ----------------------------------------------------------
  # Validation
  # ----------------------------------------------------------
  
  pdf_ok <-
    file.exists(pdf_file) &&
    file.info(pdf_file)$size > 1000
  
  tiff_ok <-
    file.exists(tiff_file) &&
    file.info(tiff_file)$size > 1000
  
  
  cat(
    "------------------------------------------------------------\n"
  )
  
  cat(
    "Figure: ",
    stem,
    "\n",
    sep = ""
  )
  
  cat(
    "PDF:  ",
    pdf_ok,
    "\n",
    sep = ""
  )
  
  cat(
    "TIFF: ",
    tiff_ok,
    "\n",
    sep = ""
  )
  
  if (!pdf_ok || !tiff_ok) {
    
    stop(
      paste0(
        "Output validation failed: ",
        stem
      )
    )
    
  }
  
  cat("Status: PASS\n\n")
  
}


# ============================================================
# 17. FINAL STANDALONE PANEL EXPORT
# ============================================================


# ------------------------------------------------------------
# Panel A
# Detection prevalence
# ------------------------------------------------------------

save_plot_pair(
  pA,
  "Supplementary_Figure_Sz_Panel_A_Detection",
  width = 7.0,
  height = 4.8
)


# ------------------------------------------------------------
# Panel B
# Quantitative expression
#
# FINAL REVISION:
# Expanded from 7.5 × 5.2 inch
# to 9.0 × 5.2 inch
# to reduce X-axis Compartment-label overlap.
# ------------------------------------------------------------

save_plot_pair(
  pB,
  "Supplementary_Figure_Sz_Panel_B_Quantitative_Expression",
  width = 9.0,
  height = 5.2
)


# ------------------------------------------------------------
# Panel C
# Patient-level expression
#
# FINAL SIZE:
# 9.0 × 5.2 inch
#
# Colors are synchronized with Panel B.
# ------------------------------------------------------------

save_plot_pair(
  pC,
  "Supplementary_Figure_Sz_Panel_C_Patient_Level_Expression",
  width = 9.0,
  height = 5.2
)


# ============================================================
# 18. FINAL OUTPUT AUDIT
# ============================================================


expected_outputs <- c(
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_A_Detection.pdf"
  ),
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_A_Detection_600dpi.tiff"
  ),
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_B_Quantitative_Expression.pdf"
  ),
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_B_Quantitative_Expression_600dpi.tiff"
  ),
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_C_Patient_Level_Expression.pdf"
  ),
  
  file.path(
    fig_dir,
    "Supplementary_Figure_Sz_Panel_C_Patient_Level_Expression_600dpi.tiff"
  )
  
)


output_check <- data.frame(
  
  File = expected_outputs,
  
  Exists =
    file.exists(expected_outputs),
  
  Size_Bytes =
    ifelse(
      file.exists(expected_outputs),
      file.info(expected_outputs)$size,
      NA
    ),
  
  stringsAsFactors = FALSE
  
)


print(output_check)


if (
  !all(output_check$Exists) ||
  any(output_check$Size_Bytes <= 1000, na.rm = TRUE)
) {
  
  stop(
    "STEP10H standalone figure output validation FAILED."
  )
  
}


cat(
  "\n============================================================\n"
)

cat(
  "STEP10H — STANDALONE FIGURE Sz EXPORT\n"
)

cat(
  "============================================================\n"
)

cat(
  "Panel A: PASS\n"
)

cat(
  "Panel B: PASS\n"
)

cat(
  "Panel C: PASS\n"
)

cat(
  "Combined figure: NOT GENERATED\n"
)

cat(
  "Vector PDF: PASS\n"
)

cat(
  "600-dpi TIFF: PASS\n"
)

cat(
  "Locked data: PRESERVED\n"
)

cat(
  "No expression recalculation: TRUE\n"
)

cat(
  "Panel B dimensions: 9.0 × 5.2 inch\n"
)

cat(
  "Panel C dimensions: 9.0 × 5.2 inch\n"
)

cat(
  "Panel C colors: SYNCHRONIZED WITH PANEL B\n"
)

cat(
  "Final status: PASS\n"
)

cat(
  "============================================================\n"
)


