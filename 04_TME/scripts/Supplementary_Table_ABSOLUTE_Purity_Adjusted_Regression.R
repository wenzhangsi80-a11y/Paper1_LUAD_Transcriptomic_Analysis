
# ============================================================
# Script:
#   Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.R
#
# Analysis module:
#   Tumor Microenvironment (TME) Analysis
#   STEP 7B6_05: Final Supplementary Table Extraction
#
# Study:
#   Integrated transcriptomic and single-cell characterization
#   of AOC3, MS4A7, and IL7R in LUAD
#
# 1. Research objective
#   To extract and verify archived ABSOLUTE purity-adjusted
#   linear regression results for the three prioritized study
#   genes (AOC3, MS4A7, and IL7R) in relation to three
#   ESTIMATE-derived tumor microenvironment scores.
#
#   This workflow generates Supplementary Table 3a and
#   supporting audit files for the ABSOLUTE purity-adjusted
#   gene-TME sensitivity analysis.
#
# 2. Data sources
#   - TCGA-LUAD primary tumor samples.
#   - ESTIMATE-derived ImmuneScore, StromalScore, and
#     ESTIMATEScore.
#   - ABSOLUTE-derived tumor purity estimates.
#   - Archived STEP7B5 and STEP7B6 analysis objects.
#
#   Original data processing and regression model fitting
#   were performed in the preceding STEP7B5 and STEP7B6
#   workflows.
#
# 3. Input files
#   - data/input_data/
#     Paper1_STEP7B5_locked_purity_adjusted_dataset.RData
#
#   - data/input_data/
#     Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData
#
#   Required archived objects:
#   STEP7B5:
#     purity_adjusted_dataset
#
#   STEP7B6:
#     purity_adjusted_models
#     purity_adjusted_results
#     purity_adjusted_results_raw
#     model_spec_7b6
#
# 4. Sample inclusion criteria
#   - Primary LUAD samples with complete analysis data.
#   - Locked STEP7B5 analysis dataset.
#   - Exactly 500 samples are expected.
#   - Unique and non-missing Sample_ID values.
#   - Complete measurements for all three prioritized
#     genes, three ESTIMATE scores, and ABSOLUTE tumor purity.
#   - All nine archived regression models must contain
#     500 observations.
#
# 5. Analysis methods
#   - Load archived RData objects into isolated environments.
#   - Verify the locked analysis dataset and model specification.
#   - Extract gene and ABSOLUTE_Purity coefficients directly
#     from the archived linear regression (lm) objects.
#   - Evaluate nine gene-TME combinations.
#   - Apply Benjamini-Hochberg (BH) correction to the nine
#     gene-term P values.
#   - Compare recalculated BH-adjusted P values with the
#     archived STEP7B6 result object.
#   - Perform dataset, model, coefficient, and final audits.
#   - Export publication-ready CSV and optional Excel files.
#
#   No regression models are refitted in this workflow.
#
# 6. Software environment
#   R version: Recorded at execution in the final audit file.
#   Key packages:
#     - here: Project-relative path management.
#     - stats: Model summaries, nobs(), and p.adjust().
#     - utils: RData loading and CSV export.
#     - openxlsx: Optional Excel workbook export.
#
#   Exact package versions should be recorded in the
#   Environment/README.md after environment verification.
#
# 7. Random seed
#   No random seed is required.
#   This workflow performs deterministic extraction,
#   statistical adjustment, and auditing of archived results.
#
# 8. Output files
#   Output directory:
#     04_TME/results/STEP7B6_05_Final_Supplementary_Table/
#
#   Main outputs:
#     Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.csv
#     Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.xlsx
#
#   Supporting audit files:
#     STEP7B6_05_Model_Coefficient_Details.csv
#     STEP7B6_05_BH_Adjustment_Audit.csv
#     STEP7B6_05_Dataset_Variable_Audit.csv
#     STEP7B6_05_Model_Audit.csv
#     STEP7B6_05_Final_Audit.txt
#
# 9. Key analytical results
#   - Three prioritized study genes:
#     AOC3, MS4A7, and IL7R.
#   - Three ESTIMATE-derived TME scores:
#     ImmuneScore, StromalScore, and ESTIMATEScore.
#   - Nine gene-TME regression models.
#   - Expected sample size: N = 500 per model.
#   - Gene-term P values are adjusted using BH correction.
#
#   Actual regression coefficients, P values, and final
#   audit status are obtained from the generated outputs.
#   No unverified numerical results are reported here.
#
# 10. Reproducibility notes
#   - Project-relative paths are managed using here::here().
#   - Original STEP7B5 and STEP7B6 RData files are retained.
#   - Archived lm objects are used without model refitting.
#   - No sample rematching or barcode remapping is performed.
#   - No ABSOLUTE purity imputation is performed.
#   - The locked analysis dataset must contain exactly
#     500 complete observations.
#   - All nine archived models must pass the predefined
#     model and coefficient checks.
#   - Recalculated BH-adjusted P values must match the
#     archived STEP7B6 results.
#   - The final audit reports PASS only when all
#     predefined conditions are satisfied.
#
# ============================================================


cat("\n============================================================\n")
cat("STEP 7B6_05: Final Supplementary Table Extraction\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

# Identify the GitHub project root.
project_dir <- here::here()

# Define input and output directories.
input_dir <- file.path(
  project_dir,
  "04_TME",
  "data",
  "input_data"
)

output_dir <- file.path(
  project_dir,
  "04_TME",
  "results",
  "STEP7B6_05_Final_Supplementary_Table"
)

# Create the output directory if it does not exist.
dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Input archived RData files.
step7b5_file <- file.path(
  input_dir,
  "Paper1_STEP7B5_locked_purity_adjusted_dataset.RData"
)

step7b6_file <- file.path(
  input_dir,
  "Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData"
)

if (!file.exists(step7b5_file)) {
  stop("STEP7B5 RData file not found: ", step7b5_file)
}

if (!file.exists(step7b6_file)) {
  stop("STEP7B6 RData file not found: ", step7b6_file)
}


# ------------------------------------------------------------
# 1. Load archived RData into isolated environments
# ------------------------------------------------------------

env_7b5 <- new.env(parent = baseenv())
env_7b6 <- new.env(parent = baseenv())

load(step7b5_file, envir = env_7b5)
load(step7b6_file, envir = env_7b6)

required_7b5 <- c(
  "purity_adjusted_dataset"
)

required_7b6 <- c(
  "purity_adjusted_models",
  "purity_adjusted_results",
  "purity_adjusted_results_raw",
  "model_spec_7b6"
)

missing_7b5 <- setdiff(
  required_7b5,
  ls(env_7b5)
)

missing_7b6 <- setdiff(
  required_7b6,
  ls(env_7b6)
)

if (length(missing_7b5) > 0) {
  stop(
    "Missing STEP7B5 objects: ",
    paste(missing_7b5, collapse = ", ")
  )
}

if (length(missing_7b6) > 0) {
  stop(
    "Missing STEP7B6 objects: ",
    paste(missing_7b6, collapse = ", ")
  )
}

purity_adjusted_dataset <-
  env_7b5$purity_adjusted_dataset

purity_adjusted_models <-
  env_7b6$purity_adjusted_models

purity_adjusted_results <-
  env_7b6$purity_adjusted_results

purity_adjusted_results_raw <-
  env_7b6$purity_adjusted_results_raw

model_spec_7b6 <-
  env_7b6$model_spec_7b6


# ------------------------------------------------------------
# 2. Verify the locked STEP7B5 dataset
# ------------------------------------------------------------

expected_genes <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

expected_traits <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)

required_dataset_columns <- c(
  "Sample_ID",
  expected_genes,
  expected_traits,
  "ABSOLUTE_Purity"
)

missing_dataset_columns <- setdiff(
  required_dataset_columns,
  colnames(purity_adjusted_dataset)
)

if (length(missing_dataset_columns) > 0) {
  stop(
    "Missing dataset columns: ",
    paste(missing_dataset_columns, collapse = ", ")
  )
}

dataset_N <- nrow(purity_adjusted_dataset)

if (dataset_N != 500) {
  stop("Unexpected STEP7B5 dataset N: ", dataset_N)
}

if (anyDuplicated(purity_adjusted_dataset$Sample_ID) > 0) {
  stop("Duplicate Sample_ID detected.")
}

if (anyNA(purity_adjusted_dataset$Sample_ID)) {
  stop("Missing Sample_ID detected.")
}

variable_audit <- data.frame(
  Variable = required_dataset_columns,
  N = vapply(
    purity_adjusted_dataset[required_dataset_columns],
    length,
    integer(1)
  ),
  Missing = vapply(
    purity_adjusted_dataset[required_dataset_columns],
    function(x) sum(is.na(x)),
    integer(1)
  ),
  Finite = vapply(
    purity_adjusted_dataset[required_dataset_columns],
    function(x) {
      if (is.numeric(x)) {
        sum(is.finite(x))
      } else {
        sum(!is.na(x))
      }
    },
    integer(1)
  ),
  stringsAsFactors = FALSE
)

if (any(variable_audit$Missing != 0)) {
  stop("Missing values detected in the locked analysis dataset.")
}

if (any(
  vapply(
    purity_adjusted_dataset[
      c(expected_genes, expected_traits, "ABSOLUTE_Purity")
    ],
    function(x) !is.numeric(x),
    logical(1)
  )
)) {
  stop("One or more analysis variables are not numeric.")
}


# ------------------------------------------------------------
# 3. Verify model specification
# ------------------------------------------------------------

if (!is.data.frame(model_spec_7b6)) {
  stop("model_spec_7b6 is not a data.frame.")
}

if (!all(c("TME_Trait", "Gene") %in%
         colnames(model_spec_7b6))) {
  stop("Unexpected model_spec_7b6 column structure.")
}

expected_combinations <- expand.grid(
  TME_Trait = expected_traits,
  Gene = expected_genes,
  stringsAsFactors = FALSE
)

spec_keys <- paste(
  model_spec_7b6$TME_Trait,
  model_spec_7b6$Gene,
  sep = "__"
)

expected_keys <- paste(
  expected_combinations$TME_Trait,
  expected_combinations$Gene,
  sep = "__"
)

if (nrow(model_spec_7b6) != 9 ||
    anyDuplicated(spec_keys) > 0 ||
    !setequal(spec_keys, expected_keys)) {
  stop("Model specification does not match the expected 3 x 3 design.")
}

if (!is.list(purity_adjusted_models) ||
    length(purity_adjusted_models) != 9) {
  stop("Expected exactly 9 archived regression models.")
}


# ------------------------------------------------------------
# 4. Extract model coefficients directly from archived lm objects
# ------------------------------------------------------------

extract_model_terms <- function(model, model_name) {
  
  if (!inherits(model, "lm")) {
    stop("Model is not an lm object: ", model_name)
  }
  
  model_formula <- formula(model)
  
  response_variable <- all.vars(
    model_formula[[2]]
  )
  
  predictor_variables <- all.vars(
    model_formula[[3]]
  )
  
  if (length(response_variable) != 1) {
    stop("Cannot identify response variable: ", model_name)
  }
  
  gene_terms <- intersect(
    predictor_variables,
    expected_genes
  )
  
  purity_terms <- intersect(
    predictor_variables,
    "ABSOLUTE_Purity"
  )
  
  if (length(gene_terms) != 1) {
    stop("Cannot uniquely identify gene term: ", model_name)
  }
  
  if (length(purity_terms) != 1) {
    stop("ABSOLUTE_Purity missing from model: ", model_name)
  }
  
  trait <- response_variable
  
  gene <- gene_terms
  
  if (!trait %in% expected_traits) {
    stop("Unexpected TME outcome: ", trait)
  }
  
  model_N <- stats::nobs(model)
  
  if (model_N != 500) {
    stop(
      "Unexpected model N for ",
      model_name,
      ": ",
      model_N
    )
  }
  
  coefficient_matrix <- summary(model)$coefficients
  
  required_terms <- c(
    gene,
    "ABSOLUTE_Purity"
  )
  
  if (!all(required_terms %in%
           rownames(coefficient_matrix))) {
    stop("Required coefficient row missing: ", model_name)
  }
  
  extract_one_term <- function(term_name, term_type) {
    
    x <- coefficient_matrix[term_name, ]
    
    data.frame(
      Model = model_name,
      TME_Trait = trait,
      Gene = gene,
      Term = term_name,
      Term_Type = term_type,
      N = model_N,
      Beta = unname(x["Estimate"]),
      SE = unname(x["Std. Error"]),
      t_value = unname(x["t value"]),
      P_value = unname(x["Pr(>|t|)"]),
      stringsAsFactors = FALSE
    )
  }
  
  rbind(
    extract_one_term(gene, "Gene"),
    extract_one_term("ABSOLUTE_Purity", "Purity_covariate")
  )
}


model_names <- names(purity_adjusted_models)

if (is.null(model_names) ||
    any(model_names == "") ||
    anyDuplicated(model_names) > 0) {
  stop("Model list must have unique, non-empty names.")
}

coefficient_details <- do.call(
  rbind,
  lapply(
    seq_along(purity_adjusted_models),
    function(i) {
      extract_model_terms(
        purity_adjusted_models[[i]],
        model_names[i]
      )
    }
  )
)

rownames(coefficient_details) <- NULL


# ------------------------------------------------------------
# 5. Separate gene terms and purity covariate terms
# ------------------------------------------------------------

gene_results <- subset(
  coefficient_details,
  Term_Type == "Gene"
)

purity_results <- subset(
  coefficient_details,
  Term_Type == "Purity_covariate"
)

if (nrow(gene_results) != 9) {
  stop("Expected 9 gene coefficient rows.")
}

if (nrow(purity_results) != 9) {
  stop("Expected 9 ABSOLUTE_Purity coefficient rows.")
}

gene_results$Key <- paste(
  gene_results$TME_Trait,
  gene_results$Gene,
  sep = "__"
)

if (anyDuplicated(gene_results$Key) > 0 ||
    !setequal(gene_results$Key, expected_keys)) {
  stop("Extracted gene-term combinations are incomplete or duplicated.")
}


# ------------------------------------------------------------
# 6. BH adjustment across the 9 gene-term P values
# ------------------------------------------------------------

gene_results$P_BH <- p.adjust(
  gene_results$P_value,
  method = "BH"
)

gene_results$Significant_FDR05 <-
  gene_results$P_BH < 0.05

gene_results$Direction <- ifelse(
  gene_results$Beta > 0,
  "Positive",
  ifelse(
    gene_results$Beta < 0,
    "Negative",
    "Zero"
  )
)


# ------------------------------------------------------------
# 7. Create publication-ready table
# ------------------------------------------------------------

trait_order <- c(
  "ImmuneScore",
  "StromalScore",
  "ESTIMATEScore"
)

gene_order <- c(
  "AOC3",
  "MS4A7",
  "IL7R"
)

gene_results$TME_Trait <- factor(
  gene_results$TME_Trait,
  levels = trait_order
)

gene_results$Gene <- factor(
  gene_results$Gene,
  levels = gene_order
)

gene_results <- gene_results[
  order(
    gene_results$TME_Trait,
    gene_results$Gene
  ),
]

publication_table <- data.frame(
  TME_score = as.character(gene_results$TME_Trait),
  Gene = as.character(gene_results$Gene),
  N = gene_results$N,
  Beta = gene_results$Beta,
  SE = gene_results$SE,
  t_value = gene_results$t_value,
  P_value = gene_results$P_value,
  BH_adjusted_P = gene_results$P_BH,
  Direction = gene_results$Direction,
  FDR_significant = gene_results$Significant_FDR05,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. Compare BH results with archived STEP7B6 result object
# ------------------------------------------------------------

saved_results <- purity_adjusted_results

saved_result_check <- "NOT_CHECKED"

if (is.data.frame(saved_results) &&
    all(c(
      "TME_Trait",
      "Gene",
      "Gene_FDR_BH"
    ) %in% colnames(saved_results))) {
  
  saved_results$Key <- paste(
    saved_results$TME_Trait,
    saved_results$Gene,
    sep = "__"
  )
  
  if (nrow(saved_results) != 9 ||
      anyDuplicated(saved_results$Key) > 0 ||
      !setequal(saved_results$Key, expected_keys)) {
    stop("Archived result object has unexpected keys.")
  }
  
  saved_aligned <- saved_results[
    match(gene_results$Key, saved_results$Key),
  ]
  
  bh_difference <- abs(
    gene_results$P_BH -
      saved_aligned$Gene_FDR_BH
  )
  
  bh_match <- all(
    bh_difference < 1e-12
  )
  
  saved_result_check <- if (bh_match) {
    "PASS"
  } else {
    "FAIL"
  }
  
  bh_audit <- data.frame(
    TME_Trait = as.character(gene_results$TME_Trait),
    Gene = as.character(gene_results$Gene),
    Model_P = gene_results$P_value,
    Recalculated_BH = gene_results$P_BH,
    Archived_BH = saved_aligned$Gene_FDR_BH,
    Absolute_Difference = bh_difference,
    Match = bh_difference < 1e-12,
    stringsAsFactors = FALSE
  )
  
  if (!bh_match) {
    stop("BH results do not match archived STEP7B6 results.")
  }
  
} else {
  stop(
    "Archived result object does not contain the required fields: ",
    "TME_Trait, Gene, Gene_FDR_BH. ",
    "Do not declare final PASS without archived-result comparison."
  )
}


# ------------------------------------------------------------
# 9. Final model and coefficient audits
# ------------------------------------------------------------

model_audit <- data.frame(
  Model = model_names,
  N = vapply(
    purity_adjusted_models,
    stats::nobs,
    numeric(1)
  ),
  Is_lm = vapply(
    purity_adjusted_models,
    inherits,
    logical(1),
    what = "lm"
  ),
  stringsAsFactors = FALSE
)

if (any(model_audit$N != 500) ||
    any(!model_audit$Is_lm)) {
  stop("Final model audit failed.")
}

if (any(!is.finite(gene_results$Beta)) ||
    any(!is.finite(gene_results$SE)) ||
    any(!is.finite(gene_results$P_value)) ||
    any(!is.finite(gene_results$P_BH))) {
  stop("Non-finite gene coefficient statistics detected.")
}

if (any(gene_results$SE <= 0) ||
    any(gene_results$P_value < 0 |
        gene_results$P_value > 1) ||
    any(gene_results$P_BH < 0 |
        gene_results$P_BH > 1)) {
  stop("Invalid coefficient statistics detected.")
}


# ------------------------------------------------------------
# 10. Export tables
# ------------------------------------------------------------

write.csv(
  publication_table,
  file.path(
    output_dir,
    "Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  coefficient_details,
  file.path(
    output_dir,
    "STEP7B6_05_Model_Coefficient_Details.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  bh_audit,
  file.path(
    output_dir,
    "STEP7B6_05_BH_Adjustment_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  variable_audit,
  file.path(
    output_dir,
    "STEP7B6_05_Dataset_Variable_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)

write.csv(
  model_audit,
  file.path(
    output_dir,
    "STEP7B6_05_Model_Audit.csv"
  ),
  row.names = FALSE,
  na = ""
)


# ------------------------------------------------------------
# 11. Export Excel workbook
# ------------------------------------------------------------

if (requireNamespace("openxlsx", quietly = TRUE)) {
  
  wb <- openxlsx::createWorkbook()
  
  openxlsx::addWorksheet(
    wb,
    "Supplementary_Table"
  )
  
  openxlsx::addWorksheet(
    wb,
    "All_Coefficients"
  )
  
  openxlsx::addWorksheet(
    wb,
    "BH_Audit"
  )
  
  openxlsx::addWorksheet(
    wb,
    "Dataset_Audit"
  )
  
  openxlsx::addWorksheet(
    wb,
    "Model_Audit"
  )
  
  openxlsx::writeData(
    wb,
    "Supplementary_Table",
    publication_table
  )
  
  openxlsx::writeData(
    wb,
    "All_Coefficients",
    coefficient_details
  )
  
  openxlsx::writeData(
    wb,
    "BH_Audit",
    bh_audit
  )
  
  openxlsx::writeData(
    wb,
    "Dataset_Audit",
    variable_audit
  )
  
  openxlsx::writeData(
    wb,
    "Model_Audit",
    model_audit
  )
  
  for (sheet in openxlsx::sheets(wb)) {
    
    openxlsx::freezePane(
      wb,
      sheet,
      firstRow = TRUE
    )
    
    openxlsx::setColWidths(
      wb,
      sheet,
      cols = 1:20,
      widths = "auto"
    )
  }
  
  openxlsx::saveWorkbook(
    wb,
    file.path(
      output_dir,
      "Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.xlsx"
    ),
    overwrite = TRUE
  )
  
} else {
  
  warning(
    "Package openxlsx is unavailable. ",
    "CSV files were exported; Excel workbook was not created."
  )
}


# ------------------------------------------------------------
# 12. Final audit summary
# ------------------------------------------------------------

overall_pass <- (
  dataset_N == 500 &&
    nrow(gene_results) == 9 &&
    nrow(purity_results) == 9 &&
    all(model_audit$N == 500) &&
    all(model_audit$Is_lm) &&
    all(gene_results$Significant_FDR05) &&
    saved_result_check == "PASS"
)

audit_lines <- c(
  "============================================================",
  "STEP 7B6_05 FINAL AUDIT",
  "============================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  "",
  paste("STEP7B5 dataset N:", dataset_N),
  paste("Unique Sample_ID:", length(unique(
    purity_adjusted_dataset$Sample_ID
  ))),
  paste("Archived models:", length(purity_adjusted_models)),
  paste("Gene coefficient rows:", nrow(gene_results)),
  paste("Purity coefficient rows:", nrow(purity_results)),
  paste("All models N=500:", all(model_audit$N == 500)),
  paste("All models are lm:", all(model_audit$Is_lm)),
  paste("All 9 gene terms BH < 0.05:",
        all(gene_results$Significant_FDR05)),
  paste("Archived BH comparison:", saved_result_check),
  "",
  paste("Overall status:",
        ifelse(overall_pass, "PASS", "FAIL")),
  "",
  "No model refitting was performed.",
  "No sample rematching was performed.",
  "No barcode remapping was performed.",
  "No purity imputation was performed.",
  "Original RData files were not modified.",
  "============================================================"
)

writeLines(
  audit_lines,
  file.path(
    output_dir,
    "STEP7B6_05_Final_Audit.txt"
  )
)

cat(
  paste(audit_lines, collapse = "\n"),
  "\n"
)

if (!overall_pass) {
  stop("STEP7B6_05 final audit did not pass.")
}

cat("\nSTEP7B6_05 COMPLETED: PASS\n")
cat("Output directory:\n", output_dir, "\n")
