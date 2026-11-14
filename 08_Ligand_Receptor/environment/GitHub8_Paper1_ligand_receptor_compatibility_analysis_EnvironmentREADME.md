# Environment and Software Dependencies

## 1. Overview

This directory documents the R environment, software dependencies, and computational requirements for the ligand-receptor compatibility analysis in the LUAD study.

The workflow integrates the predefined CCI_LR ligand-receptor resource with previously generated STEP10 single-cell pseudobulk expression data. It performs structural audits, expression-based ligand-receptor compatibility assessment, patient-level descriptive characterization, candidate-set integration, figure generation, supplementary table construction, and traceability checks.

The analysis is implemented primarily in R. Most structural and descriptive workflows use Base R, while selected steps require additional packages for efficient data manipulation, visualization, or Excel workbook generation.

## 2. Study

**Study title:**

Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD.

**Analysis module:**

Paper 1 — Workflow STEP11: Ligand-Receptor Compatibility Analysis.

**Module directory:**

`08_Ligand_Receptor/`

## 3. Software Environment

| Component | Role | Requirement |
|----|----|----|
| R | Statistical computing and analysis environment | Required |
| `here` | Project-root identification and project-relative path construction | Required by shared initialization |
| Base R | Data import, matrix operations, structural validation, set operations, descriptive calculations, and output generation | Required |
| `utils` | Recommended R utilities, including functions such as `head()` and `write.table()` | R recommended package |
| `tools` | File and package-related utility functionality where used | R recommended package |

The integrated script uses `here::here()` for project-root identification and project-relative path construction. The documented workflows do not require changing the working directory with `setwd()`.

## 4. Required R Packages

### 4.1 Project Path Management

| Package | Purpose | Workflow |
|----|----|----|
| `here` | Identifies the repository root and constructs project-relative paths | Shared script initialization |

### 4.2 Data Manipulation and Integration

| Package | Purpose | Workflow |
|----|----|----|
| `data.table` | Efficient file reading and writing, grouping, duplicate inspection, canonical-key construction, table integration, sorting, reshaping, and QC summaries | STEP11-5E-3, STEP11-5F-1, STEP11-5F-2, STEP11-5F-3, STEP11-5F-4, and related audits |
| `dplyr` | Manuscript-facing table construction and data manipulation | STEP11-5F-4D |

Important `data.table` functions used across the documented workflows include:

- `fread()`
- `fwrite()`
- `copy()`
- `uniqueN()`
- `rbindlist()`
- `setorder()`
- `setcolorder()`
- `data.table()`
- `dcast()`
- `first()`
- `merge()`

The `dplyr` package is used for constructing and organizing the final manuscript-facing supplementary table.

### 4.3 Visualization

| Package   | Purpose                                | Workflow              |
|-----------|----------------------------------------|-----------------------|
| `ggplot2` | Generation of the Figure 7 panels      | STEP11-5F-5           |
| `scales`  | Scale and label formatting for figures | STEP11-5F-5 Version 2 |

The Figure 7 plotting workflows target Arial or Arial Bold fonts. The availability of these fonts depends on the execution system.

### 4.4 Excel Workbook Generation and Inspection

| Package | Purpose | Workflow |
|----|----|----|
| `openxlsx` | Excel workbook generation, worksheet formatting, column-width configuration, and workbook export | STEP11-5F-4D |
| `openxlsx` | Conditional reading and inspection of an XLSX workbook | STEP11-6A, when the manually selected Figure 7 source or supplementary workbook is in XLSX format |

The `openxlsx` dependency is conditional for the Figure 7 traceability workflow when an XLSX input is selected. It is required for the final supplementary Excel workbook generation.

## 5. Workflow-Specific Computational Requirements

### 5.1 STEP11-0 to STEP11-3: Input and Sample-Structure Audits

These steps primarily use Base R and recommended R utilities for:

- Input-file existence and structural inspection.
- Data-frame and matrix inspection.
- Missing-value and data-type summaries.
- Patient, sample, and cell-lineage distribution checks.
- Patient-level and sample-level cell-count summaries.
- Descriptive threshold and eligibility audits.
- RDS and TSV file import and export.

Representative functions include:

`file.exists()`, `readRDS()`, `saveRDS()`, `read.delim()`, `write.table()`, `matrix()`, `table()`, `aggregate()`, `merge()`, `quantile()`, and `rowSums()`.

No additional specialized statistical or single-cell analysis package is introduced by these steps.

No random sampling or stochastic procedure is specified. A random seed is not required.

### 5.2 STEP11-4: Core-Gene Cellular Localization Integration

This step primarily uses Base R for:

- Reading and inspecting existing expression and metadata tables.
- Integrating patient-level and compartment-level information.
- Summarizing core-gene expression and detection.
- Evaluating patient-level localization patterns and structural concordance.
- Exporting descriptive summary tables and an RDS object.

Representative functions include:

`read.delim()`, `file.exists()`, `list.files()`, `grepl()`, `aggregate()`, `merge()`, `table()`, `order()`, `duplicated()`, `write.table()`, and `saveRDS()`.

The shared `here` initialization is used for project-root identification.

No random sampling or stochastic procedure is specified.

### 5.3 STEP11-5A to STEP11-5C: Resource Compatibility and Expression-Scale Audits

These steps primarily use Base R for:

- Reading the CCI_LR resource and existing STEP10 expression matrices.
- Comparing gene identifiers and compartment-specific gene sets.
- Constructing predefined sender-receiver combinations.
- Inspecting matrix dimensions and numeric integrity.
- Matching patients and compartments.
- Calculating library-size and CPM-related descriptive quantities.
- Applying the documented log2(CPM + 1) transformation where applicable.
- Evaluating expression detection and expression-scale consistency.
- Generating descriptive compatibility summaries.

Representative functions include:

`read.delim()`, `file.exists()`, `unique()`, `sort()`, `setdiff()`, `intersect()`, `expand.grid()`, `aggregate()`, `table()`, `mapply()`, `sapply()`, `as.matrix()`, `colSums()`, `sweep()`, `log2()`, `quantile()`, `median()`, `cor()`, and `write.table()`.

No additional specialized analysis package is introduced by STEP11-5A, STEP11-5B, or STEP11-5C.

These steps use existing expression data and predefined resource information. They do not perform de novo ligand-receptor inference or infer functional signaling activity.

No random sampling or stochastic procedure is specified.

### 5.4 STEP11-5D: Detection Sensitivity, Genome-Wide Screening, and Candidate-Set Definition

These steps primarily use Base R for:

- Evaluating predefined ligand and receptor detection rules.
- Summarizing patient-level compatibility across six predefined directions.
- Comparing candidate structures and exact ligand-receptor pairs.
- Auditing ligand-level, direction-level, annotation-level, and pathway-label structures.
- Defining and validating candidate sets using predefined recurrence criteria.
- Exporting descriptive audit tables and locked candidate-set records.

Representative functions include:

`read.delim()`, `file.exists()`, `dir.create()`, `unique()`, `sort()`, `intersect()`, `setdiff()`, `identical()`, `order()`, `aggregate()`, `table()`, `saveRDS()`, and `write.table()`.

No additional statistical or ligand-receptor inference package is introduced by these steps.

No random-number generation or random seed is required.

### 5.5 STEP11-5E: Candidate Overlap, Expression Magnitude, and Patient-Level Integration

STEP11-5E-1 and STEP11-5E-2 primarily use Base R for:

- Exact ligand-receptor pair identification.
- Candidate-set overlap and set comparison.
- Patient-compatibility pattern summaries.
- Direction-specific descriptive characterization.
- Candidate expression-magnitude summaries.
- Expression-field validation and descriptive distribution calculations.

Representative functions include:

`read.delim()`, `file.exists()`, `dir.create()`, `trimws()`, `unique()`, `sort()`, `intersect()`, `setdiff()`, `union()`, `rowSums()`, `table()`, `lapply()`, `sapply()`, `do.call()`, `as.matrix()`, `rowMeans()`, `median()`, `IQR()`, `min()`, `max()`, `merge()`, `aggregate()`, and `write.table()`.

STEP11-5E-3 additionally uses `data.table` for:

- Efficient input and output operations.
- Duplicate-key and canonical-key inspection.
- Candidate-set provenance assessment.
- Invariance and protected-variable conflict audits.
- Canonical patient-level integration.
- Structural QC and provenance summaries.

Representative `data.table` functions include:

`fread()`, `fwrite()`, `copy()`, `uniqueN()`, `rbindlist()`, `setorder()`, `setcolorder()`, `data.table()`, `merge()`, `first()`, `setdiff()`, and `intersect()`.

No stochastic procedure is used in these workflows. No random seed is required.

### 5.6 STEP11-5F-1 to STEP11-5F-4: Canonical Input QC and Patient-Level Context Integration

These steps use `data.table` for:

- Reading the canonical patient-level LR integration table.
- Validating input structure and candidate classes.
- Summarizing direction and patient coverage.
- Characterizing patient-level compatibility patterns.
- Integrating patient-level expression context.
- Constructing exact-pair and direction-level descriptive summaries.
- Generating the final master descriptive context table.
- Exporting TSV tables and QC summaries.

Representative functions include:

`fread()`, `fwrite()`, `uniqueN()`, `setdiff()`, `unique()`, `sort()`, `paste()`, `identical()`, `trimws()`, `setorder()`, `data.table()`, `rbindlist()`, and `dcast()`.

Base R functions are also used for file and directory operations, data validation, numeric checks, canonical-key construction, and console reporting.

The STEP11-5F-4C manuscript-table and direction-mapping audits use deterministic table inspection and integrity checks. They do not require additional statistical modeling packages.

The STEP11-5F-4D supplementary table generation uses:

- `dplyr` for manuscript-facing table construction.
- `openxlsx` for Excel workbook generation, worksheet formatting, column-width configuration, and workbook export.

No random sampling or stochastic procedure is specified in these workflows.

### 5.7 STEP11-5F-5: Figure 7 Generation

The Figure 7 plotting workflows use:

| Package      | Purpose                                      |
|--------------|----------------------------------------------|
| `data.table` | Reading and summarizing the locked LR tables |
| `ggplot2`    | Generating the Figure 7 panels               |
| `scales`     | Scale and label formatting in Version 2      |
| `here`       | Project-relative path initialization         |

The workflow generates PDF and TIFF figure files. The target font family is Arial, with Arial Bold used where specified.

Version 1 and Version 2 are retained as separate plotting outputs. Their files should not be treated as interchangeable or automatically overwritten.

The plotting steps use existing locked analytical tables and perform figure-oriented presentation and descriptive summaries. They do not introduce a stochastic procedure.

No random seed is required.

### 5.8 STEP11-6A: Figure 7 Traceability Audits

The Figure 7 traceability workflows primarily use Base R for:

- Candidate-file discovery and file-existence checks.
- Reading TSV, CSV, and RDS files.
- Inspecting LR-pair, direction, annotation, and compatibility fields.
- Comparing figure-source records with the final master table.
- Checking direction mappings and patient-level consistency.
- Generating traceability reports and audit objects.

Representative functions include:

`file.exists()`, `list.files()`, `read.delim()`, `read.csv()`, `readRDS()`, `saveRDS()`, `table()`, `sapply()`, `mapply()`, `all.equal()`, `writeLines()`, and `stopifnot()`.

The `openxlsx` package is required only when an XLSX workbook is selected for inspection.

The STEP11-6A-0 and STEP11-6A-1 source-identification and structural-audit operations primarily use Base R and do not introduce additional package dependencies.

### 5.9 STEP11-6A-2 to STEP11-6A-10: Supplementary Table Traceability Diagnostics

These diagnostics use:

- `openxlsx` to read the final Supplementary Table workbook.
- Base R to read the locked master table and perform exact-pair, string, positional, and direction-label comparisons.
- `here` through the shared project initialization for project-relative paths.

The diagnostics include:

- Exact ligand-receptor pair identity checks.
- Order-independent pair-set comparisons.
- Missing, unexpected, and duplicate pair checks.
- Normalized string and positional comparisons.
- Ligand-receptor identifier decomposition.
- Validation of the six predefined direction labels.
- ASCII-normalized direction-label comparisons, including normalization of the arrow symbol for comparison.
- Direction-frequency consistency checks.

These are deterministic, read-only traceability operations. They do not perform statistical modeling, random sampling, expression recalculation, or compatibility recalculation.

No random seed is required.

## 6. Random Seed and Stochastic Procedures

No random seed is specified or required for the documented STEP11 workflows.

The reviewed operations are primarily deterministic structural audits, set comparisons, expression-based descriptive calculations, candidate-set integration, figure generation, and traceability checks.

No random sampling, stochastic model fitting, or randomized permutation procedure is specified in the documented workflow.

If stochastic procedures are introduced in a future revision, the corresponding random seed and reproducibility settings should be documented in the relevant workflow section.

## 7. Input Data and External Resources

The analysis uses existing project data and predefined resources rather than requiring a new external analysis platform.

The main input categories include:

- The predefined CCI_LR ligand-receptor resource.
- Existing STEP10 single-cell pseudobulk expression matrices.
- Existing STEP10 cell metadata and patient/sample annotations.
- Previously generated STEP11 audit tables and candidate-set records.
- The locked canonical patient-level LR integration table.
- The final Figure 7 source tables and supplementary Excel workbook for traceability checks.

The CCI_LR resource is used as a predefined reference. The workflow assesses expression-based ligand-receptor compatibility using the existing resource and expression data; it does not constitute de novo ligand-receptor interaction inference.

No additional external software platform is identified as a required dependency in the reviewed script sections.

## 8. Reproducibility and Execution Notes

- Run the integrated script from the Git repository using the shared `here::here()` project-root initialization.
- Use project-relative paths for input and output files.
- Do not rely on a machine-specific absolute path or an implicit working-directory change.
- Keep the existing STEP10 input files and locked STEP11 intermediate outputs consistent with the documented workflow.
- Preserve the separate Figure 7 Version 1 and Version 2 output directories.
- Treat files with a `LOCKED` suffix as workflow-specific saved audit or candidate-set records; the suffix alone does not imply biological validation.
- Interpret `PASS` as satisfying the programmed structural or workflow QC checks, not as evidence of functional molecular interaction or biological validation.
- Treat pathway and annotation summaries as descriptive classifications rather than statistical pathway enrichment or pathway activation analyses.
- Record the actual R version and installed package versions from the environment used for the final analysis run.

## 9. Environment Version Record

The exact software versions should be captured from the finalized execution environment rather than inferred from the script.

| Component    | Version                                             |
|--------------|-----------------------------------------------------|
| R            | To be recorded from the final execution environment |
| `here`       | To be recorded                                      |
| `data.table` | To be recorded                                      |
| `ggplot2`    | To be recorded                                      |
| `scales`     | To be recorded                                      |
| `dplyr`      | To be recorded                                      |
| `openxlsx`   | To be recorded                                      |

Base R and recommended-package versions are associated with the installed R distribution.

A reproducible environment record can be generated during the final execution using:

``` r
R.version.string

packageVersion("here")
packageVersion("data.table")
packageVersion("ggplot2")
packageVersion("scales")
packageVersion("dplyr")
packageVersion("openxlsx")

sessionInfo()
```

The version checks should be limited to packages actually installed and required by the final script configuration.

## 10. Documentation Scope and Limitations

This README summarizes the software dependencies and computational requirements identified in the reviewed sections of `GitHub8_Paper1_ligand_receptor_compatibility_analysis.R`.

It does not claim that every package listed is required by every workflow step. Package requirements are step-specific, and `openxlsx` is conditionally required for XLSX inspection in the traceability workflow.

The analysis is based on predefined ligand-receptor annotations and existing expression data. Expression-based compatibility does not establish direct physical binding, interaction strength, communication probability, pathway activation, causal signaling, or experimental validation.

No statistical significance testing, multiple-testing correction, or expression-matched permutation-based null model is specified for the principal LR compatibility workflow.

For the complete analysis workflow and exported result inventory, refer to:

- `../scripts/README.md`
- `../results/README.md`

For the overall project organization and input-data provenance, refer to the repository-level documentation.
