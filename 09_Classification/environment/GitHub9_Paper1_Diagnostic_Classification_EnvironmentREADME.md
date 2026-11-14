# Environment and Software Requirements

## 1. Overview

This directory documents the R environment, software dependencies, and reproducibility requirements for the diagnostic classification workflows in `09_Classification`.

The workflows cover historical data and object-provenance audits, GSE32863 metadata and sample-mapping audits, model-definition and model-object checks, cross-validation, external assessment, and figure generation.

The documentation distinguishes between:

- Packages used across multiple workflow components.
- Packages required only by specific analysis or visualization steps.
- Base and recommended R functionality.
- Runtime environment information recorded during script execution.

The exact R version, platform, and package versions must be obtained from the corresponding `sessionInfo()` output generated during execution. No unverified version numbers are specified in this document.

## 2. R Environment

The classification workflows are implemented in R.

The project uses `here` to identify the repository root and construct project-relative paths. Most data inspection, object handling, statistical operations, and file management are performed using base R or recommended R packages.

The R installation must support the packages explicitly required by the individual workflow components.

### 2.1 Project-Relative Path Management

The project root is resolved using:

```r
project_root <- here::here()
```

Input, result, and workflow directories are constructed using `file.path()`.

The principal directory structure is:

```text
09_Classification/
├── data/
│   ├── input_data/
│   ├── processed_data/
│   └── sample_manifests/
├── results/
│   ├── audit_reports/
│   ├── figures/
│   ├── statistical_outputs/
│   ├── tables/
│   └── WORKFLOW13_Diagnostic_Classification/
└── Environment/
```

A typical path configuration is:

```r
project_root <- here::here()

input_dir <- file.path(
  project_root,
  "09_Classification",
  "data",
  "input_data"
)

result_dir <- file.path(
  project_root,
  "09_Classification",
  "results"
)

workflow_root <- file.path(
  result_dir,
  "WORKFLOW13_Diagnostic_Classification"
)
```

The scripts should be executed from the intended Git repository or an R project recognized by `here`. The project-relative path configuration does not require changing the working directory with `setwd()`.

## 3. R Package Dependencies

### 3.1 General Project Dependency

| Package | Purpose | Usage |
|---|---|---|
| `here` | Project-root detection and project-relative path construction | Used across multiple classification workflow components |

The `here` package is the principal shared project-path dependency.

### 3.2 Step-Specific Analysis Dependencies

The following packages are required only by the workflow components in which they are explicitly used.

| Package | Purpose | Workflow component |
|---|---|---|
| `readxl` | Reading historical Excel workbooks | STEP13B-7J-R1 |
| `openxlsx` | Excel workbook processing | STEP13B-7J-R1 |
| `pROC` | ROC curve construction, AUC calculation, and confidence-interval estimation | STEP13G and STEP13I |
| `pROC` | Included in the required-package list and loaded by the script; no ROC calculation is performed in this step | STEP13H |
| `ggplot2` | ROC and paired-expression visualization | Figures 9A–9C |
| `Cairo` | PDF and TIFF figure export | Figures 9A–9C |

These dependencies should not be interpreted as requirements for every script in the classification module.

In particular:

- `readxl` and `openxlsx` are associated with the historical Excel audit workflow.
- `pROC` is used for ROC/AUC analysis in the specified analytical workflows. Its presence in STEP13H does not mean that STEP13H independently performs ROC analysis.
- `ggplot2` and `Cairo` are specific to the figure-generation workflow.

No additional specialized statistical or bioinformatics package should be inferred unless it is explicitly loaded or called by the corresponding script.

### 3.3 Recommended R Packages

| Package | Purpose | Workflow component |
|---|---|---|
| `tools` | Standard file and file-extension utilities | STEP13C-1 |
| `utils` | General-purpose R utilities, including `head()` and `tail()` | STEP13C-4B |

These packages are part of the standard R distribution and should be distinguished from additional third-party package dependencies.

## 4. Base R Functionality

A substantial proportion of the classification workflows relies on base R functionality.

The principal operations include:

| Functionality | Representative functions | Application |
|---|---|---|
| Project and file management | `file.path()`, `file.exists()`, `dir.exists()`, `dir.create()`, `file.remove()` | Path construction and file-system checks |
| Data import and export | `read.csv()`, `write.csv()`, `readLines()`, `writeLines()`, `gzfile()` | Tabular, text, and compressed-file handling |
| R object serialization | `load()`, `saveRDS()`, `readRDS()` | Historical workspace and R object handling |
| Environment isolation | `new.env()`, `exists()`, `get()`, `ls()` | Loading and inspecting historical objects |
| Data inspection | `class()`, `dim()`, `length()`, `names()`, `vapply()` | Object and data-structure checks |
| Data comparison | `identical()`, `all.equal()`, `intersect()` | Historical object and value comparisons |
| String and metadata processing | `grep()`, `grepl()`, `sub()`, `gsub()`, `regmatches()` | Metadata extraction and identifier matching |
| Statistical analysis | `glm()`, `predict()`, `wilcox.test()`, `p.adjust()`, `complete.cases()` | Model fitting, prediction, paired comparisons, and multiple-testing correction |
| Reproducibility and reporting | `set.seed()`, `capture.output()`, `sessionInfo()` | Random-number control and runtime environment recording |

The table provides representative functions rather than an exhaustive inventory of every function used across all scripts.

## 5. Workflow-Specific Computational Requirements

### 5.1 Historical Data and Object-Provenance Audits

The historical data and object-audit workflows primarily use base R functionality, together with `here` for project-relative path management.

Depending on the workflow, these steps inspect:

- Historical `.RData` workspaces.
- Existing R objects, their classes, dimensions, and names.
- Archived GEO metadata files.
- Sample and subject identifiers.
- Historical analysis artifacts and their file structures.

Historical `.RData` workspaces may be loaded into isolated R environments using `new.env()` and `load()`.

This approach allows the scripts to inspect historical objects without directly overwriting objects in the active global environment.

These audit procedures should not be interpreted as recalculating historical expression values, refitting models, or reproducing historical ROC/AUC results unless the corresponding workflow explicitly performs those operations.

### 5.2 GSE32863 Metadata and Sample-Mapping Audits

The GSE32863 metadata-audit workflows use base R for compressed-file reading, metadata extraction, string matching, sample-code inspection, and output generation.

The principal local metadata files include:

- `GSE32863_series_matrix.txt.gz`
- `GSE32863_family.soft.gz`

The workflows use metadata fields and accession identifiers to assess sample and subject correspondence.

Some consecutive steps depend on objects generated by preceding steps. For example, the STEP13C-4E workflow requires `mapping_df` from the preceding sample-code mapping workflow.

Such dependent steps should be executed in the required sequence within the same R session, unless their required objects are restored separately.

These metadata audits do not independently reconstruct or analyze expression data unless explicitly stated in the corresponding analytical workflow.

### 5.3 Model-Definition and Historical Model-Object Audits

The model-definition and historical model-object audit workflows use base R functions for:

- Loading historical objects.
- Checking required object availability.
- Inspecting model classes and structures.
- Extracting available formulas and coefficients.
- Comparing stored values and model-related objects.
- Writing audit tables and serialized R objects.

The STEP13F model-definition lock does not perform model fitting, cross-validation, ROC/AUC calculation, or threshold optimization.

The STEP13E-3 historical object audit examines previously stored model-related objects and their available methods. It does not independently refit the historical models or recalculate the original ROC/AUC results.

### 5.4 Cross-Validation and Model Assessment

The STEP13G workflow uses `pROC` for ROC analysis and AUC estimation, together with base R functions for logistic regression, prediction, stratified fold assignment, and performance evaluation.

The recorded random seed is:

```r
set.seed(20260926)
```

The workflow uses stratified outer-fold assignment. Its implementation should be described according to the operations actually performed in the script; inner-fold stratification auditing should not be interpreted as evidence that a complete inner model-fitting and tuning procedure was executed.

The STEP13H and STEP13I workflows use the previously specified model definition and stored model-related information according to their respective workflow requirements.

STEP13I uses `pROC` for external ROC/AUC assessment. Its external assessment should not be interpreted as authorization to refit the external model, re-estimate coefficients, reselect features, or optimize the classification threshold.

Where specified in the corresponding scripts, the random seed is:

```r
20260926
```

The seed should be applied only to the workflow components that actually use random-number generation.

### 5.5 Figure Generation

The figure-generation workflow for Figures 9A–9C uses:

- `ggplot2` for graphical construction.
- `pROC` for ROC-object construction and coordinate extraction.
- `Cairo` for PDF and TIFF export.
- `here` for project-relative path resolution.

The visualization workflow reads previously generated prediction and paired-expression files.

It does not independently refit the classification model, recalculate cross-validation results, optimize thresholds, or repeat the upstream paired statistical tests and multiple-testing correction.

The plotting script specifies Arial as the font family and uses Windows-specific font registration functions such as `windowsFonts()` and `windowsFont()`.

Consequently, reproducing the figures on non-Windows operating systems may require an appropriate font-configuration adjustment. Any such adjustment should be documented separately and should not silently alter the analytical results or figure content.

## 6. Randomness and Reproducibility Settings

Random-number generation is not required by every classification workflow.

The principal explicitly documented random seed for the stochastic classification procedures is:

| Workflow | Seed | Purpose |
|---|---:|---|
| STEP13G | `20260926` | Stratified fold-generation procedure |
| STEP13H | `20260926` | Workflow-specific reproducibility setting |
| STEP13I | `20260926` | Workflow-specific reproducibility setting |

The seed should be retained in the relevant scripts where random-number generation is performed or explicitly required by the workflow.

Deterministic metadata audits, object comparisons, and file-integrity checks do not require a random seed.

The absence of a random seed in a workflow should not be interpreted as an omission when no stochastic operation is performed.

## 7. Session Information and Software Version Recording

### 7.1 Version Recording Policy

The exact R version and package versions must be obtained from the actual execution environment.

The classification scripts use `sessionInfo()` and, where appropriate, `capture.output()` to record:

- R version.
- Operating system and platform information.
- Attached packages.
- Loaded namespaces.
- Package versions available during execution.

Version information should be taken from the generated session-information files rather than manually inferred or assigned in this README.

### 7.2 Session-Information Files

The following session-information records are identified in the workflow documentation.

| Workflow component | Session-information file | Recording notes |
|---|---|---|
| Initial configuration | `WORKFLOW13_sessionInfo_CONFIG.txt` | Initial workflow environment record |
| STEP13B-7J-R1 | `STEP13B-7J-R1_sessionInfo.txt` | Historical Excel traceability workflow |
| STEP13C-1 | `STEP13C-1_sessionInfo.txt` | Workflow-specific environment record |
| STEP13C-2 | `STEP13C-2_sessionInfo.txt` | Workflow-specific environment record |
| STEP13C-3 | `STEP13C-3_sessionInfo.txt` | Session information is written by the script |
| STEP13C-4A | No separate file documented | No independent session-information export is currently documented |
| STEP13C-4F | `STEP13C-4F_sessionInfo.txt` | Subject-ID concordance audit |
| STEP13C-5 | Workflow-specific session-information file | Recorded by the workflow |
| STEP13E-1 | Workflow-specific session-information file | Recorded by the workflow |
| STEP13E-2 | Workflow-specific session-information file | Recorded by the workflow |
| STEP13E-3 | `STEP13E-3_sessionInfo.txt` | Historical model-object audit |
| STEP13F | `STEP13F_sessionInfo_FINAL.txt` / `STEP13F_sessionInfo_FAIL.txt` | Final or early-failure environment record |
| STEP13G | `STEP13G_sessionInfo_FINAL.txt` | Written when execution reaches the session-information export step |
| STEP13H | `STEP13H_sessionInfo_FINAL.txt` | Final model workflow environment record |
| STEP13I | `STEP13I_GSE32863_External_Assessment/STEP13I_sessionInfo_FINAL.txt` | External assessment environment record |

The initial configuration record is located at:

```text
results/WORKFLOW13_Diagnostic_Classification/LOGS/
WORKFLOW13_sessionInfo_CONFIG.txt
```

The STEP13H record is located at:

```text
results/WORKFLOW13_Diagnostic_Classification/
STEP13H_TCGA_Final_Model/STEP13H_sessionInfo_FINAL.txt
```

The STEP13I record is located at:

```text
results/WORKFLOW13_Diagnostic_Classification/
STEP13I_GSE32863_External_Assessment/
STEP13I_sessionInfo_FINAL.txt
```

For other workflow-specific records, the exact output location should be confirmed from the corresponding script's file-writing statements.

The Figures 9A–9C visualization code does not generate a `sessionInfo()` output or a package-version lockfile. Its actual R and package versions should therefore be recorded separately from the execution environment if formal version documentation is required.

### 7.3 Retention of Environment Records

Session-information files should be retained with the corresponding workflow outputs.

If a workflow terminates before reaching its final session-information export step, the absence of a final environment record should be documented rather than replaced with an assumed version record.

For workflows that generate both final and failure records, the appropriate record should be retained to reflect the actual execution status.

## 8. Package Installation and Environment Management

The workflow documentation describes the packages used by the scripts; it does not establish a universal package-installation or package-update procedure.

The scripts should not be assumed to install or update packages automatically.

Before executing a workflow:

1. Ensure that R is available.
2. Ensure that the packages explicitly required by that workflow are installed.
3. Confirm that the repository root is correctly recognized by `here`.
4. Ensure that the required input files and upstream workflow outputs are available.
5. Retain the generated session-information files for version traceability.

Package versions should not be changed solely to match an assumed version listed in an older analysis record. Any environment change that could affect reproducibility should be documented.

## 9. Reproducibility Notes

- Use the intended Git repository root for project-relative path resolution.
- Use `here::here()` and `file.path()` rather than hard-coded absolute paths.
- Install only the packages required by the workflow being executed.
- Distinguish third-party packages from base and recommended R packages.
- Execute sequential workflows in the required order when later steps depend on objects generated by earlier steps.
- Load historical workspaces into isolated environments when required by the audit procedure.
- Preserve historical input objects and locked upstream outputs during provenance and reproducibility audits.
- Apply random seeds only to workflows that use stochastic procedures.
- Do not interpret metadata audits or historical object inspections as model refitting or independent recalculation of previous analytical results.
- Retain the actual `sessionInfo()` output for each workflow that generates one.
- Do not infer package versions, R versions, or operating-system details that have not been verified from the execution environment.

This document describes the software environment and computational requirements of the classification module. Detailed sample provenance, input-file inventories, analysis procedures, statistical results, and output-file descriptions are maintained in the corresponding `data/README.md`, `Script/README.md`, and `Results/README.md` files.