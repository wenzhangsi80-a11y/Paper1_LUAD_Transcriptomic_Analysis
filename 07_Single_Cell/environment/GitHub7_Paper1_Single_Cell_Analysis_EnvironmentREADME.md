# Software Environment

## 1. Overview

This directory documents the software environment and R package dependencies used in the Paper 1 single-cell transcriptomic analysis.

The analysis is implemented in two R scripts:

- `GitHub7_Paper1_Single_Cell_Analysis_1.R`
- `GitHub7_Paper1_Single_Cell_Analysis_2.R`

The workflows cover single-cell object handling, cell-type visualization, core-gene expression visualization, patient and sample structure audits, raw-count matrix inspection, cell-ID reconciliation, patient-level raw UMI pseudobulk construction, descriptive expression characterization, evidence synthesis, and manuscript supplementary exports.

This README consolidates the environment information from both scripts and their associated workflow sections.

It distinguishes between:

1. Project-level R packages.
2. STEP-specific and conditional package dependencies.
3. Base R functionality.
4. Script-declared R version and random seed information.
5. Actual execution-time version records.

Package versions and operating-system details are not inferred from package-loading commands or script headers alone.

---

## 2. Programming Language and Runtime

### Programming language

- R

### R version

- R version specified in the Script 2 header: `4.5.1`.

The uploaded Script 2 sections for STEP10F-1 and STEP10F-3 document R 4.5.1.

This is a script-declared version and should not be interpreted as independent verification of the runtime used for every analysis stage.

The actual R version used for each completed workflow should be confirmed using the corresponding execution-time session information.

### Operating system

The operating system and platform-specific runtime configuration have not been fully documented in the uploaded environment records.

The Windows system version, architecture, locale, and font configuration should be recorded from the actual execution environment when available.

---

## 3. Project Path Management

The project uses the `here` package for project-relative path management.

The main project setup defines shared path variables, including:

- `project_root`
- `input_dir`
- `result_dir`

The STEP10 workflows use these shared configuration objects to construct input and output paths.

The project setup script `00_setup.R` provides the path configuration used by downstream analysis sections.

The STEP10E-0 provenance audit specifically requires `project_root` to be defined before execution.

STEP10D-D and subsequent workflows also depend on the appropriate project-level path variables being initialized.

Project-relative paths should be used instead of machine-specific absolute paths wherever supported by the project setup.

---

## 4. Core R Packages

The following packages are used across the two single-cell analysis scripts.

| Package | Main purpose | Workflow |
|---|---|---|
| `Seurat` | Construction and handling of single-cell objects and expression data | Script 1 |
| `Matrix` | Sparse matrix operations and raw-count pseudobulk aggregation | Script 1; STEP10C |
| `ggplot2` | Data visualization and publication-oriented figures | Script 1; STEP10H |
| `ggrepel` | Non-overlapping cell-type labels in UMAP visualizations | Script 1 |
| `dplyr` | Data manipulation, grouping, counting, joining, and cohort construction | Script 1; STEP10C-0E |
| `data.table` | Efficient metadata and tabular data import and processing | Script 1; STEP10C and STEP10H |
| `tidyverse` | Data preparation for Figure 9C | Script 1 |
| `scales` | Expression color-scale handling and color palette functions | Script 1; STEP10H |
| `here` | Project-root detection and project-relative path construction | Both scripts; STEP10H |
| `readr` | Tabular data import | STEP10C-0E |
| `tidyr` | Patient-by-compartment table reshaping | STEP10C-0E |
| `openxlsx` | Excel supplementary table generation | STEP10H |
| `rhdf5` | HDF5 structure inspection and selected matrix reading | STEP10C-0A, STEP10C-0B, STEP10C-0C |
| `zip` | ZIP archive inspection | STEP10C-0A |

The package inventory reflects the dependencies identified in the uploaded script sections. It does not imply that every package is required for every workflow stage.

Some packages are required only by specific steps, and some operations are conditional on input-file availability or package availability.

---

## 5. STEP-Specific Software Dependencies

### 5.1 Script 1 — Single-Cell Visualization

`GitHub7_Paper1_Single_Cell_Analysis_1.R`

The uploaded sections identify the following packages:

- `Seurat`
- `Matrix`
- `ggplot2`
- `ggrepel`
- `dplyr`
- `data.table`
- `tidyverse`
- `scales`
- `here`

Main functions include:

- Single-cell object and expression-data handling.
- UMAP visualization.
- Cell-type label placement.
- Core-gene expression visualization.
- Data preparation for Figure 9C.
- Project-relative path management.

The random seed documented for Script 1 is:

- `20250724`

The actual R and package versions should be recorded from the execution environment.

### 5.2 STEP10A — Patient / Sample Structure Audit

Required external package:

- `here`

The workflow primarily uses base R functions for metadata reading, patient and sample summaries, table generation, and reproducibility records.

Representative functions include:

- `read.delim()`
- `data.frame()`
- `table()`
- `aggregate()`
- `merge()`
- `write.table()`
- `saveRDS()`
- `sessionInfo()`

No random procedure is introduced by the uploaded STEP10A section. A random seed is therefore not required.

The workflow uses TISCH2-provided metadata and does not recalculate dimensionality reduction or cell annotations.

The generated `STEP10A_sessionInfo.txt` file records the session information associated with the audit execution.

### 5.3 STEP10B-1 and STEP10B-2 — Compartment Feasibility and Cell-Depth Audits

These workflows primarily use base R functions for:

- Metadata processing.
- Cell-count summaries.
- Data reshaping.
- Patient-level cell-depth reconciliation.
- Threshold classification.
- Descriptive summaries.
- Table export and RDS serialization.
- Session information recording.

Representative functions include:

- `read.delim()`
- `table()`
- `reshape()`
- `aggregate()`
- `combn()`
- `quantile()`
- `saveRDS()`
- `sessionInfo()`

No additional third-party package is explicitly required by the uploaded STEP10B-1 and STEP10B-2 sections beyond the shared project setup dependency on `here`.

The workflows are deterministic and do not require a random seed.

### 5.4 STEP10C — Raw-Count Matrix Inspection and Pseudobulk Construction

#### STEP10C — Sparse matrix handling

Required package:

- `Matrix`

The `Matrix` package is used for sparse count-matrix handling, including Matrix Market input and sparse row-sum aggregation.

The supported input formats identified in the uploaded sections include:

- RDS.
- Matrix Market (`.mtx`).

The original raw-count matrix must contain gene identifiers as row names and cell identifiers as column names.

Other operations use base R functions, including:

- `saveRDS()`
- `write.table()`
- `data.frame()`
- `expand.grid()`
- `mapply()`
- `rbind()`
- `colSums()`
- `sink()`
- `sessionInfo()`

No stochastic procedure is identified in the uploaded sparse matrix and pseudobulk operations.

#### STEP10C-0A and STEP10C-0B — HDF5 and ZIP inspection

Required or conditionally used packages:

- `rhdf5`
- `zip`

`rhdf5` is used for HDF5 structure inspection and selected matrix-value reading.

The STEP10C-0A section also checks ZIP archive contents.

STEP10C-0B reads selected matrix values and cell-level entries rather than loading the complete expression matrix into memory.

These steps use base R functions for file inspection, data summaries, and tabular output.

No random sampling or stochastic procedure is introduced.

#### STEP10C-0B and STEP10C-0C — HDF5 and Cell-ID reconciliation

Packages:

- `rhdf5`
- `data.table`

`rhdf5` is used to inspect HDF5 structures and read selected matrix segments.

`data.table` is used for metadata loading, cell-ID reconciliation, and tabular output.

The actual R and package versions should be documented using the generated session information files where available.

#### STEP10C-0D — GEO Raw UMI TXT File Characterization

The uploaded section uses base R functionality for:

- File and compressed-file inspection.
- GZIP reading.
- Text-line inspection.
- Cell-ID comparison.
- Metadata loading.

No additional external package dependency is explicitly identified in this section.

No random seed is required.

#### STEP10C-0D-2 — GEO Raw UMI Identifier Audit

The workflow primarily uses base R functions for:

- Reading GZIP-compressed text files.
- Extracting cell identifiers from the header.
- Scanning matrix-row identifiers.
- Comparing identifier sets.
- Exporting tabular audit outputs.

The optional H5 barcode cross-check uses `rhdf5`.

The H5 comparison is skipped if `rhdf5` is unavailable.

No random seed is required.

#### STEP10C-0E — Primary-Tumor Raw UMI Cell-ID and Compartment Reconciliation

Required packages:

- `data.table`
- `dplyr`
- `readr`
- `tidyr`

Main uses:

- `data.table`: metadata and STEP10B-1 table import.
- `dplyr`: cohort definition, grouping, counting, joins, and manifest construction.
- `readr`: tabular data import.
- `tidyr`: patient-by-compartment wide-table generation.

Base R is used for file handling, cell-ID reconciliation, integrity checks, and output generation.

The uploaded script contains an installation command for `tidyr`. Package installation should be managed consistently with the project environment and recorded package versions.

No random seed is required.

#### STEP10C — Patient-Level Raw UMI Pseudobulk Construction

Required package:

- `data.table`

Base R functions are used for:

- Reading compressed text files.
- Matching cell identifiers.
- Extracting raw UMI values.
- Aggregating counts by patient and compartment.
- Constructing pseudobulk matrices.
- Checking integer-valued raw counts.
- Generating output tables.

The complete GEO raw UMI matrix is not loaded into memory. The input file is processed sequentially, one gene row at a time.

The workflow uses deterministic identifier matching and aggregation and does not require a random seed.

### 5.5 STEP10D — Pseudobulk QC and Core-Gene Patient-Level Audit

The uploaded STEP10D sections primarily use base R functions.

No additional external package dependency is explicitly introduced in the reviewed STEP10D sections.

Main operations include:

- Reading tab-delimited pseudobulk matrices and annotations.
- Checking matrix dimensions and identifier consistency.
- Validating numeric and integer-valued raw UMI counts.
- Extracting AOC3, MS4A7, and IL7R counts.
- Calculating library sizes and detected-gene counts.
- Calculating descriptive summaries and quantiles.
- Checking patient and compartment completeness.
- Constructing patient-level detection profiles.
- Calculating patient-set intersections and overlaps.
- Generating descriptive concordance summaries.
- Writing tab-delimited audit outputs and methodological lock files.

Representative base R functions include:

- `read.delim()`
- `write.table()`
- `writeLines()`
- `file.path()`
- `dir.create()`
- `file.exists()`
- `dir.exists()`
- `list.files()`
- `exists()`
- `data.frame()`
- `aggregate()`
- `reshape()`
- `quantile()`
- `median()`
- `mean()`
- `all()`
- `setequal()`
- `anyDuplicated()`
- `order()`
- `match()`
- `ifelse()`
- `unique()`
- `table()`
- `duplicated()`
- `Reduce()`
- `intersect()`
- `combn()`
- `lapply()`
- `do.call()`

The STEP10D workflows are deterministic and do not require a random seed.

### 5.6 STEP10E-0 — Integrated Input and Provenance Lock

The uploaded section uses base R functions for:

- Project-relative path construction.
- Directory and file existence checks.
- File discovery.
- Tab-delimited file reading.
- Schema and identifier validation.
- Numeric integrity checks.
- Dataset and gene-set reconciliation.
- Current-session object availability checks.

Representative functions include:

- `file.path()`
- `list.files()`
- `file.exists()`
- `dir.exists()`
- `read.delim()`
- `exists()`

No additional package dependency is explicitly introduced.

The project setup must define `project_root` before this section is executed.

The displayed checks are deterministic and do not require a random seed.

### 5.7 STEP10E — Patient-Level Pseudobulk Expression Characterization

The uploaded STEP10E code uses base R functions for pseudobulk matrix processing, descriptive calculations, data-frame construction, and output generation.

Representative functions include:

- `read.delim()`
- `file.path()`
- `file.exists()`
- `dir.exists()`
- `dir.create()`
- `as.matrix()`
- `as.numeric()`
- `colSums()`
- `sweep()`
- `log2()`
- `quantile()`
- `mean()`
- `IQR()`
- `match()`
- `do.call()`
- `rbind()`
- `write.table()`
- `writeLines()`

No additional external R package is explicitly called in the uploaded STEP10E section.

The workflow calculates CPM and log2(CPM + 1) from raw UMI pseudobulk counts and produces descriptive expression summaries.

No random-number generation is used. A random seed is not required.

### 5.8 STEP10F-1 — Patient-Level Expression Profile

The uploaded STEP10F-1 sections use base R functionality for:

- Reading STEP10E input tables.
- Validating input schemas and data structures.
- Checking gene and cohort labels.
- Validating numeric and integer-valued raw UMI counts.
- Auditing CPM and log2(CPM + 1) values.
- Calculating descriptive expression summaries.
- Identifying complete patient sets across compartments.
- Constructing patient-level expression profiles.

Representative functions include:

- `read.delim()`
- `file.path()`
- `file.exists()`
- `dir.exists()`
- `dir.create()`
- `colnames()`
- `setdiff()`
- `unique()`
- `setequal()`
- `factor()`
- `table()`
- `paste()`
- `duplicated()`
- `as.numeric()`
- `is.finite()`
- `floor()`
- `abs()`
- `max()`
- `mean()`
- `median()`
- `quantile()`
- `min()`
- `IQR()`
- `Reduce()`
- `intersect()`
- `do.call()`
- `rbind()`
- `write.table()`

No additional external R package dependency is explicitly introduced by the uploaded STEP10F-1 sections.

The workflow uses locked STEP10E outputs and does not reread the GEO raw UMI matrix or H5 files.

The uploaded script header specifies R version 4.5.1. No random seed is required because the procedures are deterministic.

### 5.9 STEP10F-3 — Expression-Context Characterization

The uploaded STEP10F-3 sections use base R functions for:

- Patient-level descriptive summaries.
- Ranking and quantitative context comparisons.
- Data-frame splitting and merging.
- Quality-control checks.
- Final tabular output.
- Methodological lock-file generation.

Representative functions include:

- `split()`
- `interaction()`
- `lapply()`
- `do.call()`
- `merge()`
- `rank()`
- `quantile()`
- `median()`
- `mean()`
- `IQR()`
- `write.table()`
- `writeLines()`

No additional package dependency is explicitly introduced by the uploaded STEP10F-3 sections.

The uploaded environment documentation specifies R version 4.5.1. Complete session information should be obtained from the actual execution environment.

### 5.10 STEP10G — Single-Cell Evidence Synthesis

The uploaded STEP10G section uses base R functions for:

- File existence checks.
- Tab-delimited input and output.
- Data-frame subsetting and manipulation.
- Patient-level structural summaries.
- Ranking and descriptive calculations.
- Evidence-table generation.
- Quality-control checks.
- Final lock-file generation.

No additional R package is explicitly loaded in the uploaded STEP10G section.

STEP10G requires the locked output files generated by:

- STEP10F-1
- STEP10F-2
- STEP10F-3

The script checks the existence of 14 expected input files before reading them.

Project-relative paths are based on the `input_dir` variable established by the project setup script.

The uploaded STEP10G section does not independently record a package-version snapshot or a separate R version.

No random seed is required for the deterministic operations described in this section.

### 5.11 STEP10H — Final Manuscript Supplementary Export

Required and explicitly used packages:

- `data.table`
- `ggplot2`
- `openxlsx`
- `here`
- `scales`

The script explicitly checks:

- `data.table`
- `ggplot2`
- `openxlsx`

The workflow also calls:

- `here::here()`
- `scales::hue_pal()`

These packages should all be included in the consolidated dependency inventory.

Base R functions are used for:

- File and directory management.
- Structural validation.
- Data-frame construction.
- Patient-count aggregation.
- Figure output validation.

The project root is resolved using `here::here()`.

The STEP10H input and output directories remain under the historical STEP10 data directory:

`data/Paper 1_Step 10-Single-Cell/`

No random seed is specified because this workflow does not introduce a random-number generation procedure.

The exact R and package versions should be recorded from the actual execution environment.

---

## 6. Base R Functionality

In addition to the external packages listed above, the workflows rely extensively on base R and standard R functionality.

The main functional categories are summarized below.

| Functional category | Representative functions |
|---|---|
| File and directory management | `file.path()`, `file.exists()`, `dir.exists()`, `dir.create()`, `list.files()` |
| Input and output | `read.delim()`, `readLines()`, `write.table()`, `writeLines()`, `saveRDS()` |
| Data-frame construction and manipulation | `data.frame()`, `merge()`, `rbind()`, `do.call()`, `lapply()` |
| Data reshaping and aggregation | `reshape()`, `expand.grid()`, `aggregate()`, `split()` |
| Matrix operations | `as.matrix()`, `colSums()`, `sweep()`, `mapply()` |
| Numeric validation | `as.numeric()`, `is.finite()`, `floor()`, `abs()`, `max()` |
| Descriptive statistics | `mean()`, `median()`, `quantile()`, `min()`, `max()`, `IQR()` |
| Patient and identifier reconciliation | `unique()`, `setequal()`, `intersect()`, `Reduce()`, `match()`, `duplicated()` |
| Group and structure checks | `table()`, `all()`, `anyDuplicated()`, `identical()`, `setdiff()` |
| Ranking and concordance | `rank()`, `combn()`, `interaction()`, `order()` |
| Reproducibility records | `sessionInfo()` |

The list is representative rather than exhaustive. It summarizes the functions identified in the uploaded workflow sections and does not imply that every function is used in every analysis step.

---

## 7. Random Seed and Stochastic Procedures

### Script 1

The random seed documented in `GitHub7_Paper1_Single_Cell_Analysis_1.R` is:

`20250724`

This seed should be retained for the corresponding Script 1 procedures that use random operations.

### Script 2

The uploaded STEP10A–STEP10H sections primarily perform deterministic data processing, matrix inspection, identifier reconciliation, pseudobulk aggregation, descriptive summaries, quality-control checks, and supplementary output generation.

No additional random seed is specified for these uploaded Script 2 sections.

A random seed is not required for the deterministic operations described in these sections.

### Reproducibility note

The Script 1 seed should not be automatically assigned to every Script 2 step. Random seeds should be recorded only for workflows that actually use stochastic procedures.

---

## 8. Version Documentation and Session Information

### Package versions

The actual versions of the R packages have not been fully established from the uploaded README fragments.

Package versions should be recorded from the execution environment rather than inferred from:

- `library()` calls.
- `require()` calls.
- Package names in the script.
- Installation commands.
- Script headers.

The following information should be retained for a reproducible environment record:

- R version.
- R platform and operating system.
- Package versions.
- Bioconductor version, where applicable.
- Locale and relevant system configuration.
- Session information for the executed workflow.

### Session information files

The uploaded Results README fragments identify session information outputs for several audit stages, including:

- `STEP10A_sessionInfo.txt`
- `STEP10B-1_sessionInfo.txt`
- `STEP10B-2_sessionInfo.txt`
- `STEP10C-0A_sessionInfo.txt`
- `STEP10C-0B_sessionInfo.txt`
- `STEP10C-0C_sessionInfo.txt`

These files should be used to document the actual R and package versions associated with their respective executions.

Session information should be interpreted in the context of the specific workflow execution that generated it. A session record from one step should not automatically be assumed to represent every later execution if the environment changed.

A consolidated environment snapshot may be generated after the complete analysis workflow has been executed in the intended project environment.

---

## 9. System Requirements and External Dependencies

The uploaded environment documentation does not yet establish a complete operating-system or system-level dependency inventory.

The following details remain to be confirmed from the actual execution environment:

- Windows version and system architecture.
- R installation and platform details.
- Font family and font availability for publication figures.
- HDF5-related system configuration associated with `rhdf5`.
- Any additional external software or system libraries required by the installed R packages.

No additional standalone external application is explicitly confirmed as a required dependency by the uploaded script sections.

The final system requirements should be updated only when supported by the actual environment records.

---

## 10. Reproducibility and Verification Notes

### Project structure

The workflows use project-relative paths established through the shared setup configuration.

The project root and relevant input and output directories must be initialized before running downstream steps.

### Input dependencies

The single-cell workflows depend on the supplied TISCH2 metadata and the corresponding expression/count matrices.

Depending on the analysis stage, supported inputs include:

- RDS files.
- Matrix Market (`.mtx`) files.
- HDF5 (`.h5`) files.
- GZIP-compressed raw UMI text files.
- Tab-delimited metadata and intermediate result tables.

Input formats and required files should be checked against the corresponding workflow-specific documentation.

### Package availability

Package dependencies should be installed and versioned in a consistent project environment.

Conditional dependencies, such as `rhdf5` for optional H5 cross-checks, should be documented separately from packages required for the main analysis path.

### Execution and verification status

This README consolidates the software requirements documented in the uploaded script sections.

It does not claim that:

- All scripts have been executed in one identical software environment.
- All package versions have been independently verified.
- The operating system and font configuration have been fully audited.
- All conditional dependencies were available during execution.
- All generated outputs have been independently validated.

The environment record should be finalized using the actual session information and system configuration after successful execution of the complete workflows.

---

## 11. Summary of Environment Information

| Item | Current documentation status |
|---|---|
| Programming language | R |
| R version | 4.5.1 specified in Script 2 header; execution-time confirmation pending |
| Project path management | `here` and shared project setup variables |
| Core single-cell packages | `Seurat`, `Matrix`, `ggplot2`, `ggrepel`, `dplyr`, `data.table`, `tidyverse`, `scales` |
| Additional STEP-specific packages | `readr`, `tidyr`, `rhdf5`, `zip`, `openxlsx` |
| Base R | Used extensively throughout both scripts |
| Script 1 random seed | 20250724 |
| Script 2 random seed | Not required for the uploaded deterministic sections |
| Package versions | To be confirmed from execution-time session information |
| Operating system | To be confirmed |
| Font configuration | To be confirmed |
| External system dependencies | To be confirmed where applicable |
| Full environment verification | Pending execution-time review |

This README provides the consolidated software environment documentation for the Paper 1 single-cell analysis workflows. It should be updated when the final execution-time environment records become available.