# Single-Cell Analysis Scripts

## 1. Overview

This directory contains the R scripts used for the single-cell transcriptomic characterization of AOC3, MS4A7, and IL7R in lung adenocarcinoma (LUAD).

The analysis uses the TISCH2 NSCLC_GSE131907 dataset and complementary GEO raw UMI count data. The two scripts provide complementary components of the single-cell workflow:

- `GitHub7_Paper1_Single_Cell_Analysis_1.R` performs cell-level expression visualization and differential-expression dot-plot generation.
- `GitHub7_Paper1_Single_Cell_Analysis_2.R` performs patient- and compartment-level data audits, raw-UMI pseudobulk construction, core-gene expression characterization, evidence integration, and manuscript-oriented supplementary exports.

The overall workflow is organized around predefined patient eligibility criteria, cellular compartment definitions, locked pseudobulk inputs, and descriptive expression analyses.

The three prioritized study genes are:

- **AOC3**
- **MS4A7**
- **IL7R**

The analyses characterize their expression patterns across T-cell, myeloid/macrophage, and endothelial compartments. The results are intended to provide single-cell and patient-level expression context for the bulk transcriptomic findings.

This module does not establish causal relationships, exclusive cell-type localization, or cell-cell communication mechanisms.

------------------------------------------------------------------------

## 2. Scripts and Workflow Organization

### 2.1 Script 1: Cell-Level Single-Cell Visualization

**Script:** `GitHub7_Paper1_Single_Cell_Analysis_1.R`

This script focuses on cell-level visualization and expression characterization using the TISCH2 single-cell dataset.

The main components include:

1.  Loading and inspecting the single-cell dataset.
2.  Preparing the cell-level metadata and expression objects.
3.  Visualizing the major cell populations using UMAP.
4.  Displaying the expression distributions of AOC3, MS4A7, and IL7R.
5.  Exporting publication-oriented expression plots.
6.  Generating differential-expression dot plots based on the supplied dataset annotations and expression summaries.

The principal objects documented in the uploaded script sections include:

- `sce`: the main single-cell object.
- `meta`: the cell-level metadata object.
- `umap_df`: the UMAP plotting data frame.

The visualization workflow covers 12 major cell types.

The documented plotting components include:

- Figure 9A–9B: UMAP-based cell-population and target-gene expression visualization.
- Figure 9C: differential-expression dot plot.

The plotting workflow includes PDF, TIFF, PNG, and SVG export, together with publication-quality figure checks.

**Random seed:** `20250724`

The exact input and output filenames for all Script 1 sections should be taken from the final R script and its corresponding output records. They are not inferred from the Script 2 workflow.

### 2.2 Script 2: Patient-Level and Compartment-Level Expression Characterization

**Script:** `GitHub7_Paper1_Single_Cell_Analysis_2.R`

This script implements the STEP10 workflow, extending the cell-level analysis to patient-level and compartment-level characterization.

Its main components are:

1.  Patient and sample structure auditing.
2.  Primary-tumor compartment cell-count feasibility assessment.
3.  Patient-level pseudobulk input reconciliation.
4.  Raw-UMI pseudobulk construction.
5.  Pseudobulk quality control and core-gene detection auditing.
6.  Patient-level compartment concordance analysis.
7.  Patient-level pseudobulk expression characterization.
8.  Patient-level gene-by-compartment expression profiling.
9.  Expression-context and compartment-ranking analysis.
10. Integrated single-cell evidence synthesis.
11. Final manuscript supplementary table and figure export.

The workflow uses the outputs of preceding steps as locked inputs for subsequent analyses. It does not repeatedly reconstruct upstream data when a downstream step is intended only for validation or descriptive characterization.

------------------------------------------------------------------------

## 3. Data Sources and Input Files

### 3.1 TISCH2 Dataset

The principal single-cell dataset is:

**TISCH2 dataset:** `NSCLC_GSE131907`

The dataset overview reported in the original analysis documentation is:

- 203,298 cells
- 44 patients
- 58 samples

The principal TISCH2 input files are:

| Input file | Description |
|------------------------------------|------------------------------------|
| `NSCLC_GSE131907_CellMetainfo_table.tsv` | Cell-level metadata and annotation information |
| `NSCLC_GSE131907_expression.h5` | HDF5-format single-cell expression data |
| `NSCLC_GSE131907_AllDiffGenes_table.tsv` | Dataset-provided differential-expression information |
| `NSCLC_GSE131907_Expression.zip` | Additional expression-data archive inspected during input characterization |

The TISCH2-provided annotations and differential-expression information are treated as database-provided resources. The analysis does not imply that all database-level differential-expression results were independently recalculated.

### 3.2 GEO Raw UMI Count Matrix

The GEO raw UMI count matrix used for patient-level pseudobulk construction is:

`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`

This matrix is used as the source of raw UMI counts for patient-by-compartment aggregation.

The workflow performs file-format, cell-ID, gene-ID, and metadata reconciliation before using the matrix for pseudobulk construction.

The raw UMI matrix is processed in a memory-conscious manner. The documented workflow streams the matrix by gene rather than loading the complete expression matrix into memory.

The normalized expression data and database-provided differential-expression tables are not substitutes for this raw-count input.

------------------------------------------------------------------------

## 4. Patient Selection and Cellular Compartment Definitions

### 4.1 Primary-Tumor Cohort

The principal patient-level analysis uses cells satisfying both criteria:

- `Stage = Primary`
- `Source = Tumor`

The documented cohort contains:

- 65,702 cells
- 9 patients

The cohort is established from the TISCH2 cell-level metadata and reconciled with the raw UMI matrix before pseudobulk construction.

### 4.2 Cellular Compartment Definitions

The three primary cellular compartments are defined as follows:

| Compartment        | Included TISCH2 cell annotations |
|--------------------|----------------------------------|
| T-cell             | CD4Tconv, CD8T, CD8Tex           |
| Myeloid/macrophage | Mono/Macro                       |
| Endothelial        | Endothelial                      |

Dendritic cells are not included in the myeloid/macrophage compartment.

These definitions are used consistently in the patient-level cell-count audits, pseudobulk construction, expression characterization, and downstream evidence integration.

### 4.3 Locked Patient and Group Structure

The primary and sensitivity analysis groups are:

| Analysis set | Compartment        | Expected patient/group count |
|--------------|--------------------|-----------------------------:|
| Primary      | T-cell             |                            9 |
| Primary      | Myeloid/macrophage |                            9 |
| Primary      | Endothelial        |                            4 |
| Sensitivity  | Endothelial        |                            2 |

The endothelial groups are handled separately because endothelial cell availability differs across patients.

The primary endothelial analysis includes patients meeting the predefined cell-count feasibility criterion. Patients with lower endothelial cell numbers are retained in a separate sensitivity group where applicable.

The documented cell-depth categories are:

- At least 50 cells: primary endothelial eligibility.
- 20–49 cells: endothelial sensitivity eligibility.
- Fewer than 20 cells: not included in formal endothelial pseudobulk analysis.

The T-cell and myeloid/macrophage compartments are available for all nine primary-tumor patients.

The sensitivity endothelial group is not merged into the primary endothelial summaries.

### 4.4 Patient-Level Replication

The patient, rather than the individual cell, is the biological replicate for the patient-level pseudobulk and descriptive expression analyses.

Raw UMI counts are summed across cells belonging to the same patient and compartment.

Zero-count observations are retained. Patients are not excluded solely because a core gene is undetected.

------------------------------------------------------------------------

## 5. Workflow Details: Script 1

### 5.1 Cell-Level Dataset Preparation

The script prepares the main single-cell object and associated metadata for cell-population and gene-expression visualization.

The documented objects include:

- `sce`
- `meta`
- `umap_df`

The analysis uses the available dataset annotations and UMAP information to display the major cell populations.

### 5.2 Major Cell-Population UMAP

The UMAP visualization workflow includes:

1.  Preparing the plotting data.
2.  Defining the cell-population labels.
3.  Generating the major cell-type UMAP.
4.  Exporting the figure in publication-oriented formats.
5.  Performing figure-output quality checks.

The documented analysis covers 12 major cell types.

### 5.3 Core-Gene Expression Visualization

The target genes are:

- AOC3
- MS4A7
- IL7R

The workflow constructs the expression plotting data and generates target-gene expression UMAP panels.

The documented figure components include separate expression visualizations for the three genes.

### 5.4 Differential-Expression Dot Plot

The script also generates a differential-expression dot plot using the supplied dataset-level information.

The dot plot should be interpreted in the context of the underlying TISCH2 annotations and database-provided results.

It should not be interpreted as an independently performed patient-level differential-expression analysis unless supported by a separate statistical workflow.

### 5.5 Figure Export

The documented figure-export formats include:

- PDF
- TIFF
- PNG
- SVG

The script includes publication-quality output checks.

------------------------------------------------------------------------

## 6. Workflow Details: Script 2

### 6.1 STEP10A — Patient and Sample Structure Audit

**Objective:** Characterize and audit the patient-, sample-, source-, stage-, and cell-level structure of the TISCH2 NSCLC_GSE131907 dataset.

**Main procedures:**

1.  Load the TISCH2 cell-level metadata.
2.  Check required metadata columns.
3.  Audit duplicate cell IDs and missing annotations.
4.  Summarize total cells, patients, and samples.
5.  Generate patient-level sample-structure tables.
6.  Summarize cell counts per patient.
7.  Summarize cell populations per patient.
8.  Examine source and stage distributions.
9.  Cross-tabulate stage and source.
10. Classify primary-tumor and other metadata categories.
11. Document potential exclusion categories.
12. Audit the availability of AOC3, MS4A7, and IL7R.
13. Save the audit object and session information.

**Inputs:**

- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- `NSCLC_GSE131907_AllDiffGenes_table.tsv`

**Output directory:**

`results/STEP10A_Patient_Sample_Structure_Audit/`

This step is an input and metadata audit. It does not recalculate PCA, UMAP, clustering, or cell-type annotation.

### 6.2 STEP10B — Primary-Tumor Compartment Feasibility Audits

#### STEP10B-1 — Compartment Cell-Count Feasibility Audit

**Objective:** Establish the primary-tumor cell cohort and evaluate the availability of the predefined cellular compartments.

**Input:**

- `NSCLC_GSE131907_CellMetainfo_table.tsv`

**Main procedures:**

1.  Validate metadata integrity.
2.  Select the Primary and Tumor subset.
3.  Verify the expected cohort size and patient count.
4.  Audit the predefined cellular lineages.
5.  Validate cell-type-to-compartment mapping.
6.  Assign cells to the three predefined compartments.
7.  Generate patient-by-compartment cell-count tables.
8.  Examine cell-count distributions.
9.  Evaluate three-compartment representation.
10. Reconcile assigned and unmatched cells.
11. Save audit results and session information.

**Expected cohort:**

- 9 patients
- 65,702 cells

**Output directory:**

`results/STEP10B-1_Compartment_CellCount_Feasibility_Audit/`

This step establishes compartment cell counts but does not apply the minimum-cell threshold, exclude patients, analyze gene expression, or construct pseudobulk matrices.

#### STEP10B-2 — Pseudobulk Feasibility and Cell-Depth Audit

**Input:**

`STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv`

from:

`results/STEP10B-1_Compartment_CellCount_Feasibility_Audit/`

**Main procedures:**

1.  Validate the nine-patient cohort.
2.  Reconcile cell counts with STEP10B-1.
3.  Convert the patient-by-compartment table from wide to long format.
4.  Classify cell depth into `<20`, `20–49`, and `>=50` categories.
5.  Summarize cell-depth distributions.
6.  Evaluate patient-level and compartment-specific pseudobulk feasibility.
7.  Assess endothelial-specific and cross-compartment feasibility.
8.  Save audit results and session information.

**Output directory:**

`results/STEP10B-2_Pseudobulk_Feasibility_CellDepth_Audit/`

This step is a feasibility audit only. It does not perform patient exclusion, expression analysis, pseudobulk construction, or differential-expression analysis.

### 6.3 STEP10C — Patient-Level Compartment-Specific Pseudobulk Construction

**Objective:** Generate patient-level raw-UMI pseudobulk matrices for the predefined cellular compartments.

**Principal inputs:**

- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `STEP10B-1D_Primary_Tumor_CellLevel_Compartment_Mapping.tsv`
- `STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv`
- `STEP10B-2B_Patient_Compartment_CellDepth_Feasibility.tsv`

Additional intermediate input and reconciliation files include:

- `STEP10C-0E_G_PrimaryTumor_RawUMI_Pseudobulk_Input_Manifest.tsv`
- `STEP10C-0E_C_Patient_Compartment_Counts_Wide.tsv`

**Main procedures:**

1.  Define the locked patient cohort and compartment mapping.
2.  Validate metadata, patient IDs, and cell-level compartment assignments.
3.  Inspect and validate the original GEO raw UMI matrix.
4.  Reconcile metadata cell IDs with raw-matrix column names.
5.  Verify the presence of AOC3, MS4A7, and IL7R.
6.  Audit gene identifiers and duplicate gene symbols.
7.  Aggregate raw UMI counts by patient and compartment.
8.  Construct primary pseudobulk matrices for T-cell, myeloid/macrophage, and endothelial compartments.
9.  Construct the separate endothelial sensitivity pseudobulk matrix.
10. Export core-gene raw-count tables.
11. Calculate pseudobulk library sizes.
12. Audit patient-level core-gene detection.
13. Summarize core-gene detection by compartment and analysis set.
14. Reconcile expected and observed pseudobulk patient numbers.
15. Perform final quality-control checks and save the analysis objects.

The raw UMI matrix is processed gene by gene to avoid loading the entire expression matrix into memory.

**Output directory:**

`results/STEP10C_Patient_Level_Pseudobulk/`

The resulting matrices are based on summed raw UMI counts. Patient-level pseudobulk observations are used as the biological replicates for downstream descriptive characterization.

This workflow does not perform differential-expression analysis, correlation analysis, or inferential statistical testing.

#### STEP10C-0 — Expression-Matrix and Cell-ID Characterization Audits

The STEP10C-0 sections provide preparatory audits of the available expression resources and their compatibility with the TISCH2 metadata.

**STEP10C-0A and STEP10C-0B — H5 and ZIP File Inspection**

Inputs:

- `NSCLC_GSE131907_expression.h5`
- `NSCLC_GSE131907_Expression.zip`

The preliminary inspection includes file-type and internal-structure checks.

The H5 characterization includes:

- HDF5 structure inspection.
- Verification of the expected 10x-style matrix components.
- Feature and barcode dimensions.
- Selected stored expression values.
- Selected cell-level sparse entries.
- Cell-level expression-depth inspection.
- Preliminary raw-count compatibility assessment.

These are preparatory audits. They do not construct pseudobulk matrices, normalize expression, or perform statistical testing.

**STEP10C-0C — H5 Cell-ID Reconciliation Audit**

Inputs:

- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- `NSCLC_GSE131907_expression.h5`

The workflow:

1.  Audits the metadata.
2.  Extracts H5 barcodes and features.
3.  Performs global exact cell-ID matching.
4.  Identifies unmatched IDs and barcode-format differences.
5.  Reconciles the Primary and Tumor subset.
6.  Verifies the nine-patient cohort.
7.  Summarizes matching by compartment and patient.
8.  Compares compartment-level matching with STEP10B-1.
9.  Generates reconciliation summaries and audit objects.
10. Closes HDF5 handles.

This step does not perform pseudobulk construction or statistical testing.

**STEP10C-0D — GEO Raw UMI Matrix Characterization**

Input:

`GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`

The preliminary characterization includes:

- File existence and compressed-size checks.
- GZIP header inspection.
- Preliminary decompression.
- Inspection of the first five lines.
- Delimiter and header identification.
- Matrix-orientation assessment.
- Candidate cell-ID extraction.
- Comparison of candidate IDs with TISCH2 metadata.

This is a format and structure inspection, not a complete raw-count validation.

**STEP10C-0D-2 — GEO Raw UMI Cell-ID and Gene-ID Audit**

Inputs:

- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- Optional: `NSCLC_GSE131907_expression.h5`

The workflow:

1.  Verifies input-file availability.
2.  Extracts raw UMI cell IDs from the matrix header.
3.  Checks cell-ID uniqueness.
4.  Compares GEO cell IDs with TISCH2 metadata.
5.  Summarizes shared, GEO-only, and TISCH2-only cells.
6.  Optionally compares GEO IDs with H5 barcodes.
7.  Scans gene identifiers from the first field of each row.
8.  Audits duplicate gene IDs.
9.  Checks the availability of AOC3, MS4A7, and IL7R.
10. Exports lightweight audit tables and RDS objects.

The complete expression matrix is not loaded. No pseudobulk construction or statistical testing is performed.

**STEP10C-0E — Primary-Tumor Raw UMI Cell-ID and Compartment Reconciliation**

Inputs:

- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`
- `STEP10B-1F_Patient_Compartment_Cell_Counts_Wide.tsv`

The workflow:

1.  Loads and validates metadata.
2.  Checks cell-ID uniqueness.
3.  Reconstructs the predefined compartment mapping.
4.  Defines the Primary and Tumor cohort.
5.  Summarizes patient-by-compartment cell counts.
6.  Reconciles counts with STEP10B-1.
7.  Extracts GEO cell IDs from the raw-matrix header.
8.  Reconciles global cell-ID sets.
9.  Verifies Primary and Tumor cell coverage.
10. Evaluates patient- and compartment-level availability.
11. Generates the formal cell-level pseudobulk input manifest.
12. Exports reconciliation tables and the final audit summary.

The complete raw UMI matrix is not read for expression aggregation in this preparatory step.

### 6.4 STEP10D — Pseudobulk Quality Control and Core-Gene Audits

#### STEP10D-A — Pseudobulk QC and Core-Gene Patient-Level Audit

This workflow uses the locked STEP10C pseudobulk matrices and annotation files.

**Inputs:**

- `STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv`
- `STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv`
- `STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv`
- `STEP10C-D_Pseudobulk_Endothelial_Sensitivity_rawUMI.tsv`
- `STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv`
- `STEP10C-F_Sensitivity_Pseudobulk_Sample_Annotation.tsv`
- `STEP10C-G_Pseudobulk_CellCount_Audit.tsv`
- `STEP10C-H_Pseudobulk_Input_Manifest.tsv`
- `STEP10C-I_Pseudobulk_QC_Summary.tsv`

**Main procedures:**

1.  Verify the availability of STEP10C outputs.
2.  Read and validate the four pseudobulk matrices.
3.  Check locked cohort dimensions.
4.  Verify AOC3, MS4A7, and IL7R.
5.  Summarize sample-level library sizes and detected-gene counts.
6.  Extract patient-level raw UMI counts for the core genes.
7.  Generate core-gene detection and count summaries.
8.  Calculate descriptive relative abundance as core-gene raw UMI counts divided by total pseudobulk UMI counts.
9.  Validate the sample annotation files.
10. Reconcile matrix `Group_ID` values with annotation `Group_ID` values.
11. Verify locked compartment labels and expected group counts.
12. Check patient-by-compartment completeness.
13. Reconcile matrix dimensions, gene identity, and gene order.
14. Audit missing, non-finite, negative, and non-integer values.
15. Summarize library size, detected genes, zero-count genes, and detection fractions.
16. Generate a descriptive low-library flag based on the Q1 − 1.5 × IQR rule.
17. Export QC tables, summaries, and audit records.

The low-library flag is descriptive and does not automatically exclude samples.

The locked group counts are:

- Primary T-cell: 9.
- Primary Myeloid/macrophage: 9.
- Primary Endothelial: 4.
- Sensitivity Endothelial: 2.

This workflow audits existing pseudobulk matrices and does not recalculate them or perform inferential statistical testing.

#### STEP10D-B — Core-Gene Raw-Count and Detection Audit

This workflow uses the four locked STEP10C pseudobulk matrices and their corresponding Primary and Sensitivity annotation files.

**Main procedures:**

1.  Locate the locked STEP10C input files.
2.  Read the four pseudobulk matrices.
3.  Verify matrix dimensions and identifier integrity.
4.  Check raw UMI values for missing, non-finite, negative, and non-integer values.
5.  Reconcile gene identity and order across matrices.
6.  Verify the presence of AOC3, MS4A7, and IL7R.
7.  Read and validate the STEP10C annotation files.
8.  Reconcile matrix `Group_ID` values with annotation `Group_ID` values.
9.  Extract the three core genes' raw UMI counts.
10. Generate detection indicators using raw UMI count \> 0.
11. Generate descriptive detection and raw-count summaries.
12. Verify the expected group counts.
13. Export raw-count, detection, summary, and audit files.
14. Generate a final audit summary and methodological lock file.

Zero counts are retained and do not trigger sample exclusion.

#### STEP10D-C — Core-Gene Compartment-Level Expression Audit

This workflow uses the locked STEP10D-B core-gene raw-count table.

**Main procedures:**

1.  Verify required metadata and core-gene columns.
2.  Check `Group_ID` and `Dataset` integrity.
3.  Validate core-gene raw UMI count values.
4.  Verify the four locked dataset labels.
5.  Reconcile observed group counts with the predefined group-count structure.
6.  Convert the core-gene count table to long format.
7.  Inspect Primary and Sensitivity records before filtering.
8.  Retain Primary-only records for the subsequent STEP10D-D workflow.
9.  Calculate descriptive raw-count summaries by compartment and gene.
10. Generate patient-level wide-format core-gene data.
11. Calculate detection-frequency summaries.
12. Verify expected patient counts and detection-frequency arithmetic.

Primary and Sensitivity Endothelial records remain separate.

This workflow does not perform inferential statistical testing or recalculate upstream pseudobulk matrices.

#### STEP10D-D — Patient-Level Compartment Concordance Audit and Final Lock

This workflow characterizes patient-level concordance across the three Primary compartments.

**Main procedures:**

1.  Verify the required analysis objects and columns.
2.  Confirm Primary-only observations and locked compartment and gene labels.
3.  Check uniqueness of patient × compartment × gene observations.
4.  Verify the detection definition (`Raw_Count > 0`).
5.  Generate patient-level compartment availability.
6.  Calculate pairwise patient overlap.
7.  Identify the complete three-compartment patient set from the observed patient intersection.
8.  Construct patient × compartment × core-gene profiles.
9.  Generate patient-level detection-pattern summaries.
10. Summarize three-compartment detection concordance.
11. Summarize pairwise core-gene detection concordance.
12. Generate the full patient-level concordance profile.
13. Add detection-status labels for each patient, compartment, and core-gene observation.
14. Order records using the predefined patient, compartment, and gene ordering.
15. Verify expected patient counts and profile completeness.
16. Generate structural QC tables and apply the overall PASS/STOP condition.
17. Export patient-level detection and descriptive tables.
18. Export the complete three-compartment patient list.
19. Generate the final audit summary and methodological lock file.
20. Verify expected output files and inspect existing workflow outputs.

Zero-count observations are retained. Patients are not excluded based on detection status.

The complete three-compartment patient set is derived from the observed intersection of the three locked Primary compartment-specific patient sets.

**Documentation note:** The uploaded script documentation identifies two versions of the STEP10D-D final output/provenance-lock section. Their observation counts and output-directory definitions require reconciliation before the final script is considered locked.

### 6.5 STEP10E-0 — Integrated Input and Provenance Lock

**Objective:** Verify the integrity and provenance of the locked STEP10C and STEP10D inputs before downstream expression characterization.

**Main procedures:**

1.  Define the STEP10 working directory using the project root.
2.  Define locked core genes, compartments, and expected patient counts.
3.  Define a helper function for locating exactly one file.
4.  Verify the STEP10C and STEP10D workflow directories.
5.  Locate the four locked STEP10C pseudobulk matrices.
6.  Validate the STEP10C Primary annotation schema.
7.  Verify Primary eligibility and compartment counts.
8.  Read and validate the STEP10D-B core-gene raw-count table.
9.  Verify required metadata and core-gene columns.
10. Check raw-count finiteness, non-negativity, and integer-valued structure.
11. Read and validate the STEP10D-C compartment-level descriptive summary.
12. Verify its schema, dataset labels, core-gene identity, Dataset × Gene structure, and numeric integrity.
13. Inspect the STEP10D-D output directory and related files.
14. Audit the availability and dimensions of selected STEP10D-D objects in the current R session.

This workflow provides descriptive concordance and input/provenance checks. It does not recalculate pseudobulk matrices, normalize expression, or perform inferential statistical testing.

### 6.6 STEP10E — Patient-Level Pseudobulk Expression Characterization

**Objective:** Characterize AOC3, MS4A7, and IL7R expression using the locked patient-level raw-UMI pseudobulk matrices.

**Main inputs:**

- The four STEP10C pseudobulk matrices.
- The Primary and Sensitivity pseudobulk sample annotation files.

**Main procedures:**

1.  Define project-relative STEP10C input and STEP10E output directories.
2.  Check the availability of the four pseudobulk matrices and two annotation files.
3.  Define the three core genes.
4.  Read the Primary and Sensitivity pseudobulk matrices.
5.  Identify the gene identifier column in each matrix.
6.  Validate gene identity and ordering across matrices.
7.  Convert expression columns into numeric matrices.
8.  Audit raw UMI values for missing, non-finite, negative, and non-integer values.
9.  Calculate total UMI library sizes for each pseudobulk group.
10. Derive CPM values from raw UMI counts.
11. Derive log2(CPM + 1) values.
12. Match pseudobulk `Group_ID` values with sample annotations.
13. Construct patient-level expression records for each core gene and compartment.
14. Combine Primary and Sensitivity records while retaining their separate eligibility labels.
15. Generate Primary-only descriptive summaries.
16. Generate a separate Sensitivity Endothelial summary.
17. Validate expected observation, patient, and group counts.
18. Export patient-level expression tables and descriptive summaries.
19. Generate the methodological lock and final audit records.

**Locked analysis structure:**

| Analysis set | Compartment        | Groups |
|--------------|--------------------|-------:|
| Primary      | T-cell             |      9 |
| Primary      | Myeloid/macrophage |      9 |
| Primary      | Endothelial        |      4 |
| Sensitivity  | Endothelial        |      2 |

**Expression measures:**

- Raw UMI count.
- Total pseudobulk library UMI.
- CPM.
- Log2(CPM + 1).
- Detection status (`Raw UMI > 0`).

This step performs descriptive expression characterization only. It does not perform differential-expression analysis, inferential statistical testing, correlation analysis, CellChat analysis, ligand–receptor inference, or new pseudobulk aggregation.

Sensitivity Endothelial observations remain separate from Primary-only summaries.

### 6.7 STEP10F-1 — Patient-Level Core-Gene Expression Profile

**Objective:** Generate patient-level descriptive expression profiles for AOC3, MS4A7, and IL7R across the Primary cellular compartments while retaining Sensitivity Endothelial observations as a separate analysis set.

**Inputs:**

- Locked STEP10E Primary expression table.
- Locked STEP10E Sensitivity expression table.

**Main procedures:**

1.  Resolve STEP10E input and STEP10F-1 output directories using project-relative paths.
2.  Verify the presence of the Primary and Sensitivity STEP10E expression tables.
3.  Read and validate the required input columns.
4.  Verify the exact core-gene set: AOC3, MS4A7, and IL7R.
5.  Establish the predefined core-gene and compartment ordering.
6.  Audit Primary patient and compartment structure.
7.  Audit Sensitivity Endothelial patient and group structure.
8.  Verify Primary and Sensitivity eligibility labels.
9.  Verify the expected observation counts: 66 Primary and 6 Sensitivity observations.
10. Check uniqueness of patient × compartment × gene observations.
11. Validate numeric values and raw UMI integrity.
12. Recalculate CPM and compare it with the locked STEP10E values.
13. Recalculate Log2(CPM + 1) and compare it with STEP10E values.
14. Verify the detection definition (`Raw UMI > 0`).
15. Separate Primary and Sensitivity expression records.
16. Generate Primary gene-by-compartment descriptive profiles.
17. Generate a separate Sensitivity Endothelial descriptive profile.
18. Identify patients represented in all three Primary compartments.

**Descriptive profile measures:**

For each gene × compartment combination, the workflow summarizes:

- Number of patients.
- Number and percentage of detected patients.
- Raw UMI mean, median, quartiles, minimum, maximum, and interquartile range.
- CPM mean, median, quartiles, minimum, maximum, and interquartile range.
- Log2(CPM + 1) mean, median, quartiles, minimum, maximum, and interquartile range.

The workflow is descriptive and performs structural and numerical consistency checks. It does not recalculate pseudobulk expression, perform differential-expression analysis, or conduct inferential statistical testing.

### 6.8 STEP10F-2 — Gene–Compartment Expression Pattern Audit

STEP10G is documented as integrating the locked outputs of STEP10F-1, STEP10F-2, and STEP10F-3.

The STEP10F-2 output directory is identified as:

`STEP10F-2_Gene_Compartment_Expression_Pattern_Audit`

A detailed standalone description of the STEP10F-2 procedures and individual output filenames was not included in the uploaded README sections summarized here. Its scope should therefore be documented from the corresponding script section and output README rather than inferred from the STEP10G description.

### 6.9 STEP10F-3 — Expression-Context Characterization

**Objective:** Characterize patient-level gene expression and compartment-specific expression context using the locked STEP10F patient-level outputs.

The documented workflow includes:

1.  Patient-level detection and quantitative abundance summaries.
2.  Reconciliation with the locked STEP10F-1C profile.
3.  Within-gene compartment ranking by median CPM.
4.  Within-gene ranking by detection prevalence.
5.  Highest-versus-lowest median CPM contrasts.
6.  Descriptive detection–abundance classification.
7.  Manuscript-safe interpretation auditing.
8.  Separate Endothelial Sensitivity analysis.
9.  Final quality-control summary.
10. Methodological lock and final report generation.

The workflow uses patient-level expression records and previously generated profiles. It does not recalculate pseudobulk counts or reread the GEO raw UMI matrix or H5 files.

The documented analysis is descriptive. It does not perform statistical testing, differential-expression analysis, correlation analysis, or CellChat/ligand–receptor analysis.

The workflow depends on objects established in preceding script sections, including:

- `primary_df`
- `sensitivity_df`
- `primary_profile`
- `sensitivity_profile`
- `out_dir`

### 6.10 STEP10G — Single-Cell Evidence Synthesis

**Objective:** Integrate the locked outputs of STEP10F-1, STEP10F-2, and STEP10F-3 into a coherent patient-level and gene-by-compartment evidence summary.

**Input dependencies:**

STEP10G reads 14 predefined files from:

- `STEP10F-1_Patient_Level_Expression_Profile`
- `STEP10F-2_Gene_Compartment_Expression_Pattern_Audit`
- `STEP10F-3_Expression_Context_Audit`

**Main procedures:**

1.  Audit input-file availability.
2.  Read the locked STEP10F outputs.
3.  Perform basic schema validation.
4.  Check core-gene and compartment identity.
5.  Synthesize patient-level structural information.
6.  Verify the complete three-compartment patient set.
7.  Integrate gene-by-compartment quantitative evidence.
8.  Summarize within-gene quantitative and detection rankings.
9.  Generate manuscript-oriented interpretation summaries.
10. Audit unsupported claims.
11. Document the evidence hierarchy.
12. Generate manuscript evidence tables.
13. Generate methodological lock records.
14. Perform final quality-control checks.
15. Generate the final lock file.

STEP10G is intended for descriptive evidence integration and manuscript-oriented reporting.

It does not recalculate pseudobulk counts, reread GEO raw matrices or H5 files, perform differential-expression analysis, or conduct statistical testing.

### 6.11 STEP10H — Final Manuscript Supplementary Export

**Objective:** Generate manuscript-oriented supplementary tables and standalone figure panels from the locked STEP10D-B and STEP10F-1 outputs.

**Input sources:**

- Locked STEP10D-B core-gene raw-count outputs.
- Locked STEP10F-1 patient-level expression outputs.
- Locked Primary expression summary tables.

**Main procedures:**

1.  Discover and preview relevant STEP10D and STEP10F files.
2.  Check the availability of four locked input tables.
3.  Read the raw-count, patient-level expression, and Primary summary tables.
4.  Validate required columns and core-gene identity.
5.  Validate Primary patient and compartment counts.
6.  Export Supplementary Table Sz.
7.  Export Supplementary Table Saa.
8.  Prepare Figure Sz Panels A, B, and C.
9.  Export each figure panel as PDF and 600-dpi TIFF.
10. Validate the expected figure output files.

**Locked analysis scope:**

Core genes:

- AOC3
- MS4A7
- IL7R

Primary compartments:

- T-cell
- Myeloid/macrophage
- Endothelial

The Sensitivity Endothelial group is retained separately in Supplementary Table Saa.

**Path convention:**

The STEP10 data directory is resolved using:

`here::here()`

The historical STEP10 input and output placement is preserved under:

`data/Paper 1_Step 10-Single-Cell/`

STEP10H uses locked upstream expression and detection summaries. It does not perform new pseudobulk aggregation, expression normalization, differential-expression analysis, statistical testing, correlation analysis, cell-cell communication analysis, ligand–receptor analysis, or causal inference.

The script-generated PASS messages represent execution-time checks and should not be interpreted as independent verification of the underlying analysis.

------------------------------------------------------------------------

## 7. Analytical Definitions and Quality-Control Rules

### 7.1 Core Genes

The predefined core-gene set is:

- AOC3
- MS4A7
- IL7R

The same gene set is used throughout the patient-level pseudobulk, detection, expression-profile, and evidence-synthesis workflows.

### 7.2 Raw UMI Detection

Gene detection is defined as:

`Raw UMI count > 0`

A zero raw UMI count is treated as a non-detected observation.

Zero-count observations are retained in the patient-level datasets and are not used as a reason for patient exclusion.

### 7.3 Pseudobulk Aggregation

For each patient and compartment, raw UMI counts are summed across the cells assigned to that patient-compartment group.

The patient is the biological replicate.

The analysis uses the original GEO raw UMI matrix for pseudobulk aggregation, rather than normalized expression values or database-provided differential-expression statistics.

### 7.4 CPM Calculation

For a gene in a patient-level pseudobulk sample:

`CPM = (Gene raw UMI count / Total pseudobulk UMI count) × 1,000,000`

The log-transformed expression measure is:

`Log2(CPM + 1)`

These measures are used for descriptive patient-level expression characterization.

### 7.5 Primary and Sensitivity Analysis Sets

Primary and Sensitivity Endothelial observations are maintained as separate analysis sets.

Sensitivity observations are not incorporated into Primary-only descriptive summaries.

### 7.6 Patient-Level Concordance

The complete three-compartment patient set is defined by intersecting the observed patient sets in the Primary T-cell, Myeloid/macrophage, and Endothelial compartments.

The analysis does not assume that all nine patients are represented in all three compartments.

### 7.7 Quality-Control and Provenance Checks

The documented quality-control procedures include:

- Input-file existence and schema validation.
- Cell-ID and patient-ID uniqueness.
- Cell-to-compartment mapping validation.
- Gene identity and ordering checks.
- Raw UMI non-negativity and integer-value checks.
- Missing and non-finite value checks.
- Pseudobulk matrix and annotation reconciliation.
- Expected patient and group-count checks.
- Patient × compartment × gene uniqueness.
- Detection-status consistency.
- CPM and Log2(CPM + 1) consistency checks.
- Output-file existence and methodological lock checks.

A PASS result indicates that the corresponding programmed checks were satisfied during execution. It is not, by itself, independent confirmation of biological validity or statistical robustness.

------------------------------------------------------------------------

## 8. Analysis Scope and Interpretation Boundaries

The single-cell workflows provide complementary cell-level and patient-level expression evidence.

The following limitations should be maintained when interpreting or reporting the results:

1.  TISCH2 cell-type annotations and database-provided differential-expression summaries are used as supplied resources. They should not be interpreted as proof of exclusive gene localization to a particular cell type.
2.  The patient-level pseudobulk analysis is based on the available Primary and Tumor cell populations and the predefined compartment eligibility rules.
3.  The endothelial compartment has fewer eligible patients than the T-cell and myeloid/macrophage compartments. Its descriptive results should therefore be interpreted with the corresponding sample-size limitation.
4.  Endothelial Sensitivity observations are kept separate from the Primary endothelial analysis.
5.  Gene detection and quantitative abundance are distinct expression characteristics. Detection prevalence should not be treated as equivalent to expression magnitude.
6.  Patient-level descriptive expression summaries do not establish statistically significant differences between compartments.
7.  The single-cell workflow does not independently establish causal signaling relationships or cell-cell communication.
8.  The workflow does not include CellChat analysis or ligand–receptor inference.
9.  The results provide expression-context characterization rather than proof of a specific functional mechanism.
10. The single-cell evidence should be integrated with the bulk transcriptomic findings without overstating cell-type specificity or causal interpretation.

------------------------------------------------------------------------

## 9. Output Organization

The documented STEP10 workflow outputs are organized under the single-cell analysis directory.

The following workflow directories are explicitly identified in the uploaded script README sections:

| Workflow  | Output directory                                             |
|-----------|--------------------------------------------------------------|
| STEP10A   | `results/STEP10A_Patient_Sample_Structure_Audit/`            |
| STEP10B-1 | `results/STEP10B-1_Compartment_CellCount_Feasibility_Audit/` |
| STEP10B-2 | `results/STEP10B-2_Pseudobulk_Feasibility_CellDepth_Audit/`  |
| STEP10C   | `results/STEP10C_Patient_Level_Pseudobulk/`                  |
| STEP10F-1 | `STEP10F-1_Patient_Level_Expression_Profile`                 |
| STEP10F-2 | `STEP10F-2_Gene_Compartment_Expression_Pattern_Audit`        |
| STEP10F-3 | `STEP10F-3_Expression_Context_Audit`                         |

Other STEP10 workflow sections generate intermediate or final audit tables, expression profiles, supplementary tables, figures, session information, and methodological lock records.

The exact output paths for all intermediate and final files should be taken from the corresponding script definitions and execution records. They are not inferred where the uploaded README sections do not specify them.

------------------------------------------------------------------------

## 10. Manuscript-Related Outputs

The single-cell analysis contributes cell-level and patient-level expression characterization to the LUAD manuscript.

The documented outputs include:

- Major cell-population UMAP visualization.
- AOC3, MS4A7, and IL7R expression UMAP panels.
- Differential-expression dot plot.
- Patient-level raw-count and detection summaries.
- Patient-level compartment concordance summaries.
- Gene-by-compartment expression profiles.
- Expression-context and compartment-ranking summaries.
- Integrated manuscript-oriented evidence tables.
- Supplementary Table Sz.
- Supplementary Table Saa.
- Figure Sz Panels A, B, and C.

The figure and supplementary-item numbering must be checked against the final manuscript before publication.

In particular, the original high-level README and the uploaded script README sections contain different figure-number references. The original overview refers to Figure 6A–C, whereas the Script 1 fragments refer to Figure 9A–9C. The STEP10H documentation uses the provisional labels Figure Sz and Supplementary Tables Sz and Saa.

These labels are retained here pending reconciliation with the final manuscript and exported files.

------------------------------------------------------------------------

## 11. Reproducibility Notes

### 11.1 Project-Relative Paths

The STEP10 workflow uses project-relative paths, with the project root resolved through:

`here::here()`

The historical STEP10 data placement is preserved under:

`data/Paper 1_Step 10-Single-Cell/`

The scripts should be run from a valid project structure with the required input files available at their documented locations.

### 11.2 Randomness and Seeds

The documented random seed for Script 1 is:

`20250724`

No separate random seed for the STEP10 workflows is specified in the uploaded README sections. Any additional random-number settings should be documented from the final script code rather than inferred.

### 11.3 Workflow Dependencies

The two scripts are complementary but have different input and output responsibilities.

Script 1 focuses on cell-level visualization.

Script 2 depends on the TISCH2 metadata, GEO raw UMI matrix, preceding STEP10 outputs, and locked patient-level analysis objects.

The later STEP10 workflows depend on outputs generated by preceding STEP10 sections. The required upstream objects and files should therefore be available before running the corresponding downstream sections.

### 11.4 Execution and Verification Status

This README consolidates the workflow descriptions provided in the script-generated README sections.

It does not represent an independent execution or validation of the R scripts.

Script-generated PASS messages, expected group counts, and documented output checks should be interpreted as workflow-defined checks unless supported by corresponding execution logs and independently inspected output files.

Before public release, the final scripts should be checked for:

- Consistent project-relative paths.
- Consistent input and output directory definitions.
- Availability of all required upstream files and objects.
- Consistent STEP10D-D final output definitions.
- Correct figure and supplementary-table numbering.
- Consistency between the documented and actual exported filenames.

------------------------------------------------------------------------

## 12. Summary

The two single-cell analysis scripts provide a structured workflow for characterizing AOC3, MS4A7, and IL7R in the NSCLC_GSE131907 dataset.

Script 1 provides cell-level visualization and expression-display outputs.

Script 2 establishes the patient-level analysis cohort, audits cellular compartment availability, constructs raw-UMI pseudobulk matrices, evaluates core-gene detection and patient-level compartment concordance, generates descriptive expression profiles, integrates the evidence, and exports manuscript-oriented supplementary materials.

Together, these workflows provide complementary cell-level and patient-level expression evidence while preserving the predefined sample structure, raw-count provenance, and descriptive scope of the analysis.

The resulting evidence is intended to support cautious interpretation of the expression patterns of AOC3, MS4A7, and IL7R across the LUAD tumor microenvironment, without implying exclusive cell-type localization, causal mechanisms, or statistically established compartment-specific differences.
