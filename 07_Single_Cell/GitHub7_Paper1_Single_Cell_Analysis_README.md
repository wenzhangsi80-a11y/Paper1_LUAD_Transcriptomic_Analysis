# Single-Cell Transcriptomic Analysis

## 1. Overview

This directory contains the single-cell transcriptomic characterization of three prioritized study genes, AOC3, MS4A7, and IL7R, in lung adenocarcinoma (LUAD).

The analysis uses the NSCLC_GSE131907 dataset available through the Tumor Immune Single-cell Hub 2 (TISCH2) database and complementary patient-level pseudobulk profiles constructed from the original GEO raw unique molecular identifier (UMI) count matrix.

The workflow is organized into two R scripts:

- `GitHub7_Paper1_Single_Cell_Analysis_1.R`
- `GitHub7_Paper1_Single_Cell_Analysis_2.R`

Together, these scripts provide:

1.  Single-cell UMAP visualization of major annotated cellular populations.
2.  Visualization of AOC3, MS4A7, and IL7R expression across individual cells.
3.  Cell-type-associated differential-expression summaries obtained from TISCH2.
4.  Patient-level raw-UMI pseudobulk construction across predefined cellular compartments.
5.  Patient-level gene-detection and quantitative expression characterization.
6.  Structural quality-control and provenance audits.
7.  Integrated expression-context evidence synthesis.
8.  Manuscript-oriented supplementary figure and table generation.

The single-cell analyses provide downstream cellular expression characterization of the prioritized study genes. They are descriptive and observational and are not intended to establish cell-type specificity, causal mechanisms, functional relationships, or direct intercellular communication.

## 2. Directory Structure

The directory is organized into the following main components:

``` text
07_Single_Cell/
│
├── README.md
│
├── data/
│   ├── data_sources/
│   ├── input_data/
│   ├── processed_data/
│   └── sample_manifests/
│
├── environment/
│   └── README.md
│
├── results/
│   └── README.md
│
└── scripts/
    ├── GitHub7_Paper1_Single_Cell_Analysis_1.R
    ├── GitHub7_Paper1_Single_Cell_Analysis_2.R
    └── README.md
```

The directory-specific README files document the corresponding input data, software environment, scripts, generated outputs, and workflow-specific reproducibility information.

## 3. Data Source and Required Input Files

### 3.1 Data source

Single-cell data and associated annotations were obtained from the Tumor Immune Single-cell Hub 2 (TISCH2):

- Website: <https://tisch.compbio.cn>
- Dataset: `NSCLC_GSE131907`
- Original study accession: `GSE131907`

TISCH2 provides cell-level metadata, major-lineage annotations, dimensionality-reduction coordinates, and database-derived expression summaries.

### 3.2 Required input files

The principal input files include:

- `NSCLC_GSE131907_CellMetainfo_table.tsv`
- `NSCLC_GSE131907_expression.h5`
- `GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz`

The files serve different purposes:

| Input file | Main purpose |
|----|----|
| TISCH2 cell metadata | Cell annotations, patient and sample identifiers, source, stage, and major-lineage information |
| TISCH2 expression H5 file | Expression-matrix structure inspection, feature and barcode characterization, and cell-ID reconciliation |
| GEO raw UMI matrix | Original raw-count input for patient-level pseudobulk construction |

The H5 expression file and GEO raw UMI matrix are not treated as interchangeable inputs. The patient-level pseudobulk expression matrices are constructed from the GEO raw UMI matrix after cell-ID reconciliation with the locked TISCH2 target cohort.

Additional intermediate tables, manifests, annotations, and locked workflow outputs are generated or consumed by the corresponding STEP10 analysis modules.

Input-file locations must be configured using the shared project setup before execution.

## 4. Analysis Scripts

### 4.1 `GitHub7_Paper1_Single_Cell_Analysis_1.R`

**Scope: Figure 9A–C**

This script supports the single-cell visualization and database-derived expression-summary components of the study.

#### Figure 9A — Cell-Type UMAP Visualization

This analysis:

- Uses the NSCLC_GSE131907 dataset.
- Uses the previously provided UMAP coordinates.
- Displays 12 major annotated cell types.
- Applies a predefined cell-type color mapping.
- Generates publication-oriented UMAP visualizations.

The UMAP coordinates are taken from the supplied dataset annotations and are not independently recalculated in this workflow.

#### Figure 9B — Target-Gene Expression Visualization

This analysis visualizes the expression distributions of:

- AOC3
- IL7R
- MS4A7

Expression values are extracted from the Seurat RNA counts layer and normalized using the per-cell library size:

`log1p(count / library size × 10,000)`

The resulting expression values are displayed on the UMAP coordinates to characterize their distributions across individual cells.

These visualizations describe expression patterns in the analyzed dataset and do not independently establish cell-type specificity.

#### Figure 9C — Differential-Expression Dot Plot

This analysis uses the TISCH2-provided differential-expression results table.

The workflow:

- Extracts the records for AOC3, IL7R, and MS4A7.
- Organizes the results by annotated cell type.
- Summarizes the database-provided log2 fold-change values.
- Displays the corresponding proportions of expressing cells.
- Generates a cell-type-associated differential-expression dot plot.

The differential-expression statistics, log2 fold-change values, and expression percentages are obtained from the database-provided results. They are not independently recalculated by this visualization workflow.

### 4.2 `GitHub7_Paper1_Single_Cell_Analysis_2.R`

**Scope: Patient-level pseudobulk characterization, quality control, evidence synthesis, and manuscript supplementary export.**

This script implements the STEP10A–STEP10H workflow.

The main analytical stages are:

| Workflow | Main purpose |
|----|----|
| STEP10A | Patient and sample structure audit |
| STEP10B-1 | Primary-tumor compartment cell-count feasibility audit |
| STEP10B-2 | Patient-level pseudobulk feasibility and cell-depth audit |
| STEP10C-0A–0E | Expression-file characterization, Cell-ID and Gene-ID reconciliation, and pseudobulk input provenance audits |
| STEP10C | Patient-level raw UMI pseudobulk construction |
| STEP10D-A–D | Pseudobulk quality control, core-gene detection, compartment-level expression, and patient-level concordance audits |
| STEP10E | Patient-level pseudobulk quantitative expression characterization |
| STEP10F-1–3 | Patient-level expression profiling and expression-context characterization |
| STEP10G | Integrated single-cell evidence synthesis |
| STEP10H | Final manuscript supplementary figure and table export |

## 5. Dataset Overview and Cell Selection

### 5.1 TISCH2 dataset overview

The TISCH2 NSCLC_GSE131907 dataset contains 203,298 cells from 44 patients across 58 samples, according to the dataset-level information used in the manuscript.

The dataset provides the cell-level annotations used for the initial visualization and downstream patient-level cohort definition.

### 5.2 Primary-tumor target cohort

For patient-level pseudobulk analysis, the target cohort is defined using both metadata fields:

- `Stage = Primary`
- `Source = Tumor`

The locked Primary + Tumor subset comprises:

- 9 patients
- 65,702 cells

Both metadata fields are required for cohort selection. `Stage = Primary` is not assumed to be equivalent to `Source = Tumor`.

The target cohort is reconciled against the original GEO raw UMI matrix before pseudobulk construction.

### 5.3 Predefined cellular compartments

Three cellular compartments are used for patient-level pseudobulk characterization.

| Target compartment | TISCH2 major-lineage labels |
|--------------------|-----------------------------|
| T-cell             | CD4Tconv, CD8T, CD8Tex      |
| Myeloid/macrophage | Mono/Macro                  |
| Endothelial        | Endothelial                 |

Dendritic cells are not included in the locked myeloid/macrophage compartment.

The predefined mapping is applied consistently throughout the STEP10B–STEP10H workflows. No new cell annotation procedure is introduced by the patient-level pseudobulk analysis.

## 6. Patient-Level Cell-Depth Eligibility

Patient–compartment groups are classified using predefined operational cell-count criteria.

| Cell count per patient and compartment | Eligibility category | Use |
|----|----|----|
| ≥50 cells | Primary analysis | Main pseudobulk analysis |
| 20–49 cells | Sensitivity analysis | Separate endothelial sensitivity analysis |
| \<20 cells | Insufficient cell depth | Not included in formal pseudobulk expression characterization |

The resulting locked cohort structure is:

| Analysis set | Compartment        | Expected patients |
|--------------|--------------------|------------------:|
| Primary      | T-cell             |                 9 |
| Primary      | Myeloid/macrophage |                 9 |
| Primary      | Endothelial        |                 4 |
| Sensitivity  | Endothelial        |                 2 |

Three additional patients have fewer than 20 endothelial cells and are not included in the formal endothelial pseudobulk analysis.

The minimum cell-count thresholds are study-defined operational eligibility criteria. They should not be interpreted as universal biological or statistical standards for pseudobulk analysis.

No patient is excluded during the preliminary STEP10B feasibility audits. The eligibility criteria are applied when defining the downstream pseudobulk analysis groups.

## 7. Patient-Level Raw UMI Pseudobulk Construction

### 7.1 Raw-count aggregation

STEP10C constructs patient-level pseudobulk profiles from the original GEO GSE131907 raw UMI matrix.

The patient is treated as the biological replicate.

For each gene, raw UMI counts are summed across cells belonging to the same patient and predefined cellular compartment.

The workflow:

1.  Reconciles TISCH2 metadata Cell IDs with GEO raw UMI matrix column identifiers.
2.  Applies the locked Primary + Tumor cohort.
3.  Retains only the predefined target-cell columns.
4.  Matches patient and compartment assignments to the corresponding cell identifiers.
5.  Processes the raw UMI matrix sequentially, one gene row at a time.
6.  Aggregates raw counts by gene, patient, and compartment.
7.  Constructs separate patient-level pseudobulk matrices.
8.  Performs matrix-integrity and sample-count checks.

The GEO raw-only cells are not added to the TISCH2-defined target cohort.

The complete GEO raw UMI matrix is not loaded into memory at once. The matrix is processed sequentially to extract and aggregate the required target-cell columns.

### 7.2 Pseudobulk group structure

The expected primary pseudobulk structure contains:

- 9 T-cell patient-compartment groups.
- 9 Myeloid/macrophage patient-compartment groups.
- 4 Endothelial patient-compartment groups.

This produces 22 Primary patient-compartment groups.

The endothelial sensitivity dataset contains an additional 2 groups.

The total locked structure therefore contains 24 patient-compartment groups across the Primary and Sensitivity datasets.

### 7.3 Core genes

The pseudobulk analyses focus on:

- AOC3
- MS4A7
- IL7R

These genes are assessed using the locked raw UMI pseudobulk matrices and corresponding patient-level annotations.

Pseudobulk profiles represent patient-level, compartment-associated transcriptional signals. They are not interpreted as direct evidence of cell-intrinsic regulation or cell-type specificity.

## 8. Expression and Detection Measures

### 8.1 Raw UMI counts

Raw UMI counts are retained in the patient-level pseudobulk tables.

Gene detection is defined as:

`Raw UMI count > 0`

Zero-count observations are retained in the descriptive datasets.

Detection status is not used as an additional exclusion criterion.

### 8.2 Counts per million

For each patient-compartment pseudobulk sample, CPM is calculated using the total pseudobulk UMI count as the library-size denominator:

`CPM = gene raw UMI count / total pseudobulk UMI count × 1,000,000`

### 8.3 Log-transformed expression

The quantitative expression summaries also use:

`log2(CPM + 1)`

Raw UMI counts, CPM, and log2(CPM + 1) are retained or summarized as specified by the downstream workflow.

### 8.4 Descriptive summaries

The patient-level expression characterization includes:

- Number of patient-compartment groups.
- Detection count.
- Detection percentage.
- Mean expression.
- Median expression.
- First and third quartiles.
- Minimum and maximum expression.
- Interquartile range.

Primary and Sensitivity endothelial observations are summarized separately.

No inferential statistical testing is performed for the patient-level pseudobulk compartment summaries.

## 9. Quality Control and Provenance Audits

The STEP10A–STEP10D workflows establish the structural and provenance checks required before downstream expression characterization.

### 9.1 Patient and sample structure

STEP10A audits:

- Cell-level metadata integrity.
- Unique patient and sample identifiers.
- Patient-level sample structure.
- Source and stage distributions.
- Stage-by-source relationships.
- Primary-tumor patient classification.
- Cell-population distributions across patients.
- Availability of the three prioritized genes in the TISCH2 differential-expression summary table.

No PCA, UMAP, clustering, or cell annotation is recalculated in this audit step.

### 9.2 Compartment cell-count feasibility

STEP10B-1 and STEP10B-2 verify:

- The Primary + Tumor target cohort.
- The predefined compartment mapping.
- Patient-level compartment cell counts.
- Three-compartment availability.
- Patient-level cell-depth eligibility.
- Reconciliation of patient and cell counts across analysis stages.

The expected target cohort is 9 patients and 65,702 cells.

### 9.3 Expression-matrix and identifier audits

STEP10C-0A–0E characterize and reconcile the available expression inputs.

The audits include:

- File existence and file-size checks.
- ZIP and HDF5 structure inspection.
- Matrix dimensions and sparse CSC structure.
- Feature and barcode inspection.
- Sampled expression-value characterization.
- GEO raw UMI file readability and matrix-orientation assessment.
- GEO raw UMI Cell-ID uniqueness.
- Gene-ID extraction and core-gene availability.
- Exact Cell-ID matching between TISCH2 metadata and GEO raw UMI data.
- Compartment-level and patient-level Cell-ID reconciliation.
- Primary + Tumor target-cell coverage assessment.
- Construction of the cell-level pseudobulk input manifest.

Barcode-format normalization is used for diagnostic sensitivity assessment where applicable. It does not replace official exact Cell-ID matching.

H5 matrix structure and sampled integer-valued expression characteristics are not, by themselves, sufficient to establish raw-count provenance. The GEO raw UMI matrix is used as the raw-count source for pseudobulk construction.

### 9.4 Pseudobulk integrity and patient-level audits

STEP10D performs structural and descriptive quality control of the locked pseudobulk matrices.

The audits include:

- Matrix dimensions and gene identity.
- Core-gene availability.
- Raw UMI count integrity.
- Pseudobulk library-size summaries.
- Detected-gene summaries.
- Patient-compartment annotation consistency.
- Core-gene raw-count extraction.
- Gene-detection status.
- Compartment-level detection summaries.
- Patient-level compartment overlap.
- Complete three-compartment patient-set identification.
- Patient-level detection concordance.

The expected Primary compartment sizes are 9 T-cell, 9 Myeloid/macrophage, and 4 Endothelial patients.

The complete three-compartment patient set is determined from the observed intersection of the three compartment-specific patient sets. The locked downstream workflow specifies an expected size of four patients.

Low-library flags, where calculated using the Q1 − 1.5 × IQR rule, are descriptive QC indicators and do not automatically result in sample exclusion.

## 10. Integrated Expression-Context Characterization

### 10.1 STEP10E — Patient-Level Pseudobulk Expression Characterization

STEP10E uses the four locked pseudobulk expression matrices and their corresponding sample annotations.

The workflow:

- Validates the required input files.
- Verifies gene identity and ordering across matrices.
- Checks raw UMI count integrity.
- Calculates pseudobulk library sizes.
- Derives CPM and log2(CPM + 1).
- Defines gene detection using raw UMI count \> 0.
- Constructs patient-level core-gene expression tables.
- Generates separate Primary and Sensitivity descriptive summaries.
- Exports methodological and final audit records.

The Primary dataset contains 22 patient-compartment groups and 66 gene-level observations across the three core genes.

The Sensitivity endothelial dataset contains 2 patient-compartment groups and 6 gene-level observations.

### 10.2 STEP10F — Patient-Level Expression Profiling and Context Audits

STEP10F uses the locked STEP10E expression outputs.

The workflow is organized into three components:

- STEP10F-1: Patient-level core-gene expression profiling.
- STEP10F-2: Gene-compartment expression-pattern audit.
- STEP10F-3: Expression-context characterization.

The analyses integrate:

- Patient-level raw UMI counts.
- CPM.
- Log2(CPM + 1).
- Detection prevalence.
- Quantitative expression abundance.
- Compartment-level descriptive summaries.
- Patient-level expression patterns.
- Ranking and expression-context comparisons.

Primary and Sensitivity endothelial observations remain separate throughout the workflow.

The complete three-compartment patient set is identified from the intersection of the Primary T-cell, Myeloid/macrophage, and Endothelial patient sets.

The expected complete set contains four patients.

### 10.3 STEP10G — Integrated Single-Cell Evidence Synthesis

STEP10G integrates the locked outputs from:

- STEP10F-1: Patient-level expression profiles.
- STEP10F-2: Gene-compartment expression-pattern audit.
- STEP10F-3: Expression-context audit.

The synthesis focuses on:

- Patient-level detection prevalence.
- Quantitative expression abundance.
- Compartment-associated expression context.
- Manuscript-oriented evidence summaries.
- Unsupported-claim auditing.
- Evidence hierarchy documentation.

STEP10G checks the availability and selected structural properties of its expected upstream input files before generating the evidence-synthesis outputs.

This step does not recalculate pseudobulk expression or introduce new statistical inference.

## 11. Final Manuscript Supplementary Export

### STEP10H — Final Manuscript Supplementary Export

STEP10H generates manuscript-oriented supplementary tables and publication-ready standalone figure panels using locked upstream outputs.

The workflow uses four input tables:

**STEP10D-B**

- `STEP10D-B_Pseudobulk_CoreGene_RawCounts.tsv`

**STEP10F-1**

- `STEP10F-1A_Primary_Patient_Level_Expression_Profile.tsv`
- `STEP10F-1B_Sensitivity_Endothelial_Expression_Profile.tsv`
- `STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv`

### 11.1 Supplementary Table Sz

Supplementary Table Sz is derived from the locked STEP10F-1C summary table.

It contains nine gene-compartment records covering:

- AOC3
- MS4A7
- IL7R

The table reports patient counts, detection prevalence, and descriptive expression summaries.

STEP10H does not recalculate the upstream expression statistics.

### 11.2 Supplementary Table Saa

Supplementary Table Saa is derived from the locked STEP10D-B raw-count pseudobulk table.

It retains patient-level raw UMI counts for the three prioritized study genes across Primary and Sensitivity patient-compartment groups.

### 11.3 Supplementary Figure Sz

Supplementary Figure Sz is exported as three independent panels:

- Panel A: Detection prevalence.
- Panel B: Quantitative expression.
- Panel C: Patient-level expression.

The same predefined compartment color mapping is used across the panels.

Each panel is exported in:

- Vector PDF format.
- 600-dpi TIFF format.

The workflow exports three standalone panels rather than one combined figure.

### 11.4 Export quality control

STEP10H checks:

- Availability of all four locked input files.
- Required input columns.
- Core-gene identity.
- Primary compartment and patient counts.
- Expected Table Sz structure.
- Expected Table Saa structure.
- Expected Figure Sz observation count.
- Exported figure-file existence and file size.

STEP10H is an export and presentation workflow. It does not introduce new expression normalization, pseudobulk recalculation, differential-expression analysis, or statistical testing.

The figure and table identifiers Sz and Saa are workflow labels and should be replaced or confirmed against the final manuscript numbering before publication.

## 12. Main Expression-Context Findings

The single-cell analyses provide two complementary levels of expression characterization.

### 12.1 TISCH2-derived cell-population patterns

The TISCH2-derived visualizations and database-provided differential-expression summaries show distinct cell-population-associated expression patterns:

- **IL7R:** Predominantly associated with T-cell populations.
- **MS4A7:** Predominantly associated with myeloid/monocyte–macrophage populations.
- **AOC3:** Displays an endothelial-associated pattern in the database-derived annotations.

These patterns describe the database-provided cell annotations and expression summaries. They should not be interpreted as proof of exclusive cell-type localization.

### 12.2 Patient-level pseudobulk expression

Patient-level pseudobulk analysis provides complementary quantitative expression context.

The descriptive median CPM values reported in the initial analysis were:

| Gene  | T-cell | Myeloid/macrophage | Endothelial |
|-------|-------:|-------------------:|------------:|
| MS4A7 |   8.48 |             506.20 |       24.19 |
| IL7R  | 258.16 |              54.22 |      124.84 |

For AOC3, the highest median CPM was observed in the Myeloid/macrophage compartment (1.33).

AOC3 detection was observed in:

- 8/9 Myeloid/macrophage patient-compartment profiles.
- 6/9 T-cell patient-compartment profiles.
- 2/4 Primary Endothelial patient-compartment profiles.

These values are descriptive summaries from the patient-level pseudobulk analysis. They do not establish statistically significant differences between compartments.

In particular, the patient-level pseudobulk results do not provide independent evidence of endothelial enrichment for AOC3.

The TISCH2-derived summaries and patient-level pseudobulk results are complementary but are not interchangeable. They differ in data source, analytical level, and expression representation.

## 13. Manuscript Mapping

The single-cell workflows support the following manuscript components.

| Manuscript component | Supported analysis |
|----|----|
| Methods, Section 2.11 | Single-cell transcriptomic and patient-level pseudobulk analysis |
| Results, Section 3.9 | Single-cell transcriptomic characterization of the prioritized study genes |
| Figure 9A | Major annotated cellular populations |
| Figure 9B | UMAP expression distribution of AOC3, MS4A7, and IL7R |
| Figure 9C | TISCH2-derived cell-type-associated differential-expression dot plot |
| Supplementary Figure Sz | Patient-level detection and quantitative expression panels, subject to final manuscript numbering |
| Supplementary Table Sz | Gene-compartment descriptive expression summaries, subject to final manuscript numbering |
| Supplementary Table Saa | Patient-level raw UMI pseudobulk expression profiles, subject to final manuscript numbering |

The final supplementary figure and table numbering should be checked against the submitted manuscript and the corresponding supplementary materials.

Any separate ligand–receptor expression-compatibility analysis should be documented in its own workflow or explicitly mapped to its generating script. It should not be assumed to be generated by the STEP10A–STEP10H workflows described here.

## 14. Interpretation and Limitations

The single-cell analyses provide descriptive cellular expression context for AOC3, MS4A7, and IL7R in the NSCLC_GSE131907 dataset.

Several distinctions are important for interpretation:

1.  TISCH2-derived cell-type expression summaries are database-provided results and are not independently recalculated by the visualization workflow.
2.  Patient-level pseudobulk profiles are constructed by aggregating raw UMI counts within patient and predefined cellular compartment.
3.  Primary and Sensitivity endothelial observations are handled separately.
4.  Gene detection is defined as raw UMI count greater than zero.
5.  CPM and log2(CPM + 1) are used for descriptive quantitative expression characterization.
6.  No differential-expression testing or inferential statistical comparison between patient-level pseudobulk compartments is performed in the described STEP10C–STEP10H workflows.
7.  Compartment-associated expression patterns do not establish cell-type specificity, exclusivity, or cell-intrinsic regulation.
8.  Expression compatibility does not establish ligand–receptor interaction strength, communication probability, pathway activation, or experimental validation.
9.  The analyses do not establish causal or functional relationships between cell populations.

Interpretation should remain within the scope of the analyzed dataset, annotations, patient cohort, predefined compartments, and operational cell-depth eligibility criteria.

## 15. Reproducibility

To reproduce the single-cell analyses:

1.  Obtain the relevant NSCLC_GSE131907 files from TISCH2 and the original GEO raw UMI matrix.
2.  Confirm that the metadata, H5 expression file, and raw UMI matrix correspond to the intended dataset.
3.  Place the input files in the appropriate project data directories.
4.  Configure the project root and shared input/output path variables using the project setup script.
5.  Verify the required R packages and software versions.
6.  Run the scripts in the documented workflow order.
7.  Confirm the Primary + Tumor cohort and the locked patient-by-compartment cell counts.
8.  Review the Cell-ID and Gene-ID reconciliation outputs before pseudobulk construction.
9.  Confirm that Primary and Sensitivity endothelial datasets remain separate.
10. Review the generated QC, provenance, and methodological lock records.
11. Confirm the final supplementary figure and table exports against the manuscript.

The main scripts use project-relative path management and document the relevant package dependencies and reproducibility settings.

The random seed documented for Script 1 is `20250724`. The uploaded Script 2 workflows are predominantly deterministic and do not require an additional random seed.

The R version specified in the Script 2 header is 4.5.1. Actual R and package versions should be confirmed from execution-time session information.

For further details, consult:

- `scripts/README.md`
- `results/README.md`
- `environment/README.md`

These files provide the corresponding script, output, and software-environment documentation.

------------------------------------------------------------------------

**Scope of this README**

This document provides an overview of the complete single-cell analysis module and its relationship to the manuscript.

Detailed input-file inventories, intermediate outputs, package versions, execution records, and workflow-specific quality-control information are maintained in the corresponding subdirectory README files.
