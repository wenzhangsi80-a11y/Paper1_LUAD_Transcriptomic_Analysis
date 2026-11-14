# 02-WGCNA: Tumor-Only Weighted Gene Co-expression Network Analysis

## 1. Overview

This directory contains the analytical workflows and supporting documentation for the tumor-only weighted gene co-expression network analysis (WGCNA) conducted as part of the Paper 1 study:

*Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD*

The analysis investigates gene co-expression patterns in lung adenocarcinoma (LUAD), their relationships with clinical characteristics and overall survival (OS), and their associations with tumor microenvironment (TME)-related expression scores.

The tumor-only WGCNA analysis provides a complementary characterization of the three prioritized study genes, AOC3, MS4A7, and IL7R. It examines their module membership, gene significance (GS), module membership (MM), and TME associations using the primary tumor cohort. This analysis does not redefine the original candidate-gene selection or establish causal mechanisms.

The directory is organized into analysis scripts, input-data documentation, computational-environment documentation, and generated results.

## 2. Data Sources and Input Data

### 2.1 Gene expression data

The primary LUAD gene expression data were obtained from the SangerBox online platform.

- **Source platform:** SangerBox
- **Input file:** `Merge_RNA_seq_FPKM_UQ.txt`
- **Expression unit:** FPKM-UQ

This file provides the expression data used to reconstruct the tumor-only expression cohort and prepare the WGCNA input matrix.

The tumor-only analysis retains TCGA-LUAD Primary Tumor samples identified by sample type code 01.

### 2.2 Gene annotation

Gene annotation was obtained from GENCODE human release 22.

- **Source:** GENCODE
- **Annotation file:** `gencode.v22.annotation.gtf.gz`
- **Release:** GENCODE v22

The annotation reference is used for gene-identifier interpretation and gene annotation procedures where required by the analysis scripts.

### 2.3 Clinical data

Clinical information was assembled from TCGA-related clinical and sample-matching files used in the preceding study workflows.

The documented files include:

- `clinical_matched.csv`
- `clinical_matched.txt`
- `sample_group_final.csv`

A derived clinical workbook, referred to in the original workflow as `YOUR_TCGA_CLINICAL_FILE.xlsx`, was prepared for WGCNA-related cohort reconstruction and downstream clinical and survival analyses.

The exact clinical variables, sample-matching rules, inclusion criteria, and survival endpoints are defined in the corresponding R scripts and their documentation.

### 2.4 Derived and intermediate data

Additional expression matrices, sample manifests, clinical alignment files, annotation crosswalks, and intermediate R objects are generated during the workflow.

These derived files are distinguished from the original downloaded input data and are documented in the corresponding script and Results documentation.

**Data availability:** The input files are documented to describe the data sources and analytical workflow. The repository should not be assumed to contain all original downloaded or externally obtained input data. Users should obtain the relevant data from the indicated sources or use the corresponding authorized study files.

## 3. Analysis Cohort and Expression Preprocessing

The tumor-only WGCNA workflow uses **513 TCGA-LUAD Primary Tumor samples**, identified by TCGA sample type code 01.

The preprocessing workflow includes cohort reconstruction, clinical matching, expression-data quality control, gene-annotation checks, expression-prevalence filtering, transformation, and gene-variability assessment.

The principal preprocessing criteria are:

- Retain the 513 Primary Tumor samples, subject to documented sample-integrity checks.
- Preserve sample identifiers and expression-matrix sample order for cross-workflow verification.
- Match clinical records to expression samples using sample identifiers.
- Retain genes with raw FPKM-UQ values greater than 1 in at least 103 of the 513 samples.
- Transform retained expression values using `log2(FPKM-UQ + 1)`.
- Identify MAD-positive genes and prepare the expression matrix for downstream WGCNA.

PCA, hierarchical clustering, sample connectivity, and expression-distribution assessments are used as descriptive quality-control procedures. They should not be interpreted as automatic sample-exclusion criteria unless explicitly specified in the corresponding script.

The network analysis uses the prepared MAD-positive expression matrix. Subsequent integration and supplementary-analysis workflows primarily use previously generated and audited WGCNA objects.

## 4. Overall Analysis Strategy

The ten R scripts form a sequential analysis and integration workflow covering six major stages:

1.  **Cohort reconstruction and expression-data preparation:** Establish the 513-sample primary tumor cohort and prepare the expression matrix for network analysis.
2.  **WGCNA network construction and clinical analysis:** Select the soft-thresholding power, construct the co-expression network, identify modules, and evaluate module relationships with clinical characteristics and survival.
3.  **Additional survival-model assessment:** Examine proportional-hazards diagnostics and candidate time-dependent Cox models for the brown module and selected clinical variables.
4.  **TME scoring and module association analysis:** Calculate ESTIMATE-derived TME scores and assess their relationships with WGCNA module eigengenes.
5.  **Prioritized-gene and candidate-gene integration:** Integrate existing module assignments, GS, MM, and TME association results to characterize prioritized and candidate genes.
6.  **Supplementary analysis and reproducibility archiving:** Audit previously generated results, prepare Supplementary Figures 4A–4D and Supplementary Table X, and preserve relevant analysis workspaces.

The later integration and supplementary-analysis stages primarily use previously generated or locked results. They do not independently reconstruct the original WGCNA network or redefine the original candidate-gene set.

## 5. R Script Overview

The ten scripts are organized according to their analytical roles.

| Script | Main analytical role |
|----|----|
| `GitHub2_Paper1_tumor_only_WGCNA_1.R` | **Expression-data preparation and quality control.** Reconstructs the clinical cohort, extracts the 513 Primary Tumor expression samples, audits expression data and gene annotation, performs expression-prevalence filtering and transformation, assesses gene variability, and prepares the MAD-positive WGCNA input matrix. |
| `GitHub2_Paper1_tumor_only_WGCNA_2.R` | **WGCNA network construction and clinical analysis.** Recovers and validates the upstream input checkpoint, selects the soft-thresholding power, constructs the unsigned co-expression network, identifies and annotates modules, evaluates prioritized-gene module membership, and performs module–clinical and module-level survival analyses. |
| `GitHub2_Paper1_tumor_only_WGCNA_3.R` | **ME_brown survival-model assessment.** Examines proportional-hazards diagnostics and explores candidate time-dependent Cox models for ME_brown, including alternative time transformations and time-specific hazard-ratio estimates. |
| `GitHub2_Paper1_tumor_only_WGCNA_4.R` | **Gender-related time-dependent survival analysis.** Assesses the proportional-hazards assumption for Gender and compares candidate time-dependent Cox models, including time-specific Male-versus-Female hazard ratios. |
| `GitHub2_Paper1_tumor_only_WGCNA_5.R` | **Pathological-stage and integrated survival analysis.** Evaluates the proportional-hazards assumption for pathological stage, examines stage-specific contrasts, compares candidate time-dependent representations, and constructs the final integrated time-dependent Cox model. |
| `GitHub2_Paper1_tumor_only_WGCNA_6.R` | **ESTIMATE and module–TME association analysis.** Prepares an ESTIMATE-compatible expression matrix, calculates ImmuneScore, StromalScore, and ESTIMATEScore, and evaluates module eigengene associations with these TME traits using Spearman correlation and BH/FDR correction. |
| `GitHub2_Paper1_tumor_only_WGCNA_7.R` | **Gene-level GS/MM integration and prioritized-gene localization.** Integrates gene significance and module membership information with existing TME and WGCNA results, characterizes core-gene positions within their assigned modules, and summarizes the module-level context of AOC3, MS4A7, and IL7R. |
| `GitHub2_Paper1_tumor_only_WGCNA_8.R` | **Candidate-gene integration and historical-provenance audit.** Integrates the original 30-candidate gene set with previously generated tumor-only WGCNA results, including gene-identifier mapping, module assignments, module eigengene labels, MM values, ranks, percentiles, and existing TME association statistics. |
| `GitHub2_Paper1_tumor_only_WGCNA_9.R` | **Candidate-gene auditing and checkpoint preparation.** Audits existing candidate-level TME statistics and WGCNA mapping consistency, constructs a separate integrated candidate-level object, and preserves the pre-STEP19 workspace and selected locked results. |
| `GitHub2_Paper1_tumor_only_WGCNA_10.R` | **Supplementary figures, table, and workspace recovery.** Prepares and audits Supplementary Figures 4A–4D and Supplementary Table X using previously generated WGCNA results, and saves the final workspace checkpoint. |

The scripts should be interpreted in the context of their documented input requirements and upstream dependencies. Some later scripts resume from saved checkpoints rather than rerunning all preceding analyses.

For detailed step descriptions, input files, output inventories, and execution requirements, consult the corresponding documentation under `Scripts/` and `Results/`.

## 6. Main WGCNA Configuration and Module Characterization

The recorded primary WGCNA configuration includes the following parameters:

| Parameter                   | Setting                             |
|-----------------------------|-------------------------------------|
| Analysis population         | 513 TCGA-LUAD Primary Tumor samples |
| Soft-thresholding power (β) | 7                                   |
| Network type                | Unsigned                            |
| TOM type                    | Unsigned                            |
| Minimum module size         | 30                                  |
| Deep split                  | 3                                   |
| Module merge cut height     | 0.25                                |

These values describe the recorded primary network configuration. Additional preprocessing settings, correlation functions, module-detection implementation details, and any alternative analysis settings should be checked in the corresponding R scripts.

The tumor-only WGCNA analysis identified **50 co-expression modules** using the 513 Primary Tumor samples.

The module assignments of the three prioritized study genes are:

| Gene  | Tumor-only WGCNA module | Module eigengene |
|-------|-------------------------|------------------|
| AOC3  | Magenta                 | ME9              |
| MS4A7 | Blue                    | ME2              |
| IL7R  | Blue                    | ME2              |

The brown module (ME3) is retained as a survival-associated reference module in the subsequent survival-analysis workflows.

The tumor-only module colors and eigengene labels refer to this tumor-only network. In particular, the current tumor-only blue and magenta modules should not be confused with the historical blue and green modules from the earlier tumor-plus-normal analysis.

## 7. TME and Gene-Level Integration

The TME analysis evaluates relationships between WGCNA module eigengenes or individual genes and three ESTIMATE-derived traits:

- ImmuneScore
- StromalScore
- ESTIMATEScore

The module-level analysis evaluates associations between module eigengenes and the three TME traits. The gene-level integration workflows subsequently characterize the relationships between gene significance, module membership, and TME-related traits.

For the three prioritized study genes, the integration focuses on their existing module assignments and gene-level GS/MM characteristics:

- **AOC3:** Magenta module (ME9).
- **MS4A7:** Blue module (ME2).
- **IL7R:** Blue module (ME2).

The original 30-candidate gene set retains its previously established provenance. Tumor-only WGCNA is a subsequent characterization and integration analysis; it does not redefine the historical candidate set.

Existing GS, P-value, FDR, and locked MM results are inherited and audited in the relevant integration workflows where specified. Descriptive integration and consistency checks should not be interpreted as independent recalculation of the underlying statistical analyses.

The brown module (ME3) is treated separately as a survival-associated reference module and should not be conflated with the module assignments of the three prioritized study genes.

## 8. Supplementary Figures and Table

The final supplementary workflow prepares the following outputs using previously generated and locked WGCNA results:

- **Supplementary Figure 4A:** Module–trait association visualization.
- **Supplementary Figure 4B:** Module membership positioning of AOC3, MS4A7, and IL7R.
- **Supplementary Figure 4C:** Within-module associations between module membership and gene significance.
- **Supplementary Figure 4D:** Gene significance profiles of AOC3, MS4A7, and IL7R across TME traits.
- **Supplementary Table X:** WGCNA module characteristics and TME-related gene significance of the three prioritized study genes.

These supplementary outputs provide complementary descriptions of module-level and gene-level relationships. They are not independent clinical validation analyses.

## 9. Repository Navigation

The WGCNA directory separates analysis scripts, computational-environment documentation, input data, generated results, and result documentation.

| Directory or file | Purpose |
|----|----|
| `README.md` | Overall description of the tumor-only WGCNA analysis, data sources, analytical strategy, script sequence, and repository navigation. |
| `Environment/` | R version, package dependencies, software requirements, and reproducibility information. |
| `Scripts/` | The ten WGCNA analysis scripts and their detailed workflow documentation. |
| `data/` | Input expression data, clinical data, annotation references, sample manifests, and required intermediate inputs or checkpoints, as applicable. |
| `results/` | Generated analysis tables, statistical outputs, quality-control records, serialized objects, figures, and checkpoints, organized according to the actual repository structure. |
| `Results/` | Documentation and inventory of generated analysis outputs. |

Detailed documentation is provided in:

- `Environment/README.md` — computational environment and software dependencies.
- `Scripts/README.md` — script descriptions, analytical procedures, input requirements, and workflow dependencies.
- `Results/README.md` — result-directory organization and inventory of principal outputs.

Directory names and capitalization should be kept consistent with the actual repository structure. In particular, `Results/` is used for documentation, whereas `results/` stores generated analysis files.

## 10. Software Requirements

The recorded computational environment includes R 4.5.1 and packages used for network construction, survival analysis, data annotation, visualization, and statistical analysis.

The documented package families include:

- WGCNA
- dynamicTreeCut
- fastcluster
- survival
- survminer
- ggplot2
- ggpubr
- AnnotationDbi
- org.Hs.eg.db
- here
- ESTIMATE-related implementations and dependencies

The recorded ESTIMATE-related package versions include:

- `estimate`: 1.0.13
- `ESTIMATE`: 0.1.0

Some execution logs indicate that installed packages were built under R 4.5.2 or R 4.5.3, whereas the recorded R environment was R 4.5.1. These differences should be considered when reproducing the analysis.

For the complete package list, version records, and script-specific dependencies, consult `Environment/README.md` and the session information generated by the individual scripts.

Before public release, the actual R version and package versions used for the final execution should be checked against the recorded execution logs. A dependency file, such as `renv.lock`, may be added where feasible.

## 11. Running the Scripts

The general execution procedure is:

1.  Obtain the original input files described in Section 2, including the expression data from SangerBox and the relevant annotation and clinical data files.
2.  Place the input files in the designated project directories according to the repository structure and the paths specified in the scripts.
3.  Install the required R packages and verify the computational environment.
4.  Run the scripts in the documented workflow order, from Script 1 through Script 10, while respecting dependencies on intermediate objects and saved checkpoints.
5.  Verify the expected intermediate and final outputs before proceeding to downstream analyses.

Some later scripts depend on previously generated or saved R objects. They may therefore require the corresponding checkpoint files or execution of the preceding steps. The ten scripts should not be assumed to be independent, standalone scripts.

Project-relative paths are used to support reproducibility across different working environments. Users should inspect the script-specific README files for the exact input requirements and checkpoint dependencies before execution.

## 12. Reproducibility and Interpretation

The WGCNA workflows follow these general reproducibility principles:

- Document the original sources of downloaded input data.
- Preserve original input data and distinguish them from derived expression matrices.
- Use project-relative paths to locate inputs and organize outputs.
- Retain sample identifiers and sample order for cross-workflow verification.
- Record relevant R and package versions from the actual execution environment.
- Use explicit checkpoints to preserve important intermediate and final analysis objects.
- Distinguish new statistical analyses from descriptive audits, integration procedures, and workspace-recovery operations.
- Document analysis parameters and any changes to the sample or gene universe.
- Preserve the provenance of the original candidate-gene set.

The prioritized genes AOC3, MS4A7, and IL7R are exploratory study genes. Their WGCNA module membership and TME associations do not, by themselves, establish clinical diagnostic or prognostic biomarker validity.

The WGCNA and downstream analyses identify co-expression patterns and statistical associations. They should not be interpreted as evidence of direct molecular interactions, causal regulation, or clinical utility without additional supporting evidence.

The final repository documentation should be checked against the retained R scripts, actual input and output files, and runtime software-version records before public release.

## 13. Manuscript Correspondence

This directory supports the tumor-only WGCNA and related downstream analyses reported in the associated LUAD manuscript.

The main analytical components include:

- Tumor-only WGCNA and module identification.
- Module–clinical and module–TME association analyses.
- Survival modeling and time-dependent Cox assessments.
- Gene-level GS/MM integration and prioritized-gene localization.
- Supplementary Figures 4A–4D and Supplementary Table X.

The exact mapping between individual scripts, manuscript Methods and Results subsections, supplementary figure panels, and supplementary tables should be maintained in the final manuscript and repository documentation.

## 14. Citation

If you use these scripts or adapt the analytical workflow, please cite the associated manuscript when it becomes available.
