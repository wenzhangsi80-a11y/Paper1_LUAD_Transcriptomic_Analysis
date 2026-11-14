# Paper 1 — Immune Checkpoint Association Analysis

## Overview

This directory contains the immune checkpoint association and sensitivity analyses for Paper 1.

The analyses are part of the study:

**Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R across the TME in LUAD**

This module evaluates the associations between three prioritized study genes (AOC3, MS4A7, and IL7R) and 18 predefined immune checkpoint-related genes in lung adenocarcinoma (LUAD), using bulk-tumor transcriptomic data from The Cancer Genome Atlas (TCGA).

The workflow comprises three association-analysis layers and a leave-one-core-gene-out dependence audit:

1.  **Unadjusted analysis:** Spearman rank correlations.
2.  **ImmuneScore-adjusted analysis:** Partial Spearman correlations adjusted for ESTIMATE-derived ImmuneScore.
3.  **TumorPurity-adjusted sensitivity analysis:** Partial Spearman correlations adjusted for ESTIMATE-derived TumorPurity.
4.  **Leave-one-core-gene-out dependence audit:** Assessment of the sensitivity of the locked association pattern to omission of each prioritized study gene.

The main objective is to characterize how bulk-tissue expression associations vary across different adjustment strategies. The results are interpreted as statistical associations rather than evidence of direct molecular interactions, causal regulation, or functional immune activity.

------------------------------------------------------------------------

## Research Objectives

The specific objectives are to:

1.  Establish and quality-control the expression and clinical data required for immune checkpoint association analyses.
2.  Construct an analysis-ready cohort of 503 TCGA-LUAD tumor samples.
3.  Prepare expression matrices for the three prioritized study genes and the predefined 18-gene immune checkpoint-related panel.
4.  Evaluate the unadjusted associations between the three study genes and the 18 checkpoint-related genes.
5.  Assess the sensitivity of these associations to adjustment for ESTIMATE-derived ImmuneScore.
6.  Perform an independent TumorPurity-adjusted sensitivity analysis.
7.  Compare the three association-analysis layers using consistent association identifiers.
8.  Evaluate the dependence of the overall association pattern on the composition of the three-gene panel.
9.  Generate manuscript-level tables, supplementary figures, and reproducibility-oriented quality-control outputs.

------------------------------------------------------------------------

## Analysis Design

### Study Genes and Immune Checkpoint Panel

The predefined association panel consists of:

- **Three prioritized study genes:** AOC3, MS4A7, and IL7R.
- **Eighteen immune checkpoint-related genes.**
- **Fifty-four predefined gene pairs:** 3 × 18.

The 18 immune checkpoint-related genes are:

``` text
PDCD1, CD274, PDCD1LG2, CTLA4, CD80, CD86,
LAG3, HAVCR2, TIGIT, CD96, VSIR, VTCN1,
BTLA, TNFRSF4, TNFRSF9, ICOS, LGALS9, CD200
```

The complete panel is documented in Supplementary Table 12.

The three study genes are treated as prioritized genes for downstream characterization. They are not presented as a formally ranked set or as protein–protein interaction network hubs.

### Analysis Cohort

The immune checkpoint association workflow uses 503 TCGA-LUAD primary tumor samples with the required gene-expression and adjustment-variable data.

| Item                                                 | Analysis design |
|------------------------------------------------------|----------------:|
| Original primary tumor expression cohort             |     513 samples |
| Survival-complete analysis cohort                    |     503 samples |
| Excluded samples during survival-cohort construction |              10 |
| Prioritized study genes                              |               3 |
| Immune checkpoint-related genes                      |              18 |
| Predefined gene pairs                                |              54 |
| Target sample size per association                   |             503 |

All three association-analysis layers retain the same 503-sample cohort and predefined gene-pair structure.

**Cohort distinction:** The separate gene–TME purity analysis elsewhere in the manuscript uses 500 primary tumors with complete expression, TME-score, and ABSOLUTE purity data. This is a different analysis cohort and should not be confused with the 503-sample immune checkpoint association cohort.

------------------------------------------------------------------------

## Workflow Structure

The analysis is implemented in two R scripts.

### Script 1: Primary Immune Checkpoint Association Analysis

**File:** `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1.R`

**Scope:** STEP 8A through STEP 8C-8.

This script covers data preparation, primary association analyses, result extraction, supplementary table generation, and visualization.

#### STEP 8A — Data Source Inventory and Quality Control

Establishes the input data sources and evaluates the integrity of the expression and clinical datasets.

Main procedures include:

- Identification of the tumor-only WGCNA expression data.
- Identification of the corresponding clinical data.
- Expression–clinical sample ID matching.
- Assessment of sample ID uniqueness and correspondence.
- Expression matrix numerical integrity checks.
- Clinical variable completeness and validity checks.
- Survival-variable quality control.

The survival-complete cohort is constructed using non-missing, finite OS.time values greater than or equal to zero and valid OS.status values.

The original tumor expression cohort contains 513 samples. The resulting survival-complete cohort contains 503 samples, with 10 samples excluded according to the predefined eligibility criteria.

#### STEP 8B — Expression–Survival Cohort Alignment

Aligns the expression matrix with the 503-sample survival cohort using explicit sample identifiers.

The workflow checks:

- Expression and clinical sample correspondence.
- Sample ordering and uniqueness.
- Gene and sample identifier integrity.
- Missing, non-finite, or invalid expression values.
- Preservation of the original expression measurements.

The aligned expression cohort is saved as an intermediate RDS object for downstream analyses.

#### STEP 8C — Primary Immune Checkpoint Association Analysis

**STEP 8C-0 — Gene Identifier and Annotation Audit**

Assesses the gene identifier namespace in the expression matrix and checks the availability of human gene annotation resources.

Ensembl-to-gene-symbol mapping is performed using the relevant annotation packages. Mapping coverage and ambiguous mappings are assessed separately from the original expression matrix.

**STEP 8C-1 — Immune Checkpoint Panel Mapping**

Evaluates the mapping and expression availability of the predefined 18-gene immune checkpoint-related panel.

**STEP 8C-2 — Checkpoint Expression Matrix**

Constructs and quality-checks the 18-gene × 503-sample checkpoint expression matrix while preserving the predefined gene and sample order.

**STEP 8C-3 — Core-Gene Expression Matrix**

Extracts AOC3, MS4A7, and IL7R expression values using their predefined Ensembl identifiers.

The intended matrix dimensions are 3 × 503.

**STEP 8C-4 — Unadjusted Spearman Correlation**

Performs 54 two-sided Spearman correlation tests between the three prioritized study genes and the 18 immune checkpoint-related genes.

P values are adjusted using the Benjamini–Hochberg (BH) procedure across all 54 predefined comparisons.

**ESTIMATE Recovery and Alignment**

Restores the previously generated ESTIMATE object from the existing WGCNA workspace and aligns the original 513-sample ESTIMATE scores to the 503-sample analysis cohort.

The locked ESTIMATE scores are restored rather than recalculated.

**STEP 8C-5 — ImmuneScore-Adjusted Partial Spearman Analysis**

Evaluates the sensitivity of the 54 associations to adjustment for ESTIMATE-derived ImmuneScore.

The analysis uses a rank-based residual procedure:

1.  Rank the expression values of both genes and the ImmuneScore.
2.  Regress each ranked gene-expression variable on the ranked ImmuneScore.
3.  Calculate the Pearson correlation between the residual vectors.
4.  Apply BH-FDR correction across the 54 tests.

STEP 8C-5B compares the locked unadjusted and ImmuneScore-adjusted results, including correlation changes, direction transitions, significance transitions, and core-gene-level summaries.

**STEP 8C-6 — Manuscript-Level Result Extraction**

Extracts and organizes the locked results from STEP 8C-4 and STEP 8C-5.

The extraction includes:

- All 54 predefined gene pairs.
- Unadjusted and ImmuneScore-adjusted association results.
- Direction-transition classifications.
- FDR-significance transitions.
- Overall and core-gene-level summaries.
- Internal consistency checks.

No new correlation analysis or checkpoint-gene selection is performed.

**STEP 8C-7A — Manuscript Result-Definition Lock**

Establishes the manuscript-level numerical definitions using the locked STEP 8C-6 result object.

This step verifies the association structure, sample size, significance counts, direction transitions, and core-gene-level summaries.

It also identifies the positive-to-negative associations that remain significant after ImmuneScore adjustment.

No new correlation or partial correlation analysis is performed.

**STEP 8C-7B — Supplementary Table Generation**

Generates the supplementary table for the 54 predefined associations.

The locked STEP 8C-6 result object is the source of truth. The workflow verifies the predefined association structure and manuscript-level numerical definitions before exporting the supplementary results in the available formats.

No new correlation analysis, FDR correction, or significance testing is performed.

**STEP 8C-8 — Supplementary Figure Visualization**

Generates three visualization panels from the locked STEP 8C-7B results.

- **Panel A:** Comparison of unadjusted and ImmuneScore-adjusted correlation coefficients, including direction transitions and highlighted associations.
- **Panel B:** Heatmap of ImmuneScore-adjusted partial Spearman correlations for the three study genes and 18 checkpoint-related genes.
- **Panel C:** Summary of the four direction-transition categories.

The panels are exported independently in PDF and 600-dpi TIFF formats with LZW compression.

The visualization workflow also generates plotting-data files and numerical QC outputs.

No correlation coefficients, P values, or FDR values are recalculated during figure generation.

### Script 2: TumorPurity Sensitivity Analysis and Dependence Audit

**File:** `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2.R`

**Scope:** STEP 8D through STEP 8E.

This script performs the TumorPurity-adjusted sensitivity analysis, compares the three association-analysis layers, and audits the dependence of the locked association pattern on the composition of the three-gene panel.

#### STEP 8D — TumorPurity-Adjusted Sensitivity Analysis

TumorPurity is derived from the previously generated ESTIMATE score using:

$$\mathrm{TumorPurity}
=
\cos(0.6049872018 + 0.0001467884 \times \mathrm{ESTIMATEScore})$$

**STEP 8D-13 — TumorPurity-Adjusted Partial Spearman Analysis**

Performs rank-based partial Spearman correlation analysis for the same 54 predefined gene pairs, adjusting for ESTIMATE-derived TumorPurity.

The analysis retains the 503-sample cohort and applies BH-FDR correction across the 54 comparisons.

This is a separate sensitivity analysis. It does not replace or overwrite the locked STEP 8C-7B primary results.

**STEP 8D-14 — Three-Layer Paired Comparison**

Compares the following locked association results using consistent Association_ID values:

1.  Unadjusted Spearman correlations.
2.  ImmuneScore-adjusted partial Spearman correlations.
3.  TumorPurity-adjusted partial Spearman correlations.

The comparison includes:

- Overall FDR-significance counts.
- Correlation-direction patterns.
- Significance and direction transitions.
- Three-layer intersections.
- Core-gene-level summaries.
- Internal consistency checks against the locked primary results.

**STEP 8D-15 — Manuscript-Level Sensitivity Summary**

Consolidates the three-layer comparison into manuscript-level summaries and numerical QC outputs.

The summary includes global significance counts, direction-pattern counts, and core-gene-level results.

No additional statistical testing or BH-FDR recalculation is performed in this summary step.

#### STEP 8E — Leave-One-Core-Gene-Out Dependence Audit

The purpose of STEP 8E is to evaluate whether the overall 54-pair association pattern is disproportionately dependent on any one of the three prioritized study genes.

Three predefined subsets are constructed:

| Subset      | Omitted gene | Retained study genes | Number of associations |
|-------------|--------------|----------------------|-----------------------:|
| Leave_AOC3  | AOC3         | MS4A7 and IL7R       |                     36 |
| Leave_MS4A7 | MS4A7        | AOC3 and IL7R        |                     36 |
| Leave_IL7R  | IL7R         | AOC3 and MS4A7       |                     36 |

**STEP 8E-05 — Extraction of Locked Significance and Direction Results**

Extracts and standardizes the previously locked significance and direction information from the 54 associations across the three analysis layers.

The results are paired using Association_ID, and the original locked association order is restored.

No new association tests or statistical model fitting are performed.

**STEP 8E-06 — Numerical Integrity and Dependence QC**

Evaluates the numerical consistency of the locked results and summarizes the full 54-pair dataset and the three 36-pair subsets.

The QC includes:

- Global significance-count verification.
- Three-layer direction-pattern verification.
- Direction and significance arithmetic checks.
- Preservation of the locked STEP 8C-7B and STEP 8D-15 results.
- Core-gene-level detailed summaries.
- Comparison with previously locked manuscript-level values.
- Subset-size and association-identifier integrity checks.

**Final QC and Result Lock**

The final stage consolidates the STEP 8E QC results and generates the final audit object and associated output files.

The outputs include RDS, CSV, TXT, and an optional Excel workbook when `openxlsx` is available.

The entire STEP 8E workflow is a descriptive dependence and numerical integrity audit. It does not recalculate correlation coefficients, P values, or FDR values, and does not constitute independent replication.

------------------------------------------------------------------------

## Statistical Methods

### Unadjusted Spearman Correlation

Two-sided Spearman rank-correlation tests are performed for all 54 predefined gene pairs.

P values are adjusted across the 54 comparisons using the Benjamini–Hochberg procedure to control the false discovery rate (FDR).

### ImmuneScore-Adjusted Partial Spearman Correlation

Partial Spearman correlations are calculated using a rank-based residual procedure:

1.  Rank the expression values of both genes and the ImmuneScore.
2.  Regress each ranked expression variable on the ranked ImmuneScore.
3.  Calculate the Pearson correlation between the residual vectors.
4.  Apply BH correction across the 54 gene-pair comparisons.

### TumorPurity-Adjusted Partial Spearman Correlation

The same rank-based residual procedure is used with ESTIMATE-derived TumorPurity as the adjustment variable.

The same 54 gene pairs and 503 samples are retained, and P values are adjusted using BH correction across the 54 comparisons.

This analysis assesses sensitivity to an alternative adjustment variable and is not a replacement for the primary ImmuneScore-adjusted analysis.

### Leave-One-Core-Gene-Out Audit

For each omission analysis, the corresponding 36-pair subset is extracted from the previously locked results.

Significance classifications and correlation-direction patterns are summarized across the unadjusted, ImmuneScore-adjusted, and TumorPurity-adjusted layers.

No new correlation coefficients, P values, or FDR values are calculated. No additional multiple-testing procedure is introduced.

The audit therefore characterizes dependence on the composition of the predefined gene-pair set rather than providing independent validation.

------------------------------------------------------------------------

## Main Results

The following values represent the locked manuscript-level reference results documented by the analysis workflow. Their reproduction in a new execution environment should be confirmed against the corresponding QC reports and output files.

### Unadjusted Associations

- 48 of 54 associations met FDR \< 0.05.
- 53 of 54 correlations were positive.
- 1 of 54 correlations was negative.

### ImmuneScore-Adjusted Associations

- 38 of 54 associations met FDR \< 0.05.
- 21 correlations were positive.
- 33 correlations were negative.

The comparison with unadjusted results showed the following direction transitions:

| Direction transition       | Number of associations |
|----------------------------|-----------------------:|
| Positive to Positive (P→P) |                     21 |
| Negative to Negative (N→N) |                      1 |
| Positive to Negative (P→N) |                     32 |
| Negative to Positive (N→P) |                      0 |
| **Total**                  |                 **54** |

Among the 32 positive-to-negative associations, 21 remained significant after ImmuneScore adjustment at FDR \< 0.05.

### TumorPurity-Adjusted Associations

- 39 of 54 associations met FDR \< 0.05.

### Three-Layer Comparison

The locked manuscript-level sensitivity summary is:

| Analysis layer                        | FDR-significant associations |
|---------------------------------------|-----------------------------:|
| Unadjusted Spearman                   |                        48/54 |
| ImmuneScore-adjusted partial Spearman |                        38/54 |
| TumorPurity-adjusted partial Spearman |                        39/54 |
| Significant in all three layers       |                        27/54 |

The three-layer direction patterns are:

| Direction pattern (Unadjusted → ImmuneScore → TumorPurity) | Number |
|------------------------------------------------------------|-------:|
| P→N→N                                                      |     30 |
| P→P→P                                                      |     17 |
| P→P→N                                                      |      4 |
| P→N→P                                                      |      2 |
| N→N→N                                                      |      1 |
| **Total**                                                  | **54** |

These results indicate that the observed bulk-expression associations are sensitive to the adjustment strategy and are not uniformly consistent across analytical layers.

### Core-Gene-Level Three-Layer Summary

The locked numbers of FDR-significant associations for each study gene are:

| Core gene | Unadjusted | ImmuneScore-adjusted | TumorPurity-adjusted | Significant in all three |
|----|---:|---:|---:|---:|
| AOC3 | 14/18 | 13/18 | 17/18 | 8 |
| MS4A7 | 17/18 | 12/18 | 11/18 | 10 |
| IL7R | 17/18 | 13/18 | 11/18 | 9 |

### Leave-One-Core-Gene-Out Audit

The dependence audit retains the corresponding 36 associations after omitting each study gene.

The locked core-gene-level reference summaries are:

| Core gene | Unadjusted significant | ImmuneScore-adjusted significant | TumorPurity-adjusted significant | Significant in all three |
|----|---:|---:|---:|---:|
| AOC3 | 14/18 | 13/18 | 17/18 | 8 |
| MS4A7 | 17/18 | 12/18 | 11/18 | 10 |
| IL7R | 17/18 | 13/18 | 11/18 | 9 |

These values describe the original full 54-pair results by core gene. The leave-one-core-gene-out subsets are separately summarized in the STEP 8E audit outputs.

The audit is intended to characterize the dependence of the predefined association pattern on the composition of the gene-pair panel. It should not be interpreted as independent replication or as a new set of statistical tests.

------------------------------------------------------------------------

## Directory Structure

The module is organized as follows:

``` text
05_Immune_Checkpoint/
│
├── README.md
│
├── Scripts/
│   ├── README.md
│   ├── GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1.R
│   └── GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2.R
│
├── Environment/
│   └── README.md
│
└── Results/
    └── README.md
```

The R scripts generate or verify statistical outputs, supplementary tables, figures, and QC reports under the project-defined `results/` directory.

Detailed file names, output descriptions, and result-level QC information are documented in `Results/README.md`.

Software requirements and package information are documented in `Environment/README.md`.

Script organization and workflow execution details are documented in `Scripts/README.md`.

------------------------------------------------------------------------

## Key Input Data and Intermediate Objects

The primary input data and intermediate objects include:

| File or object | Purpose |
|----|----|
| `Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt` | Tumor-only WGCNA expression matrix |
| `Workflow-04-2A_clinical_aligned.txt` | Clinical data aligned with the expression cohort |
| `Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData` | Previously generated WGCNA/ESTIMATE workspace used for recovery of locked ESTIMATE scores |
| `Paper1_END_WORKSPACE_20260920_20260920_015542.RData` | Workspace checkpoint used in the STEP 8D workflow |
| `Workflow-08-STEP8B_expr_survival_503.rds` | Expression matrix aligned to the 503-sample survival cohort |
| `Workflow-08-STEP8C2_checkpoint_expr_18x503.rds` | Immune checkpoint-related gene expression matrix |
| `Workflow-08-STEP8C3_core_expr_3x503.rds` | Three prioritized study-gene expression matrix |
| `Workflow-08-STEP8C4_core_checkpoint_Spearman_54tests.rds` | Locked unadjusted association results |
| `Workflow-08-STEP8C5_core_checkpoint_partialSpearman_ImmuneScore_54tests.rds` | ImmuneScore-adjusted association results |
| `Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds` | Locked primary 54-pair supplementary association results |

The exact input paths are defined in the R scripts using project-relative path management.

------------------------------------------------------------------------

## Reproducibility and Quality Control

The workflow incorporates quality-control procedures at the data preparation, association-analysis, result-extraction, visualization, and final audit stages.

Key reproducibility principles include:

1.  **Project-relative paths:** Project paths are managed using `here::here()` and `file.path()` rather than machine-specific absolute paths.

2.  **Input preservation:** Original expression, clinical, and locked intermediate data are retained without overwriting the source files.

3.  **Sample ID alignment:** Expression, clinical, and ESTIMATE-derived variables are matched using explicit sample identifiers.

4.  **Predefined gene-pair structure:** The three prioritized study genes and 18 checkpoint-related genes define the same 54 associations across the three analysis layers.

5.  **Locked primary results:** The STEP 8C-7B supplementary association object is retained as the reference source for downstream sensitivity analyses and visualization.

6.  **Separate sensitivity analysis:** TumorPurity adjustment is performed as a separate sensitivity analysis and does not replace the primary ImmuneScore-adjusted analysis.

7.  **No redundant statistical recalculation:** Manuscript-level extraction, visualization, and the leave-one-core-gene-out audit use previously locked results rather than recalculating the original associations.

8.  **Execution-dependent QC:** Expected dimensions, significance counts, and direction patterns are reference values. Successful execution and final QC status must be confirmed using the generated output files and runtime reports.

9.  **Software documentation:** R, RStudio (if applicable), operating-system, and package versions should be recorded from the actual execution environment.

10. **Conditional Excel export:** Excel output from STEP 8E depends on the availability of `openxlsx`. The corresponding CSV, RDS, and TXT outputs are handled separately.

No random seed is specified for the deterministic procedures documented in these two scripts.

------------------------------------------------------------------------

## Interpretation and Limitations

These analyses characterize associations between bulk-tumor gene-expression measurements. They do not establish:

- Direct molecular interactions between the study genes and checkpoint-related genes.
- Direct regulation of immune checkpoint expression.
- Cell-intrinsic regulatory activity.
- Functional immune effects.
- Causal relationships.
- Therapeutic response or clinical predictive utility.

Bulk-tumor expression reflects a mixture of cellular populations and cellular states. Adjustment for ImmuneScore or estimated TumorPurity does not eliminate all potential confounding by cellular composition.

Changes in correlation magnitude, statistical significance, or direction after adjustment should therefore be interpreted as sensitivity of the statistical associations to the adjustment strategy, not as evidence of a causal or mechanistic reversal.

The leave-one-core-gene-out audit assesses dependence on the composition of the predefined gene-pair panel. It is not an independent validation study.

------------------------------------------------------------------------

## Manuscript Mapping

The analyses support the following manuscript sections and supplementary materials:

| Manuscript component | Corresponding analysis |
|----|----|
| Methods, Section 2.9 | Immune checkpoint association and sensitivity analyses |
| Results, Section 3.7 | Immune checkpoint correlations and sensitivity to bulk-tissue composition adjustment |
| Supplementary Table 4 | Three-layer sensitivity analysis of immune checkpoint-associated expression patterns |
| Supplementary Table 5 | Leave-one-core-gene-out dependence audit |
| Supplementary Table 12 | Complete list of the 18 immune checkpoint-related genes |
| Supplementary Figure | Unadjusted versus ImmuneScore-adjusted associations and adjusted association patterns |

**Figure numbering note:** The STEP 8C-8 scripts and their output filenames refer to **Supplementary Figure 4**, whereas the earlier manuscript-level summary referred to **Supplementary Figure 6**. The final figure number must be reconciled with the current manuscript and the final supplementary-material sequence before submission. The code-generated filenames should not be assumed to establish the final manuscript numbering.

The TumorPurity-adjusted results and three-layer comparisons are summarized separately from the primary visualization, while the leave-one-core-gene-out summaries are provided through the dependence-audit supplementary outputs.

------------------------------------------------------------------------

## Documentation

This module contains three supporting documentation files:

- `Scripts/README.md`: Description of the two R scripts, their workflow components, and execution organization.
- `Environment/README.md`: R environment, package dependencies, input/output directory conventions, and computational requirements.
- `Results/README.md`: Detailed statistical outputs, result summaries, QC references, supplementary tables, and figure descriptions.

This root README provides the overall research context, analytical design, workflow structure, principal findings, and manuscript mapping for the complete immune checkpoint association module.
