# Scripts — Immune Checkpoint Association and Sensitivity Analyses

## 1. Overview

This directory contains two R scripts implementing the immune checkpoint association and sensitivity analyses in lung adenocarcinoma (LUAD).

The analysis evaluates the associations between three prioritized study genes (**AOC3, MS4A7, and IL7R**) and 18 predefined immune checkpoint-related genes using bulk-tumor transcriptomic data from The Cancer Genome Atlas (TCGA).

The workflow comprises four main analytical components:

1.  **Unadjusted association analysis:** Spearman rank correlations between the three prioritized study genes and the 18 immune checkpoint-related genes.
2.  **ImmuneScore-adjusted analysis:** Partial Spearman correlations adjusted for ESTIMATE-derived ImmuneScore.
3.  **TumorPurity-adjusted sensitivity analysis:** Partial Spearman correlations adjusted for ESTIMATE-derived TumorPurity.
4.  **Leave-one-core-gene-out dependence audit:** Assessment of whether the observed association patterns depend disproportionately on any one of the three prioritized study genes.

The predefined analysis includes 54 gene-pair associations (3 core genes × 18 immune checkpoint-related genes), with 503 samples retained for each association-analysis layer.

The analyses characterize statistical associations in bulk-tumor expression data. They are not intended to establish direct molecular regulation, cell-intrinsic immune activity, causal relationships, or therapeutic effects.

------------------------------------------------------------------------

## 2. Script Organization

The workflow is implemented in two sequential R scripts.

### 2.1. GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1.R

**Scope:** STEP 8A through STEP 8C-8.

**Purpose:** Data-source inventory, sample and expression-data quality control, construction of the analysis-ready cohort, and primary immune checkpoint association analyses.

#### Covered workflow components

| Workflow | Description |
|----|----|
| STEP 8A | Data-source inventory and input-data inspection |
| STEP 8A-1 | Verification of tumor-only WGCNA source files |
| STEP 8A-2 | Strict sample ID consistency checks |
| STEP 8A-3 | Survival-variable diagnostics and clinical-data QC |
| STEP 8A-3B | Survival-variable completeness audit |
| STEP 8A-4 | Construction of the 503-sample complete survival cohort |
| STEP 8B | Expression-survival sample alignment and integrity checks |
| STEP 8B-12 | Saving the analysis-ready expression cohort |
| STEP 8C-0 to 8C-0-3 | Gene identifier and annotation-resource audits |
| STEP 8C-1 | Immune checkpoint panel mapping audit |
| STEP 8C-2 | Construction of the checkpoint expression matrix |
| STEP 8C-3 | Construction of the three-gene core expression matrix |
| STEP 8C-4 | Unadjusted Spearman association analysis |
| RECOVERY-A to RECOVERY-C | Recovery and alignment of previously saved ESTIMATE objects |
| STEP 8C-5A | ImmuneScore-adjusted partial Spearman analysis |
| STEP 8C-5B | Paired comparison of unadjusted and ImmuneScore-adjusted results |
| STEP 8C-6 | Manuscript-level result extraction and descriptive classification |
| STEP 8C-7A | Manuscript-level numerical definition and result lock |
| STEP 8C-7B | Generation of the 54-pair supplementary association table |
| STEP 8C-8 | Supplementary Figure 4 visualization and figure-panel preparation |

#### Main analytical operations

- Inspect source files, data structures, sample identifiers, and clinical variables.
- Diagnose survival-variable distributions and completeness.
- Construct the complete survival cohort and align it with the expression matrix.
- Audit gene identifier formats and annotation resources.
- Map the predefined 18-gene immune checkpoint panel.
- Extract expression values for the three prioritized study genes and 18 checkpoint-related genes.
- Verify gene identities, sample identities, ordering, and numerical integrity.
- Perform 54 unadjusted Spearman correlation tests.
- Apply Benjamini–Hochberg (BH) correction across the 54 comparisons.
- Restore previously saved ESTIMATE scores and align them with the 503-sample analysis cohort.
- Perform 54 ImmuneScore-adjusted partial Spearman correlation tests.
- Apply BH correction across the 54 adjusted comparisons.
- Compare unadjusted and ImmuneScore-adjusted association results.
- Extract and lock manuscript-level numerical summaries.
- Generate the supplementary association table and visualization data.
- Generate separate figure panels and associated QC outputs.

#### Main input files and data sources

The script uses the project-specific TCGA-LUAD expression, clinical, and previously generated WGCNA/ESTIMATE data.

Key input files and objects include:

- `Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt`
- `Workflow-04-2A_clinical_aligned.txt`
- STEP7B-7E-10 source files
- Previously saved ESTIMATE/WGCNA workspace objects
- Previously generated immune checkpoint association results, where applicable in downstream extraction and visualization steps

The exact local locations of input files are configured through the project directory structure and the corresponding script setup.

#### Main outputs

The principal outputs include:

- `Workflow-08-STEP8A4_clinical_survival_503.rds`
- Analysis-ready expression and checkpoint-gene matrices
- Unadjusted 54-pair Spearman association results
- ImmuneScore-adjusted 54-pair partial Spearman association results
- Paired association comparison and manuscript-level extraction objects
- STEP 8C-7A locked result-definition files
- STEP 8C-7B supplementary association table in RDS, CSV, and Excel formats
- STEP 8C-8 figure panels, plotting data, QC tables, and plotting objects

#### Reproducibility notes

- Project-relative paths are used.
- Previously locked association results are retained as the source of truth for downstream extraction and visualization.
- STEP 8C-5B and subsequent result-extraction procedures do not overwrite the original unadjusted or ImmuneScore-adjusted statistical results.
- STEP 8C-6, STEP 8C-7A, and STEP 8C-7B primarily perform result extraction, consistency checks, and result locking rather than introducing additional correlation analyses.
- STEP 8C-8 generates visual representations from the locked STEP 8C-7B results; it does not recalculate the underlying association statistics.
- No random seed is required for deterministic result extraction, QC, and visualization procedures. Any seed used in other analysis sections should be documented in the corresponding script.

------------------------------------------------------------------------

### 2.2. GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2.R

**Scope:** STEP 8D through STEP 8E.

**Purpose:** Independent tumor-purity sensitivity analysis and leave-one-core-gene-out dependence audit using the previously locked 54-pair association results.

#### Covered workflow components

| Workflow | Description |
|----|----|
| STEP 8D | Tumor-purity sensitivity analysis |
| STEP 8D-12 | Validation of expression, purity, and clinical-data alignment |
| STEP 8D-13 | TumorPurity-adjusted partial Spearman analysis |
| STEP 8D-14 | Three-layer paired comparison |
| STEP 8D-15 | Manuscript-ready sensitivity summary |
| STEP 8E Phase 1 | Construction and audit of leave-one-core-gene-out subsets |
| STEP 8E-05 | Extraction and standardization of locked significance and direction results |
| STEP 8E-06 | Numerical integrity QC and subset-level summary generation |
| STEP 8E-07 | Final QC evaluation and lock decision |

#### STEP 8D: TumorPurity-adjusted sensitivity analysis

The purpose of STEP 8D is to assess whether the predefined immune checkpoint association patterns are sensitive to adjustment for estimated tumor purity.

Main operations include:

- Load the previously saved Paper 1 workspace and locked STEP 8C-7B supplementary association results.
- Identify and validate the 54 predefined gene-pair associations.
- Restore the 503-sample ESTIMATE and clinical data.
- Establish the ESTIMATE sample-to-patient mapping.
- Derive the ESTIMATE-based TumorPurity covariate.
- Verify the expression data for the three core genes and 18 checkpoint-related genes.
- Confirm sample alignment, gene identity, ordering, and numerical completeness.
- Perform 54 TumorPurity-adjusted partial Spearman correlation tests.
- Apply BH correction across the 54 comparisons.
- Compare the unadjusted, ImmuneScore-adjusted, and TumorPurity-adjusted results.
- Summarize significance patterns, direction patterns, and core-gene-specific results.
- Perform numerical integrity checks and generate the manuscript-ready sensitivity summary.

The TumorPurity-adjusted analysis is a separate sensitivity analysis. It does not replace the primary unadjusted versus ImmuneScore-adjusted comparison in STEP 8C.

#### STEP 8E: Leave-one-core-gene-out dependence audit

STEP 8E evaluates the dependence of the locked 54-pair association pattern on the composition of the three predefined core genes.

Three subsets are constructed:

| Omitted core gene | Retained core genes | Number of associations |
|-------------------|---------------------|-----------------------:|
| AOC3              | MS4A7 and IL7R      |                     36 |
| MS4A7             | AOC3 and IL7R       |                     36 |
| IL7R              | AOC3 and MS4A7      |                     36 |

Main operations include:

- Validate the previously locked unadjusted, ImmuneScore-adjusted, and TumorPurity-adjusted result objects.
- Verify core-gene and checkpoint-gene identities.
- Match association records using `Association_ID`.
- Preserve the original locked association order.
- Standardize significance and correlation-direction fields.
- Summarize three-layer significance and direction patterns.
- Recheck global numerical counts and predefined direction-transition categories.
- Generate full-cohort and leave-one-core-gene-out summaries.
- Generate core-gene-level detailed audit results.
- Compare observed values with previously locked numerical definitions.
- Evaluate the final QC criteria and construct the final lock object.
- Export the consolidated audit results and final QC report.

**Analytical scope:** STEP 8E is a descriptive dependence and integrity audit of previously locked results. It does not perform new correlation tests, refit statistical models, recalculate P values, or introduce additional BH-FDR correction.

#### Main input files and objects

Key inputs include:

- `Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds`
- Previously saved Paper 1 workspace: `Paper1_END_WORKSPACE_20260920_20260920_015542.RData`
- Previously generated ESTIMATE, expression, clinical, and survival objects
- Locked STEP 8D-13 and STEP 8D-14 association results, where applicable

The locked STEP 8C-7B supplementary association results provide the reference for the downstream sensitivity and dependence audits.

#### Main outputs

STEP 8D-15 outputs:

- `STEP8D15_Manuscript_Sensitivity_Summary_LOCK.rds`
- `STEP8D15_Manuscript_Key_Findings.csv`
- `STEP8D15_QC.csv`
- `STEP8D15_FINAL_LOCK_REPORT.txt`

STEP 8E outputs:

- `STEP8E_FINAL_LOCK.rds`
- `STEP8E_54Pair_Locked_Audit.csv`
- `STEP8E_Manuscript_Key_Findings.csv`
- `STEP8E_Core_Gene_Detailed_Audit.csv`
- `STEP8E_Association_ID_QC.csv`
- `STEP8E_Global_Integrity_QC.csv`
- `STEP8E_QC.csv`
- `STEP8E_FINAL_LOCK_REPORT.txt`
- `STEP8E_Leave_One_Core_Gene_Out_Audit.xlsx` (generated when `openxlsx` is available)

#### Reproducibility notes

- Project-relative output paths are used.
- The locked STEP 8C-7B results are preserved as the reference dataset.
- The original 54-pair association definitions and sample size are retained.
- The TumorPurity-adjusted analysis is performed independently as a sensitivity analysis.
- STEP 8D-14 and STEP 8D-15 summarize and compare existing analytical layers.
- STEP 8E uses previously locked results and does not introduce new inferential tests or FDR recalculation.
- Excel export is conditional on the availability of `openxlsx`.
- Runtime execution and output-file existence should be verified separately from the final numerical lock status.

------------------------------------------------------------------------

## 3. Statistical Analysis Summary

### 3.1. Unadjusted Spearman correlation

Two-sided Spearman rank-correlation tests are performed for all 54 predefined gene pairs. P values are adjusted across the 54 comparisons using the Benjamini–Hochberg procedure.

### 3.2. ImmuneScore-adjusted partial Spearman correlation

Partial Spearman correlations are calculated using a rank-based residual procedure:

1.  Rank the expression values of each variable and the ImmuneScore.
2.  Regress each ranked expression variable on the ranked ImmuneScore.
3.  Calculate the Pearson correlation between the residuals for each gene pair.
4.  Apply BH correction across the 54 comparisons.

### 3.3. TumorPurity-adjusted partial Spearman correlation

A separate sensitivity analysis uses ESTIMATE-derived TumorPurity as the adjustment variable, applying the same rank-based residual-correlation procedure.

The same 54 predefined gene pairs and 503-sample cohort are retained, and P values are adjusted using BH correction across the 54 comparisons.

### 3.4. Leave-one-core-gene-out dependence audit

For each omitted core gene, the corresponding 36-pair subset is extracted from the previously locked results.

Significance classifications and correlation-direction patterns are summarized across the three analytical layers. No association statistics are recalculated.

------------------------------------------------------------------------

## 4. Key Analytical Results

The locked results from the three analytical layers are summarized below.

| Analytical layer or result                         | Number |
|----------------------------------------------------|-------:|
| Total predefined gene-pair associations            |     54 |
| Unadjusted associations with FDR \< 0.05           |     48 |
| ImmuneScore-adjusted associations with FDR \< 0.05 |     38 |
| TumorPurity-adjusted associations with FDR \< 0.05 |     39 |
| Associations significant in all three layers       |     27 |

The direction-transition patterns across the unadjusted, ImmuneScore-adjusted, and TumorPurity-adjusted analyses were:

| Direction pattern | Number |
|-------------------|-------:|
| P → N → N         |     30 |
| P → P → P         |     17 |
| P → P → N         |      4 |
| P → N → P         |      2 |
| N → N → N         |      1 |

These results indicate that the observed bulk-expression associations were sensitive to the adjustment strategy and were not uniformly consistent across analytical layers.

The leave-one-core-gene-out audit further characterizes how the composition of the predefined core-gene panel influences the summarized association patterns. Its results should be interpreted as a dependence audit rather than independent replication.

------------------------------------------------------------------------

## 5. Interpretation and Limitations

These analyses characterize statistical associations between bulk-tumor gene-expression measurements.

They do not establish:

- Direct molecular interactions between the studied genes.
- Direct regulation of immune checkpoint genes by AOC3, MS4A7, or IL7R.
- Cell-intrinsic regulatory activity.
- Functional immune effects or causal mechanisms.
- Clinical diagnostic performance or therapeutic response.

Bulk-tumor expression reflects a mixture of cellular populations and biological states. Adjustment for ImmuneScore or estimated TumorPurity does not eliminate all potential confounding by cellular composition.

Changes in correlation direction after adjustment should therefore be interpreted as sensitivity of the statistical associations to the adjustment strategy, rather than evidence of causal or mechanistic reversal.

The leave-one-core-gene-out analysis assesses the dependence of the summarized results on the predefined gene-pair composition. It is not an independent validation cohort or a new inferential analysis.

------------------------------------------------------------------------

## 6. Reproducibility and Execution

Both scripts document the workflow, input requirements, analytical procedures, quality-control checks, and output generation for the immune checkpoint association analysis.

To reproduce the workflow:

1.  Clone or download the project repository.
2.  Configure the required input data in the documented project directory structure.
3.  Confirm the availability and structure of the required intermediate files and saved workspace objects.
4.  Use compatible R and package versions.
5.  Execute the scripts in the intended order:
    - `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1.R`
    - `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2.R`
6.  Verify the generated output files and QC reports.

The scripts use project-relative path definitions where configured. Environment-specific input locations and required intermediate objects must be available before execution.

Exact R and package versions should be taken from the recorded execution environment rather than inferred from the script descriptions.

The workspace checkpoint generated by Script 1 is a session-preservation utility for continuing the Paper 1 revision workflow. It is not an independent statistical analysis and does not replace the reproducibility documentation of the individual analytical steps.
