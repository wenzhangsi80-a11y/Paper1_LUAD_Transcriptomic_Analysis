
# Tumor-only WGCNA: Analysis Scripts and Workflow Guide

## 1. Overview

This directory contains the R scripts used for the tumor-only weighted gene co-expression network analysis (WGCNA) in the Paper 1 LUAD transcriptomic study.

The analysis is organized into ten sequential R scripts. Together, these scripts cover:

- Clinical cohort reconstruction and expression-data preparation.
- Expression quality control and gene-variability assessment.
- WGCNA network construction and module identification.
- Module-level clinical and survival analyses.
- Time-dependent Cox regression analyses.
- ESTIMATE-based tumor microenvironment (TME) characterization.
- Gene-level module membership (MM) and gene significance (GS) analyses.
- Candidate-gene integration and mapping audits.
- Supplementary figure and table preparation.
- Workspace recovery and reproducibility archiving.

The tumor-only WGCNA analysis is based on 513 TCGA-LUAD Primary Tumor samples identified using TCGA sample type code 01.

The scripts are designed to use project-relative paths. The project root is identified using `here::here()`, and input and output paths are constructed using `file.path()`.

## 2. Script Inventory and Analytical Stages

The ten scripts are grouped according to their analytical roles rather than simply being listed as independent files.

| Analytical stage | Script | Main analytical purpose |
|---|---|---|
| 1. Cohort reconstruction and expression-data preparation | Script 1 | Reconstruct the clinical cohort, prepare the primary-tumor expression matrix, perform expression and annotation audits, apply expression-prevalence filtering, assess gene variability, and prepare the WGCNA input matrix. |
| 2. WGCNA network construction and module-level analyses | Script 2 | Construct the co-expression network, identify modules, evaluate prioritized study-gene module membership, and assess module associations with clinical traits and overall survival. |
| 3. Time-dependent survival analyses | Scripts 3–5 | Investigate proportional-hazards assumptions and evaluate time-dependent effects for the ME_brown module and selected clinical variables. |
| 4. TME association and gene-level integration | Scripts 6–7 | Calculate ESTIMATE-derived TME scores, evaluate module–TME associations, and integrate gene-level MM and GS results for prioritized study genes and other candidate genes. |
| 5. Candidate-gene integration, integrity audits, and archiving | Scripts 8–9 | Integrate the original candidate-gene set with locked WGCNA results, audit gene and module mappings, investigate historical provenance, and preserve analysis checkpoints. |
| 6. Supplementary results and workspace recovery | Script 10 | Prepare Supplementary Figures 4A–4D and Supplementary Table X, recover selected archived objects, and save the final workspace checkpoint. |

## 3. Detailed Script Descriptions

### 3.1 Script 1: Clinical Cohort Reconstruction and WGCNA Input Preparation

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_1.R`

**Purpose**

This script prepares the clinical and expression datasets required for downstream tumor-only WGCNA. It performs cohort reconstruction, expression-data quality control, annotation provenance assessment, expression-prevalence filtering, and gene-variability evaluation.

**Main analytical components**

| Workflow | Main task |
|---|---|
| Workflow-04-1 | Reconstruct the TCGA-LUAD Primary Tumor clinical cohort, clean clinical variables, audit survival information, and generate the clinical sample manifest. |
| Workflow-04-2A | Extract the 513 Primary Tumor expression samples and align clinical information with the expression matrix. |
| Workflow-04-2B | Audit expression-matrix structure, sample matching, numeric integrity, and global, sample-level, and gene-level expression distributions. |
| Workflow-04-2B-1 | Investigate extreme-expression observations and sample-level expression distributions. |
| Workflow-04-2B-2 | Provide supplementary annotation of extreme-expression features using `org.Hs.eg.db` and `AnnotationDbi`. |
| Workflow-04-2B-2A | Verify Ensembl gene identifiers and annotation provenance against GENCODE v22. |
| Workflow-04-2C | Apply the prespecified expression-prevalence filter and generate the log2-transformed expression matrix. |
| Workflow-04-2C-1 | Audit GENCODE v22 annotation provenance and cross-source annotation consistency for candidate genes. |
| Workflow-04-2C-2 | Calculate gene-level variability statistics, including variance, standard deviation, and median absolute deviation (MAD). |
| Workflow-04-2C-3 | Characterize MAD = 0 and MAD > 0 gene groups and extract the MAD-positive expression matrix. |
| Workflow-04-2D | Perform pre-analysis expression-matrix and sample-level quality control before network construction. |

**Principal input files**

- `data/Merge_RNA_seq_FPKM_UQ.txt`
- Original TCGA clinical data file stored in `data/`
- `data/reference/GENCODE_v22/gencode.v22.annotation.gtf.gz`

Additional inputs are generated by preceding workflows within the same script.

**Key processing specifications**

- Expected cohort: 513 TCGA Primary Tumor samples (sample type code 01).
- Expression-prevalence criterion: original FPKM-UQ > 1 in at least 103 of 513 samples.
- Expression transformation: `log2(FPKM-UQ + 1)`.
- Gene variability assessment: MAD calculated using `stats::mad()` with constant = 1.4826.
- The MAD-positive expression matrix is retained for subsequent network construction.

**Principal outputs**

Results are organized under:

`results/02_tumor_only_WGCNA/Workflow-04-*/`

Figures are organized under:

`figures/02_tumor_only_WGCNA/Workflow-04-2D/`

The principal downstream input is:

`Workflow-04-2D_WGCNA_input_MAD_positive_samples_by_genes.txt`

This matrix contains 513 samples in rows and the complete MAD-positive gene set in columns.

**Relationship to the next stage**

Script 2 uses the verified MAD-positive expression matrix and the corresponding upstream clinical and gene-variability information to construct the WGCNA network.

**Important scope clarification**

Workflow-04-1 through Workflow-04-2D perform data preparation, annotation audits, and pre-analysis quality control. They do not construct the WGCNA network or automatically remove samples based on descriptive QC results.

### 3.2 Script 2: WGCNA Network Construction and Module-Level Clinical and Survival Analyses

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_2.R`

**Purpose**

Construct a tumor-only WGCNA network using the 513 primary TCGA-LUAD tumor samples, identify co-expression modules, evaluate the module membership of AOC3, MS4A7, and IL7R, and investigate module-level associations with clinical characteristics and overall survival.

**Principal input files**

- `data/Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData`
- `data/Workflow-04-2D_WGCNA_input_checkpoint_DONE.txt`

The upstream checkpoint is expected to contain the expression matrix, aligned clinical data, and gene-level MAD information.

**Main analytical steps**

1. Validate upstream checkpoint files and required analysis objects.
2. Verify expression-matrix dimensions, sample identifiers, gene identifiers, sample order, and data integrity.
3. Evaluate candidate soft-thresholding powers using `WGCNA::pickSoftThreshold()`.
4. Construct the co-expression network using `WGCNA::blockwiseModules()`.
5. Identify gene modules and calculate module eigengenes.
6. Annotate and evaluate the module membership (MM) of AOC3, MS4A7, and IL7R.
7. Construct clinical traits, including age, gender, pathological stage, and TNM variables.
8. Evaluate module–clinical trait associations using the prespecified statistical tests.
9. Perform univariable module eigengene–overall survival Cox regression.
10. Assess the proportional-hazards (PH) assumption for selected modules.
11. Fit primary multivariable Cox models adjusted for age, gender, and pathological stage.
12. Prepare and audit sensitivity-analysis cohorts based on TNM variables.

**Principal analysis specifications**

| Parameter | Setting |
|---|---|
| Sample cohort | 513 Primary Tumor samples |
| Selected soft-threshold power | 7 |
| Network type | Unsigned |
| TOM type | Unsigned |
| Minimum module size | 30 |
| Deep split | 3 |
| Module merge cut height | 0.25 |
| WGCNA threading | Disabled |
| Random seed | 20261106 |

**Principal outputs**

Results and figures are organized under:

- `results/02_tumor_only_WGCNA/`
- `figures/02_tumor_only_WGCNA/`

The script generates network and module results, module–trait association tables, module-level survival results, PH diagnostics, input-validation reports, software-version records, and explicit RDS checkpoints.

The detailed output inventory is documented in `results/README.md`.

**Relationship to the next stage**

The module-level survival results and associated checkpoints provide the starting point for the additional ME_brown and clinical-variable time-dependent survival analyses in Scripts 3–5.

### 3.3 Script 3: ME_brown Time-dependent Cox Analysis

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_3.R`

**Purpose**

Continue the ME_brown survival analysis from an existing multivariable Cox and PH checkpoint. Investigate the time-dependent behavior of ME_brown and evaluate candidate time-dependent Cox models.

**Principal input**

`data/Workflow-04-2N_WGCNA_PH_multivariable_checkpoint_20260911.RData`

**Main analytical steps**

1. Inspect the existing PH test results and scaled Schoenfeld residuals.
2. Explore associations between Schoenfeld residuals and candidate time transformations.
3. Fit candidate time-dependent Cox models.
4. Compare conventional and time-dependent Cox models.
5. Estimate time-specific hazard ratios (HRs) and 95% confidence intervals.
6. Reparameterize ME_brown as a one-standard-deviation increase.
7. Evaluate coefficient and variance-covariance matrix stability.
8. Save the final model checkpoint and reporting table.

**Principal outputs**

Outputs are saved under:

`results/02_tumor_only_WGCNA/ME_brown_time_dependent_Cox/`

The detailed output filenames are documented in `results/README.md`.

**Relationship to the next stage**

The ME_brown analysis provides the context for subsequent investigation of clinical covariates, particularly the PH assumption and potential time-dependent effect of gender in Script 4.

**Important scope clarification**

This script resumes an existing analysis checkpoint. It does not reconstruct the original expression matrix, WGCNA network, or clinical model design.

### 3.4 Script 4: Gender-specific Time-dependent Cox Analysis

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_4.R`

**Purpose**

Investigate the proportional-hazards assumption for gender in the ME_brown multivariable Cox analysis and evaluate candidate time-dependent effects for male versus female patients.

**Principal input**

`data/Paper1_ME_brown_checkpoint_Step10G.RData`

**Main analytical steps**

1. Load the existing ME_brown analysis checkpoint.
2. Assess the PH assumption using `survival::cox.zph()`.
3. Inspect scaled Schoenfeld residuals for gender.
4. Verify gender-variable coding.
5. Construct a binary Male indicator, using Female as the reference group.
6. Fit linear, log1p(time), and square-root time-dependent Cox models.
7. Compare candidate models using AIC and log-likelihood.
8. Check model sample sizes and event counts.
9. Estimate time-specific HRs for male versus female at 1–5 years.
10. Evaluate coefficient and variance-covariance matrix stability.
11. Save the final model checkpoint and HR table.

**Principal outputs**

Outputs are saved under:

`results/02_tumor_only_WGCNA/Gender_time_dependent_Cox/`

**Relationship to the next stage**

The finalized gender-related survival results and checkpoint provide input to Script 5, which evaluates the PH assumption for pathological stage and defines the final integrated time-dependent Cox model.

### 3.5 Script 5: Stage-specific PH Assessment and Final Integrated Survival Model

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_5.R`

**Purpose**

Evaluate the PH assumption for pathological stage, assess a potential time-dependent Stage III effect, and summarize the final integrated survival model incorporating the brown-module eigengene and clinical covariates.

**Principal input**

`data/Paper1_ME_brown_Gender_checkpoint_Final.RData`

**Main analytical steps**

| Step | Main task |
|---|---|
| 12A–12A-3 | Examine stage distribution, missingness, survival events, and construct the complete-case stage analysis dataset. |
| 12B–12B-12 | Assess the PH assumption and perform binary stage diagnostic models using Stage I as the reference. |
| 12C | Fit and compare linear, log1p(time), and square-root time-dependent Stage III models. |
| 12D | Evaluate the incremental contribution of the Stage III time-dependent effect. |
| 13A–13C | Define and summarize the final integrated survival model. |
| 13D–13F | Calculate time-specific HRs for the brown-module eigengene, gender, and Stage III. |
| 13G | Save the final model checkpoint and selected WGCNA objects. |

**Final integrated model**

The selected model is based on the square-root time-dependent Stage III candidate model and incorporates:

- ME_brown and its time-dependent term.
- Age.
- Gender and its time-dependent term.
- Pathological stage, including a time-dependent Stage III term.

**Principal outputs**

The final checkpoint is stored under:

`results/02_tumor_only_WGCNA/Step13G_FINAL_checkpoint/`

Additional model summaries and HR tables are documented in `results/README.md`.

**Relationship to the next stage**

The finalized WGCNA and survival-analysis objects provide upstream information for the module–TME association workflow in Script 6.

**Important scope clarification**

Script 5 continues the existing survival-analysis workflow. It does not independently reconstruct the original tumor cohort or repeat WGCNA network construction.

### 3.6 Script 6: ESTIMATE-based TME Scoring and Module–TME Association Analysis

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_6.R`

**Purpose**

Calculate ESTIMATE-derived tumor microenvironment scores in primary LUAD tumor samples and evaluate their associations with tumor-only WGCNA module eigengenes.

**Principal input files**

- `Workflow-04-2C_expr_primary_tumor_log2_FPKM_UQ.txt`
- `Workflow-04-2C-1A_GENCODE_v22_complete_annotation.txt`

The script also requires previously generated WGCNA and clinical-analysis objects, including module assignments, eigengenes, expression data, and finalized clinical-analysis objects.

**Main analytical steps**

| Step | Main task |
|---|---|
| 14A | Inspect upstream WGCNA objects and perform input-data quality control. |
| 14B | Map Ensembl identifiers to gene symbols and prepare the ESTIMATE expression matrix. |
| 14C | Prepare ESTIMATE-compatible input and calculate TME scores. |
| 14D | Evaluate Spearman correlations between WGCNA module eigengenes and TME scores. |
| 14E | Verify module color-to-eigengene mapping and validate the brown module. |
| Final save | Save locked analysis objects, result tables, QC reports, and session information. |

**Main analytical methods**

- GENCODE v22-based Ensembl-to-gene-symbol mapping.
- Duplicate gene-symbol resolution based on median expression.
- ESTIMATE common-gene filtering.
- Calculation of StromalScore, ImmuneScore, and ESTIMATEScore.
- Spearman correlation analysis between module eigengenes and TME scores.
- Benjamini–Hochberg false discovery rate (FDR) correction across module–trait tests.
- Validation of module color-to-eigengene correspondence.

**Principal outputs**

Intermediate files:

`results/Workflow-03/Step14/Intermediate/`

Final locked outputs:

`results/Workflow-03/Step14/Step14_FINAL_LOCKED/`

The final directory contains RData objects, tables, QC reports, scripts, and session-information records.

**Relationship to the next stage**

The locked module–TME association results and upstream WGCNA objects are used by Script 7 for gene-level MM and GS analyses and core-gene positioning.

**Important scope clarification**

This script depends on upstream WGCNA and clinical-analysis objects. It should not be considered a fully standalone workflow unless the required objects have been loaded or generated.

### 3.7 Script 7: Gene-level MM–GS Analysis and Core-gene Integration

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_7.R`

**Purpose**

Continue the locked tumor-only WGCNA workflow by integrating module-level TME associations with gene-level module membership and gene significance results. The analysis focuses on the positioning of AOC3, MS4A7, and IL7R within their assigned modules.

**Principal input**

Previously generated and locked tumor-only WGCNA objects, including:

- Module assignments and module eigengenes.
- Module–TME association results.
- Gene-level module membership information.
- Expression data and clinical/TME scores.
- Previously generated core-gene and candidate-gene information.

**Main analytical steps**

1. Identify the module assignments of the prioritized study genes.
2. Integrate predefined module-level TME associations.
3. Verify clinical and survival cohort information.
4. Assess gene-level module membership (MM).
5. Calculate gene significance (GS) for ImmuneScore, StromalScore, and ESTIMATEScore.
6. Evaluate module-level MM–GS correlations.
7. Calculate within-module MM rankings and percentiles.
8. Integrate gene-level GS and module-level association information.
9. Prepare final core-gene and module-level summary tables.
10. Save integrated analysis objects and CSV tables for downstream visualization.

**Core-gene module assignments reported in the workflow**

| Prioritized study gene | Module assignment |
|---|---|
| AOC3 | Magenta / ME9 |
| MS4A7 | Blue / ME2 |
| IL7R | Blue / ME2 |

Brown / ME3 is retained as a survival-associated reference module.

**Principal outputs**

The script generates integrated gene-level and module-level result objects and summary tables for downstream candidate integration and supplementary visualization.

The detailed output inventory is documented in `results/README.md`.

**Relationship to the next stage**

The integrated gene-level and module-level results are used by Scripts 8 and 9 for broader candidate-gene integration, mapping-integrity audits, and preservation of locked analysis results.

**Important interpretation note**

Gene-level GS associations and module-level TME associations are distinct analytical quantities. The analyses evaluate statistical associations and do not establish causal relationships.

### 3.8 Script 8: Original Candidate-gene Integration and Historical Provenance Audit

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_8.R`

**Purpose**

Integrate the original 30 candidate genes with previously generated and locked tumor-only WGCNA results. Perform candidate-gene integrity audits, investigate historical candidate provenance, and archive the pre-STEP19 analysis workspace and relevant files.

**Principal inputs**

Previously generated WGCNA and candidate-gene objects, including:

- Original 30-candidate gene information.
- Locked core-gene MM values, rankings, and percentiles.
- Existing gene-level GS, P-value, and FDR results.
- Module assignments and eigengene information.
- Previously generated candidate integration tables.
- Existing scripts and archived RData/RDS checkpoints.

**Main analytical steps**

1. Initialize project-relative directories and load required upstream objects.
2. Integrate the original 30 candidate genes with tumor-only WGCNA module assignments.
3. Evaluate candidate-gene representation and module membership.
4. Audit gene-level TME association results using previously calculated GS, P, and FDR values.
5. Verify SYMBOL-to-ENSEMBL mapping and module color-to-eigengene correspondence.
6. Compare candidate-gene and core-gene information with locked upstream results.
7. Construct a separate integrated candidate-level object while preserving the original candidate object.
8. Investigate historical candidate provenance using project files, scripts, checkpoints, and related archived records.
9. Save a pre-STEP19 workspace checkpoint and object inventory.
10. Archive specified analysis scripts and RData/RDS files and generate archive manifests.

**Principal outputs**

- Candidate-level integration and audit tables.
- WGCNA mapping-integrity reports.
- Historical candidate-provenance audit records.
- Pre-STEP19 workspace checkpoint.
- Selected locked-object files.
- Script and RData/RDS archive manifests.

The detailed output locations and filenames are documented in `results/README.md`.

**Relationship to the next stage**

The preserved candidate-level information, mapping audits, and archived checkpoints provide the basis for the additional candidate-gene integrity and integration checks in Script 9 and subsequent supplementary result preparation.

**Important scope clarification**

This script primarily performs integration, auditing, and archiving. It does not reconstruct the WGCNA network or repeat the underlying DEG, WGCNA, GS, P-value, or FDR calculations.

### 3.9 Script 9: Candidate-gene and WGCNA Integration Integrity Audit

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_9.R`

**Purpose**

Perform additional audits of the original candidate-gene set, WGCNA module assignments, and previously generated TME association results. Preserve and verify the integrated analysis state before subsequent supplementary result preparation.

**Principal analytical components**

| Section | Main task |
|---|---|
| STEP 18A | Audit the locked 30-candidate gene-level TME association results, including existing GS, P-value, and FDR information. |
| STEP 18B | Audit WGCNA module and MM mapping integrity using existing analysis objects. |
| STEP 18C | Perform Ensembl-based module membership mapping and integrity checks, including independent recovery of the prioritized study-gene module assignments. |
| STEP 18D | Integrate existing WGCNA module and MM information with the candidate-level results while preserving previously calculated TME statistics. |
| STEP 18F | Create a pre-STEP19 workspace checkpoint, save session information, inventory important objects, and archive specified scripts and RData/RDS files. |

**Principal inputs**

Previously generated and locked analysis objects, including:

- The original 30-candidate gene set.
- Existing candidate-level TME association results.
- `moduleColors` and module membership information.
- Locked core-gene integration results.
- Existing WGCNA and TME analysis objects.
- Project scripts and archived checkpoints.

**Principal outputs**

- Candidate-level TME association audit tables.
- WGCNA module/MM mapping-integrity audit results.
- Ensembl-based module assignment audit results.
- Candidate integration and backfill objects.
- Pre-STEP19 workspace checkpoint.
- Archived script and RData/RDS files with corresponding manifests.

**Relationship to the next stage**

The locked and audited analysis objects, together with the archived STEP 19 intermediate results, are used by Script 10 to prepare supplementary figures and tables and recover the final analysis workspace.

**Important scope clarification**

- Existing GS, P-value, and FDR values are inherited rather than recalculated.
- No new WGCNA network construction or module-level statistical analysis is performed.
- The uploaded script contains STEP 18A, 18B, 18C, 18D, and 18F; no STEP 18E section is included.
- STEP 18F preserves the existing analytical state and does not recalculate analytical results.

### 3.10 Script 10: Supplementary Figures, Supplementary Table, and Final Workspace Recovery

**Script file:** `GitHub2_Paper1_tumor_only_WGCNA_10.R`

**Purpose**

Use previously generated and locked tumor-only WGCNA objects to prepare Supplementary Figures 4A–4D and Supplementary Table X. Recover selected archived analysis objects and save the final recovered workspace checkpoint.

**Principal input**

The main checkpoint directory is:

`results/checkpoints/Paper1_STEP19_PREP_20260914/`

Required archived files include:

- `R_workspace/Paper1_PRE_STEP19_WORKSPACE_20260914.RData`
- `audit/SF4C/SF4C_statistics_LOCKED_20260914.rds`
- `audit/Supplementary_Table_X/Supplementary_Table_X_TumorOnly_WGCNA_CoreGenes_FINAL.rds`

The script also uses previously generated locked WGCNA, TME, MM, GS, and core-gene integration objects.

**Main analytical and reporting components**

| Component | Main task |
|---|---|
| Supplementary Figure 4A | Prepare and visualize module–trait associations across 50 WGCNA modules and three TME traits. |
| Supplementary Figure 4B | Summarize core-gene module positioning and prepare the final module membership positioning plot. |
| Supplementary Figure 4C | Evaluate within-module MM–GS associations for ME2/blue and ME9/magenta across the three TME traits, including joint BH correction of six tests. |
| Supplementary Figure 4D | Prepare gene significance data for AOC3, MS4A7, and IL7R across ImmuneScore, StromalScore, and ESTIMATEScore and generate the final figure. |
| Supplementary Table X | Summarize module characteristics and TME-related gene significance results for the three prioritized study genes. |
| Final workspace recovery | Recover and validate selected archived objects, restore locked results, and save the final recovered workspace checkpoint. |

**Principal outputs**

Supplementary figures:

- Supplementary Figure 4A.
- Supplementary Figure 4B.
- Supplementary Figure 4C.
- Supplementary Figure 4D.

Supplementary table:

- Supplementary Table X.

Supplementary Table X is exported in CSV, XLSX, RDS, and TXT formats.

The main output directory is:

`results/checkpoints/Paper1_STEP19_PREP_20260914/`

The final recovered workspace checkpoint is generated after the required archived objects and result files have been validated.

**Relationship to the completed workflow**

Script 10 provides the final supplementary visualization, tabulation, and workspace-recovery outputs based on the previously locked tumor-only WGCNA analysis.

**Important scope clarification**

Script 10 does not reconstruct the tumor cohort, rerun the complete WGCNA network construction, or independently repeat the upstream differential expression analysis. Its role is supplementary result preparation, validation, and recovery of previously generated analysis objects.

## 4. Overall Data Flow and Script Dependencies

The main analytical sequence is summarized below.

| Order | Script | Main output or contribution | Downstream use |
|---|---|---|---|
| 1 | `GitHub2_Paper1_tumor_only_WGCNA_1.R` | Verified 513-sample MAD-positive expression matrix, aligned clinical data, and QC records | Input for network construction |
| 2 | `GitHub2_Paper1_tumor_only_WGCNA_2.R` | WGCNA network, module assignments, eigengenes, module–clinical associations, and survival-analysis checkpoints | Input for additional survival and TME analyses |
| 3 | `GitHub2_Paper1_tumor_only_WGCNA_3.R` | ME_brown time-dependent Cox analysis and checkpoint | Continued clinical survival assessment |
| 4 | `GitHub2_Paper1_tumor_only_WGCNA_4.R` | Gender-specific PH assessment and time-dependent Cox results | Input for final integrated survival model |
| 5 | `GitHub2_Paper1_tumor_only_WGCNA_5.R` | Stage-specific PH assessment and final integrated survival model | Upstream survival objects for subsequent integration |
| 6 | `GitHub2_Paper1_tumor_only_WGCNA_6.R` | ESTIMATE scores and module–TME association results | Input for gene-level MM–GS integration |
| 7 | `GitHub2_Paper1_tumor_only_WGCNA_7.R` | Integrated gene-level MM/GS results and core-gene positioning information | Input for candidate-gene integration and auditing |
| 8 | `GitHub2_Paper1_tumor_only_WGCNA_8.R` | Candidate integration, historical provenance audit, and pre-STEP19 archive | Input for further candidate integrity checks |
| 9 | `GitHub2_Paper1_tumor_only_WGCNA_9.R` | Candidate/module mapping audits and pre-STEP19 checkpoint | Input for supplementary result preparation and recovery |
| 10 | `GitHub2_Paper1_tumor_only_WGCNA_10.R` | Supplementary Figures 4A–4D, Supplementary Table X, and final recovered workspace | Final supplementary reporting outputs |

The scripts are sequential at the analytical level, but not every script independently regenerates all objects used by the next script. Several later scripts require previously generated RData/RDS checkpoints and locked in-memory analysis objects.

Before executing a downstream script, confirm that its specified upstream checkpoint and required objects are available.

## 5. General Reproducibility and Interpretation Notes

### 5.1 Project Structure and Path Management

The scripts should be executed from the root directory of the GitHub project.

Project-relative paths are constructed using:

- `here::here()` to identify the project root.
- `file.path()` to construct input and output paths.

The scripts should not depend on manually setting a computer-specific working directory using `setwd()`.

The main project directories include:

- `data/`: input datasets, reference annotations, and required checkpoints.
- `Scripts/`: analysis scripts.
- `results/`: statistical outputs, tables, QC reports, and checkpoints.
- `figures/`: generated figures.

### 5.2 Input and Checkpoint Dependencies

The first script prepares the primary-tumor expression matrix and clinical data. Subsequent scripts rely increasingly on locked upstream analysis objects and archived checkpoints.

These checkpoint files are part of the documented workflow dependencies. They should be retained and checked before executing the corresponding downstream scripts.

### 5.3 Statistical Interpretation

The WGCNA modules represent gene co-expression patterns. Module–clinical, module–TME, MM–GS, and survival-model results describe statistical associations.

These analyses do not independently establish causal relationships or confirm direct biological signaling mechanisms.

Gene-level GS associations and module-level TME associations represent distinct analytical quantities and should be interpreted according to their respective statistical definitions and multiple-testing correction scopes.

### 5.4 Reproducibility Records

Software versions, input provenance, sample manifests, QC reports, and model checkpoints should be retained with the corresponding workflow outputs.

Random seeds are workflow-specific. For example, Script 2 specifies seed 20261106, whereas several later scripts do not explicitly set a new seed. A seed specified in one script should not be assumed to apply globally to all ten scripts.

The exact R and package versions used for execution should be obtained from the corresponding runtime session information rather than inferred from the script-preparation environment.

### 5.5 Separation of Analysis and Reporting

The later scripts include substantial integration, audit, visualization, and workspace-recovery procedures.

These procedures should not be interpreted as independent repetitions of the original network construction or statistical analyses. Where previously locked GS, P-value, FDR, module membership, or survival results are reused, the corresponding scripts preserve and audit those results rather than silently recalculating them.

## 6. Related Documentation

For additional information, consult the following project documentation:

- `Environment/README.md`: consolidated R environment, package dependencies, and software-version records.
- `results/README.md`: inventory and descriptions of generated statistical outputs, tables, QC reports, and checkpoints.
- `figures/README.md`: inventory and descriptions of generated figures, where available.

## 7. Review Scope

This README consolidates the information from the ten script-specific README files for the tumor-only WGCNA analysis.

It describes the intended analytical roles, principal inputs, major procedures, outputs, and dependencies of the ten scripts. It is a workflow-level guide and does not replace the individual R scripts, their runtime validation reports, or the detailed output inventory.

The final script execution order and checkpoint dependencies should be confirmed against the actual repository files and the corresponding script configurations.
