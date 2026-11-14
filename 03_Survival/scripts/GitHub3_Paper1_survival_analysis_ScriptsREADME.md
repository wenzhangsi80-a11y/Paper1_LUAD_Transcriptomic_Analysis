# Survival Analysis Scripts

## 1. Overview

This directory contains the R scripts used for the survival analysis workflow in the Paper 1 LUAD study.

**Study:**\
Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD

The survival analysis workflow is organized into four main scripts:

| Main script | Analysis module | Primary purpose |
|----|----|----|
| `GitHub3_Paper1_survival_analysis_1.R` | STEP 6B | Survival cohort recovery, continuous univariable Cox analysis, PH diagnostics, and time-varying Cox analysis |
| `GitHub3_Paper1_survival_analysis_2.R` | STEP 6C | Clinical covariate preparation, multivariable Cox modeling, PH diagnostics, and adjusted time-varying Cox analysis |
| `GitHub3_Paper1_survival_analysis_3.R` | STEP 6D / STEP 6X | Three-gene survival reporting, Figure 4, genome-wide Cox analysis, and Figure X |
| `GitHub3_Paper1_survival_analysis_4.R` | STEP 6E | Survival cohort audit and Kaplan–Meier curve harmonization |

The workflow uses the previously established survival cohort and WGCNA-derived expression data. Subsequent scripts recover and analyze locked upstream objects rather than independently redefining the original cohort.

## 2. Analysis Workflow

### 2.1 STEP 6B — Survival Cohort Recovery and Univariable Cox Analysis

**Main script:**\
`GitHub3_Paper1_survival_analysis_1.R`

This script recovers the upstream WGCNA and survival-analysis objects, reconstructs the locked survival cohort, and performs continuous-expression univariable Cox analysis and time-varying Cox modeling for MS4A7.

Main analysis stages:

| Analysis step | Description |
|----|----|
| STEP 6B recovery audit | Identify and audit the upstream WGCNA and survival project directories |
| STEP 6B recovery | Recover upstream WGCNA checkpoints |
| Cohort reconstruction | Reconstruct and verify the locked 485-sample survival cohort |
| STEP 6B-4 recovery | Standardize core-gene expression using Z-scores |
| STEP 6B-5 | Perform continuous univariable Cox regression and Schoenfeld proportional-hazards testing |
| STEP 6B-5 checkpoint | Save the continuous Cox and PH analysis checkpoint |
| STEP 6B-6 | Evaluate time-varying Cox models for MS4A7 |
| STEP 6B-6 model comparison | Compare linear, square-root, and log1p time functions using AIC |
| STEP 6B-6 model evaluation | Compare constant and selected time-varying models using a likelihood-ratio test |
| STEP 6B-6 time-specific estimates | Estimate MS4A7 hazard ratios and 95% confidence intervals at 1–5 years |
| STEP 6B-7 | Save the MS4A7 time-varying Cox analysis checkpoint |

**Main input checkpoints:**

- `Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds`
- `Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData`

Expected input directory:

`02_WGCNA/results/statistical_outputs/`

**Main output checkpoints:**

- `Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData`
- `Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData`

Expected output directory:

`03_Survival/results/statistical_outputs/`

### 2.2 STEP 6C — Clinical Covariate Preparation and Multivariable Cox Analysis

**Main script:**\
`GitHub3_Paper1_survival_analysis_2.R`

This script prepares and audits clinical covariates, constructs the multivariable Cox design cohorts, evaluates proportional-hazards assumptions, and performs adjusted time-varying Cox analysis for MS4A7.

Main analysis stages:

| Analysis step | Description |
|----|----|
| STEP 6C-1 | Clinical covariate mapping and QC, including recovery of the locked survival cohort and patient-level matching |
| STEP 6C-1.5 | Clinical stage consistency audit, including detailed stage, M status, T category, and N category |
| STEP 6C-2 | Clinical covariate coding freeze, including age, gender, clinical stage, and TNM categories |
| STEP 6C-3 | Multivariable Cox design-cohort QC and construction of Model A and Model B complete-case cohorts |
| STEP 6C-4 | Multivariable Cox modeling, proportional-hazards diagnostics, and model QC |
| STEP 6C-5 | Adjusted time-varying Cox modeling for MS4A7, candidate time-function comparison, model selection, and time-specific HR estimation |
| STEP 6C-5 checkpoint | Preservation and QC of the adjusted time-varying Cox analysis objects |

**Model specifications:**

- **Model A:** `MS4A7_z + Age + Gender + Stage`
- **Model B:** `MS4A7_z + Age + Gender + T_stage + N_stage + M_stage`

The locked Model A cohort contains 485 patients (178 events). Model B uses 336 complete cases (137 events) because of the additional TNM covariate requirements.

The adjusted time-varying analysis evaluates three prespecified time functions: linear, square-root, and log1p. Candidate models are compared using AIC and log-likelihood, and the selected model is compared with the constant-coefficient Model A using a likelihood-ratio test.

**Main analysis objects include:**

- `clinical_covariates_locked`
- `clinical_covariates_coded`
- `modelA_cc`
- `modelB_cc`
- `cox_modelA`
- `cox_modelB`
- `ms4a7_tv_adjusted_models`
- `ms4a7_tv_adjusted_comparison`
- `best_time_function_adjusted`
- `best_tv_adjusted_model`
- `ms4a7_adjusted_LRT`
- `ms4a7_adjusted_hr_time`
- `step6c5_checkpoint_metadata`

**Main persistent output:**

`Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData`

Expected output directory:

`03_Survival/results/statistical_outputs/`

The checkpoint-preservation step verifies and saves previously fitted analysis objects. It does not refit the adjusted time-varying models.

### 2.3 STEP 6D / STEP 6X — Targeted Survival Reporting and Genome-Wide Cox Analysis

**Main script:**\
`GitHub3_Paper1_survival_analysis_3.R`

This script contains two related but distinct analysis branches:

1.  Targeted survival reporting for AOC3, MS4A7, and IL7R.
2.  Genome-wide univariable Cox analysis using the WGCNA-derived candidate gene universe.

#### 2.3.1 STEP 6D — Three-Gene Survival Reporting

The targeted reporting workflow first verifies the previously locked survival cohort and Cox model objects.

Main analysis stages:

| Analysis step | Description |
|----|----|
| STEP 6D pre-check | Verify the locked cohort, standardized survival data, and existing Cox model structures without refitting the locked models |
| STEP 6D-1 | Extract and report continuous-expression Cox results for AOC3, MS4A7, and IL7R |
| Model A reporting | Fit gene-specific adjusted Cox models using the Model A clinical covariates |
| Model B reporting | Fit gene-specific TNM sensitivity models using the Model B complete-case cohort |
| Reporting table preparation | Combine univariable, Model A, and Model B results and prepare HR, confidence interval, P-value, FDR, and PH-test summaries |
| Figure 4 data preparation | Construct and QC the three-gene forest-plot dataset |

BH-adjusted FDR values are calculated within the corresponding three-gene reporting sets.

Previously locked MS4A7 time-varying results are summarized from existing analysis objects. The reporting workflow does not recalculate the standardized expression values or refit the previously locked time-varying models.

**Main reporting objects include:**

- `univariable_table`
- `modelA_reporting`
- `modelA_fits`
- `modelA_table`
- `modelB_reporting`
- `modelB_fits`
- `modelB_table`
- `continuous_cox_table`
- `continuous_cox_table_display`
- `ms4a7_tv_summary`
- `ms4a7_tv_adjusted_summary`
- `ms4a7_time_specific_hr`
- `figure4D_data`
- `step6d1_checkpoint_metadata`

#### 2.3.2 Figure 4 — Three-Gene Cox Forest Plot

The Figure 4 section generates the compact forest plot for the three prioritized study genes.

The figure displays:

- Univariable Cox estimates.
- Model A adjusted Cox estimates.
- Model B TNM sensitivity-model estimates.
- Hazard ratios and 95% confidence intervals.
- The HR = 1 reference line.

The plotting dataset is checked for gene labels, model labels, and display ranges before export.

**Figure 4 export dimensions:** 8.8 × 4.6 inches.

Export formats include PDF, TIFF, PNG, and JPEG. Output paths are defined by the corresponding figure-export path objects in the script.

#### 2.3.3 STEP 6X-1 to STEP 6X-8B — Genome-Wide Analysis Setup

This section initializes the genome-wide survival analysis branch.

Main operations:

- Recover the locked 485-patient survival cohort.
- Read the tumor-only WGCNA candidate expression matrix.
- Align expression samples with the locked survival sample order.
- Construct and validate the survival phenotype.
- Prepare the genome-wide expression matrix.
- Define the 31,609-gene MAD-positive WGCNA survival gene universe.
- Check gene overlap and the presence of AOC3, MS4A7, and IL7R.

**Main expression input:**

`Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt`

The genome-wide workflow uses the WGCNA-derived candidate expression matrix and the previously locked survival cohort.

#### 2.3.4 STEP 6X-8C / STEP 6X-8D — Genome-Wide Univariable Cox Analysis and Audit

STEP 6X-8C fits gene-wise univariable Cox proportional-hazards models across the 31,609-gene expression universe.

Expression is standardized within the 485-patient cohort. The reported hazard ratio represents the association per one-standard-deviation increase in gene expression.

For each gene, the analysis records:

- Sample size and event count.
- Expression mean and standard deviation.
- Regression coefficient.
- Hazard ratio and 95% confidence interval.
- Nominal P value.
- Genome-wide BH-adjusted FDR.
- Model status, warnings, and errors.

**Main result object:**

`genome_cox_results_6X`

STEP 6X-8D audits the genome-wide results by checking model status, numerical validity, FDR distribution, association direction, top-ranked genes, and the results for AOC3, MS4A7, and IL7R.

The genome-wide analysis does not perform additional genome-wide PH testing in this workflow.

#### 2.3.5 STEP 6X-8E — Checkpoint Recovery and Gene-Universe Validation

These sections recover and validate the genome-wide Cox results without refitting the models.

Main operations:

- Recover the STEP 6X-8C result and gene-universe objects.
- Check result dimensions and required columns.
- Verify the survival cohort size and event count.
- Compare the Cox results with the defined gene universe.
- Validate Ensembl identifiers, gene counts, and uniqueness.
- Verify the presence of the three prioritized study genes.
- Inspect the availability of the STEP 6X-8D audit output.

STEP 6X-8E-final prepares and freezes the Figure X dataset from `genome_cox_results_6X`. It preserves the original statistical source values and records plotting metadata and core-gene annotations.

**Important checkpoint note:** The STEP 6X-8E-final workflow includes a subsequent core-gene annotation correction and QC sequence. The final saved checkpoint should be inspected to confirm that it contains the corrected object, because the described correction sequence is not followed by another save operation.

#### 2.3.6 STEP 6X-9A to STEP 6X-9C — Figure X Data Audit and Visualization

| Analysis step | Description |
|----|----|
| STEP 6X-9A | Audit the frozen Figure X dataset, including dimensions, missing values, statistical distributions, significant-gene directions, reference lines, and core-gene annotations |
| STEP 6X-9B | Prepare the Figure X plotting dataset, preserve source values and gene order, define plotting classes and coordinates, and save the plotting-preparation checkpoint |
| STEP 6X-9C-R1/R2 | Inspect and recover the STEP 6X-9B plotting checkpoint and perform object-level QC |
| STEP 6X-9C-R4 | Refine core-gene label placement, add leader lines and the FDR threshold annotation, and export the final Figure X |

The STEP 6X-9A audit object, `audit_6X9A`, is an in-memory object and is not written as a separate audit file in the supplied workflow.

The STEP 6X-9B checkpoint is saved by the described code even if `step6X9B_pass` is FALSE. Downstream use therefore requires inspection of this QC flag.

The STEP 6X-9C-R4 section performs visual refinement and figure export. It does not refit genome-wide Cox models or recalculate gene-level statistical results.

### 2.4 STEP 6E — Survival Cohort Audit and Kaplan–Meier Curve Harmonization

**Main script:**\
`GitHub3_Paper1_survival_analysis_4.R`

This script recovers and audits the locked survival cohort and performs the prespecified median-split Kaplan–Meier analysis for AOC3, MS4A7, and IL7R.

Main analysis stages:

| Analysis stage | Description |
|----|----|
| Setup | Detect the project root and define input and output directories |
| STEP 6E Phase 1 | Recover the survival cohort and inspect checkpoint provenance |
| Phase 1 cohort audit | Validate sample identifiers, survival fields, event coding, and gene expression values |
| Zero-time audit | Inspect the four observations with `OS.time = 0` and retain them |
| Gene expression audit | Verify that the three gene expression values are already standardized |
| Phase 2A | Audit median-split thresholds and group sizes before Kaplan–Meier estimation |
| Phase 2B | Fit Kaplan–Meier curves, calculate log-rank tests, and export plots and tables |
| Phase 2B/2C | Re-run the locked KM workflow with grouped Cox HR annotations and revised plot formatting |
| Output audit | Inspect exported plots, result tables, and analysis objects |

**Authoritative Phase 1 cohort source:**

`survival_gene_primary_z` in the STEP 6B5 checkpoint.

**Main input checkpoints:**

- `Paper1_Step6D1_Survival_Reporting_checkpoint.RData`
- `Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData`
- STEP 6B6 checkpoint files, where available
- `Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData`

The script also refers to the previous KM workflow file:

`Paper 1 – STEP_6E_KM_curve_harmonization.R`

**Main analysis objects include:**

- `dat`
- `zero_dat`
- `gene_summary`
- `split_summary`
- `zero_time_split`
- `group_membership`
- `km_fits`
- `km_plots`
- `logrank_results`
- `grouped_cox_results`
- `file_audit`

**Median-split grouping rule:**

For each gene, the median is calculated in the locked 485-case survival cohort:

- Low: expression ≤ cohort median.
- High: expression \> cohort median.

The four original zero-time observations are retained.

**Important distinction:** The grouped Cox hazard ratios from STEP 6E compare the median-defined High and Low expression groups. They are not the continuous-expression Cox hazard ratios reported in the preceding survival-analysis workflow.

The uploaded script contains repeated audit and rerun sections. These sections are retained under the path-only editing scope and should not be interpreted as separate independent cohorts.

## 3. Key Cohort and Analysis Definitions

The four scripts share the previously locked survival cohort and related upstream analysis objects.

| Analysis component | Definition |
|----|----|
| Primary survival cohort | 485 patients |
| Observed events | 178 |
| Censored observations | 307 |
| Model A | 485 complete cases; 178 events |
| Model B | 336 complete cases; 137 events |
| Prioritized study genes | AOC3, MS4A7, and IL7R |
| Genome-wide gene universe | 31,609 MAD-positive WGCNA candidate genes |
| Genome-wide Cox effect scale | HR per one-standard-deviation increase in expression |
| KM grouping | Gene-specific median split within the locked 485-case cohort |
| Zero-time observations | Four observations with `OS.time = 0`, retained in STEP 6E |

The cohort definitions and sample identifiers should remain consistent across the four scripts. Any change to the locked cohort, clinical coding, expression values, or gene universe requires a separately documented analytical decision.

## 4. Main Checkpoints and Output Locations

### 4.1 Upstream WGCNA Inputs

Expected upstream input directory:

`02_WGCNA/results/statistical_outputs/`

Important upstream files include:

- `Workflow-04-2N_WGCNA_OS_multivariable_cohort_QC.rds`
- `Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData`
- `Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt`

### 4.2 Survival Analysis Checkpoints

Expected survival checkpoint directory:

`03_Survival/results/statistical_outputs/`

Important checkpoint files include:

- `Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData`
- `Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6D1_Survival_Reporting_checkpoint.RData`

Additional STEP 6X genome-wide and Figure X checkpoints are generated by the corresponding sections of Script 3. Their exact filenames and locations should be taken from the script and verified against the actual output directory.

### 4.3 Figures, Tables, and Other Outputs

The scripts generate or prepare:

- Continuous and adjusted Cox reporting tables.
- MS4A7 time-varying Cox summaries and time-specific HR estimates.
- Figure 4 three-gene Cox forest plot.
- Genome-wide Cox result tables and audit outputs.
- Figure X genome-wide survival visualization.
- Kaplan–Meier plots for AOC3, MS4A7, and IL7R.
- Log-rank and grouped Cox result tables.
- Checkpoint and output-file audit information.

Persistent output files should be documented in the corresponding results README after their filenames and locations have been confirmed from the scripts or actual execution outputs.

## 5. Execution and Reproducibility Notes

1.  Execute the scripts within the configured project structure using project-relative paths managed through `here::here()`.
2.  Ensure that the required upstream WGCNA inputs and survival checkpoints are available before running downstream analysis sections.
3.  Preserve the locked 485-patient survival cohort, sample identifiers, event coding, and original standardized gene expression values.
4.  Do not silently modify clinical covariate coding, model reference categories, or complete-case cohort definitions.
5.  Do not refit previously locked models during reporting or checkpoint-preservation steps unless a separately documented analysis step explicitly requires model fitting.
6.  Preserve the distinction between continuous-expression Cox analyses and median-grouped Cox analyses.
7.  Retain the four observations with `OS.time = 0` in the STEP 6E workflow.
8.  Do not interpret the presence of a checkpoint file as evidence that its QC passed.
9.  Inspect `step6X9B_pass` before using the STEP 6X-9B plotting checkpoint downstream.
10. Inspect the final STEP 6X-8E checkpoint to confirm that the subsequent core-gene annotation correction is represented in the saved object.
11. Treat STEP 6X-9A audit results as in-memory objects unless a separate output file is confirmed.
12. Record actual execution details, R version, package versions, and any section-specific random seeds in the corresponding Environment documentation.
13. Verify generated checkpoint files, figures, and tables after execution.

## 6. Related Documentation

- `../README.md` — Survival analysis module overview.
- `../Environment/README.md` — R environment, package dependencies, and computational settings.
- `../results/README.md` — Survival analysis results, checkpoints, figures, and tables.
