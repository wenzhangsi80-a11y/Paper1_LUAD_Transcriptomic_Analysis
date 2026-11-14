# 03_Survival: Survival Analysis

## 1. Study Overview

**Study title:**

Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD

The `03_Survival/` directory contains the survival analysis workflow for the lung adenocarcinoma (LUAD) study.

This module evaluates the associations between the expression of three previously prioritized study genes — **AOC3, MS4A7, and IL7R** — and overall survival (OS) in LUAD. It combines targeted gene-level survival analyses, clinical covariate-adjusted Cox regression, genome-wide univariable Cox analysis, and exploratory Kaplan–Meier (KM) analyses.

The workflow uses previously processed expression data, clinical information, and upstream WGCNA and survival-analysis checkpoints. The primary survival cohort is retained consistently across the relevant analysis stages.

The three study genes are treated as prioritized genes for downstream characterization. The analyses describe observational transcriptomic associations and do not establish causal relationships, clinical diagnostic validity, or clinical utility.

## 2. Analysis Objectives

The survival analysis module has five main objectives:

1. Recover, reconstruct, and verify the previously established primary LUAD survival cohort.
2. Evaluate the associations between AOC3, MS4A7, IL7R expression and overall survival using continuous-expression Cox regression.
3. Assess proportional-hazards assumptions and investigate the time-dependent association of MS4A7 expression with overall survival.
4. Evaluate the effects of clinical covariates through multivariable Cox regression and summarize the three-gene survival results.
5. Conduct an independent genome-wide univariable Cox analysis and prepare complementary survival visualizations, including KM curves.

## 3. Directory Structure

The module is organized as follows:

```text
03_Survival/
├── Environment/
│   └── README.md
├── scripts/
│   ├── README.md
│   ├── GitHub3_Paper1_survival_analysis_1.R
│   ├── GitHub3_Paper1_survival_analysis_2.R
│   ├── GitHub3_Paper1_survival_analysis_3.R
│   └── GitHub3_Paper1_survival_analysis_4.R
└── results/
    ├── audit_reports/
    ├── figures/
    ├── statistical_outputs/
    ├── tables/
    └── README.md
```

The four analysis scripts are organized by analytical stage. Their detailed execution sequences, input dependencies, and output files are documented in `scripts/README.md` and `results/README.md`.

## 4. Data Sources and Upstream Dependencies

The survival analysis uses previously processed LUAD expression data, clinical information, and analysis checkpoints.

### 4.1 Upstream WGCNA inputs

The tumor-only WGCNA workflow provides the expression and candidate-gene information required for downstream survival analyses.

The principal upstream directory is:

```text
02_WGCNA/results/statistical_outputs/
```

Important upstream data include:

- WGCNA candidate-gene expression data.
- Clinical information aligned with the tumor expression cohort.
- `Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData`

The genome-wide survival analysis uses the WGCNA-derived MAD-positive gene universe.

### 4.2 Survival-analysis checkpoints

Previously generated survival-analysis checkpoints are reused at different stages of the workflow.

Important checkpoint files include:

- `Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData`
- `Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6D1_Survival_Reporting_checkpoint.RData`

These checkpoints preserve previously generated analysis objects and results. Their roles and authoritative data sources are defined by the individual scripts.

The Step 6E workflow uses `survival_gene_primary_z` from the Step 6B5 continuous-Cox PH checkpoint as its authoritative Phase 1 cohort source. Other candidate checkpoints are inspected for provenance and consistency rather than treated as interchangeable sources.

## 5. Locked Survival Cohort

The primary survival analysis uses the previously established cohort of **485 LUAD cases**.

| Cohort characteristic | Number |
|---|---:|
| Total cases | 485 |
| Death events | 178 |
| Censored observations | 307 |

This cohort is used in the principal continuous-expression Cox analyses, the primary clinical-adjustment model, the genome-wide univariable Cox analysis, and the KM workflow.

Some analyses apply model-specific complete-case filtering. In particular, the TNM-adjusted sensitivity model uses a smaller cohort.

| Analysis cohort | Patients | Events | Censored |
|---|---:|---:|---:|
| Primary survival cohort | 485 | 178 | 307 |
| Model A: age, gender, and clinical stage | 485 | 178 | 307 |
| Model B: age, gender, and TNM variables | 336 | 137 | 199 |
| Genome-wide univariable Cox analysis | 485 | 178 | 307 |

The Model B cohort is a model-specific complete-case subset. It does not redefine the original 485-case survival cohort.

The Step 6E KM workflow also identifies four observations with `OS.time = 0`. These observations are retained, and their original survival times are not modified.

## 6. Main Analysis Workflow

The survival analysis is implemented through four main R scripts.

### 6.1 Script 1: Initial Survival Analysis and Continuous Cox Modeling

**Script:**

`GitHub3_Paper1_survival_analysis_1.R`

This script establishes and evaluates the primary survival-analysis framework.

The main procedures include:

1. Recovery and audit of upstream WGCNA and survival-analysis checkpoints.
2. Reconstruction and verification of the locked primary survival cohort.
3. Clinical and expression sample matching.
4. Extraction of AOC3, MS4A7, and IL7R expression.
5. Z-score standardization of the study-gene expression values.
6. Continuous univariable Cox proportional-hazards regression.
7. Schoenfeld residual-based proportional-hazards testing.
8. Evaluation of the time-dependent association between MS4A7 expression and overall survival.
9. Comparison of prespecified time functions using AIC.
10. Likelihood-ratio comparison of constant-coefficient and time-varying models.
11. Estimation of time-specific hazard ratios.

This stage establishes the principal survival cohort and the initial gene-level survival results used by subsequent analyses.

### 6.2 Script 2: Clinical Covariate Reconstruction and Multivariable Cox Modeling

**Script:**

`GitHub3_Paper1_survival_analysis_2.R`

This script evaluates the association between MS4A7 expression and overall survival after adjustment for clinical covariates.

The workflow includes:

**STEP 6C-1: Clinical Covariate Mapping and QC**

- Recovery of the locked survival cohort and corresponding clinical information.
- Matching of clinical records using `Sample`.
- Patient-level source auditing using `patient_id`.
- Inspection of age, gender, clinical stage, and TNM variables.
- Evaluation of missingness and stage consistency.
- Construction of `clinical_covariates_locked`.

**STEP 6C-1.5: Clinical Stage Consistency Audit**

- Cross-tabulation of detailed clinical stage.
- Inspection of M, T, and N categories.
- Review of unknown categories, including TX, NX, and MX.
- Verification of the original 485-case cohort.

**STEP 6C-2: Clinical Covariate Coding Freeze**

- Construction of `clinical_covariates_coded`.
- Definition of age, gender, clinical stage, and TNM categories.
- Preservation of the original survival cohort.

**STEP 6C-3 and STEP 6C-4: Multivariable Cox Modeling and PH Assessment**

Two clinical-adjustment designs are evaluated:

- **Model A:** `MS4A7_z + Age + Gender + Stage`
- **Model B:** `MS4A7_z + Age + Gender + T_stage + N_stage + M_stage`

Model A uses 485 cases and 178 events. Model B uses 336 complete cases and 137 events.

The workflow uses `survival::coxph()` and `survival::cox.zph()` to estimate adjusted hazard ratios and evaluate proportional-hazards assumptions.

**STEP 6C-5: Adjusted Time-Varying Cox Model for MS4A7**

This stage evaluates prespecified time functions for the MS4A7 expression effect:

1. Linear: `MS4A7_z × t`
2. Square-root: `MS4A7_z × sqrt(t)`
3. Logarithmic: `MS4A7_z × log1p(t)`

Candidate models are compared using AIC. The selected model is compared with the constant-coefficient Model A using likelihood-based statistics and AIC. Time-specific adjusted hazard ratios and 95% confidence intervals are estimated at 1, 2, 3, 4, and 5 years.

A checkpoint is used to preserve and verify the selected model, numerical estimates, and relevant PH results.

### 6.3 Script 3: Three-Gene Survival Reporting and Genome-Wide Cox Analysis

**Script:**

`GitHub3_Paper1_survival_analysis_3.R`

This script contains two complementary analytical branches.

**A. STEP 6D: Three-Gene Cox Survival Reporting**

The reporting workflow summarizes the survival associations of AOC3, MS4A7, and IL7R.

The main procedures include:

- Verification of previously generated Cox model objects.
- Extraction and summarization of continuous-expression univariable Cox results.
- Gene-specific adjusted Cox analyses using Model A and Model B.
- Benjamini–Hochberg false discovery rate (BH-FDR) adjustment within the corresponding reporting models.
- Summarization of previously generated MS4A7 time-varying Cox results.
- Preparation of time-specific hazard-ratio summaries.
- Preparation of the Figure 4D forest-plot dataset.
- Quality-control checks of cohort sizes, model structures, and reporting datasets.

The Figure 4 workflow visualizes the three-gene continuous-expression Cox results, including hazard ratios, 95% confidence intervals, P values, and adjusted FDR values for the univariable and adjusted models.

The previously generated MS4A7 time-varying Cox results are summarized from locked objects rather than refitted in the reporting section.

**B. STEP 6X: Genome-Wide Univariable Cox Analysis**

This is a separate genome-wide analysis branch using the locked 485-case survival cohort and the 31,609-gene MAD-positive expression universe derived from the tumor-only WGCNA workflow.

The workflow includes:

1. Recovery and validation of the locked survival cohort.
2. Recovery and alignment of the genome-wide expression matrix.
3. Construction of the survival phenotype.
4. Definition and validation of the WGCNA-based gene universe.
5. Gene-wise univariable Cox proportional-hazards regression.
6. Genome-wide BH-FDR adjustment.
7. Numerical, warning, and model-status auditing.
8. Recovery and validation of the genome-wide Cox results.
9. Preparation and freezing of the Figure X dataset.
10. Plotting-data distribution auditing and preparation.
11. Recovery of plotting checkpoints and final Figure X visualization.

Gene expression is standardized within the 485-patient analysis cohort. Genome-wide hazard ratios therefore represent the association per one-standard-deviation increase in gene expression.

The genome-wide branch is distinct from the targeted three-gene reporting workflow. It does not intentionally modify the previously locked Step 6B, Step 6C, or Step 6D-1 reporting results.

### 6.4 Script 4: Survival Cohort Audit and Kaplan–Meier Curve Harmonization

**Script:**

`GitHub3_Paper1_survival_analysis_4.R`

This script implements STEP 6E and provides complementary exploratory survival analyses based on median-defined expression groups.

The workflow contains three main components:

**Phase 1: Survival Cohort and Input Audit**

- Recovery and verification of the locked 485-case survival cohort.
- Inspection of checkpoint provenance.
- Validation of survival variables and sample identifiers.
- Audit of AOC3, MS4A7, and IL7R expression values.
- Verification that the input expression values are already standardized.
- Identification and retention of four zero-time observations.

No automatic re-standardization or survival-time modification is performed.

**Phase 2A: Median-Split Audit**

For each study gene, the median expression value is calculated within the locked cohort.

The grouping rule is:

- Expression less than or equal to the median: Low.
- Expression greater than the median: High.

Group sizes, median thresholds, and zero-time observations by group are audited before survival estimation.

**Phase 2B and Phase 2B/2C: KM Analysis and Plot Harmonization**

The workflow includes:

- Kaplan–Meier survival estimation.
- Confidence interval display.
- Number-at-risk table generation.
- Two-sided log-rank tests.
- Grouped Cox regression comparing High versus Low expression.
- Estimation of grouped hazard ratios and 95% confidence intervals.
- Harmonization of plot typography and statistical annotations.
- Export of TIFF and PDF figures.
- Audit of exported plot files and analysis results.

The grouped Cox hazard ratios compare median-defined expression groups and are distinct from the continuous-expression Cox estimates reported in the preceding analyses.

## 7. Statistical and Reporting Principles

The following principles are applied across the survival analysis module:

- The previously locked primary survival cohort is retained throughout the relevant analysis stages.
- Model-specific complete-case filtering is applied only where required by the corresponding model design.
- Continuous-expression Cox analyses and median-group survival analyses are reported as distinct analytical approaches.
- Model A and Model B results are reported separately.
- The TNM-adjusted Model B is interpreted as a sensitivity analysis based on its smaller complete-case cohort.
- Proportional-hazards diagnostics are performed for applicable Cox models.
- The MS4A7 time-varying Cox analysis evaluates prespecified time functions and preserves the selected model and its estimates.
- BH-FDR adjustment is applied within the corresponding three-gene reporting models and across the genome-wide tested gene universe, respectively.
- Genome-wide effect sizes are expressed as hazard ratios per one-standard-deviation increase in expression.
- Figure preparation and plotting do not constitute additional statistical model fitting unless explicitly specified by the corresponding analysis stage.
- The four zero-time observations are retained in the STEP 6E workflow.
- No causal interpretation is made from the observed survival associations.

## 8. Output Organization

All survival-analysis outputs are stored under:

```text
03_Survival/results/
```

The main output categories are:

| Directory | Contents |
|---|---|
| `audit_reports/` | Cohort, clinical-variable, expression, checkpoint, and model-quality audits |
| `figures/` | Survival figures and figure exports |
| `statistical_outputs/` | Statistical results, model objects, and analysis checkpoints |
| `tables/` | Summarized survival results and reporting tables |

The scripts also generate workflow-specific subdirectories for their analysis outputs.

Important output groups include:

- **Initial survival analysis:** cohort reconstruction, continuous Cox results, PH diagnostics, and MS4A7 time-varying Cox checkpoints.
- **Clinical-adjustment analysis:** clinical covariate audits, Model A and Model B results, PH diagnostics, and the adjusted MS4A7 time-varying Cox checkpoint.
- **Three-gene reporting:** STEP 6D reporting tables, Figure 4D plotting data, and reporting checkpoints.
- **Genome-wide survival analysis:** STEP 6X Cox results, audit records, data-freeze and plotting checkpoints, and Figure X exports.
- **KM analysis:** STEP 6E cohort and expression audits, median-split summaries, KM fit and plot objects, grouped Cox results, and harmonized KM figures.

For the detailed output-file inventory and exact output locations, see `results/README.md`.

## 9. Software Environment and Reproducibility

The survival analysis is implemented in R and uses packages for survival modeling, data processing, statistical reporting, and visualization.

The project uses project-root-relative paths constructed with `here::here()` where specified by the individual scripts.

The following reproducibility principles apply:

- Upstream WGCNA and survival-analysis checkpoints are loaded from their designated project directories.
- Previously locked cohorts and model objects are reused where specified.
- Clinical covariate coding is preserved after the STEP 6C-2 coding freeze.
- Model-specific complete-case cohorts are not used to redefine the original primary survival cohort.
- Previously fitted models are not unnecessarily refitted during reporting or checkpoint-recovery stages.
- The continuous-expression Cox estimates are kept distinct from the median-group KM and grouped Cox results.
- Checkpoint existence alone does not establish successful quality control; the corresponding audit results and PASS/FAIL status must also be inspected.
- R and package versions should be documented in `Environment/README.md` based on the finalized local execution environment.
- Actual output files, numerical results, and QC outcomes should be confirmed against the corresponding execution records.

## 10. Interpretation and Scope

The survival analysis module provides complementary assessments of the associations between the prioritized study genes and overall survival in LUAD.

The targeted Cox analyses evaluate continuous-expression associations and clinical-adjusted effects. The genome-wide analysis provides a separate assessment across the WGCNA-derived gene universe. The KM workflow offers an exploratory comparison of survival between median-defined expression groups.

These analyses are observational and are intended to characterize survival associations within the analyzed cohort. They do not establish causal effects, clinical diagnostic validity, or clinical utility.

## 11. Related Documentation

The following documentation provides additional details for the `03_Survival/` module:

- `Environment/README.md`: R version, package dependencies, and software environment.
- `scripts/README.md`: individual script descriptions, execution stages, and input dependencies.
- `results/README.md`: detailed output-directory structure, checkpoint descriptions, and result-file inventory.