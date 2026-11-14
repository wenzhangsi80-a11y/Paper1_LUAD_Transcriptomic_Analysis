# Diagnostic Classification

## Study

**Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R across the TME in LUAD**

## Research Objective

This module contains the diagnostic classification workflow for the LUAD study. It evaluates the exploratory ability of three prioritized study genes — **AOC3, MS4A7, and IL7R** — to distinguish tumor from normal lung tissue using transcriptomic data.

The workflow integrates historical data provenance audits, sample and subject-level metadata assessments, expression-scale and model-input traceability analyses, model-definition locking, internal cross-validation, final model fitting, external dataset assessment, and manuscript figure generation.

The classification analyses are exploratory tissue-discrimination assessments. They do not establish clinical diagnostic utility, clinical readiness, or a clinically validated diagnostic biomarker.

## Directory Structure

``` text
09_Classification/
├── data/
│   ├── data_sources/
│   ├── input_data/
│   └── processed_data/
├── Environment/
│   └── README.md
├── scripts/
│   ├── README.md
│   └── GitHub9_Paper1_Diagnostic_Classification.R
├── results/
│   ├── README.md
│   └── WORKFLOW13_Diagnostic_Classification/
└── README.md
```

The detailed input-data inventory, script organization, software environment, and generated output inventory are documented in the corresponding subdirectory README files.

## Data Sources and Cohort Definitions

### 1. Historical TCGA ROC Branch

The historical TCGA ROC analysis branch is based on the archived Xiangtao-derived data and historical RData workspace.

The historical branch is documented with the following sample composition:

| Sample group | Number of samples |
|--------------|------------------:|
| Tumor        |               516 |
| Normal       |                59 |
| Total        |               575 |

The historical workspace is located under:

``` text
data/input_data/Paper1_Figure2_Figure12_TCGA_LUAD_ROC_etal/
```

The historical expression source is documented as TCGA-LUAD STAR RNA-seq TPM, with the recorded transformation:

``` r
log2(value + 1)
```

Historical sample-count and expression-object provenance are audited separately from the harmonized TCGA dataset used in the broader LUAD analyses.

Importantly, the historical 575-sample branch must not be conflated with the harmonized TCGA expression dataset described below. The historical count-provenance audit also distinguishes literal occurrences of the number 516 from evidence establishing a valid sample or patient count.

### 2. Harmonized TCGA Dataset

The separately harmonized TCGA expression dataset contains:

| Sample group        | Number of samples |
|---------------------|------------------:|
| Primary Tumor       |               513 |
| Solid Tissue Normal |                59 |
| Metastatic          |                 2 |
| Total               |               574 |

The tumor-versus-normal subset contains 572 eligible samples after excluding the two metastatic samples.

These counts belong to the harmonized TCGA data branch and are not interchangeable with the historical 575-sample ROC branch.

### 3. External GSE32863 Dataset

The external dataset is GSE32863, based on GPL6884.

The audited sample structure is:

| Characteristic                      | Number |
|-------------------------------------|-------:|
| Total samples                       |    116 |
| Tumor samples                       |     58 |
| Normal samples                      |     58 |
| Unique subjects                     |     59 |
| Complete Normal–Tumor subject pairs |     57 |
| Tumor-only subjects                 |      1 |
| Normal-only subjects                |      1 |

The two singleton subjects are retained separately:

- Tumor-only subject: 3023
- Normal-only subject: 3035

Subject identifiers are obtained from explicit GEO metadata, while sample codes and tumor/normal labels are derived from GEO sample titles. The pairing workflow evaluates their concordance using explicit metadata rather than sample order or sample counts.

The two singleton subjects are not artificially paired.

The archived GEO input files include:

``` text
GSE32863_series_matrix.txt.gz
GSE32863_family.soft.gz
```

## Core Model Specification

The prespecified classification model is a multivariable logistic regression model:

``` r
TumorStatus ~ AOC3 + MS4A7 + IL7R
```

The model specification is:

| Component             | Definition                        |
|-----------------------|-----------------------------------|
| Model type            | Multivariable logistic regression |
| Predictors            | AOC3, MS4A7, IL7R                 |
| Outcome               | TumorStatus                       |
| Outcome levels        | Normal, Tumor                     |
| Reference level       | Normal                            |
| Positive class        | Tumor                             |
| Probability threshold | 0.50                              |

The three genes are treated as prioritized study genes in the classification workflow. Their inclusion in the same model does not imply a causal relationship, a mechanistic signaling axis, or a validated biological interaction.

## Main Workflow Components

All workflow-specific outputs are organized under:

``` text
results/WORKFLOW13_Diagnostic_Classification/
```

### STEP13A — Clinical Claim De-escalation and Interpretation Lock

Establishes conservative interpretation boundaries for the diagnostic classification analysis.

The workflow distinguishes exploratory tumor-versus-normal tissue discrimination from clinical diagnostic validity and restricts unsupported claims regarding clinical utility, diagnostic biomarkers, and clinical implementation.

### STEP13B — Historical TCGA Sample-Source and Count-Provenance Audits

Audits the provenance and consistency of historical TCGA sample counts and archived classification-related objects.

The workflow includes historical count comparisons, sample-set comparisons, metastatic-sample provenance checks, historical RData inspections, and traceability assessments of the Xiangtao-derived ROC Excel branch.

The historical Xiangtao Excel audit evaluates the expected sample-gene structure and its correspondence with historical ROC objects. Structural concordance is not equivalent to independent validation of the original ROC analysis or proof that different TCGA data branches are interchangeable.

### STEP13C — GSE32863 Input and Sample-Identity Lock

Audits the provenance and structure of the external GSE32863 dataset before downstream classification analyses.

The workflow includes:

- Historical expression-object inventory and traceability.
- Expression-matrix provenance assessment.
- Sample annotation and Normal/Tumor group concordance.
- Historical metadata and sample-ID source audits.
- GEO Series Matrix and Family SOFT metadata inspection.
- Title-derived sample-code mapping.
- Candidate sample-code pairing assessment.
- Explicit GEO subject-ID concordance assessment.
- Final subject-level matched-pair structure locking.

The final audited structure contains 57 complete Normal–Tumor subject pairs and two singleton subjects.

A title-derived sample code is not automatically treated as a patient identifier. Patient-level pairing conclusions rely on explicit subject metadata and the corresponding concordance checks.

### STEP13E — Historical Expression-Scale and Model-Input Traceability

Audits the historical TCGA and GSE32863 expression objects and the compatibility of the archived model inputs.

The workflow examines:

- Historical TCGA ROC expression-scale provenance.
- Historical GSE32863 expression-object structure and dimensions.
- Sample and gene identifiers.
- Core-gene availability.
- Sample-level expression-value correspondence.
- Historical prediction-vector and ROC-object traceability.
- Cross-platform model-input compatibility.

These steps document the provenance and structure of historical analytical objects. They do not independently reconstruct expression matrices, repeat probe selection, or recalculate historical ROC/AUC statistics.

### STE P13F — Model Definition Lock

Records and audits the prespecified classification model, outcome coding, reference level, positive class, cross-validation policy, and restrictions on downstream analysis.

The model definition is locked before the subsequent model assessment and external evaluation steps.

The lock also documents restrictions on external model refitting, coefficient estimation, feature selection, threshold optimization, cutpoint searching, and model tuning.

### STEP13G — Internal Classification Assessment

Evaluates the prespecified model using stratified 10-fold outer cross-validation in the historical TCGA ROC branch.

The workflow:

- Generates stratified outer folds.
- Fits the model using each outer training set.
- Obtains predictions for the corresponding held-out outer test set.
- Combines the held-out predictions into an out-of-fold prediction dataset.
- Generates inner-fold assignments within outer training sets for stratification auditing.

**Methodological clarification:** The uploaded implementation uses 10-fold stratified outer cross-validation with inner-fold stratification auditing. The inner-fold assignments are not used for inner-loop model fitting, validation-based model selection, feature selection, or hyperparameter tuning. Therefore, this implementation should not be described as a fully executed conventional nested cross-validation procedure.

The out-of-fold predictions provide an internal resampling-based assessment and do not constitute independent external validation.

### STEP13H — Final TCGA Logistic Regression Model

Fits the prespecified multivariable logistic regression model using the historical TCGA ROC branch.

The model is:

``` r
TumorStatus ~ AOC3 + MS4A7 + IL7R
```

The workflow audits complete-case retention and extracts the model coefficients and associated statistics, including:

- Regression coefficients.
- Standard errors.
- z statistics.
- P values.
- Odds ratios.
- 95% Wald confidence intervals.

The fitted model object and coefficient information are retained for traceability and subsequent external assessment.

The STEP13H model coefficients are transferred unchanged to the external GSE32863 assessment.

### STEP13I — GSE32863 External Model Assessment and Paired Expression Analysis

Applies the final TCGA model from STEP13H to the historical GSE32863 expression dataset.

The external assessment includes:

- Model-definition and model-status verification.
- Locked coefficient-transfer auditing.
- Historical GSE32863 object and sample-count checks.
- Sample-ID and group concordance assessment.
- Core-gene availability and expression-value checks.
- External model-input completeness assessment.
- Predicted-probability generation.
- External ROC/AUC analysis and confidence-interval estimation.
- Auditing of prohibited external-analysis operations.

The TCGA model coefficients are applied without external model refitting or coefficient re-estimation. External feature selection, threshold optimization, cutpoint searching, and model tuning are not performed.

The same workflow also includes a separate paired-expression analysis of AOC3, MS4A7, and IL7R in the 57 complete Normal–Tumor subject pairs.

Paired Wilcoxon signed-rank tests are used for the three gene-level comparisons, with Benjamini–Hochberg adjustment applied to the three P values.

The paired-expression analysis is separate from the transferred-model assessment. It does not modify the locked model, coefficients, predictions, AUC, or classification threshold.

### STEP13J — Static Final Manuscript Result Lock and Figure Preparation

Organizes and documents the finalized classification results used for manuscript presentation.

The figure-generation workflow uses previously generated and locked analytical results rather than repeating model fitting or upstream statistical analyses.

The documented manuscript figure components include:

- **Figure 9A:** TCGA out-of-fold ROC curve.
- **Figure 9B:** GSE32863 external ROC curve.
- **Figure 9C:** GSE32863 paired-expression visualization.

The plotting workflow uses the corresponding STEP13G and STEP13I prediction and traceability outputs. Statistical annotations are transferred from the finalized upstream results rather than recalculated during figure generation.

The figures are exported as PDF and 600-dpi TIFF files under:

``` text
results/WORKFLOW13_Diagnostic_Classification/FIGURES/
```

## Key Classification Results and Interpretation

### Internal TCGA Assessment

The locked out-of-fold ROC result documented for Figure 9A is:

| Metric                  |             Value |
|-------------------------|------------------:|
| AUC                     |          0.987945 |
| 95% confidence interval | 0.977145–0.998745 |

This result represents internal resampling-based discrimination in the historical TCGA branch.

### External GSE32863 Assessment

The locked external ROC result documented for Figure 9B is:

| Metric                  |               Value |
|-------------------------|--------------------:|
| AUC                     |           0.9973246 |
| 95% confidence interval | 0.9929578–1.0000000 |

The external AUC describes rank-based discrimination from the transferred TCGA model in the GSE32863 dataset.

An important interpretation limitation is that all external predicted probabilities were below the prespecified 0.50 classification threshold, despite the high ROC AUC and the higher ranking of tumor samples relative to normal samples. At the fixed 0.50 threshold, the external classification sensitivity was 0 and specificity was 1.

Accordingly, the high external AUC should not be interpreted as evidence of satisfactory classification at the prespecified threshold or as proof of clinical diagnostic performance. Cross-dataset expression-scale compatibility is also relevant when interpreting transferred-model performance.

### Paired Expression Assessment

Figure 9C presents the paired expression patterns of AOC3, MS4A7, and IL7R in the 57 complete GSE32863 subject pairs.

The visualization uses the finalized paired-expression traceability data and previously generated statistical annotations. It does not independently establish the validity of the upstream statistical analysis.

## Reproducibility

The project uses project-relative paths resolved through `here::here()` and `file.path()`.

The workflow-specific random seed documented for the model assessment and subsequent classification steps is:

``` r
20260926
```

R session information is exported by the relevant workflow steps.

For detailed reproducibility information, refer to:

- `data/README.md` — input data sources, input files, and data organization.
- `scripts/README.md` — script structure, workflow implementation, and execution details.
- `results/README.md` — generated output directories, reports, tables, model objects, and figures.
- `Environment/README.md` — R version, package dependencies, and software environment.

## Scope Restrictions and Interpretation Boundaries

The classification workflow is designed as an exploratory transcriptomic tissue-discrimination analysis.

The following restrictions apply:

- No unsupported clinical diagnostic claims.
- No claim that the three genes constitute a causal or mechanistic signaling axis.
- No external model refitting or coefficient re-estimation in GSE32863.
- No external feature selection or gene selection.
- No external threshold optimization or cutpoint searching.
- No external model tuning.
- No artificial Normal–Tumor subject pairing.
- No inference of patient identity from sample order or sample counts.
- No interpretation of an audit PASS status as independent validation of every upstream analytical decision.

The internal cross-validation and external dataset assessment provide complementary exploratory evaluations. Neither establishes clinical diagnostic utility, prospective diagnostic performance, or clinical implementation readiness.

The results should therefore be interpreted as exploratory evidence of tumor-versus-normal tissue discrimination for the three prioritized study genes in the analyzed datasets.
