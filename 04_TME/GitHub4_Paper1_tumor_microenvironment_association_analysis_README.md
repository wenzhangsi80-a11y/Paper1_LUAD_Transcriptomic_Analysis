# Tumor Microenvironment (TME) Analysis

**Study:** Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD

## 1. Overview

This directory contains the analytical workflows, input-data references, computational environment information, and results associated with tumor microenvironment (TME) characterization in the LUAD study.

The analyses investigate the relationships between three prioritized study genes — AOC3, MS4A7, and IL7R — and ESTIMATE-derived TME scores in primary TCGA-LUAD tumor samples. The workflow combines gene–TME association analysis, ABSOLUTE tumor-purity-adjusted regression, historical ESTIMATE source and signature audits, ESTIMATE score reconstruction, and leave-one-gene-out (LOO) sensitivity analyses.

The TME module comprises three complementary R scripts:

1.  **Gene–TME association and purity-adjusted analysis:** evaluates associations between the three prioritized study genes and ESTIMATE-derived TME scores, examines associations after adjustment for ABSOLUTE-derived tumor purity, and audits historical ESTIMATE inputs and signature membership.
2.  **ESTIMATE reconstruction and LOO sensitivity analysis:** recovers and validates the original ESTIMATE signatures, expression inputs, and reference scores; reconstructs the original scores; and evaluates the sensitivity of the relevant ESTIMATE scores to the removal of AOC3 or IL7R from their corresponding signatures.
3.  **Final ABSOLUTE regression-table extraction and audit:** extracts coefficient statistics from previously fitted and archived STEP7B6 regression models, independently checks the extracted results against the archived results, and generates Supplementary Table 3a and supporting audit files. This script does not refit the regression models.

The three study genes were carried forward from preceding exploratory bulk-transcriptomic analyses. They should not be interpreted as the initial DEG-ranking results or as the PPI hub genes.

The analyses characterize statistical associations and computational score sensitivity in bulk-tissue transcriptomic data. They do not establish causal relationships, direct cellular interactions, or a mechanistic signaling axis.

## 2. Study Cohort and Data Sources

### 2.1 Primary TCGA-LUAD Cohort

The principal gene–TME association and ESTIMATE reconstruction/LOO analyses use 513 primary TCGA-LUAD tumor samples.

Gene-expression data, clinical/sample information, and ESTIMATE-related reference objects are obtained from previously established Paper 1 analyses. The workflows use locked upstream WGCNA and ESTIMATE objects where specified.

Sample identifiers and sample ordering are checked when expression matrices, TME scores, and tumor-purity information are aligned.

The ABSOLUTE purity-adjusted regression workflow uses a separate locked analysis dataset. The STEP7B5 dataset used by the archived STEP7B6 regression workflow is expected to contain 500 samples. This sample count is distinct from the 513-sample primary cohort used in the original gene–TME association analysis.

### 2.2 Prioritized Study Genes

The three genes investigated throughout the TME module are:

- AOC3
- MS4A7
- IL7R

Their roles differ between the association and signature-removal analyses. MS4A7 is retained in the broader gene–TME association and purity-adjusted regression analyses but is not a target of signature-specific LOO analysis because it is absent from both original ESTIMATE signatures according to the documented signature-membership audit.

### 2.3 ESTIMATE-Derived TME Measures

The gene–TME association and purity-adjusted regression analyses evaluate three ESTIMATE-derived measures:

- **StromalScore:** a transcriptome-derived measure related to stromal content.
- **ImmuneScore:** a transcriptome-derived measure related to immune content.
- **ESTIMATEScore:** a composite measure reflecting immune- and stromal-related signals.

ESTIMATEScore is a composite measure and should not be interpreted as a separate cellular compartment. These scores are computational estimates derived from bulk gene-expression data, not direct measurements of individual cell populations.

### 2.4 ABSOLUTE Tumor-Purity Information

The purity-adjusted analysis uses tumor-purity information from the TCGA PanCanAtlas ABSOLUTE source.

The original ABSOLUTE workflow includes source-file discovery, source-data inspection, identifier matching, barcode reconciliation, and purity-value integrity checks. The workflow constructs a locked purity-adjusted dataset from the eligible matched samples.

The STEP7B5 locked dataset used by the archived STEP7B6 regression models is expected to contain 500 samples. The actual dataset size and model observation counts must be confirmed using the saved dataset, model objects, and execution records.

Samples without an eligible ABSOLUTE purity value are handled according to the predefined upstream matching and exclusion procedures. The final extraction script does not rematch samples, remap barcodes, or impute missing purity values.

### 2.5 Historical ESTIMATE Reference Data

The ESTIMATE reconstruction and LOO workflow uses previously generated and locked reference data, including:

- Original ESTIMATE expression matrix: 58,387 genes × 513 samples.
- Original locked ESTIMATE score matrix: 3 × 513, comprising StromalScore, ImmuneScore, and ESTIMATEScore.
- Original StromalSignature and ImmuneSignature definitions.
- Official ESTIMATE common-gene expression input in GCT format.

The principal locked reference files include:

`WORKFLOW-03_Step14A-14E_FINAL_LOCKED.RData`

`STEP10B3NC_official_filterCommonGenes_output.gct`

`Paper1_STEP10B3ND2_official_SI_geneset_substitution_validation_LOCKED.RData`

These files provide the reference expression, score, signature, and validation objects required by the reconstruction and LOO workflow.

### 2.6 Archived ABSOLUTE Regression Inputs for Final Table Extraction

The standalone final extraction script uses two previously generated RData files:

- `Paper1_STEP7B5_locked_purity_adjusted_dataset.RData`
- `Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData`

The expected archived objects include:

- `purity_adjusted_dataset`
- `purity_adjusted_models`
- `purity_adjusted_results`
- `purity_adjusted_results_raw`
- `model_spec_7b6`

These archived objects provide the locked analysis dataset, fitted regression models, model specifications, and previously calculated regression results.

The final extraction script loads and checks these objects. It does not repeat the upstream ABSOLUTE sample matching or refit the archived regression models.

## 3. Analysis Scripts

The TME module contains three complementary R scripts.

| Script | Main analytical scope |
|----|----|
| `GitHub4_Paper1_tumor_microenvironment_association_analysis_1.R` | Gene–TME associations, manuscript table and figure preparation, ABSOLUTE source validation, purity-adjusted regression, and historical ESTIMATE source/signature audits |
| `GitHub4_Paper1_tumor_microenvironment_association_analysis_2.R` | ESTIMATE signature recovery, input and implementation audits, original-score reconstruction and validation, and formal AOC3-LOO and IL7R-LOO sensitivity analyses |
| `Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.R` | Extraction and independent auditing of coefficient statistics from archived STEP7B6 regression models; generation of Supplementary Table 3a and supporting audit files |

The first two scripts perform the original TME analytical workflows. The third script is a separate downstream extraction and audit workflow that uses the archived STEP7B5 dataset and STEP7B6 regression checkpoint.

The third script does not replace the original regression analysis and does not independently repeat model fitting.

Detailed script organization, workflow stages, and input/output relationships are documented in `scripts/README.md`.

## 4. Analysis Workflow

### 4.1 Gene–TME Association Analysis

The first script evaluates the associations between the expression of AOC3, MS4A7, and IL7R and the three ESTIMATE-derived TME measures.

The nine prespecified comparisons are:

| Study gene | TME measures                             |
|------------|------------------------------------------|
| AOC3       | StromalScore, ImmuneScore, ESTIMATEScore |
| MS4A7      | StromalScore, ImmuneScore, ESTIMATEScore |
| IL7R       | StromalScore, ImmuneScore, ESTIMATEScore |

The statistical procedures include:

- Two-sided Spearman rank-correlation tests.
- `exact = FALSE` in `cor.test()`.
- Complete-observation analysis for the aligned expression and TME variables.
- Validation of the expected 513-sample analytical cohort.
- Benjamini–Hochberg (BH) adjustment across the nine prespecified tests.

The nominal significance threshold is P \< 0.05, and the multiple-testing significance threshold is BH-FDR \< 0.05.

The principal analysis objects include:

- `core_tme_alignment`
- `gene_tme_results`
- `gene_tme_manuscript_table`

The original manuscript-ready correlation table is derived from these association results. Formatting and structural validation do not constitute additional statistical analyses.

### 4.2 Figure Generation

The first script constructs a three-by-three visualization of the nine gene–TME associations.

The intended panel organization is:

- Rows: AOC3, MS4A7, and IL7R.
- Columns: ImmuneScore, StromalScore, and ESTIMATEScore.
- Samples per panel: 513 primary LUAD tumors, subject to the corresponding data-integrity checks.

The visualization uses previously generated association results and their statistical annotations. The figure-generation workflow does not independently recalculate the nine gene–TME correlation tests.

The script documents the following figure outputs:

- `Figure6_CoreGenes_TME_Associations_600dpi_STEP6F.tiff`
- `Figure6_CoreGenes_TME_Associations_STEP6F.pdf`

The TIFF output is specified at 600 dpi with LZW compression, and the PDF output uses the `cairo_pdf` device.

**Manuscript figure-number verification is required.** The script-level documentation identifies this output as Figure 6, whereas earlier manuscript-level documentation identified the original gene–TME association figure as Figure 5. The final figure number should be checked against the actual figure files and the current manuscript before publication. The output filenames above are retained as documented by the script.

### 4.3 ABSOLUTE Source Validation and Purity-Adjusted Regression

The first script investigates the availability and integrity of ABSOLUTE-derived tumor-purity information for the 513 primary LUAD samples.

The source-validation workflow includes:

- Discovery and inspection of the relevant ABSOLUTE source file.
- Verification of required columns and data types.
- Inspection of missing values, duplicated identifiers, and purity-value integrity.
- Examination of call-status distributions.
- Matching ABSOLUTE records to primary TCGA-LUAD sample identifiers.
- Identification of unmatched samples and unavailable purity values.
- Investigation of potential barcode-level relationships without automatic remapping.

After sample matching, the gene–TME alignment data are combined with usable ABSOLUTE purity measurements to construct the purity-adjusted analysis dataset.

Nine multiple linear regression models are specified:

$$\text{TME trait} \sim \text{Study gene} + \text{ABSOLUTE\_Purity}$$

The models cover all three study genes and all three ESTIMATE-derived TME measures.

The original regression workflow includes:

- Extraction of gene-expression and tumor-purity coefficients.
- Reporting of model-level R-squared and adjusted R-squared.
- BH adjustment across the nine gene-level P values.
- Model-formula, sample-size, and result-integrity checks.
- Archiving of the fitted models and regression results.

The purity-adjusted analysis uses a separate locked dataset expected to contain 500 samples. The actual model sample sizes must be taken from the archived fitted models and execution records.

These regression models describe purity-adjusted associations. They do not establish causal effects.

### 4.4 Historical ESTIMATE Source and Signature-Membership Audits

The first script also examines historical ESTIMATE-related objects, source information, and direct membership of the prioritized study genes in the ESTIMATE signatures.

The audit includes:

- Inspection of available historical ESTIMATE input and intermediate objects.
- Recovery of historical objects in isolated R environments where appropriate.
- Inspection of the installed `estimate` package and available source information.
- Extraction and inspection of StromalSignature and ImmuneSignature gene sets.
- Direct membership checks for AOC3, MS4A7, and IL7R.
- Explicit distinction between MS4A6A and MS4A7.
- Comparison of signature membership with available historical input genes.
- Checks of signature-vector structure, gene counts, and duplicated entries.

The signature-membership workflow uses the `tidyestimate` package and its `gene_sets` data object when available.

The corresponding audit object is:

`step7be9a_audit`

The documented locked audit file is:

`Paper1_STEP7B7E9A_ESTIMATE_signature_membership_audit_LOCKED.RData`

The audit includes file-existence, file-size, and read-back checks where specified by the script.

Direct signature membership is relevant to interpreting ESTIMATE-derived scores, but membership alone does not quantify the magnitude or direction of an individual gene's contribution to a score. These historical-source and signature-membership audits do not themselves recalculate ESTIMATE scores or perform LOO analyses.

### 4.5 ESTIMATE Signature Recovery and Reference-Object Validation

The second script focuses on the recovery and validation of the original ESTIMATE scoring framework before conducting formal LOO analyses.

The workflow recovers the original signature definitions and locked Step14 reference objects, including:

- `expr_estimate`
- `estimate_scores_wgcna`
- Original StromalSignature and ImmuneSignature definitions.

The original signature definitions are also recovered and checked against the validated D-2 checkpoint:

`Paper1_STEP10B3ND2_official_SI_geneset_substitution_validation_LOCKED.RData`

The workflow inspects the official ESTIMATE implementation, including:

- `filterCommonGenes()`
- `estimateScore()`
- Package-level signature data and `SI_geneset` objects.

The specified ESTIMATE implementation reference is version 1.0.13. The installed package version and successful recovery or installation must be verified from the actual execution records.

### 4.6 Original ESTIMATE Score Reconstruction

The second script reconstructs the original StromalScore, ImmuneScore, and ESTIMATEScore using the recovered original signatures and filtered common-gene expression input.

The reconstructed scores are compared with the locked Step14 reference scores.

The validation includes:

- Expression and score matrix dimensions.
- Sample identity and ordering.
- Finite-value and matrix-integrity checks.
- Pearson and Spearman correlations.
- Maximum absolute difference.
- Mean absolute difference.
- Root mean squared error (RMSE).

The documented numerical reconstruction criterion includes a maximum absolute difference below $1 \times 10^{-8}$.

This stage is a reference-reconstruction and validation procedure. It is distinct from the formal LOO sensitivity analyses.

### 4.7 Formal Leave-One-Gene-Out Sensitivity Analysis

The formal LOO analyses evaluate the sensitivity of the predefined ESTIMATE signatures to the removal of two study genes.

The validated original signatures and filtered common-gene expression input are used as the reference for these analyses.

**AOC3-LOO**

AOC3 is removed from the original StromalSignature while the original ImmuneSignature is retained.

The workflow:

- Constructs and validates the reduced StromalSignature.
- Calculates ESTIMATE scores using the official `estimateScore()` implementation.
- Parses and validates the resulting score matrix.
- Checks the 513-sample identity and ordering.
- Compares the AOC3-LOO scores with the original locked Step14 reference scores.
- Evaluates numerical score differences and the number of samples with changed scores.
- Assesses the association between AOC3 expression and the AOC3-LOO StromalScore using Spearman correlation.

**IL7R-LOO**

IL7R is removed from the original ImmuneSignature while the original StromalSignature is retained.

The workflow:

- Constructs and validates the reduced ImmuneSignature.
- Calculates ESTIMATE scores using the official `estimateScore()` implementation.
- Parses and validates the resulting score matrix.
- Checks the 513-sample identity and ordering.
- Compares the IL7R-LOO scores with the original locked Step14 reference scores.
- Evaluates numerical score differences and the number of samples with changed scores.
- Assesses the association between IL7R expression and the IL7R-LOO ImmuneScore using Spearman correlation.

MS4A7 is not subjected to signature-specific LOO analysis because the validated signature-membership audit indicates that it is absent from both the original StromalSignature and ImmuneSignature.

The original ESTIMATE signature and reference objects are restored and audited after the LOO experiments to check for unintended changes.

The LOO procedures assess computational sensitivity to removing a gene from its corresponding predefined signature. They do not establish causal effects of that gene on the tumor microenvironment.

### 4.8 Final Extraction and Audit of Archived ABSOLUTE Purity-Adjusted Regression Results

The standalone script `Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.R` performs the final extraction and audit of the previously fitted STEP7B6 ABSOLUTE purity-adjusted regression results.

The workflow is identified as:

`STEP7B6_05`

It uses the archived input files:

- `Paper1_STEP7B5_locked_purity_adjusted_dataset.RData`
- `Paper1_STEP7B6_purity_adjusted_regression_checkpoint.RData`

The script loads the archived objects in isolated environments and checks the expected dataset, fitted models, archived regression results, and model specifications.

The expected analysis includes nine models covering:

- Study genes: AOC3, MS4A7, and IL7R.
- TME traits: ImmuneScore, StromalScore, and ESTIMATEScore.
- Adjustment covariate: ABSOLUTE_Purity.

The STEP7B5 locked dataset is expected to contain 500 samples. The extraction script verifies the sample identifiers, required variables, data completeness, and expected sample count. Each of the nine archived models is expected to contain 500 observations; the actual model observation counts and audit status must be confirmed from the generated audit outputs.

The script performs the following tasks:

1.  Loads and validates the archived STEP7B5 dataset.
2.  Loads and validates the archived STEP7B6 regression checkpoint.
3.  Checks the availability and structure of the nine archived fitted models.
4.  Extracts the gene-expression and ABSOLUTE-purity coefficients from the existing models.
5.  Extracts the corresponding standard errors, test statistics, and P values.
6.  Recalculates BH-adjusted P values across the nine gene-term P values.
7.  Compares the recalculated adjusted P values with the archived `Gene_FDR_BH` results.
8.  Audits dataset variables, model observation counts, and coefficient extraction.
9.  Exports Supplementary Table 3a and supporting audit files.

**This script extracts and audits existing fitted models; it does not refit the regression models.** It also does not repeat ABSOLUTE source discovery, sample matching, barcode reconciliation, or purity imputation. The original archived RData files are retained as inputs and are not modified by this extraction workflow.

The script-defined final audit includes a stringent condition requiring all nine gene-term BH-adjusted P values to be below 0.05, together with successful comparison against the archived BH-adjusted results. This is an audit criterion, not a statement that all nine associations have been confirmed as statistically significant. The actual PASS/FAIL status and statistical results must be taken from the generated output files.

The main publication output is:

`Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.csv`

An Excel version is also generated when the optional `openxlsx` package is available.

This supplementary table represents the ABSOLUTE purity-adjusted regression sensitivity analysis. It is distinct from the original ESTIMATE-based gene–TME correlation table and does not replace the original 513-sample Spearman correlation analysis.

## 5. Main Analysis Objects and Outputs

The following table summarizes the principal objects and output categories documented across the three scripts. Detailed filenames, directory locations, and checkpoint contents are maintained in `results/README.md`.

| Category | Main objects or outputs | Purpose |
|----|----|----|
| Locked WGCNA and ESTIMATE inputs | Upstream expression and ESTIMATE-related objects | Reference data for gene–TME association analysis |
| Gene–TME alignment | `core_tme_alignment` | Aligned gene-expression and TME data |
| Correlation results | `gene_tme_results` | Nine prespecified Spearman associations |
| Original manuscript correlation table | `gene_tme_manuscript_table` and associated exported table | Manuscript-ready ESTIMATE gene–TME correlation results |
| Figure output | TIFF and PDF figure files | Visualization of the nine gene–TME associations |
| ABSOLUTE source audit | Source-validation and sample-matching audit objects | Source integrity and sample-identifier reconciliation |
| Original purity-adjusted analysis | `purity_adjusted_dataset`, `purity_adjusted_models`, and regression-result objects | Nine gene–TME regression models adjusted for ABSOLUTE purity |
| Historical signature audit | `step7be9a_audit` | Historical ESTIMATE signature membership and provenance |
| Locked signature audit | `Paper1_STEP7B7E9A_ESTIMATE_signature_membership_audit_LOCKED.RData` | Serialized signature-membership audit |
| Original-score reconstruction | Reconstructed ESTIMATE score objects and comparison results | Validation against locked Step14 reference scores |
| AOC3-LOO output | `STEP10B3NE_AOC3_LOO_estimateScore.gct` | ESTIMATE scores after AOC3 removal from StromalSignature |
| IL7R-LOO output | `STEP10B3NE_IL7R_LOO_estimateScore.gct` | ESTIMATE scores after IL7R removal from ImmuneSignature |
| Final LOO checkpoint | `Paper1_STEP10B3NE_formal_LOO_sensitivity_analysis_LOCKED.RData` | Locked LOO scores, comparisons, associations, and integrity-audit objects |
| Final ABSOLUTE extraction | `Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.csv` | Supplementary Table 3a generated from archived STEP7B6 models |
| Final extraction audit | STEP7B6_05 audit CSV and TXT files | Verification of dataset variables, model structure, coefficient extraction, and BH adjustment |

The documented LOO output locations are:

``` text
10B-3N-E/
├── AOC3_LOO/
│   └── STEP10B3NE_AOC3_LOO_estimateScore.gct
├── IL7R_LOO/
│   └── STEP10B3NE_IL7R_LOO_estimateScore.gct
└── Audit/
    └── Paper1_STEP10B3NE_formal_LOO_sensitivity_analysis_LOCKED.RData
```

The final ABSOLUTE extraction workflow writes its outputs to:

``` text
04_TME/
└── results/
    └── STEP7B6_05_Final_Supplementary_Table/
        ├── Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.csv
        ├── Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.xlsx
        ├── STEP7B6_05_Model_Coefficient_Details.csv
        ├── STEP7B6_05_BH_Adjustment_Audit.csv
        ├── STEP7B6_05_Dataset_Variable_Audit.csv
        ├── STEP7B6_05_Model_Audit.csv
        └── STEP7B6_05_Final_Audit.txt
```

The Excel workbook, when generated, contains the following worksheets:

- `Supplementary_Table`
- `All_Coefficients`
- `BH_Audit`
- `Dataset_Audit`
- `Model_Audit`

These are script-documented output names and locations. Successful generation, file availability, and integrity should be confirmed from the execution logs and saved audit outputs.

## 6. Quality Control and Reproducibility

The three scripts incorporate structural, numerical, and source-integrity checks at different stages of the TME workflow.

The principal quality-control procedures include:

1.  Verification of upstream WGCNA and ESTIMATE object availability and structure.
2.  Validation of gene-expression and TME-score alignment.
3.  Confirmation of the expected 513-sample primary tumor cohort where applicable.
4.  Verification of the nine prespecified gene–TME comparisons and BH adjustment.
5.  Inspection of ABSOLUTE source columns, purity values, and sample identifiers.
6.  Validation of the original purity-adjusted regression dataset and fitted model results.
7.  Inspection of historical ESTIMATE objects, package information, and signature membership.
8.  Verification of the original ESTIMATE signature definitions and filtered common-gene input.
9.  Comparison of reconstructed original scores with locked Step14 reference scores.
10. Validation of reduced LOO signatures, score matrices, sample identity, and ordering.
11. Restoration and integrity checks of original signature and reference objects.
12. Read-back and validation of designated locked checkpoints.
13. Validation of the STEP7B5 dataset and the nine archived STEP7B6 regression models.
14. Recalculation and comparison of BH-adjusted gene-term P values in the final extraction workflow.
15. Verification of the final supplementary table and supporting audit files.

Project-relative paths are used for project input and output management through the project setup and `here::here()`.

The specified ESTIMATE implementation reference for the reconstruction and LOO workflow is version 1.0.13. The first script also uses packages including `openxlsx`, `ggplot2`, `dplyr`, `tidyr`, and `tidyestimate` for the documented export, visualization, and signature-membership tasks. The final extraction script requires `here` and uses the R `stats` and `utils` packages; `openxlsx` is optional for Excel workbook generation.

Exact R and package versions must be recorded from the actual execution environment. The second script header declares R 4.5.1; this declaration should not be treated as confirmation of the runtime version without the corresponding execution record. Similarly, package versions observed during later source or provenance audits should not automatically be interpreted as the versions used in historical analyses.

No random seed is specified for the deterministic correlation, regression, ESTIMATE score reconstruction, LOO scoring, or final archived-model extraction procedures documented in these workflows.

For detailed software requirements and version-recording information, see `environment/README.md`.

## 7. Interpretation and Analytical Limitations

The TME analyses should be interpreted within the following scope:

- **Observational associations:** Spearman correlations and purity-adjusted regression coefficients describe statistical associations between gene expression and ESTIMATE-derived scores. They do not establish causality.
- **Bulk-tissue estimates:** ESTIMATE scores represent transcriptome-derived immune- and stromal-related signals. They are not direct measurements of specific immune or stromal cell populations.
- **Tumor-purity adjustment:** ABSOLUTE-adjusted models account for the specified tumor-purity covariate but do not eliminate all potential confounding or dependence associated with bulk-tissue composition.
- **Separate analysis cohorts:** The original gene–TME association and ESTIMATE reconstruction/LOO workflows use the 513-sample primary cohort. The ABSOLUTE purity-adjusted regression workflow uses a separate locked dataset expected to contain 500 samples. These sample counts should not be treated as interchangeable.
- **Signature membership:** Direct membership of a gene in an ESTIMATE signature does not independently quantify its contribution to the corresponding score.
- **LOO sensitivity:** Removing AOC3 or IL7R from its corresponding signature evaluates the computational sensitivity of the predefined scoring procedure. It does not remove all possible sources of biological or statistical dependence.
- **MS4A7:** MS4A7 remains a prioritized study gene in the broader TME association and purity-adjusted regression analyses but is not a signature-specific LOO target because it is absent from both original ESTIMATE signatures according to the documented membership audit.
- **Final table extraction:** Supplementary Table 3a is generated by extracting and auditing coefficients from previously fitted STEP7B6 models. It is not the result of a new regression-fitting procedure.
- **Statistical significance:** The final extraction script's requirement that all nine gene-term BH-adjusted P values be below 0.05 is a predefined audit criterion. Actual significance must be assessed from the generated statistical results and audit report.
- **Historical provenance:** Package and source inspections describe the available audit environment and recovered materials. They do not independently prove the exact software environment used in every historical scoring run.

The module does not claim independent clinical validation, direct cellular interactions, or causal regulation of the tumor microenvironment.

Numerical findings, validation outcomes, and execution status should be reported only when supported by the corresponding generated outputs, audit objects, and execution records.

## 8. Directory Structure

The TME analysis module is organized as follows:

``` text
04_TME/
├── README.md
│
├── data/
│   ├── input_data/
│   └── sample_manifests/
│
├── environment/
│   └── README.md
│
├── results/
│   ├── STEP7B6_05_Final_Supplementary_Table/
│   └── README.md
│
└── scripts/
    ├── GitHub4_Paper1_tumor_microenvironment_association_analysis_1.R
    ├── GitHub4_Paper1_tumor_microenvironment_association_analysis_2.R
    ├── Supplementary_Table_ABSOLUTE_Purity_Adjusted_Regression.R
    └── README.md
```

The input files required by the final extraction workflow are expected under:

`data/input_data/`

The final extraction outputs are written under:

`results/STEP7B6_05_Final_Supplementary_Table/`

The corresponding scripts and their detailed workflow descriptions are maintained under `scripts/`.

Generated outputs from the original TME analyses may also be retained in their designated project output locations. Detailed output inventories, checkpoint descriptions, and result interpretation are provided in `results/README.md`.

## 9. Manuscript Correspondence

This module supports the manuscript analyses concerning ESTIMATE-derived TME associations, ABSOLUTE tumor-purity-adjusted associations, and signature-removal sensitivity assessment.

The original gene–TME correlation analysis and its manuscript-ready table should remain distinct from the additional ABSOLUTE purity-adjusted sensitivity analysis.

The new extraction workflow supports the following manuscript item:

**Methods 2.9. ABSOLUTE purity-adjusted gene–TME sensitivity analysis.**

**Results 3.6.** Report the actual purity-adjusted regression findings using the estimates and statistical results in the generated supplementary table. The wording of this Results subsection should be determined by the verified results rather than assuming that all nine associations are significant or consistent.

**Supplementary Table 3a. ABSOLUTE purity-adjusted linear regression analyses of AOC3, MS4A7, and IL7R in relation to ESTIMATE-derived tumor microenvironment scores in primary LUAD samples.**

Supplementary Table 3a is an additional sensitivity-analysis output derived from the archived STEP7B6 regression models. It should not be confused with the original ESTIMATE gene–TME correlation table.

Earlier manuscript-level documentation associated the original gene–TME correlation analysis with Figure 5 and Supplementary Table 3, and the signature-removal sensitivity analysis with Supplementary Table 5. However, the script-level documentation identifies the figure output as Figure 6. The final figure and supplementary-table numbering must therefore be reconciled with the current manuscript and final exported files before repository release. This README does not treat the initial manuscript numbering as independently verified.

The relevant manuscript Methods and Results sections should be updated to match the final analysis outputs and manuscript organization.

## 10. Related Documentation

The following README files provide complementary documentation for this module:

- `environment/README.md` — R environment, package dependencies, computational functions, version recording, and reproducibility requirements.
- `scripts/README.md` — Script organization, workflow stages, input dependencies, and analysis procedures.
- `results/README.md` — Generated result files, figures, tables, statistical outputs, audit reports, and locked checkpoints.

The root project README provides the overall analytical strategy, data scope, major workflow components, and interpretation limits. The three subordinate README files provide the corresponding environment, script, and output details.

## 11. Citation

If you use or adapt these analysis scripts, please cite the associated Paper 1 manuscript when it becomes available.
