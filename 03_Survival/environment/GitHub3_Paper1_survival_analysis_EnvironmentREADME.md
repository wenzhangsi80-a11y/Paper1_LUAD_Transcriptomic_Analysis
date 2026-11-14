# Environment and Reproducibility

## 1. Analysis Module

**Directory:** `03_Survival/`

**Study title:**

Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD

This document summarizes the computational environment, software dependencies, statistical settings, and reproducibility requirements for the four survival-analysis scripts in the `03_Survival/` module.

The module includes:

- Survival cohort reconstruction and quality-control auditing.
- Continuous-expression univariable Cox regression.
- Clinical covariate-adjusted multivariable Cox regression.
- Proportional-hazards diagnostics.
- Time-varying Cox regression for MS4A7.
- Three-gene survival reporting.
- Genome-wide univariable Cox regression.
- Kaplan–Meier estimation, log-rank testing, and grouped Cox regression.
- Statistical reporting and figure generation.

The analysis is implemented in R. Exact software versions must be recorded from the finalized local execution environment rather than inferred from individual script sections.

## 2. R Environment

| Component | Information |
|---|---|
| Programming language | R |
| R version | To be verified from the local R session |
| Operating system | To be recorded from the execution environment |
| Project path management | `here::here()` |
| Input and output paths | Project-relative paths |
| Session information | To be recorded after final execution |
| Random seed | No single module-wide random seed is specified in the reviewed scripts |

The project uses `here` for project-root detection and project-relative path construction.

The common directory objects include:

```r
project_root <- here::here()

input_dir <- file.path(project_root, "data")
result_dir <- file.path(project_root, "results")
figure_dir <- file.path(project_root, "figures")
```

Individual scripts may define additional workflow-specific directories under these project locations.

Local absolute computer paths should not be embedded in committed analysis scripts.

## 3. R Package Dependencies

The following packages are identified across the four survival-analysis scripts.

| Package | Main functions or functionality | Purpose |
|---|---|---|
| `here` | `here::here()` | Project-root detection and project-relative path management |
| `survival` | `Surv()`, `coxph()`, `cox.zph()`, `survfit()`, `survdiff()` | Survival objects, Cox regression, PH diagnostics, Kaplan–Meier estimation, and log-rank testing |
| `survminer` | Kaplan–Meier visualization functions | KM curves and number-at-risk tables |
| `ggplot2` | Plot construction and formatting | Forest plots, genome-wide survival visualization, and KM plot formatting |
| `dplyr` | Data manipulation | Clinical data preparation and result processing |
| `tidyr` | Data reshaping | Reporting-table preparation |
| `tibble` | Structured data frames | Statistical reporting tables |
| `readr` | Tabular data input and output | Data import and export |
| `stringr` | String manipulation | Text and identifier handling |
| `openxlsx` | Excel workbook generation | Excel-based result reporting |
| `grid` | Graphics-related functions | Plot layout, graphical annotations, and visualization support |

The table consolidates dependencies identified in the supplied script documentation. The actual package loading requirements should be confirmed against the complete scripts and finalized execution environment.

### 3.1 Base R and the `stats` Package

The workflows also use base R and functions from the recommended `stats` package.

Examples include:

**File and directory operations**

- `file.path()`
- `file.exists()`
- `file.info()`
- `dir.create()`
- `dir.exists()`
- `normalizePath()`

**Data handling and validation**

- `read.table()`
- `as.matrix()`
- `match()`
- `intersect()`
- `setdiff()`
- `setequal()`
- `identical()`
- `stopifnot()`
- `summary()`
- `quantile()`
- `table()`

**Statistical calculations**

- `stats::coef()`
- `stats::vcov()`
- `stats::AIC()`
- `stats::logLik()`
- `stats::sd()`
- `stats::p.adjust()`
- `stats::pnorm()`
- `stats::pchisq()`
- `stats::format.pval()`

**Checkpoint and reporting operations**

- `load()`
- `save()`
- `write.csv()`
- `sprintf()`

These functions are provided by base R or the recommended `stats` package and do not require separate installation as third-party packages.

## 4. Statistical Methods and Computational Settings

### 4.1 Primary Survival Cohort

The formally locked primary LUAD survival cohort contains:

| Characteristic | Number |
|---|---:|
| Total cases | 485 |
| Death events | 178 |
| Censored observations | 307 |

The original cohort is retained throughout the relevant analysis stages.

Model-specific complete-case filtering is applied only when required by the corresponding statistical model.

The principal clinical-adjustment cohorts are:

| Model | Covariates | Patients | Events |
|---|---|---:|---:|
| Model A | MS4A7 expression, age, gender, and clinical stage | 485 | 178 |
| Model B | MS4A7 expression, age, gender, and T/N/M stage | 336 | 137 |

The Model B cohort is a model-specific complete-case subset and does not redefine the original 485-case cohort.

### 4.2 Continuous-Expression Cox Regression

The targeted survival analyses evaluate the associations between AOC3, MS4A7, and IL7R expression and overall survival.

The main statistical procedures include:

- Z-score standardization of study-gene expression.
- Continuous-expression univariable Cox regression.
- Clinical covariate-adjusted Cox regression.
- Hazard ratio and 95% confidence interval estimation.
- Efron handling of tied event times where specified.
- Benjamini–Hochberg false discovery rate adjustment for the corresponding reporting analyses.

The continuous-expression hazard ratios represent the association with survival per unit increase in the standardized expression value.

### 4.3 Proportional-Hazards Diagnostics

The proportional-hazards assumption is evaluated using Schoenfeld residual-based diagnostics.

The principal function is:

```r
survival::cox.zph()
```

The clinical-adjustment workflow uses the Kaplan–Meier transformation:

```r
transform = "km"
```

These diagnostics are applied to the relevant constant-coefficient Cox models.

The selected time-varying Cox model containing `tt()` terms is not subjected to `cox.zph()` in the described workflow.

### 4.4 Adjusted Time-Varying Cox Regression

The MS4A7 time-varying analysis evaluates three prespecified time functions:

1. Linear: `MS4A7_z × t`
2. Square-root: `MS4A7_z × sqrt(t)`
3. Logarithmic: `MS4A7_z × log1p(t)`

The models adjust for age, gender, and clinical stage.

The computational procedures include:

- Efron handling of tied event times.
- AIC-based comparison of candidate time functions.
- Likelihood-ratio comparison with the constant-coefficient model.
- Estimation of time-specific adjusted hazard ratios.
- Calculation of corresponding 95% confidence intervals.

Time-specific estimates are evaluated at 1, 2, 3, 4, and 5 years.

Previously generated model objects, comparison statistics, and PH diagnostic results are preserved in the corresponding checkpoint.

### 4.5 Genome-Wide Univariable Cox Analysis

The genome-wide survival analysis uses:

| Setting | Specification |
|---|---|
| Survival cohort | 485 cases |
| Death events | 178 |
| Censored observations | 307 |
| Gene universe | 31,609 MAD-positive WGCNA candidate genes |
| Expression input | Previously prepared log2-transformed expression values |
| Standardization | Gene-wise standardization within the 485-patient cohort |
| Statistical model | Gene-wise univariable Cox proportional-hazards regression |
| Multiple testing | Genome-wide BH-FDR adjustment |

The hazard ratio represents the association per one-standard-deviation increase in gene expression.

Genome-wide proportional-hazards testing is not described as part of the STEP 6X-8C model-fitting procedure.

Subsequent checkpoint recovery, result auditing, data freezing, and visualization stages do not refit the original genome-wide Cox models unless explicitly specified by the corresponding script.

### 4.6 Kaplan–Meier and Grouped Cox Analysis

The STEP 6E workflow performs exploratory survival analyses using median-defined expression groups for AOC3, MS4A7, and IL7R.

The grouping rule is:

- Low: expression less than or equal to the gene-specific median.
- High: expression greater than the gene-specific median.

The procedures include:

- Kaplan–Meier survival estimation.
- Confidence interval display.
- Number-at-risk table generation.
- Two-sided log-rank testing.
- Grouped Cox regression comparing High versus Low expression.
- Hazard ratio and 95% confidence interval estimation.

Four observations with `OS.time = 0` are retained. Their original survival times are not modified.

The grouped Cox hazard ratios are distinct from the continuous-expression Cox estimates.

No automatic re-standardization of the previously standardized gene-expression values is performed in STEP 6E.

## 5. Multiple Testing and Reporting

The survival analysis uses Benjamini–Hochberg false discovery rate adjustment in the relevant reporting branches.

The adjustment procedures are applied separately according to the analysis design:

- Three-gene targeted Cox reporting: BH-FDR adjustment within each corresponding reporting model.
- Genome-wide univariable Cox analysis: BH-FDR adjustment across the tested gene universe.

The three-gene reporting results and genome-wide results are maintained as distinct statistical analyses.

Plotting and figure refinement are kept separate from the underlying statistical estimates. Visual adjustments should not be interpreted as additional model fitting or changes to the original gene-level results.

## 6. Graphics and Export Settings

The reporting workflows generate survival figures and statistical visualizations.

### 6.1 Figure 4 Compact Forest Plot

The Figure 4 forest-plot workflow uses the export settings specified in its analysis script.

| Parameter | Setting |
|---|---|
| PDF | Cairo PDF |
| TIFF | 600 dpi, LZW compression |
| PNG | 600 dpi |
| JPEG | 600 dpi, quality 100 |
| Figure dimensions | 8.8 × 4.6 inches |

### 6.2 Figure X Genome-Wide Survival Visualization

The final Figure X visualization uses:

| Parameter | Setting |
|---|---|
| PDF | Cairo PDF |
| TIFF | 600 dpi, LZW compression |
| PNG | 600 dpi |
| JPEG | 600 dpi, quality 100 |
| Figure dimensions | 8.8 × 5.8 inches |
| FDR reference threshold | 0.05 |
| log2HR reference | 0 |
| Corresponding −log10(FDR) reference | −log10(0.05) |
| STEP 6X-9B x-axis limits | −0.75 to 0.75 |
| STEP 6X-9B y-axis limits | 0 to 4.60 |
| Core-gene display order | AOC3, MS4A7, IL7R |

The final visual refinement changes label positions and visual annotations without refitting the Cox models or intentionally changing the original gene-level statistical estimates.

### 6.3 Kaplan–Meier Figures

The STEP 6E workflow exports gene-specific KM figures in TIFF and PDF formats.

The plots include confidence intervals, number-at-risk tables, and statistical annotations as defined by the corresponding plotting workflow.

Exact filenames and output paths are documented in `results/README.md`.

## 7. Random Seed and Stochastic Procedures

No single random seed is specified across all four survival-analysis scripts.

The reviewed STEP 6E procedures do not explicitly invoke randomized resampling.

For other analysis sections, a seed should be documented only when random sampling, resampling, or another stochastic procedure is explicitly implemented.

A random seed should not be inferred or retrospectively assigned to analysis sections that do not specify one.

Where applicable, section-specific seeds should be recorded in the corresponding script and execution record.

## 8. Checkpoint and Data Management

The survival-analysis workflow reuses previously generated WGCNA and survival-analysis checkpoints.

Important upstream and intermediate checkpoints include:

- `Workflow-04-2D_final_WGCNA_candidate_checkpoint.RData`
- `Paper1_Step6B5_continuous_Cox_PH_checkpoint.RData`
- `Paper1_Step6B6_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6C5_adjusted_MS4A7_time_varying_Cox_checkpoint.RData`
- `Paper1_Step6D1_Survival_Reporting_checkpoint.RData`

Checkpoint files are stored in their designated project directories, including:

```text
02_WGCNA/results/statistical_outputs/
03_Survival/results/statistical_outputs/
```

The scripts use checkpoints to recover previously fitted models, preserve selected model specifications, verify cohort consistency, and prepare downstream reporting objects.

Checkpoint recovery and reporting stages should not be interpreted as new model fitting unless explicitly stated in the corresponding analysis script.

The existence of a checkpoint file alone does not establish successful quality control. The associated audit results and PASS/FAIL status should also be inspected.

## 9. Package and Environment Verification

Before reproducing the survival analyses, verify the following:

1. The required R packages are installed.
2. The project root is correctly detected by `here::here()`.
3. The required upstream WGCNA and survival-analysis checkpoints are available.
4. The locked 485-case survival cohort can be recovered.
5. The clinical covariates and expression data can be aligned using the specified sample identifiers.
6. The genome-wide expression matrix can be aligned to the primary survival cohort.
7. The required output directories can be created.
8. The plotting devices required for PDF and raster exports are available.
9. The relevant model and checkpoint quality-control checks pass.

The following commands can be used to record the actual software environment:

```r
R.version.string

packageVersion("here")
packageVersion("survival")
packageVersion("survminer")
packageVersion("ggplot2")
packageVersion("dplyr")
packageVersion("tidyr")
packageVersion("tibble")
packageVersion("readr")
packageVersion("stringr")
packageVersion("openxlsx")

sessionInfo()
```

Package versions should be recorded only for packages actually installed and used in the finalized execution environment.

The operating system, execution date, relevant input checkpoint versions, and any section-specific random seeds should also be documented.

## 10. Reproducibility Notes

- All analysis paths should be constructed relative to the repository root.
- Absolute local computer paths should not be included in committed scripts.
- The original 485-case survival cohort should be retained unless a model explicitly requires a smaller complete-case subset.
- Model-specific filtering should not redefine the original locked cohort.
- Previously standardized expression values should not be standardized again unless explicitly required by a separate analysis.
- The four zero-time observations should remain unchanged in the STEP 6E workflow.
- Continuous-expression Cox estimates and median-group survival estimates should be reported separately.
- Previously fitted statistical models should not be unnecessarily refitted during checkpoint recovery or figure preparation.
- Actual software versions, checkpoint locations, output files, and QC outcomes should be verified from the local execution environment.
- No unverified package versions, random seeds, or execution outcomes are reported in this document.

## 11. Related Documentation

- `03_Survival/README.md`: Overall study module description, analysis objectives, cohort definitions, and workflow organization.
- `03_Survival/scripts/README.md`: Individual script descriptions, execution stages, and input dependencies.
- `03_Survival/results/README.md`: Detailed output-directory structure, checkpoint descriptions, and result-file inventory.