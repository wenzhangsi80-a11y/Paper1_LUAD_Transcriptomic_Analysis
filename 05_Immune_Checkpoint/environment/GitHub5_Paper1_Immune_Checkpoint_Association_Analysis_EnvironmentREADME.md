# Environment

## 1. Software Environment

### Programming Language

- Programming language: R
- R version: To be recorded from the actual execution environment.
- RStudio: Optional; record the version if used.
- Operating system: To be recorded from the actual execution environment.
- Project root: Identified using `here::here()`.
- Path management: Project-relative paths using `here::here()` and `file.path()`.

The exact software environment should be documented using the output of `sessionInfo()` after the complete workflow has been executed.

### Recommended Environment Record

The following information should be added after execution:

| Item                    | Environment information                    |
|-------------------------|--------------------------------------------|
| R version               | To be confirmed                            |
| RStudio version         | To be confirmed, if applicable             |
| Operating system        | To be confirmed                            |
| Platform                | To be confirmed from `sessionInfo()`       |
| Package versions        | To be recorded from `sessionInfo()`        |
| Project path management | `here`                                     |
| Random seed             | Not required for the documented procedures |

No software or package version is inferred from the script files alone.

------------------------------------------------------------------------

## 2. Required R Packages and Dependencies

The two analysis scripts use a combination of base R, standard R statistical functions, and additional packages for path management, gene annotation, data manipulation, statistical result handling, and visualization.

### 2.1 Core and Statistical Dependencies

| Package | Requirement | Purpose |
|----|----|----|
| Base R | Required | Data import, matrix extraction, sample matching, data manipulation, quality control, summaries, and file operations |
| stats | R recommended package | Spearman correlation tests, linear regression, residual extraction, ranking, and BH-FDR adjustment |
| grDevices | R recommended package | PDF and TIFF graphics devices |
| here | Required | Project-root identification and project-relative path management |

Functions used in the workflow include `cor.test()`, `rank()`, `lm()`, `residuals()`, `p.adjust()`, `merge()`, `match()`, `complete.cases()`, `write.csv()`, `saveRDS()`, `file.path()`, and `file.exists()`.

### 2.2 Gene Annotation Dependencies

| Package | Requirement | Purpose |
|----|----|----|
| AnnotationDbi | Required for annotation mapping | Ensembl-to-gene-symbol annotation and identifier mapping |
| org.Hs.eg.db | Required for human gene annotation | Human gene annotation database |
| biomaRt | Optional / availability check | Alternative annotation resource availability check |

Gene annotation is performed separately from the original expression matrix. The original expression data are retained without modification.

### 2.3 Data Manipulation and Statistical Result Handling

| Package | Requirement | Purpose |
|----|----|----|
| dplyr | Required by the relevant workflow sections | Data-frame selection, filtering, transformation, and summary preparation |
| tidyr | Required by the visualization/data-preparation sections that use it | Data reshaping and preparation |
| readr | Required by the relevant workflow sections | Tabular data import and export |
| openxlsx | Conditional dependency | Excel workbook generation |

The final STEP 8E Excel workbook is generated only when `openxlsx` is available. CSV and RDS outputs are generated using base R functions.

### 2.4 Visualization Dependencies

| Package | Requirement | Purpose |
|----|----|----|
| ggplot2 | Required for figure generation | Scatter plots, heatmaps, graphical customization, and figure construction |
| cowplot | Required for the relevant figure-layout sections | Combining plot panels, legends, and figure components |

The visualization workflow uses the locked STEP 8C-7B results. It does not introduce new correlation calculations or statistical tests.

### Package Version Documentation

Package versions should be recorded from the actual execution environment rather than inferred from package names or script contents.

Recommended command:

``` r
sessionInfo()
```

For a focused package-version record:

``` r
packages_used <- c(
  "here",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "biomaRt",
  "dplyr",
  "tidyr",
  "readr",
  "openxlsx",
  "ggplot2",
  "cowplot"
)

package_versions <- data.frame(
  Package = packages_used,
  Version = vapply(
    packages_used,
    function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(packageVersion(pkg))
      } else {
        NA_character_
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(package_versions)
```

The resulting package-version table should be updated with the versions actually used for the final analysis.

------------------------------------------------------------------------

## 3. Project Directory Structure

The analysis uses project-relative paths. The project root is identified with `here::here()`, and file paths are constructed using `file.path()`.

### Input Directory

``` text
data/
└── processed_data/
    └── Paper1_tumor-only WGCNA/
```

The principal input files and intermediate objects include:

| Input file | Purpose |
|----|----|
| `Workflow-04-2C_WGCNA_candidate_log2_FPKM_UQ.txt` | Gene-expression matrix used for the immune checkpoint association workflow |
| `Workflow-04-2A_clinical_aligned.txt` | Clinical information aligned with the expression data |
| `Paper1_WORKFLOW03_14F_4D_2E_FINAL.RData` | Existing workspace containing the locked ESTIMATE-related object used in STEP 8C recovery |
| `Paper1_END_WORKSPACE_20260920_20260920_015542.RData` | Workspace checkpoint used by the STEP 8D sensitivity workflow |
| `Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds` | Locked primary 54-pair association results used in downstream sensitivity analyses and visualization |

The workspace files serve different workflow stages and should not be treated as interchangeable.

### Output Directory

The main statistical output directory is:

``` text
results/
└── statistical_outputs/
    └── Paper1_tumor-only WGCNA/
```

The workflow uses subdirectories under this project-defined output location, including:

- STEP 8A–8C data preparation, association analysis, and result extraction outputs.
- `Workflow-08-STEP8C8_Supplementary_Figure4/` for Supplementary Figure 4 visualization outputs.
- `Workflow-08-STEP8E_Leave-One-Core-Gene-Out_Audit/` for the final leave-one-core-gene-out audit outputs.

The exact output paths are defined by the project-level path variables in the scripts. The scripts should be executed from the configured project structure rather than relying on machine-specific absolute paths.

------------------------------------------------------------------------

## 4. Analysis Cohort and Computational Design

The immune checkpoint association workflow uses 503 TCGA-LUAD primary tumor samples with the required expression and adjustment-variable data.

The predefined association panel consists of:

- Three prioritized study genes: `AOC3`, `MS4A7`, and `IL7R`.
- Eighteen immune checkpoint-related genes.
- Fifty-four predefined gene pairs (3 × 18).
- A target sample size of 503 for each association pair.

The 18 immune checkpoint-related genes are:

``` text
PDCD1, CD274, PDCD1LG2, CTLA4, CD80, CD86,
LAG3, HAVCR2, TIGIT, CD96, VSIR, VTCN1,
BTLA, TNFRSF4, TNFRSF9, ICOS, LGALS9, CD200
```

The same predefined gene-pair structure and sample identifiers are retained across the three association-analysis layers.

The separate 500-primary-tumor cohort used in the gene–TME purity analysis elsewhere in the project is not the cohort used in this 503-sample immune checkpoint workflow.

------------------------------------------------------------------------

## 5. Statistical Computing Environment

### 5.1 Unadjusted Spearman Correlation

The primary unadjusted analysis evaluates the 54 predefined gene pairs using two-sided Spearman rank-correlation tests.

- Correlation method: Spearman.
- R implementation: `stats::cor.test()`.
- Exact calculation setting: `exact = FALSE`.
- Multiple-testing correction: Benjamini–Hochberg (BH).
- Correction scope: all 54 predefined gene-pair tests.

### 5.2 ImmuneScore-Adjusted Partial Spearman Correlation

The primary adjustment analysis uses a rank-based residual procedure:

1.  Rank the expression values of the two genes and the ImmuneScore.
2.  Regress each ranked expression variable on the ranked ImmuneScore.
3.  Calculate the Pearson correlation between the two residual vectors.
4.  Apply BH-FDR correction across the 54 predefined tests.

This procedure uses standard R statistical functions, including `rank()`, `lm()`, `residuals()`, `cor.test()`, and `p.adjust()`.

### 5.3 TumorPurity-Adjusted Sensitivity Analysis

The separate TumorPurity-adjusted sensitivity analysis uses the same rank-based residual-correlation procedure, with ESTIMATE-derived TumorPurity as the adjustment variable.

TumorPurity is calculated using:

``` r
TumorPurity <- cos(
  0.6049872018 + 0.0001467884 * ESTIMATEScore
)
```

The analysis retains the predefined 54 gene pairs and the 503-sample cohort. BH-FDR correction is applied across the 54 tests.

This sensitivity analysis does not replace the primary unadjusted or ImmuneScore-adjusted analyses.

### 5.4 Leave-One-Core-Gene-Out Dependence Audit

The STEP 8E audit uses previously locked association results to assess the dependence of the overall association pattern on the composition of the three-core-gene panel.

The three predefined subsets are:

| Subset      | Omitted gene | Retained core genes | Number of pairs |
|-------------|--------------|---------------------|----------------:|
| Leave_AOC3  | AOC3         | MS4A7, IL7R         |              36 |
| Leave_MS4A7 | MS4A7        | AOC3, IL7R          |              36 |
| Leave_IL7R  | IL7R         | AOC3, MS4A7         |              36 |

This is a descriptive dependence audit of existing results.

- No new correlation coefficients are calculated.
- No new P values are calculated.
- No new FDR values are calculated.
- No additional multiple-testing correction is introduced.
- The locked primary association results are retained.

The audit should not be interpreted as independent replication or a new inferential analysis.

------------------------------------------------------------------------

## 6. Randomness and Reproducibility Settings

No random seed is required for the procedures documented in these two scripts.

The documented operations consist of deterministic data processing, sample alignment, correlation testing, regression-based residual calculations, BH-FDR adjustment, descriptive summaries, quality-control checks, and figure export.

No random sampling or stochastic model-fitting procedure is specified in these workflow sections.

The patient-level bootstrap or other stochastic procedures, if used in separate project analyses, are outside the scope of this Environment README.

------------------------------------------------------------------------

## 7. Graphics and Figure Export Environment

The STEP 8C-8 visualization workflow generates figure panels from the locked STEP 8C-7B supplementary association results.

### Plotting Framework

- `ggplot2`: plot construction and graphical customization.
- `dplyr`: data preparation and filtering.
- `cowplot`: figure composition and layout support.

### Graphics Settings

| Setting | Specification |
|----|----|
| Font family | Arial |
| PDF graphics device | `grDevices::cairo_pdf()` |
| TIFF graphics device | `grDevices::tiff()` |
| TIFF resolution | 600 dpi |
| TIFF compression | LZW |
| Figure dimensions | Explicitly specified in the corresponding script section |

Panel-specific settings include:

| Figure component | Width (inches) | Height (inches) |
|----|---:|---:|
| Panel A | 7.2 | 6.0 |
| Panel B | 9.0 | 5.8 |
| Panel C | To be confirmed from the corresponding script section |  |

The Panel B heatmap palette is defined as:

- Low: `#82B5EC`
- Midpoint: `white`
- High: `#F57C6E`

The visualization modules use the locked STEP 8C-7B association table:

``` text
Workflow-08-STEP8C7B_Supplementary_Table_Immune_Checkpoint_54.rds
```

No new correlations, P values, or FDR values are calculated during figure construction. The previously locked statistical results are retained.

------------------------------------------------------------------------

## 8. Workspace Checkpoint and Session Documentation

The workflow includes a workspace checkpoint procedure to support continuation and reproducibility.

The checkpoint documentation may include:

- Complete R workspace (`.RData`).
- R session information (`.txt`).
- Workspace object inventory (`.txt`).
- Workflow continuation README (`.txt`).

These files document the computational state and support subsequent workflow execution. They are not additional statistical analyses.

The workspace checkpoint should be distinguished from the final analysis outputs and should not be treated as evidence that all downstream quality-control checks have passed.

------------------------------------------------------------------------

## 9. Reproducibility Notes

1.  **Project-relative paths:** The project root is identified using `here::here()`, and paths are constructed using `file.path()`. Machine-specific absolute paths should not be required.

2.  **Input preservation:** Original expression, clinical, and previously locked result files are retained. The workflow does not overwrite the original input data.

3.  **Sample alignment:** Sample identifiers are explicitly matched across expression matrices, clinical information, and ESTIMATE-derived variables. The intended analysis cohort contains 503 samples.

4.  **Predefined association panel:** The three prioritized study genes and 18 immune checkpoint-related genes define 54 gene pairs. The predefined panel is retained throughout the analysis.

5.  **Locked primary results:** The STEP 8C-7B results serve as the reference for downstream sensitivity analyses and visualization. The tumor-purity sensitivity analysis is separate from the primary ImmuneScore-adjusted analysis.

6.  **No redundant statistical recalculation:** STEP 8C-6 through STEP 8C-8 primarily extract, summarize, verify, or visualize previously generated results. STEP 8E extracts and audits locked results without recalculating correlations, P values, or FDR values.

7.  **Conditional dependencies:** The availability of `openxlsx` determines whether the optional Excel workbook can be generated. CSV and RDS exports use base R functions.

8.  **Software version transparency:** R, RStudio (if applicable), operating-system, and package versions should be recorded from the actual execution environment. No version should be inferred from the script header or an unrelated project workflow.

9.  **Execution-dependent verification:** Package availability, input-file availability, output-file generation, file integrity, and final QC status must be checked after execution. A predefined QC target or an expected output filename does not by itself demonstrate successful execution.

10. **Interpretation boundary:** The workflow characterizes bulk-tumor expression associations and their sensitivity to adjustment variables. It does not establish direct molecular interactions, causal regulation, cell-intrinsic checkpoint activity, functional immune effects, or therapeutic response.

------------------------------------------------------------------------

## 10. Scope of This Environment README

This document covers the computational environment and dependencies for the two immune checkpoint association scripts:

| Script | Workflow coverage | Main computational components |
|----|----|----|
| `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_1.R` | STEP 8A through STEP 8C-8 | Input and cohort QC, gene annotation, expression extraction, unadjusted and ImmuneScore-adjusted associations, result extraction, supplementary table preparation, and figure generation |
| `GitHub5_Paper1_Immune_Checkpoint_Association_Analysis_2.R` | STEP 8D through STEP 8E | TumorPurity-adjusted sensitivity analysis, three-layer comparison, and leave-one-core-gene-out dependence audit |

This README documents the software and computational requirements of the immune checkpoint association module. Detailed analytical procedures and generated results are described separately in the module-level `Scripts/README.md` and `Results/README.md` files.
