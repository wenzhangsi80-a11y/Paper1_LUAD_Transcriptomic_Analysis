# Software Environment

## GSE32863 Figure 2D Expression Assessment

This document describes the software environment required to reproduce the GSE32863 paired expression assessment and sample-manifest generation.

## 1. R Version

The analysis was developed and tested using:

- R 4.5.1
- A compatible R 4.5.x environment

The R version and platform information should be retained with the generated session information.

## 2. Required R Packages

### 2.1 Figure 2D Expression-Analysis Script

The following four packages are required by:

`scripts/01_DEG/GitHub1_Paper1_differential_expression_GSE32863.R`

| Package | Purpose |
|---|---|
| `here` | Project-root detection and reproducible relative paths |
| `GEOquery` | Loading the GEO Series Matrix |
| `Biobase` | Accessing ExpressionSet expression and phenotype data |
| `ggplot2` | Statistical visualization and Figure 2D generation |

### 2.2 Sample-Manifest Export Script

The following packages are required by:

`scripts/01_DEG/GitHub1A_GSE32863_Sample_Manifest_Export.R`

| Package | Purpose |
|---|---|
| `here` | Project-root detection and relative paths |
| `GEOquery` | Loading the GSE32863 Series Matrix |
| `Biobase` | Accessing sample-level phenotype data |

The sample-manifest export script does not require `ggplot2`.

No additional plotting or data-manipulation packages are required by these two scripts as currently specified.

## 3. Package Installation

Install CRAN packages:

    install.packages(c(
      "here",
      "ggplot2"
    ))

Install Bioconductor packages:

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    BiocManager::install(c(
      "GEOquery",
      "Biobase"
    ))

Package installation is intentionally separated from the analysis scripts. The scripts check whether their required packages are available and stop if dependencies are missing. They do not install or update packages automatically.

## 4. Environment Recording

The analysis workflow records software information to support reproducibility.

Generated environment and package records include:

- `environment/GSE32863_sessionInfo.txt`
- `results/01_DEG/statistical_outputs/GSE32863_package_versions.csv`

The session information should document the R version, platform, and loaded package information. The package-version table records the versions of the principal packages used in the analysis.

The sample-manifest and cross-file audit reports are stored in:

`results/01_DEG/audit_reports/`

These records should be retained with the corresponding analysis outputs.

## 5. Reproducibility Recommendations

For reproducible execution:

1. Use R 4.5.1 or a compatible R 4.5.x environment.
2. Install the required packages before running the scripts.
3. Preserve the original GSE32863 Series Matrix and GPL6884 annotation files in `data/input_data/`.
4. Run the scripts from the root directory of the GitHub repository.
5. Generate and check the sample manifests before running the paired expression assessment.
6. Retain the generated session information, package-version records, and audit reports with the analysis outputs.

Package versions are recorded at execution time rather than hard-coded in advance.

A formal `renv.lock` file should only be generated after the final analysis environment has been installed, tested, and confirmed.