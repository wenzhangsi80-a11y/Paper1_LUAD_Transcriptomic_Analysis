------------------------------------------------------------------------

editor_options: markdown: wrap: 72 ---

# Differential Expression Assessment of AOC3, MS4A7, and IL7R in GSE32863

## Overview

This module performs an independent paired expression assessment of AOC3, MS4A7, and IL7R using the Gene Expression Omnibus (GEO) dataset GSE32863 and its GPL6884 platform annotation.

The analysis supports the independent expression assessment presented in Figure 2D and the corresponding description in Section 2.2, Differential Expression Analysis, of the associated manuscript.

The analysis focuses on the three prioritized study genes and evaluates tumor--normal expression differences using matched lung tumor and adjacent non-tumor lung tissue samples from the same subjects.

Two scripts are provided:

| Script | Purpose |
|:---|:---|
| `GitHub1_Paper1_differential_expression_GSE32863.R` | Probe annotation, target-gene probe selection, paired expression analysis, statistical testing, and Figure 2D generation |
| `GitHub1A_GSE32863_Sample_Manifest_Export.R` | Sample-manifest generation, complete-pair identification, singleton documentation, and sample ID consistency audits |

These scripts have separate responsibilities. The sample-manifest export script does not perform expression analysis, statistical testing, or figure generation.

This analysis does not perform diagnostic classification, ROC analysis, prediction-model fitting, or external validation of a prediction model.

## Data Sources

The expression assessment requires the following files downloaded from GEO:

| File | Description |
|:---|:---|
| `GSE32863_series_matrix.txt.gz` | GSE32863 series matrix containing expression data and sample-level metadata |
| `GPL6884-11607.txt` | GPL6884 platform annotation used for probe-to-gene mapping |

GEO records:

- GSE32863: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32863>
- GPL6884: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6884>

Place both files in `data/input_data/` within the project directory, retaining the filenames shown above.

## Cohort and Sample Structure

The audited GSE32863 cohort comprises:

- 116 samples in total;
- 58 lung tumor samples;
- 58 adjacent non-tumor lung samples;
- 59 unique subjects;
- 57 subjects with complete tumor--normal pairs;
- one tumor-only subject;
- one normal-only subject.

The two singleton samples are retained in the complete sample manifest but excluded from paired statistical testing.

Accordingly, the paired expression analysis uses 57 complete pairs (114 samples).

## Analysis Workflow

### 1. Sample-manifest generation and audit

Run:

``` r
scripts/01_DEG/GitHub1A_GSE32863_Sample_Manifest_Export.R
```

This script:

1.  Loads the GSE32863 Series Matrix and extracts sample-level metadata.
2.  Uses GEO accession identifiers, subject identifiers, and the original tissue annotations to construct the complete sample manifest.
3.  Identifies the 57 subjects with exactly one tumor and one adjacent-normal sample.
4.  Identifies the two singleton samples that cannot be included in paired testing.
5.  Exports the complete sample manifest, complete-pair manifest, and singleton exclusion manifest.
6.  Checks cohort counts and the consistency of GEO sample identifiers across the exported manifests.

The script generates sample-manifest audit reports in `results/audit_reports/`.

### 2. Probe annotation and target-gene expression assessment

Run:

``` r
scripts/01_DEG/GitHub1_Paper1_differential_expression_GSE32863.R
```

This script:

1.  Loads the GSE32863 expression matrix and sample-level metadata.
2.  Loads the GPL6884 platform annotation and maps expression probes to gene symbols.
3.  Aligns sample metadata to the expression matrix using GEO accession identifiers.
4.  Identifies complete tumor--normal pairs using subject identifiers.
5.  For each target gene, retains the probe with the highest mean expression across all 116 samples.
6.  Extracts the selected-probe expression values for the 57 complete pairs.
7.  Evaluates tumor--normal differences using two-sided paired Wilcoxon signed-rank tests.
8.  Generates Figure 2D and exports the corresponding expression, probe-selection, and statistical output files.

Probe selection is performed across all 116 samples, whereas paired expression comparisons and statistical tests use only the 57 complete pairs.

Unmatched expression probes and probe-annotation matching information are documented in the corresponding audit output. Probes without valid gene-symbol annotations are not used for target-gene selection.

## Statistical Analysis

Tumor and matched adjacent non-tumor expression values are compared using two-sided paired Wilcoxon signed-rank tests.

The statistical output records the number of complete pairs, group-specific expression summaries, paired differences, and the corresponding test results.

Figure annotations display unadjusted P values. If Benjamini--Hochberg-adjusted P values are included in the final statistical output, they represent an additional multiple-testing summary across the three target-gene tests and should be interpreted separately from the P values displayed in the figure.

## Figure

The principal figure generated by this module is:

- `Figure2D_GSE32863.pdf`
- `Figure2D_GSE32863.tiff`

The figure presents the expression distributions of AOC3, MS4A7, and IL7R in the 57 matched tumor--normal pairs.

The figure corresponds to panel D of Figure 2 in the manuscript.

## Software Requirements

The analysis was developed for R 4.5.1 and is intended to run in a compatible R 4.5.x environment.

Required packages for the Figure 2D expression-analysis script:

- `here`
- `GEOquery`
- `Biobase`
- `ggplot2`

The sample-manifest export script uses `here`, `GEOquery`, and `Biobase`.

All required packages must be installed before running the corresponding script. The analysis scripts check package availability but do not install or update packages automatically.

See `environment/README.md` for software environment information.

## Running the Scripts

Run the scripts from the root directory of the GitHub repository.

First, generate and audit the sample manifests:

``` r
source(
  "scripts/01_DEG/GitHub1A_GSE32863_Sample_Manifest_Export.R"
)
```

Then, perform the paired expression assessment and generate Figure 2D:

``` r
source(
  "scripts/01_DEG/GitHub1_Paper1_differential_expression_GSE32863.R"
)
```

Both scripts use `here::here()` to identify the project root and construct input and output paths relative to that directory.

## Output Organization

The principal outputs are organized as follows:

| Directory | Contents |
|:---|:---|
| `data/sample_manifests/` | Complete sample manifest, 57 complete-pair manifest, and singleton exclusion manifest |
| `results/figures/` | Figure 2D in PDF and TIFF formats |
| `results/tables/` | Probe-selection records and expression data tables |
| `results/statistical_outputs/` | Paired statistical results, analysis summaries, and package-version records |
| `results/audit_reports/` | Sample-manifest summary and cross-file sample ID audit reports |
| `environment/` | R session information and software environment documentation |

See `data/README.md` and `results/README.md` for detailed file inventories.

## Reproducibility Notes

The expression assessment uses the input expression values provided in the GSE32863 Series Matrix. The script does not introduce additional normalization, scaling, centering, z-score transformation, or cross-platform correction.

The figure y-axis is labelled neutrally as "Expression"; the analysis does not assume a specific normalized log2 scale beyond the values supplied in the Series Matrix.

Probe selection, sample matching, and paired testing are explicitly implemented in the analysis script. The sample manifests and audit reports provide traceability for the cohort composition and the 57 complete pairs used in the statistical analysis.

The reported expression differences describe associations within this dataset and should not be interpreted as evidence of causal mechanisms or clinical diagnostic performance.
