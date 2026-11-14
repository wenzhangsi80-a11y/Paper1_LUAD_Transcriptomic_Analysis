# Functional Enrichment Analysis

## 1. Overview

This directory contains the functional enrichment analysis workflows for the study:

**Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD**

The purpose of this module is to characterize biological processes and pathways associated with the transcriptomic findings in lung adenocarcinoma (LUAD), with particular attention to the three prioritized study genes:

- **AOC3**
- **MS4A7**
- **IL7R**

The functional enrichment analysis integrates differential-expression-based over-representation analysis (ORA), genome-wide preranked gene set enrichment analysis (GSEA), gene identifier and background-universe audits, pathway redundancy assessment, and descriptive cross-database evidence integration.

The module is implemented in two sequential R scripts:

1.  `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`
2.  `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

Together, these scripts provide a structured workflow from differential-expression-based enrichment analysis to genome-wide pathway analysis and final descriptive evidence review.

The analyses are designed to preserve the provenance of input gene lists, identifier-mapping decisions, ranking definitions, enrichment parameters, and downstream audit records.

## 2. Directory Organization

The Functional Enrichment module is organized as follows:

``` text
06_Functional_Enrichment/
│
├── README.md
│
├── data/
│   ├── data_sources/
│   ├── input_data/
│   ├── processed_data/
│   └── sample_manifests/
│
├── environment/
│   └── README.md
│
├── scripts/
│   ├── README.md
│   └── GitHub6_Paper1_Functional_Enrichment_Analysis_1.R
│   └── GitHub6_Paper1_Functional_Enrichment_Analysis_2.R
│
└── results/
    └── README.md
```

The directory structure above describes the intended organization of the module. Individual workflow outputs are stored in their corresponding result subdirectories.

For detailed file inventories and workflow-specific output locations, see:

- `scripts/README.md` — script inventory, analysis stages, major inputs, and workflow dependencies.
- `results/README.md` — result-directory organization and key output files.
- `environment/README.md` — R environment, package requirements, annotation resources, and reproducibility information.

## 3. Overall Analysis Workflow

The functional enrichment analysis is organized into two connected stages.

### 3.1. Script 1 — Differential-Expression-Based Functional Enrichment and Integration

**Script:** `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`

This script implements the initial functional enrichment workflow, covering STEP9A–STEP9C and the associated STEP9D-1 input-audit procedures.

The major analytical components include:

1.  Preparation and verification of differential-expression-derived gene sets.
2.  Gene identifier mapping and annotation assessment.
3.  Definition and auditing of enrichment input sets and background universes.
4.  KEGG over-representation analysis.
5.  GO Biological Process (GO-BP) over-representation analysis.
6.  GO semantic similarity and redundancy assessment.
7.  Integration and comparison of GO-BP and KEGG enrichment findings.
8.  Descriptive biological-theme assessment for the prioritized study genes.
9.  Genome-wide input and background-universe provenance checks.
10. Preparation and locking of inputs for subsequent genome-wide preranked GSEA.

This stage establishes the input definitions, identifier conventions, enrichment results, and provenance information required by the downstream analysis.

The workflow distinguishes between gene-list-based enrichment analysis and subsequent genome-wide ranking-based enrichment analysis. These analytical approaches are complementary but are not interchangeable.

### 3.2. Script 2 — Genome-Wide Preranked GSEA and Final Evidence Integration

**Script:** `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

This script continues from the locked inputs and audit records established by Script 1. It covers STEP9D-1I through STEP9D-5C.

The major analytical components include:

1.  Verification of the differential-expression source and genome-wide gene universe.
2.  Inspection and locking of the GSEA ranking metric.
3.  Construction and auditing of the genome-wide preranked vector.
4.  Ensembl-to-Entrez mapping diagnostics and KEGG collision-policy assessment.
5.  Formal genome-wide KEGG GSEA.
6.  KEGG GSEA result and reproducibility audits.
7.  GO-BP annotation coverage and input-policy assessment.
8.  Formal genome-wide GO-BP GSEA.
9.  GO-BP output and provenance verification, including designated output repairs.
10. Comparison of historical Reactome reference observations with the locked KEGG and GO-BP results.
11. Strict direct and contextual pathway-theme evidence assessment.
12. Final evidence-consistency, metadata, and concordance-lock procedures.

The script distinguishes formal enrichment execution from read-only audits, output repairs, metadata repairs, and descriptive evidence integration.

The final Reactome-related workflows organize existing evidence and evaluate biological-theme correspondence. They do not establish pathway equivalence, direct molecular interactions, or causal relationships.

## 4. Analytical Framework

### 4.1. Over-Representation Analysis

The initial enrichment workflow uses differential-expression-derived gene sets to identify over-represented biological pathways.

The principal enrichment resources include:

- KEGG pathways.
- Gene Ontology Biological Process (GO-BP).

The analyses use explicitly defined input gene sets and background universes. Identifier conversion, annotation coverage, and input provenance are assessed as part of the workflow.

The KEGG and GO-BP analyses use their respective identifier requirements and annotation resources. Their input identifiers and mapping procedures should not be treated as interchangeable.

### 4.2. Genome-Wide Preranked GSEA

The subsequent genome-wide GSEA workflow uses a locked ranking vector derived from the differential-expression results.

The ranking statistic is defined as:

$$R_i =
\operatorname{sign}(\log_2FC_i)
\times
\left[-\log_{10}\left(\max(P_i,10^{-300})\right)\right]$$

where:

- $\log_2FC_i$ is the gene-level log2 fold change.
- $P_i$ is the nominal differential-expression P value.
- $10^{-300}$ is the lower numerical floor applied to the P value.

The ranking metric uses nominal P values rather than adjusted P values.

The KEGG and GO-BP GSEA workflows apply their respective identifier and annotation requirements:

- **KEGG GSEA:** Entrez-compatible identifiers, with an explicitly documented collision-resolution policy.
- **GO-BP GSEA:** Ensembl identifiers with human annotation resources.

The ranking definition and upstream input provenance are retained for downstream reproducibility audits.

### 4.3. Identifier Mapping and Collision Resolution

Identifier mapping is an important component of the genome-wide enrichment workflow.

The analysis distinguishes:

- One-to-one mappings.
- One-to-many mappings.
- Many-to-one mappings.
- Unmapped identifiers.

For KEGG GSEA, the documented final collision-resolution policy is **Policy B — Maximum Absolute Ranking**.

Within an Entrez collision group, the Ensembl gene with the largest absolute value of the previously locked ranking statistic is retained.

The original genome-wide Ensembl ranking vector is not recalculated as part of this collision-resolution procedure.

GO-BP and KEGG identifier handling are documented separately because their annotation requirements differ.

### 4.4. Pathway Redundancy and Cross-Database Integration

GO semantic similarity and related integration procedures are used to assess overlapping or redundant biological annotations.

These analyses help organize enrichment findings into interpretable biological themes.

The integration of KEGG and GO-BP results is descriptive and should not be interpreted as independent replication, because the analyses share the underlying TCGA-LUAD biological context and related transcriptomic inputs.

### 4.5. Reactome Theme Evidence Assessment

The final workflow compares historical Reactome reference observations with the locked genome-wide KEGG and GO-BP results.

The assessment focuses on:

- IL7R
- MS4A7
- AOC3

The workflow distinguishes between broad theme concordance and stricter direct/contextual theme evidence.

The strict theme-evidence procedure uses predefined gene-specific keyword dictionaries applied to pathway descriptions. Leading-edge and core-enrichment information may be retained for supporting inspection but do not determine the strict description-based theme classification.

The final evidence review is descriptive. It does not generate a new enrichment test, biological-strength ranking, composite gene score, or causal model.

## 5. Key Analysis Definitions and Audit References

The functional enrichment workflows retain explicit input definitions, annotation policies, and audit reference values.

The following values are documented as workflow-specific reference criteria:

| Analysis component                   | Documented reference |
|--------------------------------------|---------------------:|
| Genome-wide Ensembl ranking universe |         57,828 genes |
| KEGG GSEA result pathways            |                  353 |
| KEGG pathways with FDR \< 0.05       |                  178 |
| GO-BP GSEA result pathways           |                6,316 |
| GO-BP pathways with FDR \< 0.05      |                2,003 |

These values are used in the corresponding workflow audits and consistency checks.

**They should be interpreted as documented reference or audit values, not as independent confirmation that every workflow has successfully executed or passed its final audit.** Actual execution status must be established from the corresponding generated reports, saved objects, and output-verification records.

Other workflow-specific values, including input-set sizes, annotation coverage, and gene-level evidence counts, are documented in the relevant result and audit files.

## 6. Principal Outputs

The Functional Enrichment module generates several categories of outputs.

### 6.1. Differential-Expression-Based Enrichment Results

These include:

- KEGG ORA result tables.
- GO-BP ORA result tables.
- Enrichment significance summaries.
- Gene identifier and background-universe audit records.
- GO semantic similarity and redundancy-assessment outputs.
- Integrated GO-BP and KEGG evidence tables.

### 6.2. Genome-Wide GSEA Results

These include:

- Locked genome-wide ranking vectors.
- KEGG-compatible ranked input tables.
- GO-BP ranked input tables.
- KEGG GSEA result objects and tables.
- GO-BP GSEA result objects and tables.
- FDR-significant pathway tables.
- GSEA parameter and provenance records.
- Result-integrity and reproducibility audit reports.

### 6.3. Reactome and Theme-Evidence Outputs

These include:

- Original Reactome reference information.
- Reactome–KEGG–GO-BP theme-concordance summaries.
- Strict direct/contextual theme-evidence tables.
- Gene-level descriptive evidence summaries.
- STEP9D-5B evidence-audit objects.
- STEP9D-5C final theme-evidence and concordance-lock records.

The detailed filenames and their corresponding result directories are indexed in `results/README.md`.

## 7. Reproducibility and Quality-Control Principles

The following principles apply throughout the Functional Enrichment module.

### 7.1. Input and Provenance Preservation

- Preserve the provenance of differential-expression input files.
- Maintain explicit gene identifier and background-universe definitions.
- Retain locked ranking vectors and their associated metric definitions.
- Document the annotation resources and mapping policies used by each enrichment procedure.
- Keep formal analysis outputs distinguishable from subsequent audit and repair outputs.

### 7.2. Analytical Immutability

Downstream audit and integration steps should not silently alter upstream analytical definitions.

In particular:

- Do not recalculate differential-expression results during downstream GSEA audits.
- Do not redefine the locked genome-wide ranking metric.
- Do not introduce artificial jitter or arbitrary tie-breaking.
- Do not silently change identifier-mapping or collision-resolution policies.
- Do not modify locked source results during read-only audits.
- Document designated file-repair and metadata-repair operations separately from formal enrichment execution.

### 7.3. Audit Status and Execution Verification

The scripts contain multiple levels of structural, numerical, provenance, and reproducibility checks.

These checks may include:

- Required input-file availability.
- Input and output dimensions.
- Gene identifier uniqueness and mapping consistency.
- Ranking-vector integrity.
- Enrichment-result structure and numerical validity.
- RDS–CSV concordance.
- Provenance consistency.
- Method-immutability checks.
- Final audit and output-reload verification.

A predefined audit criterion or expected output filename does not establish a successful execution.

Final PASS/FAIL status must be taken from the actual execution output and the corresponding saved audit records.

## 8. Interpretation and Limitations

The functional enrichment results provide pathway-level and biological-process-level context for the transcriptomic findings in LUAD.

Several interpretive limitations should be considered:

1.  ORA and GSEA address different enrichment questions and should be interpreted according to their respective input definitions and statistical frameworks.
2.  KEGG and GO-BP use different annotation systems and identifier requirements.
3.  Cross-database pathway correspondence does not establish identical biological mechanisms.
4.  Keyword-based theme matching depends on pathway-description terminology and predefined dictionaries.
5.  Historical Reactome observations and the revised KEGG/GO-BP analyses share the underlying TCGA-LUAD biological context and should not be described as independent external validation.
6.  Enrichment results indicate statistical representation of gene sets or pathways; they do not establish causal regulation or direct molecular interactions.
7.  These analyses do not independently establish clinical diagnostic or prognostic biomarker validity.

The prioritized genes AOC3, MS4A7, and IL7R are evaluated within the broader transcriptomic and pathway context. Their enrichment associations should not be interpreted as proof of gene-specific causal mechanisms.

## 9. Software and Environment

The analysis scripts were developed with a target R environment of R 4.5.1.

The workflows use R packages for:

- Gene identifier mapping and annotation.
- KEGG and GO enrichment analysis.
- Preranked GSEA.
- GO semantic similarity.
- Data manipulation and result processing.
- Visualization and audit reporting.

Important annotation resources include:

- GENCODE v22.
- `org.Hs.eg.db`.
- `AnnotationDbi`.
- Ensembl-based gene annotation resources.

The actual R version, installed package versions, and runtime environment should be recorded from the execution environment rather than inferred from the target version.

For detailed package requirements and environment documentation, see `environment/README.md`.

## 10. Running the Workflows

Run the scripts from the root directory of the GitHub repository.

The recommended order is:

**Step 1 — Initial functional enrichment and input preparation**

``` r
source(
  "06_Functional_Enrichment/scripts/GitHub6_Paper1_Functional_Enrichment_Analysis_1.R"
)
```

**Step 2 — Genome-wide GSEA and final evidence integration**

``` r
source(
  "06_Functional_Enrichment/scripts/GitHub6_Paper1_Functional_Enrichment_Analysis_2.R"
)
```

Script 2 depends on the relevant locked inputs and outputs generated or established by the preceding workflow. It should therefore be run only after the required upstream files and objects have been prepared.

The scripts use project-relative path management. Run them within the expected repository structure and ensure that all required input files and R packages are available.

Consult `scripts/README.md` for the detailed workflow sequence and dependencies before execution.

## 11. Documentation Guide

The four README files serve complementary purposes:

| Documentation file | Main purpose |
|----|----|
| `README.md` | Overall study context, analytical framework, two-script workflow, major results, and interpretation |
| `scripts/README.md` | Script inventory, workflow stages, inputs, outputs, and execution dependencies |
| `results/README.md` | Result-directory organization, principal result files, figures, tables, and audit outputs |
| `environment/README.md` | R version, package requirements, annotation resources, and reproducibility environment |

This separation avoids duplicating detailed script descriptions, file inventories, and package lists in the module-level README.

## 12. Summary

The Functional Enrichment module integrates differential-expression-based ORA, genome-wide preranked GSEA, gene identifier and annotation audits, pathway redundancy assessment, and descriptive Reactome–KEGG–GO-BP evidence integration.

The two scripts form a sequential workflow:

- **Script 1** establishes and audits the differential-expression-based enrichment results, gene identifiers, background definitions, and downstream input provenance.
- **Script 2** performs genome-wide ranking-based enrichment analysis, verifies the KEGG and GO-BP results, and organizes the pathway-theme evidence for the prioritized study genes.

Together, these workflows provide a reproducible framework for describing the functional and pathway-level context of AOC3, MS4A7, and IL7R in LUAD, while maintaining explicit distinctions between statistical enrichment, descriptive evidence integration, and biological interpretation.
