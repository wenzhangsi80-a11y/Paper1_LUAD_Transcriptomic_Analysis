# Functional Enrichment Analysis — Environment and Software Requirements

## 1. Overview

This directory documents the computational environment, software dependencies, annotation resources, and reproducibility requirements for the Functional Enrichment Analysis module of Paper 1.

**Study title:**  
*Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD*

The Functional Enrichment Analysis module integrates two R scripts:

- `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`
- `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

Together, these scripts implement the following analytical components:

1. Input-data inventory and quality-control assessment.
2. Differentially expressed gene (DEG) result review and input preparation.
3. Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation analysis (ORA).
4. GO–KEGG functional integration and semantic redundancy assessment.
5. Genome-wide preranked KEGG and GO Biological Process (GO-BP) GSEA.
6. Annotation coverage, identifier-mapping, and ranking-integrity audits.
7. Historical Reactome reference and genome-wide enrichment evidence comparison.
8. Descriptive theme-evidence integration and final reproducibility audits.

This README consolidates the environment information documented across both analysis scripts. It describes the intended software environment and the dependencies used by the workflows; it does not claim that every package is required by every individual analysis step or that all documented analyses have been independently rerun in a newly established environment.

---

## 2. Project Structure and Path Management

The module follows the project-level directory structure:

```text
06_Functional_Enrichment/
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
├── results/
│   ├── STEP9A_Input_Inventory_and_QC/
│   ├── STEP9B_DEG_Analysis/
│   ├── STEP9C_Overrepresentation_Analysis/
│   └── STEP9D_GenomeWide_Preranked_GSEA/
│
└── scripts/
    ├── GitHub6_Paper1_Functional_Enrichment_Analysis_1.R
    └── GitHub6_Paper1_Functional_Enrichment_Analysis_2.R
```

### 2.1 Project-relative paths

The project uses `here` for project-root identification and project-relative path construction.

The preferred approach is:

```r
library(here)

project_root <- here::here()

functional_enrichment_dir <- file.path(
  project_root,
  "06_Functional_Enrichment"
)
```

Subdirectories should be constructed using `file.path()` rather than machine-specific absolute paths.

The scripts may also use shared directory variables, including `step9d_root`, `step9d_dir`, and step-specific directory variables. These variables must be initialized by the preceding workflow setup before dependent sections are executed.

### 2.2 Path reproducibility

The following conventions should be maintained:

- Avoid machine-specific absolute paths in the formal GitHub version.
- Use `here::here()` to identify the project root.
- Use `file.path()` to construct input and output paths.
- Preserve the established directory and output-file names.
- Do not assume that a script section relying on previously initialized objects or path variables can be executed independently.
- Keep external source-data locations and their provenance documented separately from the project-relative output structure.

---

## 3. R Environment

### 3.1 Target R version

| Item | Specification |
|---|---|
| Target/project R version | R 4.5.1 |
| Actual execution version | To be confirmed from runtime records |
| Bioconductor version | To be confirmed from runtime records |
| Operating system | Not specified in the available environment records |
| Project path management | `here` |
| File and directory operations | Base R, `utils`, and `tools` |

R 4.5.1 is the target version documented for this project. It should not be interpreted as proof that every analysis section was executed under this exact runtime version.

The actual R version should be recorded using:

```r
R.version.string
```

The installed R and package environment should be documented using:

```r
sessionInfo()
```

Where available, the Bioconductor version should also be recorded.

### 3.2 Runtime-version documentation

Package versions should be obtained from the actual R session rather than inferred from package names, script headers, or the current target R version.

For important packages, the following commands can be used:

```r
packageVersion("here")
packageVersion("data.table")
packageVersion("clusterProfiler")
packageVersion("org.Hs.eg.db")
packageVersion("AnnotationDbi")
packageVersion("GOSemSim")
```

Additional package versions should be recorded when the corresponding packages are installed and used.

The `packageVersion()` function should only be called for packages available in the current R library.

No fixed package-version numbers are asserted in this README because a complete, verified package-version inventory from the actual execution environment has not been established for all workflow sections.

---

## 4. R Package Dependencies

The packages are grouped according to their documented roles in the two analysis scripts. Their inclusion in this table does not mean that every package is required for every workflow step.

### 4.1 Core project and data-handling packages

| Package | Main role | Relevant workflow components |
|---|---|---|
| `here` | Project-root identification and project-relative path management | Both scripts; STEP9A–STEP9D |
| `data.table` | Efficient import and handling of large tabular datasets | Input inventory and data preparation |
| Base R | General data processing, object manipulation, validation, and file operations | Both scripts |
| `utils` | Reading and writing tables, package-version inspection, and other utility functions | Data handling and audit outputs |
| `tools` | File and path-related utility operations | File inventory and validation |

### 4.2 Functional enrichment and annotation packages

| Package | Main role | Relevant workflow components |
|---|---|---|
| `clusterProfiler` | GO/KEGG ORA, identifier conversion, and preranked GSEA | STEP9C and STEP9D |
| `org.Hs.eg.db` | Human gene annotation and identifier-to-annotation relationships | GO/KEGG annotation and mapping |
| `AnnotationDbi` | Access to annotation databases and identifier-mapping operations | Ensembl, Entrez, SYMBOL, and GO annotation audits |
| `GOSemSim` | GO semantic similarity calculations | GO-BP semantic redundancy analysis |

The principal documented enrichment functions include:

- `clusterProfiler::bitr()`
- `clusterProfiler::enrichGO()`
- `clusterProfiler::enrichKEGG()`
- `clusterProfiler::GSEA()`
- `clusterProfiler::gseGO()`
- `GOSemSim::godata()`
- `GOSemSim::mgoSim()`

The exact functions required depend on the specific analysis step.

### 4.3 Visualization packages

| Package | Main role | Relevant workflow components |
|---|---|---|
| `enrichplot` | Visualization of functional enrichment results | Enrichment visualization |
| `ggplot2` | General graphical presentation of enrichment results | Functional enrichment figures |

These packages support visualization and should not be interpreted as necessary for every data-processing or audit-only section.

### 4.4 Theme-evidence and concordance packages

| Package | Main role | Relevant workflow components |
|---|---|---|
| `dplyr` | Data-frame manipulation and evidence-table organization | STEP9D-5 |
| `stringr` | Text and keyword matching in pathway descriptions | STEP9D-5 and STEP9D-5B |
| `readr` | Tabular input and output operations | STEP9D-5 |
| `tibble` | Structured tabular data handling | STEP9D-5 |
| `purrr` | Functional iteration and data-processing operations | STEP9D-5 |

These packages are associated particularly with the historical Reactome reference comparison and the direct/contextual theme-evidence workflow.

### 4.5 Optional or diagnostic packages

| Package | Documented role | Important qualification |
|---|---|---|
| `biomaRt` | Optional Ensembl gene-biotype and annotation diagnostics | Depends on external Ensembl service availability and release information |
| `fgsea` | Availability and implementation inspection in STEP9D-2A; may also be relevant to the installed GSEA backend | Its inspection in STEP9D-2A is not itself a GSEA execution |

The presence of a package in a software-availability check does not establish that it is directly called by every formal analysis step.

In particular:

- `biomaRt` is an optional diagnostic resource and is not used to redefine the locked genome-wide ranking.
- STEP9D-2A inspects the availability and implementation of `clusterProfiler` and `fgsea`; this inspection must not be described as an enrichment analysis.
- The actual package dependencies of a formal GSEA execution should be interpreted in the context of the installed `clusterProfiler` version and the functions used by that execution.

---

## 5. Annotation Resources and Identifier Systems

### 5.1 Human gene annotation

The module uses human gene annotation resources, including:

- GENCODE v22 Ensembl-to-GeneSymbol annotation.
- Existing Paper 1 gene-annotation files generated during the preceding WGCNA workflow.
- `org.Hs.eg.db`.
- `AnnotationDbi`.
- Optional Ensembl annotation information accessed through `biomaRt`.

An existing annotation file used in the upstream workflow is:

```text
Workflow-04-2B-2_top100_annotated.txt
```

This file is part of the previously established Paper 1 annotation workflow and should not be interpreted as a newly generated annotation resource for every functional enrichment step.

### 5.2 Identifier systems

Different enrichment workflows use different identifier systems.

| Analysis component | Identifier system | Annotation context |
|---|---|---|
| Upstream genome-wide ranking | Ensembl gene identifiers | Locked DEG-derived ranking |
| GO-BP GSEA | ENSEMBL | `org.Hs.eg.db` |
| KEGG ORA | Entrez Gene IDs | KEGG-compatible identifier mapping |
| KEGG GSEA | Entrez Gene IDs | Locked Ensembl-to-Entrez mapping and collision-resolution policy |
| Gene-level evidence review | Gene symbols | AOC3, MS4A7, and IL7R |

Identifier conversion and mapping policies are part of the analytical provenance. They should not be silently reconstructed or replaced when downstream steps use previously locked input objects.

### 5.3 Annotation database provenance

The following information should be recorded whenever available:

- R version.
- Bioconductor version.
- `org.Hs.eg.db` version.
- `AnnotationDbi` version.
- `clusterProfiler` version.
- Ensembl release or annotation provenance.
- `biomaRt` availability and relevant Ensembl service information when the optional diagnostic is used.
- Identifier-mapping policy and the corresponding upstream lock-file provenance.

Annotation results may differ between database releases. Therefore, annotation package versions and resource provenance are relevant components of reproducibility.

---

## 6. Software Dependencies by Workflow

The following table summarizes the principal environment requirements across the module.

| Workflow | Main software dependencies | Main computational purpose |
|---|---|---|
| STEP9A — Input Inventory and QC | `here`, `data.table`, Base R, `utils`, `tools` | Input-file inventory, structure inspection, and QC |
| STEP9B — DEG Analysis and Input Preparation | `here`, Base R and relevant upstream data-handling packages | DEG result preparation and ranking-input provenance |
| STEP9C — GO/KEGG ORA | `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `here` | GO/KEGG over-representation analysis and identifier mapping |
| STEP9C — GO Semantic Redundancy | `GOSemSim`, `org.Hs.eg.db`, `AnnotationDbi`, Base R | GO semantic similarity and descriptive clustering |
| STEP9C — GO–KEGG Integration | `here`, Base R, `utils` | Integration and audit of previously generated enrichment results |
| STEP9D-1 — Input and Ranking Provenance Audits | `here`, Base R, `utils`, `tools` | Input inspection, lock verification, and provenance checks |
| STEP9D-2 — Ranking Metric and Software Audits | `here`, Base R; `clusterProfiler` and `fgsea` for software inspection | Ranking methodology and implementation inspection |
| STEP9D-3 — KEGG GSEA | `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `here` | Entrez-based KEGG preranked GSEA and result auditing |
| STEP9D-4 — GO-BP GSEA | `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `here` | Ensembl-based GO-BP preranked GSEA and result auditing |
| STEP9D-5 — Reactome Theme Concordance | `dplyr`, `stringr`, `readr`, `tibble`, `purrr`, Base R | Descriptive pathway-theme comparison and evidence integration |

This table describes the documented package roles and workflow dependencies. The precise package set required to execute an individual script section should be determined from its actual function calls and execution environment.

---

## 7. Randomness and Seed Management

### 7.1 Random-seed policy

The reviewed Environment README sections do not specify a global random seed for the formal Functional Enrichment workflows.

The documented operations primarily involve:

- Deterministic file and object inspection.
- Identifier mapping and annotation.
- Enrichment-result processing.
- Ranking-vector validation.
- GO semantic similarity calculations.
- Descriptive keyword matching.
- CSV/RDS serialization and reload verification.

No artificial jitter or artificial tie-breaking is specified for the formal KEGG or GO-BP ranking procedures.

### 7.2 Ranking ties

The formal genome-wide ranking workflows preserve ranking ties rather than introducing artificial numerical perturbations.

The documented restrictions include:

- No jitter.
- No artificial tie-breaking.
- No modification of locked ranking values.
- No recalculation of the original ranking metric during downstream result audits.

For the formal KEGG and GO-BP GSEA workflows, the ranking inputs are inherited from the upstream locked ranking and identifier-mapping procedures.

### 7.3 Reproducibility interpretation

The absence of an explicitly specified seed in the reviewed sections should not be interpreted as proof that every possible operation in every package version is independent of stochastic implementation details.

If a future workflow introduces a stochastic operation, its seed and relevant software versions should be documented in that workflow's own provenance record.

---

## 8. Locked Analytical Inputs and Environment Dependencies

The Functional Enrichment module relies on previously generated and locked upstream inputs.

Important examples include:

- DEG result tables and the associated source-definition records.
- STEP9D-1 input and provenance locks.
- STEP9D-2B ranking metric lock.
- STEP9D-2C genome-wide ranking vector.
- STEP9D-3A KEGG identifier-mapping and input-policy records.
- STEP9D-4A GO-BP input and annotation-policy lock.
- Existing KEGG and GO-BP GSEA result objects and provenance records.
- Historical Reactome reference information used in STEP9D-5.

These inputs should be treated as authoritative according to their respective workflow definitions.

### 8.1 Ranking and annotation provenance

The following distinctions are important:

- The genome-wide ranking is inherited from the upstream locked analysis.
- KEGG GSEA uses the locked Entrez-based ranking representation.
- GO-BP GSEA uses the locked Ensembl-based ranking representation.
- Annotation and identifier mapping are governed by their corresponding upstream policy locks.
- Downstream audits should verify the stored inputs rather than silently reconstructing them.

### 8.2 Expected values versus verified results

Some workflow sections contain predefined expected values used for audit and integrity checks.

Examples include:

- KEGG final ranking size: 34,120.
- KEGG pathways returned: 353.
- KEGG pathways with FDR < 0.05: 178.
- GO-BP input ranking size: 18,008.
- GO-BP pathways returned: 6,316.
- GO-BP pathways with FDR < 0.05: 2,003.
- Genome-wide ranking size: 57,828.

These values are documented as workflow-specific audit references. Their presence in a script or README does not independently establish that the corresponding execution passed.

The final verified values and audit status should be obtained from the actual saved results and execution records.

---

## 9. Execution Order and Session Dependencies

The two scripts contain multiple sequential workflow sections. Some later sections rely on objects, path variables, and audit records created earlier in the same R session.

Examples include:

- Shared project and workflow directory variables.
- Locked input and ranking objects.
- GO–KEGG integration objects.
- GSEA result and provenance objects.
- STEP9D-5B evidence objects and metadata.
- STEP9D-5C final evidence-review objects.

Therefore:

1. Execute the workflow in the intended order.
2. Initialize project-relative directories before running dependent sections.
3. Preserve required intermediate RDS objects and lock files.
4. Do not assume that an isolated downstream section can run without its upstream dependencies.
5. Confirm that the expected input files and objects exist before executing a dependent audit or repair operation.

### 9.1 Audit and repair operations

The module includes several types of operations that should be distinguished in reproducibility records:

- **Formal analysis:** performs the designated enrichment calculation.
- **Read-only audit:** checks existing inputs, results, metadata, and provenance without recalculating the analysis.
- **Diagnostic inspection:** examines software implementation, object structure, or annotation coverage.
- **Output repair:** reconstructs designated serialized result files or metadata from an authoritative existing object.
- **Metadata repair:** updates specifically designated metadata fields while preserving substantive evidence objects.
- **Final evidence lock:** consolidates and verifies previously generated evidence.

An audit should not be described as a new enrichment analysis. Conversely, an output or metadata repair that replaces files or modifies an RDS object should not be described as strictly read-only.

Where an official output is repaired, the original file should be backed up and the repaired output should be reloaded and verified.

---

## 10. Recommended Environment Recording

To support reproducibility, each formal execution should retain the actual software environment information.

The following commands provide a basic record:

```r
# Record the R version
R.version.string

# Record the current R session
sessionInfo()

# Record selected package versions
packages_to_check <- c(
  "here",
  "data.table",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "GOSemSim",
  "enrichplot",
  "ggplot2",
  "dplyr",
  "stringr",
  "readr",
  "tibble",
  "purrr",
  "biomaRt",
  "fgsea"
)

installed_packages <- rownames(installed.packages())

package_versions <- data.frame(
  Package = packages_to_check,
  Version = vapply(
    packages_to_check,
    function(pkg) {
      if (pkg %in% installed_packages) {
        as.character(utils::packageVersion(pkg))
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

The resulting environment record should be associated with the corresponding analysis outputs or provenance files where appropriate.

For formal reproducibility documentation, the following information should be retained:

- Actual R version.
- Bioconductor version, when applicable.
- Installed versions of relevant R packages.
- Annotation database versions.
- Input-file provenance.
- Ranking and identifier-mapping policy locks.
- Relevant script version or Git commit identifier.
- Execution date and audit status.
- Any package or external-service availability limitations.

A package being unavailable during an optional diagnostic should be recorded as such; it should not automatically be interpreted as failure of an unrelated formal analysis step.

---

## 11. Reproducibility and Interpretation Notes

The following principles apply across the Functional Enrichment module:

1. Preserve the locked upstream DEG and ranking definitions.
2. Preserve the established Ensembl and Entrez identifier policies.
3. Record annotation database and package versions from the actual execution environment.
4. Do not infer installed software versions from the target R version.
5. Do not introduce jitter or artificial tie-breaking into locked ranking values.
6. Distinguish formal enrichment calculations from result audits and diagnostic inspections.
7. Preserve original outputs before any authorized repair or metadata modification.
8. Verify repaired outputs through file reload and object-concordance checks.
9. Distinguish predefined audit expectations from successfully verified execution results.
10. Preserve the descriptive nature of the Reactome, KEGG, and GO-BP theme-comparison workflows.

The historical Reactome references and the subsequent genome-wide KEGG and GO-BP analyses provide descriptive evidence comparisons within the same TCGA-LUAD study context. They should not be described as independent external validation.

Keyword-based theme matching and pathway-level correspondence do not establish identical biological mechanisms, direct molecular interactions, causal relationships, or clinical diagnostic validity.

---

## 12. Environment Documentation Status

This README consolidates the environment and software information documented in:

- `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`
- `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

It is intended to serve as the central environment reference for the `06_Functional_Enrichment` GitHub module.

**Environment documentation status:**

- R target version documented: R 4.5.1.
- Package roles and workflow dependencies consolidated.
- Annotation resources and identifier systems documented.
- Random-seed and ranking-tie policies summarized.
- Project-relative path conventions documented.
- Runtime package versions: to be confirmed from actual execution records.
- Actual operating-system and Bioconductor versions: to be confirmed from runtime records.
- Successful execution of all workflows: to be established from the corresponding analysis outputs and audit records.

This README documents the intended computational environment and reproducibility requirements. It does not independently certify that all analyses have been executed successfully in an identical software environment.