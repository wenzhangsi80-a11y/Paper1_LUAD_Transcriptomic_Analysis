# Functional Enrichment Analysis — Script Documentation

## 1. Overview

This directory contains the R scripts used for the functional enrichment analysis in Paper 1, entitled:

**Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R across the TME in LUAD**

The functional enrichment workflow is implemented in two main R scripts:

- `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`
- `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

The scripts are organized sequentially to support input provenance assessment, analytical definition recovery, over-representation analysis (ORA), genome-wide preranked gene set enrichment analysis (GSEA), result validation, and descriptive integration of pathway-level evidence.

The workflow distinguishes two complementary analytical components:

1. **STEP9C — DEG-based over-representation analysis (ORA):** evaluates GO and KEGG enrichment separately for tumor-upregulated and tumor-downregulated differentially expressed genes.
2. **STEP9D — Genome-wide preranked GSEA:** evaluates pathway-level enrichment using a genome-wide ranking statistic derived from differential expression results.

The workflow also includes semantic redundancy assessment, GO–KEGG evidence integration, gene-level theme comparison, and reproducibility audits.

The scripts use project-relative paths through the `here` package. They are intended to be executed from the root directory of the GitHub project.

---

## 2. Script Organization and Workflow Summary

| Analysis stage | Script | Main purpose |
|---|---|---|
| Input inventory and analytical definition recovery | Script 1 | Inventory upstream files, recover the DEG definition, and establish the analytical input framework |
| DEG-based GO and KEGG ORA | Script 1 | Perform enrichment analyses for tumor-upregulated and tumor-downregulated DEGs |
| GO semantic redundancy assessment | Script 1 | Identify semantically related GO-BP terms without changing the original enrichment results |
| GO–KEGG integration and descriptive evidence synthesis | Script 1 | Integrate pathway-level evidence and generate gene-associated theme summaries |
| Genome-wide GSEA input recovery and locking | Scripts 1–2 | Audit DEG source files, reconcile gene universes, and establish the final GSEA input |
| Ranking statistic construction | Script 2 | Construct and validate the genome-wide preranked gene vector |
| Identifier mapping and collision-policy assessment | Script 2 | Establish the Ensembl-to-Entrez mapping strategy for KEGG GSEA |
| Formal KEGG preranked GSEA | Script 2 | Perform genome-wide KEGG pathway enrichment |
| Formal GO-BP preranked GSEA | Script 2 | Perform genome-wide GO Biological Process enrichment |
| GSEA result audits and output repair | Script 2 | Validate result files, provenance, and reproducibility; repair designated output files where required |
| Reactome reference comparison | Script 2 | Compare historical gene-specific Reactome themes with the locked KEGG and GO-BP results |
| Final descriptive evidence and concordance locking | Script 2 | Produce gene-level theme evidence summaries and final audit records |

---

## 3. Script 1: `GitHub6_Paper1_Functional_Enrichment_Analysis_1.R`

### 3.1 Purpose

This script establishes the upstream analytical definitions and performs the DEG-based functional enrichment workflow.

Its main objectives are to:

- Inventory and assess the provenance of upstream WGCNA, DEG, and candidate-gene files.
- Recover and document the final DEG testing definition.
- Perform GO and KEGG over-representation analyses.
- Assess redundancy among enriched GO-BP terms.
- Integrate GO and KEGG evidence descriptively.
- Prepare and audit the source information required for subsequent genome-wide preranked GSEA.

The script includes both analytical procedures and supporting input audits. Not every workflow step performs a new statistical analysis.

### 3.2 Main workflow components

| Workflow | Main operations |
|---|---|
| STEP9A — Input Inventory and Definition Lock | Inventory upstream expression matrices, WGCNA module assignments, module–TME results, candidate-gene results, and DEG-related files; document input availability, provenance, and unresolved definitions |
| STEP9B — DEG Definition Recovery | Recover the final DEG testing universe, statistical test, effect-size definition, multiple-testing correction, and DEG selection criteria |
| STEP9C — DEG-based Functional Enrichment | Perform GO and KEGG ORA separately for tumor-upregulated and tumor-downregulated DEGs |
| STEP9C-5 — KEGG ORA recovery and locking | Reconcile the KEGG input vectors and background definition, assess the effective enrichment universe, and recover the designated results and analytical definitions |
| STEP9C-6 — GO-BP ORA and result recovery | Generate or recover GO-BP ORA results and establish the corresponding result definitions |
| STEP9C-7A — GO-BP Semantic Redundancy Assessment | Calculate semantic similarities among enriched GO-BP terms and identify descriptive clusters |
| STEP9C-7B — GO–KEGG Integration | Integrate same-direction GO-BP and KEGG evidence, examine cross-direction relationships, and assess direct study-gene overlap |
| STEP9C-7C — Theme-level Evidence Aggregation | Summarize related enrichment findings into descriptive biological themes |
| STEP9C-7D — Final Theme Evidence Lock | Generate the final descriptive theme-level evidence tables and save the designated locked integration object |
| STEP9D-1A to STEP9D-1I — GSEA Input Audit and Lock | Identify candidate input files, assess source and universe concordance, recover DEG definitions, and establish the final genome-wide GSEA input lock when the required criteria are satisfied |

### 3.3 Locked DEG definition

The DEG definition recovered in STEP9B is based on the following criteria:

- The testing universe consists of genes with finite expression values across the samples and non-zero combined expression variance.
- Differential expression is assessed using the Wilcoxon rank-sum test.
- The effect size is defined as the difference between the mean log-transformed expression values in tumor and normal samples:

  `log2FC = mean[Tumor: log2(FPKM-UQ + 1)] - mean[Normal: log2(FPKM-UQ + 1)]`

- Multiple-testing correction is performed using the Benjamini–Hochberg procedure.
- DEGs are defined using:
  - BH-adjusted FDR < 0.05
  - Absolute log2FC ≥ 1
- Direction is determined by the sign of log2FC.

The log2FC definition represents a difference between group means on the log-transformed expression scale. It is not a ratio of raw FPKM-UQ values.

### 3.4 Main analytical universes and reference counts

The following quantities refer to different analytical stages and must not be treated as interchangeable.

| Analytical component | Reference quantity | Definition |
|---|---:|---|
| DEG testing universe | 57,828 genes | Genes eligible for the recovered genome-wide DEG testing procedure |
| Tumor-upregulated DEGs | 9,311 genes | Genes meeting the locked tumor-upregulated DEG criteria |
| Tumor-downregulated DEGs | 6,754 genes | Genes meeting the locked tumor-downregulated DEG criteria |
| KEGG ORA background | 35,365 unique Entrez IDs | Locked Entrez-based background used for KEGG ORA |
| KEGG ORA tumor-up input | 5,194 unique Entrez IDs | Mapped tumor-upregulated DEG input |
| KEGG ORA tumor-down input | 4,463 unique Entrez IDs | Mapped tumor-downregulated DEG input |

The corresponding Ensembl-to-Entrez mapping audit recorded:

- 57,828 genes in the testing universe.
- 34,631 genes with at least one Entrez mapping.
- 23,197 genes without an Entrez mapping.
- 34,138 genes mapping to exactly one Entrez ID.
- 493 genes mapping to more than one Entrez ID.

These values describe the mapping audit and should be interpreted separately from the final KEGG ORA background.

### 3.5 GO-BP semantic redundancy and GO–KEGG integration

The semantic redundancy analysis is performed on the locked GO-BP enrichment results.

The workflow uses:

- Wang semantic similarity.
- Complete-linkage hierarchical clustering.
- Similarity thresholds of 0.50, 0.60, 0.70, and 0.80.
- A threshold of 0.70 as the primary descriptive clustering reference.

The purpose is to identify groups of semantically related GO-BP terms and reduce redundant descriptions when summarizing biological themes.

The complete enrichment results remain unchanged. Semantic clustering does not replace the original statistical results or establish a new significance criterion.

The subsequent GO–KEGG integration:

- Examines same-direction GO-BP and KEGG evidence.
- Audits cross-direction relationships.
- Records pathway and term cluster assignments.
- Examines direct overlap with the three study genes.
- Aggregates related findings into descriptive themes.

The integration workflow distinguishes the complete integration universe from the subset containing non-zero evidence. It does not create a new statistical ranking, composite score, or top-N selection.

The final theme evidence object is saved as:

`STEP9C7D_FINAL_LOCK.rds`

### 3.6 STEP9D-1: Genome-wide GSEA input audit

The final part of Script 1 prepares the source and definition framework for genome-wide preranked GSEA.

The workflow includes:

- Candidate input-file discovery.
- Structural and provenance assessment.
- Source-to-universe gene-set concordance checks.
- Recovery of the DEG definition.
- Comparison of existing and reconstructed DEG classifications.
- Examination of boundary disagreements and gene-level classification.
- Final input-lock assessment.

The STEP9D-1I procedure is an input-locking and validation step. It does not construct the final ranking vector or perform GSEA.

The formal lock is established only when the required input and concordance criteria are satisfied. The existence of a lock file alone should not be interpreted as evidence that all criteria passed.

---

## 4. Script 2: `GitHub6_Paper1_Functional_Enrichment_Analysis_2.R`

### 4.1 Purpose

Script 2 continues the genome-wide preranked GSEA workflow using the source definitions and audit information established in the preceding steps.

Its main objectives are to:

- Reassess and finalize the GSEA input.
- Construct the genome-wide ranking vector using the locked ranking formula.
- Establish identifier-mapping policies for KEGG and GO-BP analysis.
- Perform formal KEGG and GO-BP preranked GSEA.
- Audit the resulting enrichment outputs and provenance.
- Compare historical Reactome themes with the genome-wide GSEA results.
- Generate and lock the final descriptive gene-level concordance evidence.

This script is a continuation of the preceding workflow. Some steps depend on objects and path variables created in earlier sections of the same R session.

### 4.2 Main workflow components

| Workflow | Main operations |
|---|---|
| STEP9D-1I — Final GSEA Input Lock | Reassess the final source, DEG universe, and gene-level concordance requirements |
| STEP9D-2A — Ranking Implementation Audit | Inspect the ranking source, statistical fields, software availability, and implementation requirements |
| STEP9D-2B — Ranking Metric Lock | Establish and validate the ranking statistic definition |
| STEP9D-2C — Ranking Vector Construction | Construct the named genome-wide ranking vector and validate its integrity |
| STEP9D-3A — Identifier Mapping and Policy Preparation | Audit Ensembl-to-Entrez mapping, annotation coverage, and candidate collision-handling strategies |
| STEP9D-3A-6 — Formal KEGG Mapping Policy Lock | Lock the selected Entrez collision-handling policy |
| STEP9D-3B — Formal KEGG GSEA | Execute genome-wide preranked KEGG enrichment |
| STEP9D-3C — KEGG Result and Reproducibility Audit | Validate pathway results, ranking integrity, provenance, and exported files |
| STEP9D-4A — GO-BP Input and Annotation Policy | Establish the GO-BP identifier and annotation framework |
| STEP9D-4B — Formal GO-BP GSEA | Execute genome-wide preranked GO Biological Process enrichment |
| STEP9D-4B — Output and Provenance Repair | Repair designated exported result and provenance files using the authoritative result object |
| STEP9D-4C — GO-BP Result and Reproducibility Audit | Validate GO-BP results, exported files, ranking consistency, and provenance |
| STEP9D-5 — Initial Reactome Concordance Audit | Compare historical Reactome themes with the genome-wide KEGG and GO-BP results |
| STEP9D-5B — Strict Theme Evidence Audit | Apply predefined direct and contextual theme dictionaries to pathway descriptions |
| STEP9D-5C — Final Evidence and Concordance Lock | Review, reconcile, and lock the final descriptive gene-level theme evidence |

### 4.3 Genome-wide ranking statistic

The ranking statistic is defined as:

`sign(log2FC) × [-log10(max(P_value, 1e-300))]`

This statistic combines the direction of differential expression with the magnitude of nominal statistical evidence.

The ranking construction workflow:

1. Reads the source data from the established GSEA input lock.
2. Checks the availability and validity of the required statistical fields.
3. Applies the locked ranking formula.
4. Applies the predefined lower bound to P values to avoid infinite ranking statistics.
5. Examines zero values and direction consistency.
6. Constructs a named ranking vector.
7. Saves and reloads the ranking object for integrity assessment.

The ranking metric is not recalculated independently in downstream GSEA steps. No artificial jitter or arbitrary tie-breaking is introduced to alter the ranking.

### 4.4 Identifier mapping and KEGG collision policy

KEGG GSEA requires Entrez Gene identifiers, whereas the original genome-wide ranking is based on Ensembl gene identifiers.

The mapping workflow assesses:

- Ensembl-to-Entrez annotation coverage.
- Unmapped and multiply mapped genes.
- Duplicate Entrez identifiers.
- Directional consistency among genes sharing an Entrez identifier.
- Retention of the three study genes.

Three candidate collision-handling strategies are evaluated:

- **Policy A:** Strict one-to-one mapping.
- **Policy B:** Retain the Ensembl gene with the largest absolute value of the already locked ranking statistic for each Entrez identifier.
- **Policy C:** Select a deterministic representative based on ranking order.

The formally locked KEGG strategy is Policy B.

This policy selects the gene with the largest absolute value of the existing ranking statistic within each Entrez collision group. It does not recalculate the ranking statistic, introduce artificial tie-breaking, or use arbitrary first-match selection.

### 4.5 Formal KEGG preranked GSEA

The formal KEGG GSEA is implemented using `clusterProfiler` and the associated annotation packages.

The main analytical settings include:

- Organism: *Homo sapiens*.
- Identifier type: Entrez Gene ID.
- Gene-set size: 10–500.
- Enrichment method: preranked GSEA.
- Exponent: 1.
- Multiple-testing correction: Benjamini–Hochberg.

The formal execution uses the locked genome-wide ranking vector and the established Ensembl-to-Entrez mapping policy.

The resulting KEGG analysis is distinct from the DEG-based KEGG ORA performed in STEP9C.

### 4.6 Formal GO-BP preranked GSEA

GO-BP GSEA is performed using Ensembl identifiers directly, rather than applying the KEGG Entrez mapping policy.

The main settings include:

- Ontology: Biological Process (BP).
- Identifier type: ENSEMBL.
- Annotation database: `org.Hs.eg.db`.
- Gene-set size: 10–500.
- Exponent: 1.
- P-value cutoff: 1 for result retention.
- Multiple-testing correction: Benjamini–Hochberg.

The formal analysis is implemented using `clusterProfiler::gseGO()`.

The GO-BP analysis uses the locked genome-wide ranking source and its designated GO-BP annotation policy. It is not a repetition of the GO-BP ORA performed in STEP9C.

### 4.7 GSEA result validation and output repair

The KEGG and GO-BP result audits assess:

- Result-object structure and required columns.
- Pathway or GO-term counts.
- Gene identifiers and ranking-vector consistency.
- NES, nominal P values, and adjusted P values.
- The number of results meeting the FDR criterion.
- Concordance between RDS and CSV outputs.
- Core-enrichment and leading-edge information.
- Ranking order and tie-related characteristics.
- Provenance records and analytical parameters.
- Reload integrity of the saved audit outputs.

The GO-BP output repair step uses the authoritative `gseaResult` RDS object to reconstruct designated CSV and provenance files.

Existing contaminated outputs are backed up before replacement. The authoritative GSEA result object is not modified, and `gseGO()` is not rerun during this repair procedure.

The repair operation should therefore be understood as an output-file and provenance correction, not as a new enrichment analysis.

### 4.8 Reactome reference and theme concordance

The final stages compare historical gene-specific Reactome GSEA observations with the locked genome-wide KEGG and GO-BP results.

The historical Reactome reference themes are:

| Study gene | Historical Reactome reference |
|---|---|
| IL7R | `REACTOME_FCGR_ACTIVATION` |
| MS4A7 | `REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS` |
| AOC3 | `REACTOME_CD22_MEDIATED_BCR_REGULATION` |

The comparison is descriptive and uses predefined theme-related terms to identify potentially related evidence in the KEGG and GO-BP results.

The workflow distinguishes two levels of assessment:

**Initial STEP9D-5 audit**

- Examines pathway identifiers, descriptions, and supporting gene-level fields.
- Summarizes matching pathways, significant hits, maximum absolute NES, and best adjusted P values.
- Uses preliminary descriptive evidence categories.

**STEP9D-5B strict theme evidence audit**

- Uses predefined direct and contextual keyword dictionaries.
- Bases primary theme matching on the pathway or GO-term `Description` field.
- Retains `leading_edge` and `core_enrichment` as supporting evidence rather than using them as the primary keyword-matching fields.
- Distinguishes five evidence groups:
  - `Direct_significant`
  - `Contextual_significant`
  - `Direct_non_significant`
  - `Contextual_non_significant`
  - `No_theme_match`

**STEP9D-5C final evidence and concordance lock**

- Reviews the authoritative STEP9D-5B evidence objects and gene-level summaries.
- Checks expected dimensions, arithmetic consistency, file concordance, and Reactome reference inheritance.
- Generates the final descriptive evidence tables and audit records.
- Preserves the distinction between direct and contextual evidence.
- Does not assign biological strength rankings or composite scores.

The Reactome reference and genome-wide GSEA results share the same general TCGA-LUAD study context. Their comparison is therefore a descriptive cross-analysis comparison, not independent external validation.

---

## 5. Key Analytical Definitions and Reference Values

The following reference values are documented to facilitate reproducibility and output inspection. They should be checked against the actual executed result objects and audit reports before being reported as verified final results.

| Analysis component | Reference value |
|---|---:|
| DEG testing universe | 57,828 genes |
| Tumor-upregulated DEGs | 9,311 genes |
| Tumor-downregulated DEGs | 6,754 genes |
| KEGG ORA background | 35,365 unique Entrez IDs |
| KEGG ORA tumor-up input | 5,194 unique Entrez IDs |
| KEGG ORA tumor-down input | 4,463 unique Entrez IDs |
| Genome-wide GSEA ranking source | 57,828 genes |
| Non-zero genome-wide ranking values | 57,815 |
| KEGG GSEA ranking vector | 34,120 unique Entrez IDs |
| GO-BP annotated gene universe | 18,008 genes |
| KEGG pathways | 353 |
| KEGG pathways with adjusted P < 0.05 | 178 |
| GO-BP terms | 6,316 |
| GO-BP terms with adjusted P < 0.05 | 2,003 |

The KEGG ORA background of 35,365 unique Entrez identifiers and the KEGG GSEA ranking vector of 34,120 unique Entrez identifiers arise from different procedures and must remain separate in reporting.

Likewise, the GO-BP annotated gene universe is not equivalent to the complete genome-wide DEG testing universe.

### 5.1 Reference gene-level theme summary

The STEP9D-5C documentation includes the following predefined gene-by-database summary values:

| Study gene | Database | Total terms/pathways | Theme matches | Significant theme matches |
|---|---|---:|---:|---:|
| IL7R | KEGG | 353 | 0 | 0 |
| IL7R | GO-BP | 6,316 | 41 | 15 |
| MS4A7 | KEGG | 353 | 7 | 3 |
| MS4A7 | GO-BP | 6,316 | 251 | 83 |
| AOC3 | KEGG | 353 | 1 | 1 |
| AOC3 | GO-BP | 6,316 | 58 | 9 |

These values are reference expectations recorded in the workflow documentation. They should not be interpreted as independently verified results unless confirmed by the corresponding final RDS and CSV outputs.

---

## 6. Input Dependencies

The functional enrichment scripts depend on upstream outputs generated during the Paper 1 transcriptomic analysis.

The principal input categories include:

1. Tumor-only WGCNA expression data and gene-module assignments.
2. Locked module–tumor microenvironment association results.
3. Candidate-gene and core-gene audit results.
4. Differential expression source files and recovered DEG definitions.
5. Gene annotation resources, including GENCODE v22 and `org.Hs.eg.db`.
6. Locked STEP9C ORA results and integration objects.
7. STEP9D input-lock and ranking-vector objects.
8. Existing Reactome reference results used for descriptive theme comparison.

The exact file names and locations are defined in the corresponding workflow sections and output inventories.

Script 2 is a sequential continuation of Script 1. Some sections rely on objects, directory variables, and previously loaded checkpoints retained in the R session. Therefore, execution should follow the documented workflow order rather than assuming that every section can be run independently in a fresh R session.

---

## 7. Output Organization

The scripts organize outputs by analytical stage under the functional enrichment results directory.

The principal output categories include:

| Output category | Main contents |
|---|---|
| Input inventory and locked definitions | Input inventories, DEG definition records, source and universe concordance audits |
| DEG-based ORA | GO and KEGG enrichment result objects and exported result tables |
| Semantic redundancy assessment | GO-BP similarity matrices, cluster assignments, and representative-term summaries |
| GO–KEGG integration | Integrated evidence tables, theme-level summaries, and locked integration objects |
| GSEA input and ranking | Input-lock records, ranking metric definitions, ranking vectors, and gene-level audits |
| Identifier mapping | Annotation coverage, Entrez collision audits, and formal mapping-policy records |
| Formal KEGG GSEA | KEGG GSEA result objects, CSV exports, and provenance records |
| Formal GO-BP GSEA | GO-BP GSEA result objects, CSV exports, and provenance records |
| Reproducibility audits | Structural, statistical, file-concordance, and provenance audit reports |
| Reactome concordance | Gene-specific theme evidence tables and comparison summaries |
| Final evidence lock | Final gene-level concordance tables, evidence summaries, and audit records |

Important intermediate and final records include:

- `STEP9C7D_FINAL_LOCK.rds`
- `STEP9D1I_FINAL_LOCK.rds`
- The STEP9D-2C ranking-vector RDS and CSV outputs.
- The STEP9D-3A-6 KEGG mapping-policy lock.
- The formal KEGG and GO-BP GSEA result RDS objects.
- The STEP9D-5B theme evidence objects.
- The STEP9D-5C final evidence and concordance lock files.

The precise output filenames and subdirectory locations should be taken from the executed scripts and their corresponding result inventories.

---

## 8. Relationship to the Manuscript

The functional enrichment workflow supports the transcriptomic functional characterization described in the Paper 1 manuscript.

Its main contributions are:

- DEG-based GO and KEGG enrichment for tumor-upregulated and tumor-downregulated genes.
- Genome-wide preranked KEGG and GO-BP enrichment.
- Descriptive assessment of biological themes associated with AOC3, MS4A7, and IL7R.
- Comparison of historical gene-specific Reactome observations with the genome-wide enrichment results.

The results should be interpreted as complementary transcriptomic evidence. The enrichment analyses do not establish causal signaling relationships, cell-specific mechanisms, or independent clinical biomarker performance.

The Reactome comparison is a descriptive correspondence analysis rather than an independent validation cohort.

---

## 9. Interpretation and Reporting Considerations

The following distinctions should be maintained when interpreting and reporting the results:

1. **ORA and preranked GSEA are different analytical approaches.** ORA evaluates predefined DEG subsets against a background universe, whereas preranked GSEA evaluates pathway-level enrichment across a ranked genome-wide gene list.

2. **Analytical universes are procedure-specific.** The DEG testing universe, KEGG ORA background, KEGG GSEA ranking vector, and GO-BP annotated universe must not be interchanged.

3. **Semantic clustering is descriptive.** GO-BP semantic similarity is used to identify redundant or related terms; it does not replace the original enrichment statistics.

4. **GO–KEGG integration is evidence synthesis.** It does not generate a new statistical ranking, composite score, or causal network.

5. **Study-gene theme matching is descriptive.** Direct and contextual keyword matches should not be interpreted as proof that a gene regulates a pathway.

6. **Output repair is not reanalysis.** Designated file and metadata repairs preserve the authoritative analytical result objects and should be documented separately from formal enrichment execution.

7. **Reference values are not automatically verified results.** Final numerical claims should be checked against the saved result objects and execution logs.

8. **Reactome concordance is not independent validation.** The historical Reactome observations and the genome-wide GSEA analyses should be presented as complementary analyses within the study context.

9. **Bulk transcriptomic enrichment does not establish cell-intrinsic activity.** Enrichment patterns may reflect differences in cellular composition as well as transcriptional changes within individual cell populations.

---

## 10. Reproducibility and Execution Notes

- Open the GitHub project in RStudio and execute the scripts from the project root.
- Use `here::here()` and project-relative paths to resolve input and output locations.
- Avoid manual working-directory changes with `setwd()`.
- Confirm that all upstream input files, locked objects, and required checkpoints are available before execution.
- Install and load the required R packages, including `here`, `clusterProfiler`, `fgsea`, `org.Hs.eg.db`, `AnnotationDbi`, and the supporting data-manipulation and visualization packages used in the scripts.
- Follow the sequential workflow dependencies, especially the STEP9D input-lock, ranking-vector, and identifier-mapping stages.
- Preserve the authoritative RDS result objects and designated backup files when performing output repairs.
- Inspect the generated audit reports and execution logs before treating any expected value or PASS status as verified.
- Record the R version, package versions, and session information in the corresponding reproducibility outputs.
- Use the final saved result objects and audit reports as the source of truth for manuscript-level numerical reporting.

This README summarizes the responsibilities and dependencies of the two functional enrichment scripts. Detailed intermediate variables and step-specific audit records are intentionally retained in the corresponding output files rather than duplicated in this directory-level script guide.