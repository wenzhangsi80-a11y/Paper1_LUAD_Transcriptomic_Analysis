# Analysis Scripts

## 1. Overview

This directory contains the R script used to perform patient-level ligand–receptor (LR) expression-compatibility analysis in the lung adenocarcinoma (LUAD) study.

**Main script:**

`GitHub8_Paper1_ligand_receptor_compatibility_analysis.R`

**Analysis module:** Ligand–Receptor Compatibility Analysis

**Study:**

*Integrated transcriptomic and single-cell characterization of AOC3, MS4A7, and IL7R in LUAD.*

The workflow evaluates the expression compatibility of predefined ligand–receptor records across three cellular compartments: endothelial cells, myeloid/macrophage cells, and T cells. It integrates previously generated single-cell expression and pseudobulk outputs with a predefined ligand–receptor reference resource.

The analysis is organized into sequential input audits, expression-compatibility screening, candidate-set characterization, patient-level integration, descriptive summarization, figure generation, and manuscript-table traceability checks.

The script uses project-relative paths and is intended to be executed within the GitHub project context.

## 2. Project Setup and Input Data

### 2.1 Project-relative directory initialization

The script initializes the project and module directories using `here::here()`.

The main module directories are:

- `data/input_data/`
- `data/sample_manifests/`
- `environment/`
- `results/`
- `scripts/`

Machine-specific absolute paths and manual working-directory changes are not required by the revised workflow.

### 2.2 Primary input files

The main external input files include:

1.  `NSCLC_GSE131907_CCI_LR.tsv`
2.  `NSCLC_GSE131907_CellMetainfo_table.tsv`

These files are expected under:

`08_Ligand_Receptor/data/input_data/`

The CCI_LR file provides the predefined ligand–receptor reference resource. Its records define the candidate interaction universe used by downstream compatibility analyses; LR pairs are not inferred de novo from the analyzed expression data.

### 2.3 Previously generated STEP10 inputs

The workflow also reuses previously generated outputs from:

`Paper 1_Step 10-Single-Cell/`

Important inputs include:

- `STEP10E-1_Primary_Patient_Level_CoreGene_Expression.tsv`
- `STEP10F-1C_Primary_Gene_Compartment_Expression_Profile.tsv`
- `STEP10F-2A_Primary_Gene_Compartment_Pattern_Profile.tsv`
- `STEP10F-2D_Complete_Patient_Within_Patient_Ranks.tsv`
- `STEP10F-2E_Complete_ThreeCompartment_Pattern_Summary.tsv`
- `STEP10F-2F_Complete_Patient_Dominant_Compartment.tsv`
- `STEP10F-1E_Complete_ThreeCompartment_Patient_Profile.tsv`
- `STEP10C-E_Primary_Pseudobulk_Sample_Annotation.tsv`

Three existing compartment-specific raw-UMI pseudobulk matrices are particularly important:

- `STEP10C-A_Pseudobulk_Tcell_Primary_rawUMI.tsv`
- `STEP10C-B_Pseudobulk_MyeloidMacrophage_Primary_rawUMI.tsv`
- `STEP10C-C_Pseudobulk_Endothelial_Primary_rawUMI.tsv`

The three pseudobulk files are located by their exact filenames. Their specific upstream subdirectory locations are not assumed by the file-discovery steps.

### 2.4 Downstream input dependencies

Subsequent workflows read previously generated outputs from the module's `results/` directory. These include:

- STEP11-3 patient and sample eligibility audit outputs.
- STEP11-5D-1 genome-wide LR compatibility candidate tables.
- STEP11-5D-2 candidate structural audit tables.
- STEP11-5E-2 candidate expression summaries.
- STEP11-5E-3 canonical patient-level integration tables.
- STEP11-5F-4 final master descriptive context table.
- STEP11-5F-5 Figure 7 plotting outputs.

The downstream workflows are designed to reuse these locked or previously generated results rather than repeatedly processing the original expression matrices.

## 3. Cellular Compartments and Direction Definitions

Three predefined cellular compartments are considered:

| Compartment        | Included cell populations  | Abbreviation |
|--------------------|----------------------------|--------------|
| Endothelial        | Endothelial cells          | E            |
| Myeloid/macrophage | Mono/Macro                 | M            |
| T-cell             | CD4Tconv, CD8T, and CD8Tex | T            |

Six sender–receiver directions are evaluated:

1.  Endothelial → Myeloid/macrophage (E→M)
2.  Myeloid/macrophage → Endothelial (M→E)
3.  Endothelial → T-cell (E→T)
4.  T-cell → Endothelial (T→E)
5.  Myeloid/macrophage → T-cell (M→T)
6.  T-cell → Myeloid/macrophage (T→M)

The matched patient-level compatibility analyses use the four patients represented in all three predefined compartments.

## 4. Workflow Organization

### 4.1 Initial input and metadata audits

| Workflow | Main purpose |
|----|----|
| `00_setup` | Initialize the project root and module directories. |
| `STEP11-0` | Audit the structure and semantic fields of the CCI_LR reference resource. |
| `STEP11-1` | Audit the structure of the cell metadata table. |
| `STEP11-1A` | Perform detailed metadata checks, including patient, sample, cell-type, source, and stage fields. |
| `STEP11-2` | Audit Primary-stage patient-level target-cell count distributions and candidate thresholds. |
| `STEP11-3` | Audit Primary-stage sample-level compartment coexistence and patient-level eligibility. |
| `STEP11-4` | Integrate STEP10 expression and localization outputs with STEP11-3 patient-level eligibility information. |

#### STEP11-2 and STEP11-3

These workflows characterize the distribution of target-cell populations across patients and samples.

The main operations include:

- Patient-by-cell-type and sample-by-cell-type count construction.
- Patient-level coverage and cell-count distribution summaries.
- Candidate cell-count threshold assessment.
- Combined T-cell population assessment.
- Same-sample coexistence of the three predefined compartments.
- Tumor-source-specific and multiple-sample patient audits.
- Patient-level eligibility and structural QC.

These steps do not perform CellChat inference, LR interaction inference, gene-expression analysis, or final patient eligibility selection.

#### STEP11-4: Core-Gene Cellular Localization Integration

This workflow integrates previously generated STEP10 localization and pseudobulk outputs with the STEP11-3 eligibility results.

The main components include:

1.  STEP10 input-file inventory and structural audits.
2.  STEP11-3 output-file inventory and structural audits.
3.  Required-column and core-gene checks.
4.  Compartment-level expression and localization integration.
5.  Patient-level localization and within-patient rank summaries.
6.  Dominant-compartment and expression-detection summaries.
7.  Patient eligibility feasibility assessment.
8.  Multiple-sample patient summaries.
9.  Structural concordance checks.
10. TSV and RDS output generation.

This step reuses existing upstream expression results. It does not perform new raw-UMI processing, expression recalculation, cell-level reanalysis, or patient-level CellChat inference.

### 4.2 CCI_LR resource and pseudobulk compatibility audits

| Workflow | Main purpose |
|----|----|
| `STEP11-5A` | Audit compatibility between the CCI_LR reference resource and STEP10 gene-expression profiles. |
| `STEP11-5B-0` | Discover and structurally audit candidate STEP10 genome-wide expression files. |
| `STEP11-5B-3A` | Locate the three existing compartment-specific STEP10C pseudobulk files. |
| `STEP11-5B-3B` | Audit the structure of the existing pseudobulk matrices. |
| `STEP11-5B-3C` | Decode patient and compartment information and inspect the matched-patient design. |
| `STEP11-5C-1` | Audit raw-UMI matrix integrity and expression-scale compatibility. |
| `STEP11-5C-2` | Assess normalized-expression stability and matched-patient scale consistency. |

#### STEP11-5A: CCI_LR × STEP10 Gene-Expression Compatibility Audit

This step evaluates the structural compatibility between the predefined CCI_LR resource and the previously generated STEP10 expression profiles.

The main operations include:

- Input-file existence and structural checks.
- Required-column verification.
- CCI_LR resource integrity assessment.
- STEP10 expression-table integrity assessment.
- Global ligand/receptor gene-universe compatibility.
- Focused-compartment gene compatibility.
- Construction of the six predefined sender–receiver combinations.
- Direction-specific candidate interaction summaries.
- Pathway and annotation summaries.
- Core-gene compartment-presence and CCI_LR resource-membership checks.
- Interpretation-status documentation.

The resulting candidates represent structural expression compatibility for downstream screening. They are not validated cell–cell communication events.

#### STEP11-5B: Pseudobulk File Discovery and Structural Audits

These steps locate and inspect the three existing STEP10C pseudobulk matrices.

The main checks include:

- File availability and matrix dimensions.
- Gene identifiers, ordering, and uniqueness.
- Patient identifiers and compartment labels.
- Numeric raw-UMI matrix structure.
- Matched-patient representation across compartments.
- Gene overlap with the CCI_LR ligand/receptor resource.

The matrices are reused as existing inputs. The workflow does not re-aggregate raw UMI counts or recalculate pseudobulk profiles.

#### STEP11-5C-1: Expression-Scale and Normalization Compatibility Audit

This step checks the integrity of the three existing raw-UMI pseudobulk matrices.

The main operations include:

- Gene identity, ordering, and uniqueness checks.
- Identification of the common patient set.
- Numeric, finite, non-missing, and non-negative raw-UMI checks.
- Library-size and detected-gene summaries.
- Descriptive library-depth comparisons across compartments.
- Calculation of CPM and log2(CPM + 1) for audit purposes.
- Transformation-integrity checks.
- Descriptive inspection of AOC3, MS4A7, and IL7R expression and detection.

#### STEP11-5C-2: Normalized Expression Stability and Matched-Patient Scale Audit

This step examines normalized-expression distributions and scale stability across the matched patient set.

The workflow includes:

- Reuse of the three existing STEP10C pseudobulk matrices.
- Verification of gene identity and matched-patient representation.
- Library-size and normalized-expression summaries.
- LR-gene overlap and detection summaries.
- Core-gene expression inspection.
- Rank preservation between raw UMI and normalized expression.
- Patient-wise scale-stability assessment.

The step requires exactly four matched patients. Its intermediate objects and audit summaries are generated in memory; no standalone output files are exported.

### 4.3 Genome-wide LR compatibility screening and candidate-set locking

| Workflow | Main purpose |
|----|----|
| `STEP11-5D-0` | Evaluate LR detection-rule and patient-consistency sensitivity. |
| `STEP11-5D-1` | Perform genome-wide patient-level LR compatibility screening. |
| `STEP11-5D-2` | Audit the structural properties of genome-wide LR candidate records. |
| `STEP11-5D-3` | Lock the candidate-set definitions and inherited direction mapping. |

#### STEP11-5D-0: LR Detection and Patient-Consistency Sensitivity Audit

This step evaluates three predefined detection rules:

- UMI \> 0
- CPM \> 0
- CPM ≥ 1

The workflow:

1.  Verifies the existing pseudobulk matrices and gene identities.
2.  Identifies the four patients shared across all three compartments.
3.  Constructs matched-patient expression matrices.
4.  Derives CPM and log2(CPM + 1).
5.  Identifies LR genes represented in the expression matrices.
6.  Evaluates ligand and receptor detection across six directions.
7.  Assesses same-patient ligand–receptor compatibility.
8.  Summarizes patient consistency at the 2/4, 3/4, and 4/4 levels.
9.  Compares results across the three detection rules.
10. Exports four sensitivity-audit tables.

This step evaluates alternative detection conditions. It does not establish a preferred biological threshold or lock the final LR candidate set.

#### STEP11-5D-1: Genome-wide Patient-Level LR Compatibility Screening

This workflow screens the predefined CCI_LR resource using the existing matched-patient pseudobulk matrices.

The main operations include:

- Gene identity and numeric-integrity checks.
- Matched-patient identification and matrix reordering.
- CPM and log2(CPM + 1) calculation.
- CCI_LR ligand/receptor gene-universe assessment.
- Core-gene separation checks.
- Evaluation of ligand detection in the sender compartment and receptor detection in the receiver compartment.
- Compatibility summaries at the ≥2/4, ≥3/4, and 4/4 patient-consistency levels.
- Direction-level, global, pathway, and annotation summaries.
- Output record-count and logical consistency checks.

The primary compatibility assessment uses the predefined CPM \> 0 detection rule and requires ligand and receptor detection in the same patient.

No LR statistical testing, P-value calculation, FDR correction, pathway enrichment, or biological ranking is performed.

#### STEP11-5D-2: Genome-wide LR Candidate Structural Audit

This step uses previously generated STEP11-5D-1 candidate tables.

The main operations include:

- Candidate-table availability and schema checks.
- Construction and verification of exact LR-pair and directional identifiers.
- Repeated interaction-record characterization.
- Ligand-level receptor expansion summaries.
- Directional representation and symmetry audits.
- Annotation and pathway structure summaries.
- Core-gene separation confirmation.
- Structural comparison of the ≥3/4 and 4/4 candidate sets.
- Export of structural audit tables and QC results.

The step does not recalculate expression, pseudobulk matrices, LR detection, or CellChat results. Annotation and pathway summaries are descriptive structural summaries, not enrichment analyses.

#### STEP11-5D-3: Genome-wide LR Candidate Set Definition Lock

This step formally locks the previously audited candidate universe and inherited Direction_ID mapping.

The workflow:

- Verifies the required candidate and structural-audit input tables.
- Checks inherited exact LR-pair identifiers.
- Evaluates predefined candidate counts and the nested relationship between candidate sets.
- Verifies the six existing Direction_ID categories.
- Compares sender/receiver mappings across candidate sets.
- Confirms the inherited separation of AOC3, MS4A7, and IL7R.
- Exports the locked RDS object and summary/QC tables.

The expected record and exact-pair counts are predefined QC conditions. They are not newly inferred biological results.

### 4.4 Candidate overlap, expression characterization, and canonical integration

| Workflow | Main purpose |
|----|----|
| `STEP11-5E-1` | Characterize candidate overlap and patient-consistency patterns. |
| `STEP11-5E-2` | Summarize inherited patient-level ligand and receptor expression magnitudes. |
| `STEP11-5E-3` | Audit duplicate structure, canonical keys, and candidate provenance. |
| `STEP11-5E-3B` | Generate the formal canonical patient-level LR integration table. |

#### STEP11-5E-1: Candidate Overlap and Patient-Consistency Characterization

This step compares the previously generated ≥3/4 and 4/4 candidate sets.

The main operations include:

- Input-table availability and required-column checks.
- Exact LR-pair identifier verification.
- Identification of shared, ≥3/4-only, and 4/4-only pairs.
- Candidate-set nesting and overlap arithmetic.
- Calculation of the 4/4 retention proportion within the ≥3/4 set.
- Inspection of inherited patient-level compatibility patterns.
- Characterization of non-compatible patient records.
- Verification of the six sender–receiver direction mappings.
- Generation of pair-level, patient-pattern, directional, dropout, and QC summaries.

Patient-level compatibility values are inherited from STEP11-5D-1 and are not recalculated.

#### STEP11-5E-2: Candidate Expression-Magnitude Characterization

This step summarizes existing expression values in the previously generated candidate tables.

The workflow includes:

- Candidate identity and schema checks.
- Exact LR-pair identity verification.
- Candidate record counts and nesting checks.
- Numeric, finite, and non-negative checks of existing log2(CPM + 1) fields.
- Inspection of inherited patient-compatibility values.
- Descriptive candidate consistency classification.
- Construction of patient-level ligand and receptor expression records.
- Candidate-level and direction-specific expression summaries.
- Sender–receiver mapping and core-gene separation checks.
- Export of expression summaries and QC results.

Expression values are inherited from STEP11-5D-1. Raw UMI counts, pseudobulk expression, CPM, log2(CPM + 1), LR detection, and patient compatibility are not recalculated.

#### STEP11-5E-3: Patient-Level Candidate Integration and Provenance Audits

This workflow contains consecutive duplicate-structure, canonical-key, and provenance audits.

The main components include:

1.  **Duplicate-structure pre-audit**
    - Identifies repeated `Exact_LR_Pair × Direction_ID × Patient_Index` keys.
    - Characterizes duplicate multiplicity.
    - Inspects variation in interaction names, pathways, annotations, expression values, and compatibility.
2.  **Pre-collapse and canonical-key audit**
    - Checks the previously generated patient-level expression table.
    - Verifies exact LR-pair identity.
    - Defines the canonical key as `Exact_LR_Pair × Direction_ID × Patient_Index`.
    - Evaluates invariance of protected variables within each canonical key.
    - Retains interaction-name multiplicity as provenance.
    - Assesses key uniqueness, row-count conservation, and conflicts.
3.  **Candidate-set provenance audit**
    - Evaluates the expected `ge3` and `ge3 | ge4` provenance patterns.
    - Compares the exact LR-pair universes and candidate-set nesting.
    - Distinguishes provenance fields from canonical biological identifiers.
    - Audits interaction-name multiplicity.
4.  **Formal canonical integration**
    - Generates the final canonical patient-level representation.
    - Retains `Candidate_Set` and `interaction_name` provenance explicitly.
    - Exports the canonical integration table and structural QC summaries.

The canonical key is used to define the patient-level LR record identity. Candidate-set and interaction-name fields are retained as provenance rather than treated as biological invariants.

The audit does not delete or collapse records in the original STEP11-5E-2 input table. The formal canonical output is a structured integration of inherited records, not a new expression or compatibility analysis.

### 4.5 Canonical QC and patient-level descriptive context

| Workflow | Main purpose |
|----|----|
| `STEP11-5F-1` | Audit the final canonical input table. |
| `STEP11-5F-2` | Summarize patient-level compatibility consistency. |
| `STEP11-5F-3` | Integrate patient-level expression context. |
| `STEP11-5F-4` | Generate the final LR-pair descriptive master table and directional summaries. |

#### STEP11-5F-1: Canonical Input QC

This step validates the formal canonical patient-level integration table.

The main checks include:

- Expected input dimensions and required schema.
- Exact LR-pair identity and canonical-key uniqueness.
- Candidate consistency classes and provenance fields.
- Patient identifiers and six locked direction mappings.
- Direction-level record arithmetic.
- Inherited core-gene separation.
- Numeric, finite, and non-negative properties of existing expression fields.
- Completeness and type of the inherited compatibility field.
- Candidate-class, direction, patient-coverage, and final QC summaries.

The canonical key remains `Exact_LR_Pair × Direction_ID × Patient_Index`.

#### STEP11-5F-2: Patient-Level Compatibility Consistency Analysis

This step summarizes the inherited patient-level compatibility values.

The workflow:

- Verifies the canonical schema and predefined global structure.
- Checks canonical-key identity and uniqueness.
- Audits candidate-class consistency and direction mappings.
- Verifies inherited compatibility and provenance fields.
- Reshapes compatibility values into patient-level patterns.
- Generates descriptive summaries by LR pair, direction, patient, candidate class, and compatibility pattern.
- Exports descriptive summaries and final QC results.

Compatibility fractions are descriptive summaries. They are not interaction-strength estimates or communication probabilities.

#### STEP11-5F-3: Patient-Level Expression Context Integration

This step integrates inherited expression values into patient-level LR context summaries.

The main operations include:

- Exact input-schema and canonical-key checks.
- Exact LR-pair and direction-mapping verification.
- Candidate-class and provenance checks.
- Core-gene separation confirmation.
- Expression-field integrity and mean/difference definition checks.
- Inherited compatibility-field validation.
- Patient-level expression-context construction.
- LR-pair, direction, candidate-class, compatibility, annotation, and pathway summaries.
- Global expression-context summary and final QC export.

The expression and compatibility values are inherited from the locked canonical input. They are not recalculated.

#### STEP11-5F-4: Final Patient-Level LR Context Summary and Directional Consistency Audit

This step generates the final descriptive summaries of the locked genome-wide patient-level LR candidate universe.

The workflow:

1.  Validates the canonical input schema and predefined global structure.
2.  Audits canonical keys and exact LR-pair identity.
3.  Verifies the six sender–receiver mappings.
4.  Checks candidate consistency classes and provenance fields.
5.  Confirms inherited core-gene separation.
6.  Checks expression-field integrity and arithmetic definitions.
7.  Validates compatibility, annotation, and pathway fields.
8.  Generates patient-level compatibility patterns.
9.  Generates directional consistency summaries.
10. Constructs a final master table with one row per exact LR pair.
11. Generates candidate-class, annotation, pathway, and global context summaries.
12. Exports the final QC table.

The expected 2,405 patient-level LR-direction pattern records are a predefined structural QC criterion rather than a newly estimated biological result.

The final master table is a descriptive integration of inherited values. It does not represent new statistical inference or functional validation.

### 4.6 Figure 7 plotting

#### STEP11-5F-5: Final Figure 7 Plotting

The Figure 7 workflow generates four independent panels from the locked upstream LR results.

The two plotting implementations are retained separately.

**Version 1 — Basic plotting version**

Version 1 is the original basic plotting implementation.

Input files:

- `STEP11-5E-3B_FINAL_Canonical_Patient_Level_LR_Integration.tsv`
- `STEP11-5F-4C_FINAL_Master_LR_Descriptive_Context_Table.tsv`

Panels:

- Figure 7A: Patient-level candidate consistency.
- Figure 7B: Six sender–receiver directions.
- Figure 7C: Patient-by-direction compatibility heatmap.
- Figure 7D: Resource-level annotation summary.

Each panel is exported independently as PDF and TIFF. The plotting QC summary is printed to the console.

**Version 2 — Revised plotting version**

Version 2 is the subsequent revised plotting implementation. It uses the same locked canonical patient-level table and final master descriptive table.

Panels:

- Figure 7A: Patient-level consistency.
- Figure 7B: Six sender–receiver directions with explicit 50:50 sender/receiver visual segmentation.
- Figure 7C: Patient-by-direction compatibility percentage heatmap.
- Figure 7D: Resource-level annotation summary.

Figure 7C calculates descriptive compatibility percentages from the inherited `Compatible` field. Figure 7D summarizes unique exact LR pairs by the standardized resource annotation categories.

The workflow includes structural checks for the expected patient, direction, canonical-record, and exact LR-pair counts. Direction-specific expected counts are retained as predefined QC values.

The two versions are stored separately to preserve the original plotting implementation and prevent output overwriting.

Both versions are presentation-only. They do not recalculate upstream expression or compatibility values and do not perform statistical testing or biological ranking.

### 4.7 Master-table and supplementary-table preparation

#### STEP11-5F-4C: Master Table Audits

The script performs several read-only audits of the locked final master table.

**Master descriptive table audit**

The audit inspects:

- Table dimensions and column names.
- Column classes and missing values.
- Exact duplicate rows.
- Exact LR-pair frequencies.
- Sender, receiver, and direction distributions.
- Patient and candidate-consistency fields.
- Pathway and resource-annotation fields.
- Expression-context fields.
- Manuscript-facing candidate columns.

The audit report is exported as:

`STEP11-5F-4C_Master_LR_Descriptive_Context_Table_Audit.txt`

**Manuscript table mapping audit**

This audit inspects the potential mapping between the locked master-table fields and manuscript-facing Supplementary Table fields.

It covers:

- Exact LR-pair identity.
- Direction and direction-label information.
- Direction-specific compatibility fractions.
- Patient-level compatibility summaries.
- Candidate consistency and provenance.
- Pathway and resource annotation.
- Core-gene involvement.
- Ligand and receptor expression context.

It also constructs a provisional mapping between manuscript-facing fields and source-table columns. This mapping audit does not generate an additional Supplementary Table or Excel workbook.

**Direction mapping audit**

This audit examines the relationship between the six direction-specific compatibility fraction fields and the observed sender–receiver direction labels.

The main checks include:

- Direction-label integrity.
- Missingness and range summaries.
- Direction-field and label overlap.
- Single-direction observation summaries.
- Candidate direction mapping and ambiguity.
- Overall compatibility fraction consistency.

The output is:

`STEP11-5F-4C_Direction_Mapping_Audit.txt`

These audits do not recalculate LR compatibility or modify the locked master table.

#### STEP11-5F-4D: Final Supplementary Table Structure Lock and Generation

This step constructs a manuscript-facing Supplementary Table from the locked master table.

The expected input structure is:

- 585 rows.
- 36 columns.
- One row per exact LR pair.

The manuscript-facing table contains 16 columns.

The main operations include:

- Input dimension and required-field checks.
- Exact LR-pair uniqueness verification.
- Locked six-direction mapping.
- Direction-specific compatibility formatting.
- Manuscript-facing field selection and renaming.
- Patient-compatibility terminology checks.
- Ligand and receptor expression-summary formatting.
- Data-transfer and output-structure checks.
- Excel workbook generation.
- Final text audit report generation.

The Excel workbook contains two worksheets:

- `Supplementary Table`
- `Table Note`

Main outputs:

- `Supplementary_Table_Sx_Patient_Level_LR_Descriptive_Context_FINAL.xlsx`
- `STEP11-5F-4D_FINAL_Supplementary_Table_Lock_Audit.txt`

The workflow transfers and formats existing locked values. It does not recalculate LR compatibility, patient counts, expression summaries, or statistical results.

### 4.8 Figure and supplementary-table traceability audits

| Workflow | Main purpose | Output type |
|----|----|----|
| `STEP11-6A-0` | Identify candidate Figure 7 source files. | Console diagnostics |
| `STEP11-6A-1` | Inspect the structural characteristics of candidate Figure 7 panel sources. | Console diagnostics |
| `STEP11-6A` | Perform formal Figure 7 source traceability checks. | RDS, conditional TSVs, TXT |
| `STEP11-6A-2` | Compare exact LR-pair identity and ordering between the master and supplementary tables. | Console diagnostics |
| `STEP11-6A-3` | Diagnose exact pair-set discrepancies and normalized identifier differences. | Console diagnostics |
| `STEP11-6A-6` | Localize pair-set, duplicate, and ordering discrepancies. | Console diagnostics |
| `STEP11-6A-10` | Audit ASCII-normalized direction-label traceability. | Console diagnostics |

#### STEP11-6A-0: Figure 7 Source Identification Audit

This step checks eight predefined candidate source files from the STEP11-5E and STEP11-5F result directories.

The audit includes:

- Candidate-file existence checks.
- Dimensions and column-name inspection.
- Exact LR-pair field discovery.
- Direction-field discovery.
- Compatibility-field discovery.
- Detailed inspection of candidate Figure 7B and Figure 7C sources.
- Detailed inspection of the final master table.

This step provides console-based source inspection only. It does not automatically select a Figure 7 source or replace the formal traceability audit.

#### STEP11-6A-1: Figure 7 Panel Source Structural Audit

This step inspects the structural characteristics of candidate Figure 7A–7D source files.

The main inputs include:

- The STEP11-5F-4C locked master LR descriptive context table.
- The Figure 7 Version 2 directional LR-pair count summary.
- The STEP11-5F-4A patient-direction compatibility pattern.
- Relevant STEP11-5F-2 and STEP11-5F-3 candidate source summaries.

The script reports file availability, dimensions, column names, selected records, and detected LR-pair and direction-related fields.

This step produces console output only.

#### STEP11-6A: Formal Figure 7 Traceability Audit

This step evaluates the traceability of the generated Figure 7 panels to their locked source tables.

The required in-memory objects include:

- `master`
- `supp_table`

The actual Figure 7 plotting-data file must be manually confirmed by setting `figure7_source_file`. The script does not automatically select the first candidate source file.

The audit includes:

- Locked master and supplementary-table dimension checks.
- Required master-field checks.
- Exact LR-pair identity checks.
- Direction-field identification and normalization.
- Direction-to-master traceability checks.
- Pair-level master-field transfer comparisons.
- Figure 7A consistency-class traceability.
- Figure 7B direction and count traceability.
- Figure 7C LR-pair, direction, compatibility-fraction, and patient-pattern traceability.
- Figure 7D resource-annotation source checks.
- Core-gene/LR-node separation checks.
- Final traceability status reporting.

Potential outputs include:

- `STEP11-6A_Figure7_Source_Audit_Object.rds`
- `STEP11-6A_Figure7_Direction_Summary.tsv`, when available.
- `STEP11-6A_Figure7_Field_Traceability.tsv`, when available.
- `STEP11-6A_LR_Pairs_Missing_From_Master.tsv`, when missing pairs are identified.
- `STEP11-6A_FINAL_Traceability_Audit.txt`

The workflow does not recalculate LR compatibility, patient-level expression, or statistical results. It does not modify the locked master or supplementary table.

#### STEP11-6A-2, STEP11-6A-3, and STEP11-6A-6: Supplementary Table Exact LR-Pair Diagnostics

These read-only diagnostic steps compare the exact LR-pair identifiers in the locked master table and the generated Supplementary Table.

The checks include:

- Master and supplementary-table row counts.
- Unique exact LR-pair counts.
- Exact vector identity.
- Order-independent set identity.
- Missing and unexpected LR pairs.
- Duplicate LR-pair identifiers.
- Positional ordering differences.
- Whitespace and non-breaking-space normalization.
- Diagnostic ligand/receptor identifier decomposition.
- Classification of exact identity, ordering differences, and pair-set discrepancies.

String normalization is used only for diagnosis. It does not modify the locked source tables or regenerate the supplementary data.

#### STEP11-6A-10: ASCII Direction Traceability Audit

This step compares the direction labels in the locked master table with those in the Supplementary Table.

The main operations include:

- Reading the master `Direction_Labels` field.
- Reading the manuscript-facing direction field from the Supplementary Table.
- Converting Unicode arrows to ASCII arrows for comparison.
- Normalizing whitespace and direction separators.
- Comparing direction strings by row position.
- Validating tokens against the six predefined sender–receiver directions.
- Comparing direction-label frequencies.
- Reporting direction completeness and traceability status.

The ASCII normalization is diagnostic only. The source master table and Excel workbook are not modified.

STEP11-6A-2, STEP11-6A-3, STEP11-6A-6, and STEP11-6A-10 produce console-based diagnostics only. They do not independently export additional tables, figures, or audit files.

## 5. Key Candidate Definitions and Data Structures

### 5.1 Matched-patient analysis

The patient-level compatibility analysis is restricted to the four patients represented in all three predefined cellular compartments.

The patient-consistency categories include:

- Compatible in at least three of four patients.
- Compatible in all four patients.

The ≥3/4 and 4/4 candidate sets are compared and retained as distinct candidate definitions.

### 5.2 Exact LR-pair identity

The canonical exact LR-pair identifier is:

`Exact_LR_Pair`

It represents the ligand–receptor gene pair and is distinguished from resource-level interaction identifiers such as `interaction_name`.

The canonical patient-level key is:

`Exact_LR_Pair × Direction_ID × Patient_Index`

This key defines the patient-level LR record used in the downstream integration and descriptive summaries.

### 5.3 Provenance fields

`Candidate_Set` and `interaction_name` are retained as provenance information.

They are not treated as canonical biological identifiers or required to be invariant within a canonical patient-level key.

### 5.4 Locked candidate structure

The downstream workflow checks the predefined candidate structure:

- 585 exact LR pairs.
- 492 pairs in the 4/4 consistency class.
- 93 pairs in the 3/4-only consistency class.
- Four matched patients.
- Six predefined sender–receiver directions.

These values are used as structural QC expectations in the relevant downstream steps.

### 5.5 Core-gene separation

AOC3, MS4A7, and IL7R are the prioritized study genes characterized in the broader LUAD analysis.

Their separation from the genome-wide LR candidate records is audited explicitly. The LR compatibility workflow does not establish a signaling chain mediated by these three genes.

## 6. Output Organization

The workflow writes its results under:

`08_Ligand_Receptor/results/`

The main output groups include:

| Output directory or workflow | Main contents |
|----|----|
| `STEP11-4_CoreGene_Cellular_Localization_Integration/` | Core-gene localization integration tables and locked RDS object |
| `STEP11-5A_CCI_LR_STEP10_Compatibility_Audit/` | CCI_LR and STEP10 compatibility summaries and audit object |
| `STEP11-5D-0_LR_Detection_Consistency_Sensitivity_Audit/` | Detection-rule sensitivity audit tables |
| `STEP11-5D-2_Genomewide_LR_Candidate_Structural_Audit/` | Candidate structural audit tables |
| `STEP11-5D-3_Candidate_Set_Definition_Lock/` | Locked candidate-set object and QC summaries |
| `STEP11-5E-2_Candidate_Expression_Magnitude/` | Candidate expression summaries and patient-level expression records |
| `STEP11-5E-3_Patient_Level_Candidate_Integration/` | Canonical patient-level integration and provenance summaries |
| `STEP11-5F-1_Canonical_Input_QC/` | Canonical input QC tables |
| `STEP11-5F-2_Patient_Level_Compatibility_Consistency/` | Patient-level compatibility summaries |
| `STEP11-5F-3_Patient_Level_Expression_Context/` | Patient-level expression-context summaries |
| `STEP11-5F-4_Patient_Level_LR_Final_Context_Summary/` | Final master table, patient-direction patterns, descriptive summaries, and audit reports |
| `STEP11-5F-4D_FINAL_Supplementary_Table/` | Final Supplementary Table Excel workbook and lock audit |
| `STEP11-5F-5_Figure7_Final_Plotting/Version1/` | Original basic Figure 7 plotting outputs |
| `STEP11-5F-5_Figure7_Final_Plotting/Version2/` | Revised Figure 7 plotting outputs |

Some audit steps generate only in-memory objects or console output. These are not independent exported result files and are therefore not listed as file outputs.

## 7. Analytical Boundaries and Interpretation

The workflow characterizes recurrent expression compatibility between predefined ligand and receptor genes across selected cellular compartments.

The following distinctions are important:

- Expression compatibility does not establish direct molecular interaction.
- Detection of a ligand and receptor does not demonstrate that the corresponding proteins interact.
- Compatibility fractions are not estimates of interaction specificity, interaction strength, or communication probability.
- Expression magnitude does not measure interaction strength.
- Resource annotations describe reference-resource classifications and do not establish functional pathway activation.
- Candidate-set recurrence across four matched patients is not population-level validation.
- The CPM ≥ 1 sensitivity analysis evaluates robustness to a stricter detection threshold; it is not independent validation of functional communication.
- No expression-matched permutation test or other background null-model analysis is performed.
- No new CellChat communication-probability calculation, pathway enrichment, statistical significance testing, or causal signaling inference is performed in the downstream descriptive integration steps.

The candidate LR tables should therefore be interpreted as expression-based compatibility results within the predefined reference resource and matched-patient subset.

## 8. Reproducibility Notes

### 8.1 Software and packages

The workflow uses R and the packages required by the individual script sections.

Important package dependencies include:

- `here` for project-relative path initialization.
- `data.table` for table processing and data manipulation.
- `ggplot2` for Figure 7 plotting.
- `scales` for plotting color palettes in the revised Figure 7 implementation.
- `openxlsx` for Supplementary Table generation and Excel-file reading.
- `readr`, `dplyr`, and `stringr` for the supplementary-table preparation workflow.

Additional packages may be checked or loaded by individual workflow sections.

The exact R version and package versions should be documented in the accompanying environment records rather than inferred from the script documentation.

### 8.2 Random seed

No single global random seed is specified for the complete workflow. Any seed requirement should be documented at the individual analysis step where a stochastic operation is actually performed.

### 8.3 Input and output paths

All module input and output paths are constructed relative to the GitHub project root.

The script does not require a machine-specific absolute path or manual `setwd()` call.

### 8.4 Reproducibility and interpretation safeguards

- Existing STEP10 pseudobulk matrices are reused rather than re-aggregated.
- Locked upstream expression and compatibility values are audited and summarized rather than recalculated in downstream integration steps.
- Candidate-set definitions and Direction_ID mappings are inherited from the preceding locked outputs.
- Provenance fields are retained and distinguished from canonical biological identifiers.
- Figure and supplementary-table workflows preserve the locked analytical values.
- Traceability audits are designed to identify discrepancies rather than silently alter source data.

The script revision focuses on project-relative path management and standardized documentation while preserving the original analytical procedures and parameters.
