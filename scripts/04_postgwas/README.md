# Post-GWAS cTWAS and SMR Workflow

This module contains the post-GWAS follow-up workflow used for ALPS latent factors (`F1`, `F2`) and neurodegenerative disease GWAS (`AD`, `PD`, `LBD`). It covers:

- cross-trait LDSC/MiXeR/pleioFDR follow-up inputs
- locus definition for coloc and PWCoCo
- trait-wise multigroup cTWAS
- SMR using bulk brain and cell-type eQTL references
- integration of SNP-level and gene-level evidence into publication-ready candidate tables

The workflow is organized for reproducibility rather than as a single monolithic launcher. Scripts are grouped by analysis stage and can be rerun with updated GWAS inputs as long as the required file formats are preserved.

## What This Module Produces

- trait-wise multigroup cTWAS results for `F1`, `F2`, `AD`, `PD`, and optionally `LBD`
- bulk-brain SMR results using BrainMeta v2 cortex and GTEx v8 brain tissues
- Bryois 2022 cell-type SMR results across 8 major brain cell types
- candidate gene integration tables combining `conjFDR`, `coloc`, `PWCoCo`, `cTWAS`, and `SMR`

## Important Terminology

- The cTWAS stage implemented here should be described as **trait-wise multigroup cTWAS**.
- The factor-disease comparison stage should be described as **cTWAS-based cross-trait comparison / gene prioritization**, not as a full pairwise multigroup cTWAS joint model.
- `F1` and `F2` SMR runs are performed with `--disable-freq-ck` because the standardized factor GWAS exports retain `MAF` rather than effect-allele frequency in the `FREQ` column.

## Recommended Execution Order

### A. Locus definition and coloc/PWCoCo preparation

1. `14_define_coloc_pwcoco_loci.py`
2. `15_run_coloc_factor_ndd.R`
3. `16_prepare_pwcoco_r17_inputs.R`
4. `17_summarize_pwcoco_r17.py`
5. `18_prepare_pwcoco_region_inputs.R`
6. `19_summarize_pwcoco_regions.py`

### B. cTWAS

1. `20_prepare_mctwas_hg19_refs.R`
   - downloads Broad/FUSION GTEx brain weights
   - builds b37 LD block references for cTWAS
2. `21_run_mctwas_single_trait.R`
   - runs one trait at a time
   - supports initialization from successful traits when needed
3. `22_summarize_mctwas_f1f2_ndd.py` or `22_summarize_mctwas_f1f2_ndd.R`
4. `23_extract_ctwas_rich_results.R`
5. `24_prepare_ctwas_cross_trait_publishable_tables.py`

### C. SMR

1. `25_prepare_smr_sumstats.py`
   - converts standardized GWAS files into SMR input format
2. `26_run_smr_single_trait.py`
   - single-trait bulk SMR runner for BrainMeta or GTEx
3. `27_launch_smr_brainmeta_batch.ps1`
4. `29_resume_smr_brainmeta_missing.ps1` or `29b_resume_smr_brainmeta_missing.py`
5. `33_launch_smr_gtex_queue.py`

### D. Bryois 2022 cell-type SMR

1. `32_prepare_bryois_celltype_fastqtl.py`
   - filters fastQTL nominal files and generates `.epi/.esi` update tables
2. `34_prepare_bryois_reference.py`
   - builds and refreshes BESD references for one cell type and chromosome
3. `35_run_smr_bryois_single.py`
4. `38_run_bryois_all_until_done.py`
   - continuous controller that fills missing cell-type SMR outputs until all tasks finish

### E. Summary and publication-ready tables

1. `30_summarize_smr_and_integrate.py`
   - legacy broad summary builder
2. `37_refresh_smr_supplementary.py`
   - publication-oriented SMR table curation
3. `39_build_candidate_master.py`
   - integrates `SMR`, `cTWAS`, `coloc`, and `PWCoCo` into gene-level candidate tables

## Core Scripts vs Helper Scripts

### Core scripts to retain for future reruns

- `14_define_coloc_pwcoco_loci.py`
- `15_run_coloc_factor_ndd.R`
- `18_prepare_pwcoco_region_inputs.R`
- `19_summarize_pwcoco_regions.py`
- `20_prepare_mctwas_hg19_refs.R`
- `21_run_mctwas_single_trait.R`
- `24_prepare_ctwas_cross_trait_publishable_tables.py`
- `25_prepare_smr_sumstats.py`
- `26_run_smr_single_trait.py`
- `32_prepare_bryois_celltype_fastqtl.py`
- `34_prepare_bryois_reference.py`
- `35_run_smr_bryois_single.py`
- `37_refresh_smr_supplementary.py`
- `38_run_bryois_all_until_done.py`
- `39_build_candidate_master.py`

### Helper / queue / recovery scripts

- `27_launch_smr_brainmeta_batch.ps1`
- `28_write_initial_brainmeta_launch_status.py`
- `29_resume_smr_brainmeta_missing.ps1`
- `29b_resume_smr_brainmeta_missing.py`
- `33_launch_smr_gtex_queue.py`
- `36_launch_bryois_celltype_queue.py`

These helper scripts are retained because they document how the full analysis was operationalized on a desktop environment, but the core rerun logic should rely on the scripts listed above.

## Input Expectations

Standardized GWAS files should contain:

- `SNP`
- `CHR`
- `BP`
- `A1`
- `A2`
- `FREQ`
- `BETA`
- `SE`
- `P`
- `N`

Additional notes:

- `AD`, `PD`, and `LBD` use effect-allele frequency for SMR.
- `F1` and `F2` standardized exports retain `MAF`; therefore, factor-side SMR uses relaxed frequency-check settings and should be interpreted more cautiously than disease-side SMR.

## Publication-Facing Table Strategy

The final supplementary workbook should emphasize:

- Bonferroni-significant SMR hits with HEIDI support
- FDR-significant SMR hits with HEIDI support
- `p_SMR < 1e-4` exploratory SMR hits as secondary material
- cTWAS primary genes defined by posterior support
- a candidate master table integrating SNP-level and gene-level evidence

Avoid exporting every raw per-chromosome file into the publication workbook. Keep full raw outputs on disk and export only curated tables into the workbook.
