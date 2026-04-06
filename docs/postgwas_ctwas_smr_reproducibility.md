# Post-GWAS cTWAS/SMR Reproducibility Notes

This document describes how to rerun the ALPS factor vs NDD post-GWAS analysis branch used for cTWAS and SMR follow-up.

## Scope

The branch covers:

- trait-wise multigroup cTWAS for `F1`, `F2`, `AD`, `PD`, and `LBD`
- bulk-brain SMR using BrainMeta v2 cortex and GTEx v8 brain tissues
- Bryois 2022 cell-type SMR across 8 major brain cell types
- integration of `conjFDR`, `coloc`, `PWCoCo`, `cTWAS`, and `SMR` into candidate gene tables

## External Dependencies

### cTWAS

- R with `ctwas`
- FUSION-format brain expression weights
- PLINK executable
- b37/hg19 LD reference panel

### SMR

- SMR executable from Yang Lab
- 1000G EUR LD reference
- BrainMeta v2 cortex BESD
- GTEx v8 brain BESD
- Bryois 2022 cell-type eQTL summary files for conversion into BESD

## Key Methodological Decisions

### cTWAS naming

The implemented workflow is **trait-wise multigroup cTWAS**, followed by cross-trait comparison. It is not a full pairwise factor-disease multigroup cTWAS joint model.

### SMR frequency handling

`AD`, `PD`, and `LBD` retain effect-allele frequency and can be run with standard frequency checks.

`F1` and `F2` standardized factor GWAS exports retain `MAF` in the `FREQ` column. Because SMR expects effect-allele frequency, factor-side SMR is run using `--disable-freq-ck`. These results are still useful for prioritization, but they should be described as more limited than the disease-side SMR analyses.

## Minimal Rerun Order

1. Prepare coloc/PWCoCo regions
2. Prepare cTWAS references
3. Run one cTWAS trait at a time
4. Prepare SMR sumstats
5. Run BrainMeta bulk SMR
6. Run GTEx bulk SMR
7. Build Bryois BESD references and run cell-type SMR
8. Refresh curated supplementary tables
9. Build candidate master integration tables

## Files Intended For Publication

The following outputs are the preferred publication-facing products:

- curated SMR tables from `summary_curated/`
- cTWAS primary and secondary gene tables
- candidate master and shortlist tables

The following outputs should generally remain as working files rather than direct supplementary workbook sheets:

- raw per-chromosome `.smr` files
- raw Bryois BESD preparation intermediates
- queue status files
- transient smoke-test outputs

## Recommendation For Future Data Updates

If new GWAS data are substituted:

1. regenerate standardized GWAS input files first
2. rerun `25_prepare_smr_sumstats.py`
3. rerun cTWAS and SMR from their reference-preparation entry points
4. regenerate curated supplementary tables and candidate master integration

Do not reuse old candidate tables after changing GWAS inputs.
