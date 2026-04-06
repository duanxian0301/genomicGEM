# ALPS GenomicSEM Workflow

This repository contains the cleaned, reproducible workflow used to analyze ALPS traits with GenomicSEM, including LDSC quality control, EFA/CFA model building, native WSL factor GWAS reruns, and LDSC validation for the final latent factors.

## Overview

The project was organized around two latent factors, `F1` and `F2`, derived from five ALPS-related GWAS traits:

- `aALPS`
- `Left_ALPS`
- `Right_ALPS`
- `mALPS`
- `pALPS`

The final workflow includes:

1. harmonizing original GWAS summary statistics
2. running multivariable LDSC
3. exploring and refining factor structure with EFA and CFA
4. generating native WSL factor GWAS results for `F1` and `F2`
5. validating the resulting factor GWAS with LDSC
6. extracting and clumping `Q_SNP` signals for factor-specific heterogeneity follow-up

## Repository layout

- `scripts/01_ldsc_efa_cfa/`
  LDSC preprocessing, exploratory factor analysis, and CFA refinement.
- `scripts/02_factor_gwas/`
  Factor GWAS input preparation, native WSL execution, chunk merging, and standardized export.
- `scripts/03_validation/`
  LDSC validation scripts for the native factor GWAS outputs.
- `scripts/04_postgwas/`
  Post-GWAS follow-up workflow for MiXeR, pleioFDR, coloc, PWCoCo, cTWAS, and SMR.
- `docs/`
  Environment notes and execution details.
- `results/`
  Suggested location for shareable summary outputs, tables, and figure-ready exports.
- `manuscript/`
  Suggested location for result text, supplementary material drafts, and figure/table legends.

## Recommended execution order

1. `scripts/01_ldsc_efa_cfa/step1_ldsc_alps.R`
   Munges summary statistics and runs initial LDSC QC for the five ALPS traits.
2. `scripts/01_ldsc_efa_cfa/step2_efa_alps.R`
   Runs EFA and factor-number selection.
3. `scripts/01_ldsc_efa_cfa/step3_cfa_alps.R`
   Runs the initial exploratory CFA evaluation.
4. `scripts/01_ldsc_efa_cfa/step3b_cfa_model_search.R`
   Explores alternative CFA structures.
5. `scripts/01_ldsc_efa_cfa/step3c_cfa_constrained_search.R`
   Searches constrained CFA variants.
6. `scripts/01_ldsc_efa_cfa/step3d_cfa_refined_final.R`
   Finalizes the exploratory refined CFA model.
7. `scripts/01_ldsc_efa_cfa/step3f_cfa_usermodel_native_wsl.R`
   Re-runs formal CFA/model comparison in native WSL via `GenomicSEM::usermodel()`.
8. `scripts/02_factor_gwas/step4_prepare_factor_gwas_inputs.R`
   Prepares chunked factor GWAS inputs.
9. `scripts/02_factor_gwas/step4b_ols_sensitivity_mini_compare.R`
   Runs a mini sensitivity comparison of `OLS = FALSE` vs `OLS = TRUE` factor-GWAS input processing.
10. `docs/step11_wsl_native_environment.md`
   Documents the working native WSL GenomicSEM environment.
11. `scripts/02_factor_gwas/step11_run_factor_gwas_native_wsl.R`
   Runs the native WSL factor GWAS.
12. `scripts/02_factor_gwas/step11b_merge_native_factor_gwas_results.R`
    Merges native chunked factor GWAS outputs.
13. `scripts/02_factor_gwas/step12_export_native_standard_factor_gwas.R`
    Exports standardized `txt` files for `F1` and `F2`.
14. `scripts/03_validation/step13_ldsc_native_factor_validation.R`
    Runs LDSC validation on the native rerun outputs.
15. `scripts/03_validation/step26_extract_qsnp_and_overlap.py`
    Extracts factor-specific `Q_SNP` result tables and overlap summaries.
16. `scripts/03_validation/step27_qsnp_clump_and_plots.py`
    Performs PLINK clumping for factor and `Q_SNP` signals and generates Manhattan/QQ plots.
17. `scripts/03_validation/step28_rebuild_supplement_with_qsnp.py`
    Rebuilds supplementary tables with the `Q_SNP` summaries and clumped lead loci.
18. `scripts/04_postgwas/README.md`
    Documents the post-GWAS follow-up workflow, including cTWAS and SMR.

## Utility launchers

- `scripts/02_factor_gwas/resume_native_factor_gwas_lowload.ps1`
  Resumes the native factor GWAS with lower parallelism.
- `scripts/02_factor_gwas/resume_native_factor_gwas_3x8.ps1`
  Resumes the native factor GWAS with three parallel groups.

## Input and output conventions

Expected raw GWAS columns:

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

Standardized factor GWAS exports also follow the same layout for downstream tools.

## Notes

- This repository tracks scripts and lightweight documentation only.
- Large generated artifacts, temporary debug scripts, WSL scratch folders, and failed patched attempts are excluded by `.gitignore`.
- Native WSL reruns were retained because they produced LDSC-valid factor GWAS outputs, whereas earlier patched Windows attempts did not.
- For formal CFA/model-fit reporting, prefer the `usermodel()` workflow in `step3f_cfa_usermodel_native_wsl.R`.
- The earlier `lavaan::sem(sample.nobs = 200)` scripts should be treated as exploratory model-search helpers, not as the primary confirmatory fit framework.
- ALPS factor-GWAS input preparation should use `OLS = TRUE` because the source ALPS GWAS were continuous traits analyzed with linear models.
