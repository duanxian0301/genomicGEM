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

## Repository layout

- `scripts/01_ldsc_efa_cfa/`
  LDSC preprocessing, exploratory factor analysis, and CFA refinement.
- `scripts/02_factor_gwas/`
  Factor GWAS input preparation, native WSL execution, chunk merging, and standardized export.
- `scripts/03_validation/`
  LDSC validation scripts for the native factor GWAS outputs.
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
   Runs the initial CFA evaluation.
4. `scripts/01_ldsc_efa_cfa/step3b_cfa_model_search.R`
   Explores alternative CFA structures.
5. `scripts/01_ldsc_efa_cfa/step3c_cfa_constrained_search.R`
   Searches constrained CFA variants.
6. `scripts/01_ldsc_efa_cfa/step3d_cfa_refined_final.R`
   Finalizes the refined CFA model.
7. `scripts/02_factor_gwas/step4_prepare_factor_gwas_inputs.R`
   Prepares chunked factor GWAS inputs.
8. `docs/step11_wsl_native_environment.md`
   Documents the working native WSL GenomicSEM environment.
9. `scripts/02_factor_gwas/step11_run_factor_gwas_native_wsl.R`
   Runs the native WSL factor GWAS.
10. `scripts/02_factor_gwas/step11b_merge_native_factor_gwas_results.R`
    Merges native chunked factor GWAS outputs.
11. `scripts/02_factor_gwas/step12_export_native_standard_factor_gwas.R`
    Exports standardized `txt` files for `F1` and `F2`.
12. `scripts/03_validation/step13_ldsc_native_factor_validation.R`
    Runs LDSC validation on the native rerun outputs.

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
