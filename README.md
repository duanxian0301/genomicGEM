# ALPS GenomicSEM Workflow

This repository contains the cleaned and reproducible analysis scripts for the ALPS GenomicSEM workflow, from LDSC QC through native WSL factor GWAS reruns and LDSC validation of the final latent factors.

## Repository layout

- `scripts/01_ldsc_efa_cfa/`
  LDSC preprocessing, exploratory factor analysis, and CFA model refinement.
- `scripts/02_factor_gwas/`
  Factor GWAS input preparation, native WSL execution, chunk merging, and standardized export.
- `scripts/03_validation/`
  LDSC validation for the native rerun factor GWAS outputs.
- `docs/`
  Environment notes and execution details for the native WSL setup.

## Analysis order

1. `scripts/01_ldsc_efa_cfa/step1_ldsc_alps.R`
   Build munged summary statistics and run initial LDSC QC for the five ALPS traits.
2. `scripts/01_ldsc_efa_cfa/step2_efa_alps.R`
   Run EFA and factor-number selection.
3. `scripts/01_ldsc_efa_cfa/step3_cfa_alps.R`
   Run initial CFA evaluation.
4. `scripts/01_ldsc_efa_cfa/step3b_cfa_model_search.R`
   Explore alternative CFA structures.
5. `scripts/01_ldsc_efa_cfa/step3c_cfa_constrained_search.R`
   Search constrained CFA variants.
6. `scripts/01_ldsc_efa_cfa/step3d_cfa_refined_final.R`
   Finalize the refined CFA model.
7. `scripts/02_factor_gwas/step4_prepare_factor_gwas_inputs.R`
   Prepare factor GWAS chunked inputs.
8. `docs/step11_wsl_native_environment.md`
   Record the working native WSL GenomicSEM environment.
9. `scripts/02_factor_gwas/step11_run_factor_gwas_native_wsl.R`
   Run the native WSL factor GWAS.
10. `scripts/02_factor_gwas/step11b_merge_native_factor_gwas_results.R`
    Merge native chunked factor GWAS outputs.
11. `scripts/02_factor_gwas/step12_export_native_standard_factor_gwas.R`
    Export standardized `txt` files for `F1` and `F2`.
12. `scripts/03_validation/step13_ldsc_native_factor_validation.R`
    Run LDSC validation on the native rerun outputs.

## Utility launchers

- `scripts/02_factor_gwas/resume_native_factor_gwas_lowload.ps1`
  Resume the native factor GWAS with low parallelism.
- `scripts/02_factor_gwas/resume_native_factor_gwas_3x8.ps1`
  Resume the native factor GWAS with three parallel groups.

## Kept outputs

The repository tracks scripts and lightweight documentation only. Large generated artifacts, temporary debug scripts, WSL scratch folders, and patched failed attempts are excluded by `.gitignore`.
