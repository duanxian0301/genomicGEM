# GenomicSEM Workflow

This repository contains the cleaned analysis scripts for the ALPS GenomicSEM workflow.

## Final workflow

1. `step1_ldsc_alps.R`
   Build munged summary statistics and run initial LDSC QC for the five ALPS traits.
2. `step2_efa_alps.R`
   Run EFA and factor-number selection.
3. `step3_cfa_alps.R`
   Run initial CFA evaluation.
4. `step3b_cfa_model_search.R`
   Explore alternative CFA structures.
5. `step3c_cfa_constrained_search.R`
   Search constrained CFA variants.
6. `step3d_cfa_refined_final.R`
   Finalize the refined CFA model.
7. `step4_prepare_factor_gwas_inputs.R`
   Prepare factor-GWAS chunked inputs.
8. `step11_wsl_native_environment.md`
   Notes for the working native WSL GenomicSEM environment.
9. `step11_run_factor_gwas_native_wsl.R`
   Native WSL factor-GWAS run script.
10. `step11b_merge_native_factor_gwas_results.R`
    Merge native chunked factor-GWAS outputs.
11. `step12_export_native_standard_factor_gwas.R`
    Export standardized `txt` files for `F1` and `F2`.
12. `step13_ldsc_native_factor_validation.R`
    Run LDSC validation on the native rerun outputs.

## Utility launchers

- `resume_native_factor_gwas_lowload.ps1`
  Resume the native factor GWAS with low parallelism.
- `resume_native_factor_gwas_3x8.ps1`
  Resume the native factor GWAS with three parallel groups.

## Notes

- Earlier patched Windows factor-GWAS attempts, smoketests, and temporary debug scripts were removed.
- Generated intermediate artifacts are ignored by `.gitignore`.
