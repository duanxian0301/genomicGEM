# step01_step13_genomicgem_core

This folder preserves the original `genomicGEM` pipeline before the post-GWAS extension.

## Contents

- `scripts/`: the main stepwise scripts from LDSC through factor-GWAS validation
- `docs/`: environment notes and workflow documentation
- `launchers/`: PowerShell helper scripts for resuming native factor-GWAS runs

## Core workflow

1. `step1_ldsc_alps.R`
2. `step2_efa_alps.R`
3. `step3_cfa_alps.R`
4. `step3b_cfa_model_search.R`
5. `step3c_cfa_constrained_search.R`
6. `step3d_cfa_refined_final.R`
7. `step4_prepare_factor_gwas_inputs.R`
8. `step11_run_factor_gwas_native_wsl.R`
9. `step11b_merge_native_factor_gwas_results.R`
10. `step12_export_native_standard_factor_gwas.R`
11. `step13_ldsc_native_factor_validation.R`

## Notes

- File names were preserved during repository cleanup to keep script references stable.
- Only the directory layout was reorganized for readability; the original workflow code was not removed.
