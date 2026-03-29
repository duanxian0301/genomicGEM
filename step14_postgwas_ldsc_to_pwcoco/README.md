# step14_postgwas_ldsc_to_pwcoco

This folder extends the `genomicGEM` pipeline from the factor-GWAS stage to the post-GWAS cross-trait analyses between glymphatic system latent factors and neurodegenerative diseases.

It contains the code used for the following analyses:

1. LDSC for `F1/F2 x AD/PD/LBD`
2. Extended LDSC for `F1/F2 + ALPS-family traits x AD/PD/LBD`
3. MiXeR input preparation and downstream result export
4. pleioFDR input preparation, batch configuration, and summary generation
5. candidate locus definition for coloc and PWCoCo
6. classical coloc for factor-level candidate regions
7. PWCoCo input preparation, sensitivity analyses, and summary generation
8. supplementary table export and integrated figure generation

## Structure

- `scripts/`: ordered scripts used in this stage
- `manuscript/`: Chinese draft of Methods and Results for the LDSC-to-PWCoCo section

## Notes

- Large intermediate outputs and summary statistics are not versioned here.
- The working analysis outputs were generated under `D:\\文章\\GS\\postgwas`.
- Scripts in this folder preserve the original numbering used during the analysis workflow.
