# genomicGEM

This repository contains the cleaned script collection for the glymphatic-system GenomicSEM workflow and its downstream post-GWAS extensions.

The repository is organized into sequential analysis stages so that the original `genomicGEM` workflow is preserved and the later cross-trait work can continue without mixing files at the top level.

## Repository layout

### `step01_step13_genomicgem_core`

The original GenomicSEM workflow, including:

1. LDSC preprocessing for the ALPS traits
2. EFA and CFA model development
3. factor-GWAS input preparation
4. native factor-GWAS rerun in WSL
5. standardized factor-GWAS export
6. LDSC validation of the native rerun

### `step14_postgwas_ldsc_to_pwcoco`

The downstream post-GWAS extension of the factor workflow, including:

1. LDSC for `F1/F2 x AD/PD/LBD`
2. extended LDSC for ALPS-family traits
3. MiXeR input preparation and result export
4. pleioFDR / conjFDR workflows
5. candidate locus definition
6. classical coloc
7. PWCoCo and sensitivity analyses
8. supplementary table and integrated figure generation
9. Chinese manuscript draft for the LDSC-to-PWCoCo section

## Notes

- Large intermediate outputs and summary statistics are not versioned in this repository.
- The working analysis outputs were generated outside this repository, mainly under `D:\\文章\\GS\\postgwas`.
- Repository organization was updated to keep the original `genomicGEM` code intact while adding the downstream post-GWAS extension in a separate stage folder.
