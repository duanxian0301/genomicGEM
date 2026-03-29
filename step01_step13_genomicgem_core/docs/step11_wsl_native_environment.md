## Native WSL GenomicSEM environment for factor GWAS

This workflow avoids the earlier Windows patch that replaced internal
`lavaan::sem()` DWLS reorder calls. The WSL route uses a source install of the
current `GenomicSEM` master tarball plus an isolated R library.

### Why this environment

- The copied legacy package in WSL still failed with `ReorderModel` /
  `ReorderModelnoSNP not found`.
- A source install from `GenomicSEM_master.tar.gz` in WSL resolves the
  `usermodel()` and `userGWAS()` failures.
- Native `DWLS` + `fix_measurement=TRUE` now runs without the old monkey patch.

### WSL library

- R executable: `/usr/bin/Rscript`
- R version: `4.5.0`
- Isolated library: `/home/shenjing/R/genomicsem_fix_lib`

### Installed core packages

- `GenomicSEM`
- `lavaan`
- `data.table`
- `plyr`
- `e1071`
- `gdata`
- `doParallel`
- `foreach`
- `splitstackshape`
- `R.utils`
- `iterators`
- `mgsub`
- `readr`
- `stringr`
- `Matrix`
- `dplyr`
- `Rcpp`

### Path aliases in WSL

To avoid Unicode path issues, WSL symlinks are used:

- GWAS root: `/home/shenjing/gs_paths/gwas`
- External trait root: `/home/shenjing/gs_paths/ext`

### Minimal validation that passed

- `usermodel()` on the ALPS 2-factor model
- `userGWAS()` on a small SNP subset with:
  - `estimation = "DWLS"`
  - `fix_measurement = TRUE`
  - refined 2-factor model

### Recommended execution pattern

1. Run factor GWAS in WSL with `step11_run_factor_gwas_native_wsl.R`.
2. Merge chunks into a separate native-results directory.
3. Export standard summary statistics from the native outputs.
4. Re-run single-trait LDSC QC on native `F1` and `F2`.
