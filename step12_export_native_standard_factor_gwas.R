library(data.table)

root_dir <- "D:/文章/GS/GWAS"
native_dir <- file.path(root_dir, "step11_factor_gwas_native_results")
merged_dir <- file.path(native_dir, "merged")
std_dir <- file.path(merged_dir, "standard_txt")

dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(std_dir, showWarnings = FALSE, recursive = TRUE)

f1_file <- file.path(merged_dir, "ALPS_F1_factorGWAS_native_merged.tsv.gz")
f2_file <- file.path(merged_dir, "ALPS_F2_factorGWAS_native_merged.tsv.gz")

if (!file.exists(f1_file) || !file.exists(f2_file)) {
  stop("Merged native factor GWAS files are required before export.")
}

proxy_n <- fread(file.path(root_dir, "step10_ldsc_validation", "factor_N_proxy.tsv"))
n_lookup <- setNames(proxy_n$N_proxy, proxy_n$factor)

export_one <- function(infile, outfile, n_value) {
  dt <- fread(infile)
  out <- data.table(
    SNP = dt$SNP,
    CHR = dt$CHR,
    BP = dt$BP,
    A1 = dt$A1,
    A2 = dt$A2,
    FREQ = dt$MAF,
    BETA = dt$est,
    SE = dt$SE,
    P = dt$Pval_Estimate,
    N = as.integer(round(n_value))
  )
  fwrite(out, outfile, sep = "\t")
}

export_one(
  f1_file,
  file.path(std_dir, "ALPS_F1_factorGWAS_native_standard.txt"),
  n_lookup[["F1"]]
)

export_one(
  f2_file,
  file.path(std_dir, "ALPS_F2_factorGWAS_native_standard.txt"),
  n_lookup[["F2"]]
)

meta <- data.table(
  factor = c("F1", "F2"),
  standard_file = c(
    file.path(std_dir, "ALPS_F1_factorGWAS_native_standard.txt"),
    file.path(std_dir, "ALPS_F2_factorGWAS_native_standard.txt")
  ),
  n_proxy = c(n_lookup[["F1"]], n_lookup[["F2"]]),
  notes = "Native WSL GenomicSEM rerun; FREQ uses MAF; N uses weighted proxy N."
)

fwrite(meta, file.path(std_dir, "native_standard_txt_metadata.tsv"), sep = "\t")
