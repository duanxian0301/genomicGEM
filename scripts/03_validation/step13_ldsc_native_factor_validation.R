library(GenomicSEM)
library(data.table)

root_dir <- "D:/文章/GS/GWAS"
std_dir <- file.path(root_dir, "step11_factor_gwas_native_results", "merged", "standard_txt")
ascii_work_dir <- "D:/codex/GenomicSEM/step13_native_ldsc_work"
out_dir <- "D:/codex/GenomicSEM/step13_native_ldsc_validation"

dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

source_files <- c(
  file.path(std_dir, "ALPS_F1_factorGWAS_native_standard.txt"),
  file.path(std_dir, "ALPS_F2_factorGWAS_native_standard.txt")
)
trait.names <- c("F1_native", "F2_native")

if (!all(file.exists(source_files))) {
  stop("Missing native standard files in: ", std_dir)
}

files <- file.path(ascii_work_dir, paste0(trait.names, ".txt"))
for (i in seq_along(source_files)) {
  file.copy(source_files[i], files[i], overwrite = TRUE)
}

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

proxy_n <- fread(file.path(root_dir, "step10_ldsc_validation", "factor_N_proxy.tsv"))
N <- as.numeric(proxy_n$N_proxy[match(c("F1", "F2"), proxy_n$factor)])

munge(
  files = files,
  hm3 = hm3,
  trait.names = trait.names,
  N = N
)

traits <- file.path(ascii_work_dir, paste0(trait.names, ".sumstats.gz"))

ldsc_out <- ldsc(
  traits = traits,
  sample.prev = rep(NA, length(traits)),
  population.prev = rep(NA, length(traits)),
  ld = ld,
  wld = wld,
  trait.names = trait.names,
  ldsc.log = file.path(out_dir, "native_factor_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "native_factor_ldsc.rds"))
saveRDS(names(ldsc_out), file.path(out_dir, "native_factor_ldsc_names.rds"))

summary_dt <- data.table(
  trait = colnames(ldsc_out$S),
  h2 = diag(ldsc_out$S)
)

fwrite(summary_dt, file.path(out_dir, "native_factor_ldsc_summary.tsv"), sep = "\t")
