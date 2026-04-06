library(GenomicSEM)
library(data.table)

factor_std_dir <- "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt"
ndd_gwas_dir <- "D:/文章/4NDD/NDDGWAS"

ascii_work_dir <- "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work"
out_dir <- "D:/文章/GS/postgwas/01_ldsc_f1f2_ndd"

dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

factor_meta <- fread(file.path(factor_std_dir, "native_standard_txt_metadata.tsv"))
n_proxy <- setNames(as.numeric(factor_meta$n_proxy), factor_meta$factor)

trait_inputs <- data.table(
  trait = c("F1", "F2", "AD", "PD", "LBD"),
  source_file = c(
    file.path(factor_std_dir, "ALPS_F1_factorGWAS_native_standard.txt"),
    file.path(factor_std_dir, "ALPS_F2_factorGWAS_native_standard.txt"),
    file.path(ndd_gwas_dir, "AD.txt"),
    file.path(ndd_gwas_dir, "PD.txt"),
    file.path(ndd_gwas_dir, "LBD.txt")
  ),
  N = c(n_proxy["F1"], n_proxy["F2"], NA_real_, NA_real_, NA_real_)
)

if (!all(file.exists(trait_inputs$source_file))) {
  missing_files <- trait_inputs[!file.exists(source_file)]
  print(missing_files)
  stop("Some LDSC input files are missing.")
}

trait_inputs[, work_txt := file.path(ascii_work_dir, paste0(trait, ".txt"))]
trait_inputs[, sumstats := file.path(ascii_work_dir, paste0(trait, ".sumstats.gz"))]

for (i in seq_len(nrow(trait_inputs))) {
  ok <- file.copy(trait_inputs$source_file[i], trait_inputs$work_txt[i], overwrite = TRUE)
  if (!ok) {
    stop("Failed to copy input file for trait: ", trait_inputs$trait[i])
  }
}

need_munge <- !file.exists(trait_inputs$sumstats)
if (any(need_munge)) {
  munge(
    files = trait_inputs$work_txt[need_munge],
    hm3 = hm3,
    trait.names = trait_inputs$trait[need_munge],
    N = trait_inputs$N[need_munge]
  )
}

ldsc_out <- ldsc(
  traits = trait_inputs$sumstats,
  sample.prev = rep(NA, nrow(trait_inputs)),
  population.prev = rep(NA, nrow(trait_inputs)),
  ld = ld,
  wld = wld,
  trait.names = trait_inputs$trait,
  ldsc.log = file.path(out_dir, "f1f2_vs_ndd_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "f1f2_vs_ndd_ldsc.rds"))

S <- as.matrix(ldsc_out$S)
rownames(S) <- colnames(S)

ratio <- tcrossprod(1 / sqrt(diag(S)))
R <- S * ratio
rownames(R) <- rownames(S)
colnames(R) <- colnames(S)

se_cov <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S))
se_cov[lower.tri(se_cov, diag = TRUE)] <- sqrt(diag(ldsc_out$V))
se_cov[upper.tri(se_cov)] <- t(se_cov)[upper.tri(se_cov)]

ratio_vec <- ratio[lower.tri(ratio, diag = TRUE)]
V_std <- ldsc_out$V * tcrossprod(ratio_vec)

se_rg <- matrix(0, nrow(R), ncol(R), dimnames = dimnames(R))
se_rg[lower.tri(se_rg, diag = TRUE)] <- sqrt(diag(V_std))
se_rg[upper.tri(se_rg)] <- t(se_rg)[upper.tri(se_rg)]

diag_dt <- data.table(
  trait = colnames(S),
  h2 = diag(S),
  h2_se = diag(se_cov),
  h2_z = diag(S) / diag(se_cov),
  intercept = diag(ldsc_out$I)
)

pair_grid <- as.data.table(expand.grid(
  trait1 = c("F1", "F2"),
  trait2 = c("AD", "PD", "LBD"),
  stringsAsFactors = FALSE
))

pair_grid[, `:=`(
  covariance = mapply(function(x, y) S[x, y], trait1, trait2),
  covariance_se = mapply(function(x, y) se_cov[x, y], trait1, trait2),
  rg = mapply(function(x, y) R[x, y], trait1, trait2),
  rg_se = mapply(function(x, y) se_rg[x, y], trait1, trait2)
)]

pair_grid[, z_cov := covariance / covariance_se]
pair_grid[, p_cov := 2 * pnorm(abs(z_cov), lower.tail = FALSE)]
pair_grid[, z_rg := rg / rg_se]
pair_grid[, p_rg := 2 * pnorm(abs(z_rg), lower.tail = FALSE)]
pair_grid[, fdr_rg := p.adjust(p_rg, method = "BH")]
pair_grid[, fdr_cov := p.adjust(p_cov, method = "BH")]

input_manifest <- trait_inputs[, .(trait, source_file, copied_txt = work_txt, sumstats, N)]

fwrite(input_manifest, file.path(out_dir, "input_manifest.tsv"), sep = "\t")
fwrite(diag_dt, file.path(out_dir, "f1f2_vs_ndd_ldsc_h2.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"), file.path(out_dir, "f1f2_vs_ndd_ldsc_rg_matrix.tsv"), sep = "\t")
fwrite(pair_grid, file.path(out_dir, "f1f2_vs_ndd_requested_pairs.tsv"), sep = "\t")
