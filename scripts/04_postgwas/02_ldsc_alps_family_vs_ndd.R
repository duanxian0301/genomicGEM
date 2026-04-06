library(GenomicSEM)
library(data.table)

out_dir <- "D:/文章/GS/postgwas/01b_ldsc_alps_family_vs_ndd"
ascii_work_dir <- "D:/codex/GenomicSEM/postgwas_ldsc_alps_family_vs_ndd_work"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

trait_map <- data.table(
  trait = c(
    "F1", "F2",
    "aALPS", "Left_ALPS", "mALPS", "pALPS", "Right_ALPS",
    "Mean_ALPS", "tALPS",
    "AD", "PD", "LBD"
  ),
  sumstats = c(
    "D:/文章/GS/postgwas/01_ldsc_f1f2_ndd/../..",
    "D:/文章/GS/postgwas/01_ldsc_f1f2_ndd/../..",
    "D:/文章/GS/GWAS/aALPS.sumstats.gz",
    "D:/文章/GS/GWAS/Left_ALPS.sumstats.gz",
    "D:/文章/GS/GWAS/mALPS.sumstats.gz",
    "D:/文章/GS/GWAS/pALPS.sumstats.gz",
    "D:/文章/GS/GWAS/Right_ALPS.sumstats.gz",
    "D:/文章/GS/GWAS/step10_ldsc_validation/sumstats_inputs/Mean_ALPS.sumstats.gz",
    "D:/文章/GS/GWAS/step10_ldsc_validation/sumstats_inputs/tALPS.sumstats.gz",
    "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work/AD.sumstats.gz",
    "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work/PD.sumstats.gz",
    "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work/LBD.sumstats.gz"
  )
)

trait_map[trait == "F1", sumstats := "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work/F1.sumstats.gz"]
trait_map[trait == "F2", sumstats := "D:/codex/GenomicSEM/postgwas_ldsc_f1f2_ndd_work/F2.sumstats.gz"]

if (!all(file.exists(trait_map$sumstats))) {
  missing_files <- trait_map[!file.exists(sumstats)]
  print(missing_files)
  stop("Some required .sumstats.gz files are missing.")
}

ldsc_out <- ldsc(
  traits = trait_map$sumstats,
  sample.prev = rep(NA, nrow(trait_map)),
  population.prev = rep(NA, nrow(trait_map)),
  ld = ld,
  wld = wld,
  trait.names = trait_map$trait,
  ldsc.log = file.path(out_dir, "alps_family_vs_ndd_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "alps_family_vs_ndd_ldsc.rds"))

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

alps_family <- c("F1", "F2", "aALPS", "Left_ALPS", "mALPS", "pALPS", "Right_ALPS", "Mean_ALPS", "tALPS")
ndd_traits <- c("AD", "PD", "LBD")

pair_grid <- as.data.table(expand.grid(
  trait1 = alps_family,
  trait2 = ndd_traits,
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

pair_grid[, panel := fifelse(
  trait1 %in% c("F1", "F2"), "factor",
  fifelse(trait1 %in% c("Mean_ALPS", "tALPS"), "composite", "original")
)]

pair_grid[, fdr_rg_global := p.adjust(p_rg, method = "BH")]
pair_grid[, fdr_cov_global := p.adjust(p_cov, method = "BH")]
pair_grid[, fdr_rg_within_panel := p.adjust(p_rg, method = "BH"), by = panel]
pair_grid[, fdr_cov_within_panel := p.adjust(p_cov, method = "BH"), by = panel]

comparison_table <- pair_grid[, .(
  trait1, trait2, panel, rg, rg_se, p_rg, fdr_rg_global, fdr_rg_within_panel
)]

ad_rank <- comparison_table[trait2 == "AD"][order(p_rg)]
pd_rank <- comparison_table[trait2 == "PD"][order(p_rg)]
lbd_rank <- comparison_table[trait2 == "LBD"][order(p_rg)]

fwrite(trait_map, file.path(out_dir, "input_manifest.tsv"), sep = "\t")
fwrite(diag_dt, file.path(out_dir, "alps_family_vs_ndd_ldsc_h2.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"), file.path(out_dir, "alps_family_vs_ndd_ldsc_rg_matrix.tsv"), sep = "\t")
fwrite(pair_grid, file.path(out_dir, "alps_family_vs_ndd_requested_pairs.tsv"), sep = "\t")
fwrite(comparison_table, file.path(out_dir, "alps_family_vs_ndd_comparison_table.tsv"), sep = "\t")
fwrite(ad_rank, file.path(out_dir, "AD_ranked_rg.tsv"), sep = "\t")
fwrite(pd_rank, file.path(out_dir, "PD_ranked_rg.tsv"), sep = "\t")
fwrite(lbd_rank, file.path(out_dir, "LBD_ranked_rg.tsv"), sep = "\t")
