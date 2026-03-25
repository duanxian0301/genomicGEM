.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/home/shenjing/gs_paths/gwas"
work_dir <- file.path(root_dir, "ols_sensitivity_mini_compare")
sample_dir <- file.path(work_dir, "sampled_raw")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

source_files <- c(
  aALPS = file.path(root_dir, "aALPS.txt"),
  Left_ALPS = file.path(root_dir, "Left_ALPS.txt"),
  mALPS = file.path(root_dir, "mALPS.txt"),
  pALPS = file.path(root_dir, "pALPS.txt"),
  Right_ALPS = file.path(root_dir, "Right_ALPS.txt")
)

ref_file <- "/mnt/d/文章/Mibiogen and AD/多变量LDSC/reference.1000G.maf.0.005.txt"
ldsc_file <- file.path(root_dir, "step3_cfa_results", "ALPS5_ALL_CFA_ldsc.rds")
stopifnot(all(file.exists(source_files)), file.exists(ref_file), file.exists(ldsc_file))

sample_n <- 120000L
set.seed(20260325)

message("Creating mini sampled raw files...")
for (nm in names(source_files)) {
  dt <- fread(source_files[[nm]], nrows = sample_n)
  fwrite(dt, file.path(sample_dir, paste0(nm, ".txt")), sep = "\t")
}

run_sumstats <- function(ols_flag, out_prefix) {
  files <- file.path(sample_dir, paste0(names(source_files), ".txt"))
  outdir <- file.path(work_dir, out_prefix)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  oldwd <- getwd()
  setwd(outdir)
  on.exit(setwd(oldwd), add = TRUE)

  ss <- sumstats(
    files = files,
    ref = ref_file,
    trait.names = names(source_files),
    se.logit = rep(FALSE, length(files)),
    OLS = rep(ols_flag, length(files)),
    linprob = rep(FALSE, length(files)),
    N = rep(NA, length(files)),
    betas = NULL,
    info.filter = 0.6,
    maf.filter = 0.01,
    keep.indel = FALSE,
    parallel = FALSE,
    cores = 1
  )
  saveRDS(ss, file.path(outdir, paste0(out_prefix, "_sumstats.rds")))
  fwrite(as.data.frame(ss), file.path(outdir, paste0(out_prefix, "_sumstats.tsv.gz")), sep = "\t")
  ss
}

ss_false <- run_sumstats(FALSE, "ols_false")
ss_true <- run_sumstats(TRUE, "ols_true")

common_cols <- intersect(names(ss_false), names(ss_true))
join_cols <- intersect(c("SNP", "CHR", "BP", "A1", "A2"), common_cols)
cmp <- merge(
  as.data.table(ss_false)[, ..common_cols],
  as.data.table(ss_true)[, ..common_cols],
  by = join_cols,
  suffixes = c("_false", "_true")
)

fwrite(cmp, file.path(work_dir, "sumstats_false_vs_true_merged.tsv.gz"), sep = "\t")

LDSCoutput <- readRDS(ldsc_file)
model <- "
F1 =~ aALPS + Left_ALPS + Right_ALPS
F2 =~ aALPS + mALPS + pALPS
F1 ~~ F2
Left_ALPS ~~ lvar*Left_ALPS
mALPS ~~ mvar*mALPS
lvar > 0.0001
mvar > 0.0001
F1 ~ SNP
F2 ~ SNP
"

pick_chunk <- function(ss_obj, n = 5000L) {
  dt <- as.data.table(ss_obj)
  dt[seq_len(min(nrow(dt), n))]
}

message("Running mini userGWAS with OLS=FALSE input...")
gw_false <- userGWAS(
  covstruc = LDSCoutput,
  SNPs = pick_chunk(ss_false),
  estimation = "DWLS",
  model = model,
  printwarn = TRUE,
  sub = c("F1~SNP", "F2~SNP"),
  cores = 1,
  toler = 1e-50,
  parallel = FALSE,
  GC = "none",
  MPI = FALSE,
  smooth_check = TRUE,
  fix_measurement = TRUE,
  Q_SNP = FALSE
)

message("Running mini userGWAS with OLS=TRUE input...")
gw_true <- userGWAS(
  covstruc = LDSCoutput,
  SNPs = pick_chunk(ss_true),
  estimation = "DWLS",
  model = model,
  printwarn = TRUE,
  sub = c("F1~SNP", "F2~SNP"),
  cores = 1,
  toler = 1e-50,
  parallel = FALSE,
  GC = "none",
  MPI = FALSE,
  smooth_check = TRUE,
  fix_measurement = TRUE,
  Q_SNP = FALSE
)

f1_false <- as.data.table(gw_false[[1]])
f2_false <- as.data.table(gw_false[[2]])
f1_true <- as.data.table(gw_true[[1]])
f2_true <- as.data.table(gw_true[[2]])

fwrite(f1_false, file.path(work_dir, "mini_F1_ols_false.tsv"), sep = "\t")
fwrite(f2_false, file.path(work_dir, "mini_F2_ols_false.tsv"), sep = "\t")
fwrite(f1_true, file.path(work_dir, "mini_F1_ols_true.tsv"), sep = "\t")
fwrite(f2_true, file.path(work_dir, "mini_F2_ols_true.tsv"), sep = "\t")

compare_one <- function(dt_false, dt_true, factor_name) {
  m <- merge(
    dt_false[, .(SNP, est_false = est, se_false = SE, p_false = Pval_Estimate)],
    dt_true[, .(SNP, est_true = est, se_true = SE, p_true = Pval_Estimate)],
    by = "SNP"
  )
  data.table(
    factor = factor_name,
    n_snps = nrow(m),
    cor_est = cor(m$est_false, m$est_true, use = "pairwise.complete.obs"),
    cor_se = cor(m$se_false, m$se_true, use = "pairwise.complete.obs"),
    cor_logp = cor(-log10(m$p_false), -log10(m$p_true), use = "pairwise.complete.obs"),
    mean_abs_delta_est = mean(abs(m$est_false - m$est_true), na.rm = TRUE),
    median_abs_delta_est = median(abs(m$est_false - m$est_true), na.rm = TRUE),
    mean_abs_ratio_se = mean(abs(m$se_false / m$se_true - 1), na.rm = TRUE)
  )
}

summary_dt <- rbindlist(list(
  compare_one(f1_false, f1_true, "F1"),
  compare_one(f2_false, f2_true, "F2")
))

fwrite(summary_dt, file.path(work_dir, "mini_ols_sensitivity_summary.tsv"), sep = "\t")

