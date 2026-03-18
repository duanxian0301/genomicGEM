library(GenomicSEM)
library(lavaan)
library(data.table)

gwas_dir <- "D:/文章/GS/GWAS"
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "step3_cfa_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(gwas_dir)

trait_names <- c("aALPS", "Left_ALPS", "mALPS", "pALPS", "Right_ALPS")
traits <- file.path(gwas_dir, paste0(trait_names, ".sumstats.gz"))

if (!all(file.exists(traits))) {
  stop("Missing .sumstats.gz files for one or more traits.")
}

# Stable 2-factor CFA model derived from the EFA step.
cfa_model <- "
F1 =~ aALPS + Left_ALPS + Right_ALPS
F2 =~ aALPS + mALPS + pALPS
F1 ~~ F2
"
writeLines(cfa_model, file.path(output_dir, "ALPS_2factor_CFA_model.txt"))

run_ldsc <- function(select_mode = NULL, log_name) {
  args <- list(
    traits = traits,
    sample.prev = rep(NA, length(traits)),
    population.prev = rep(NA, length(traits)),
    ld = ld_ref_dir,
    wld = ld_ref_dir,
    trait.names = trait_names,
    ldsc.log = log_name
  )
  if (!is.null(select_mode)) {
    args$select <- select_mode
  }
  do.call(ldsc, args)
}

fit_lavaan_cfa <- function(covstruc, label) {
  S <- as.matrix(covstruc[[2]])
  rownames(S) <- colnames(S)

  fit <- sem(
    model = cfa_model,
    sample.cov = S,
    estimator = "ML",
    sample.nobs = 200,
    std.lv = TRUE,
    sample.cov.rescale = FALSE
  )

  saveRDS(fit, file.path(output_dir, paste0(label, "_fit.rds")))

  pe <- parameterEstimates(fit, standardized = FALSE)
  ss <- standardizedSolution(fit)
  merge_key <- c("lhs", "op", "rhs")
  pe2 <- merge(pe, ss[, c(merge_key, "est.std")], by = merge_key, all.x = TRUE)
  fwrite(pe2, file.path(output_dir, paste0(label, "_parameters.tsv")), sep = "\t")

  fit_summary <- data.frame(
    analysis = label,
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    aic = fitMeasures(fit, "aic")
  )
  fwrite(fit_summary, file.path(output_dir, paste0(label, "_fit.tsv")), sep = "\t")

  list(fit = fit, fit_summary = fit_summary)
}

message("Running EVEN-chromosome LDSC for CFA validation...")
LDSCoutput_even <- run_ldsc(
  select_mode = "EVEN",
  log_name = file.path(output_dir, "ALPS5_EVEN_CFA")
)
saveRDS(LDSCoutput_even, file.path(output_dir, "ALPS5_EVEN_CFA_ldsc.rds"))

message("Running genome-wide LDSC for final CFA...")
LDSCoutput_all <- run_ldsc(
  select_mode = NULL,
  log_name = file.path(output_dir, "ALPS5_ALL_CFA")
)
saveRDS(LDSCoutput_all, file.path(output_dir, "ALPS5_ALL_CFA_ldsc.rds"))

message("Running lavaan ML CFA on EVEN chromosomes...")
even_out <- fit_lavaan_cfa(LDSCoutput_even, "ALPS5_EVEN_2factor_CFA")

message("Running lavaan ML CFA on all autosomes...")
all_out <- fit_lavaan_cfa(LDSCoutput_all, "ALPS5_ALL_2factor_CFA")

summary_table <- rbind(even_out$fit_summary, all_out$fit_summary)
fwrite(summary_table, file.path(output_dir, "ALPS5_2factor_CFA_fit_summary.tsv"), sep = "\t")

message("Step 3 CFA finished successfully.")
