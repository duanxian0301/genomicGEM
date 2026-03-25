library(GenomicSEM)
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

# Best-fitting unconstrained 2-factor model from the EFA/model-search stage.
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

fit_usermodel_cfa <- function(covstruc, label) {
  fit <- usermodel(
    covstruc = covstruc,
    estimation = "DWLS",
    model = cfa_model,
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = FALSE
  )

  saveRDS(fit, file.path(output_dir, paste0(label, "_fit.rds")))
  fwrite(as.data.frame(fit$results), file.path(output_dir, paste0(label, "_results.tsv")), sep = "\t")

  mf <- fit$modelfit
  fit_summary <- data.frame(
    analysis = label,
    chisq = if ("chisq" %in% names(mf)) as.numeric(mf[["chisq"]]) else NA_real_,
    df = if ("df" %in% names(mf)) as.numeric(mf[["df"]]) else NA_real_,
    p_chisq = if ("p_chisq" %in% names(mf)) as.numeric(mf[["p_chisq"]]) else NA_real_,
    cfi = if ("CFI" %in% names(mf)) as.numeric(mf[["CFI"]]) else NA_real_,
    srmr = if ("SRMR" %in% names(mf)) as.numeric(mf[["SRMR"]]) else NA_real_
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

message("Running GenomicSEM usermodel CFA on EVEN chromosomes...")
even_out <- fit_usermodel_cfa(LDSCoutput_even, "ALPS5_EVEN_2factor_CFA")

message("Running GenomicSEM usermodel CFA on all autosomes...")
all_out <- fit_usermodel_cfa(LDSCoutput_all, "ALPS5_ALL_2factor_CFA")

summary_table <- rbind(even_out$fit_summary, all_out$fit_summary)
fwrite(summary_table, file.path(output_dir, "ALPS5_2factor_CFA_fit_summary.tsv"), sep = "\t")

message("Step 3 CFA finished successfully.")
