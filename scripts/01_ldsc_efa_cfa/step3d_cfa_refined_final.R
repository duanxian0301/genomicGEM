library(GenomicSEM)
library(data.table)

output_dir <- "D:/文章/GS/GWAS/step3_cfa_results"
final_dir <- file.path(output_dir, "refined_final")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

even_ldsc <- readRDS(file.path(output_dir, "ALPS5_EVEN_CFA_ldsc.rds"))
all_ldsc <- readRDS(file.path(output_dir, "ALPS5_ALL_CFA_ldsc.rds"))

refined_model <- "
F1 =~ aALPS + Left_ALPS + Right_ALPS
F2 =~ aALPS + mALPS + pALPS
F1 ~~ F2
Left_ALPS ~~ lvar*Left_ALPS
mALPS ~~ mvar*mALPS
lvar > 0.0001
mvar > 0.0001
"

writeLines(refined_model, file.path(final_dir, "ALPS_2factor_refined_model.txt"))

fit_refined <- function(covstruc, label) {
  fit <- usermodel(
    covstruc = covstruc,
    estimation = "DWLS",
    model = refined_model,
    CFIcalc = TRUE,
    std.lv = TRUE,
    imp_cov = FALSE
  )

  fwrite(as.data.frame(fit$results), file.path(final_dir, paste0(label, "_results.tsv")), sep = "\t")
  saveRDS(fit, file.path(final_dir, paste0(label, "_fit.rds")))

  mf <- fit$modelfit
  fit_table <- data.frame(
    analysis = label,
    chisq = if ("chisq" %in% names(mf)) as.numeric(mf[["chisq"]]) else NA_real_,
    df = if ("df" %in% names(mf)) as.numeric(mf[["df"]]) else NA_real_,
    p_chisq = if ("p_chisq" %in% names(mf)) as.numeric(mf[["p_chisq"]]) else NA_real_,
    cfi = if ("CFI" %in% names(mf)) as.numeric(mf[["CFI"]]) else NA_real_,
    srmr = if ("SRMR" %in% names(mf)) as.numeric(mf[["SRMR"]]) else NA_real_
  )
  fwrite(fit_table, file.path(final_dir, paste0(label, "_fit.tsv")), sep = "\t")

  fit_table
}

even_fit <- fit_refined(even_ldsc, "EVEN_refined")
all_fit <- fit_refined(all_ldsc, "ALL_refined")

summary_table <- rbind(even_fit, all_fit)
fwrite(summary_table, file.path(final_dir, "ALPS_refined_CFA_fit_summary.tsv"), sep = "\t")
