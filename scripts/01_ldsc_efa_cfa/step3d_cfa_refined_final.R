library(lavaan)
library(data.table)

output_dir <- "D:/文章/GS/GWAS/step3_cfa_results"
final_dir <- file.path(output_dir, "refined_final")
dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)

even_ldsc <- readRDS(file.path(output_dir, "ALPS5_EVEN_CFA_ldsc.rds"))
all_ldsc <- readRDS(file.path(output_dir, "ALPS5_ALL_CFA_ldsc.rds"))

get_S <- function(covstruc) {
  S <- as.matrix(covstruc[[2]])
  rownames(S) <- colnames(S)
  S
}

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

fit_refined <- function(S, label) {
  fit <- sem(
    model = refined_model,
    sample.cov = S,
    estimator = "ML",
    sample.nobs = 200,
    std.lv = TRUE,
    sample.cov.rescale = FALSE
  )

  pe <- parameterEstimates(fit, standardized = TRUE)
  fwrite(pe, file.path(final_dir, paste0(label, "_parameters.tsv")), sep = "\t")

  theta <- lavInspect(fit, "theta")
  write.csv(theta, file.path(final_dir, paste0(label, "_theta.csv")), row.names = TRUE)

  fit_table <- data.frame(
    analysis = label,
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    aic = fitMeasures(fit, "aic")
  )
  fwrite(fit_table, file.path(final_dir, paste0(label, "_fit.tsv")), sep = "\t")
  saveRDS(fit, file.path(final_dir, paste0(label, "_fit.rds")))

  fit_table
}

even_fit <- fit_refined(get_S(even_ldsc), "EVEN_refined")
all_fit <- fit_refined(get_S(all_ldsc), "ALL_refined")

summary_table <- rbind(even_fit, all_fit)
fwrite(summary_table, file.path(final_dir, "ALPS_refined_CFA_fit_summary.tsv"), sep = "\t")
