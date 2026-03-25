.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/home/shenjing/gs_paths/gwas"
output_dir <- file.path(root_dir, "step3_cfa_results", "native_usermodel")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

even_ldsc <- readRDS(file.path(root_dir, "step3_cfa_results", "ALPS5_EVEN_CFA_ldsc.rds"))
all_ldsc <- readRDS(file.path(root_dir, "step3_cfa_results", "ALPS5_ALL_CFA_ldsc.rds"))

models <- list(
  common_factor = "
    F1 =~ aALPS + Left_ALPS + mALPS + pALPS + Right_ALPS
  ",
  two_factor = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
  ",
  refined_two_factor = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
    Left_ALPS ~~ lvar*Left_ALPS
    mALPS ~~ mvar*mALPS
    lvar > .0001
    mvar > .0001
  ",
  bifactor = "
    G =~ aALPS + Left_ALPS + mALPS + pALPS + Right_ALPS
    F1 =~ Left_ALPS + Right_ALPS
    F2 =~ mALPS + pALPS
    G ~~ 0*F1
    G ~~ 0*F2
    F1 ~~ 0*F2
  "
)

fit_one <- function(covstruc, dataset_label, model_name, model_syntax) {
  fit <- tryCatch(
    usermodel(
      covstruc = covstruc,
      estimation = "DWLS",
      model = model_syntax,
      CFIcalc = TRUE,
      std.lv = TRUE,
      imp_cov = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      dataset = dataset_label,
      model = model_name,
      converged = FALSE,
      chisq = NA_real_,
      df = NA_real_,
      p_chisq = NA_real_,
      cfi = NA_real_,
      srmr = NA_real_,
      aic = NA_real_,
      status = paste("error:", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  saveRDS(fit, file.path(output_dir, paste0(dataset_label, "__", model_name, "__fit.rds")))
  fwrite(as.data.frame(fit$results), file.path(output_dir, paste0(dataset_label, "__", model_name, "__results.tsv")), sep = "\t")

  summ <- fit$modelfit
  data.frame(
    dataset = dataset_label,
    model = model_name,
    converged = TRUE,
    chisq = if ("chisq" %in% names(summ)) as.numeric(summ[["chisq"]]) else NA_real_,
    df = if ("df" %in% names(summ)) as.numeric(summ[["df"]]) else NA_real_,
    p_chisq = if ("p_chisq" %in% names(summ)) as.numeric(summ[["p_chisq"]]) else NA_real_,
    cfi = if ("CFI" %in% names(summ)) as.numeric(summ[["CFI"]]) else NA_real_,
    srmr = if ("SRMR" %in% names(summ)) as.numeric(summ[["SRMR"]]) else NA_real_,
    aic = NA_real_,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

out <- list()
for (dataset_label in c("EVEN", "ALL")) {
  covstruc <- if (dataset_label == "EVEN") even_ldsc else all_ldsc
  for (model_name in names(models)) {
    message("Running ", dataset_label, " / ", model_name)
    out[[paste(dataset_label, model_name, sep = "_")]] <- fit_one(covstruc, dataset_label, model_name, models[[model_name]])
  }
}

summary_df <- rbindlist(out, fill = TRUE)
fwrite(summary_df, file.path(output_dir, "ALPS_CFA_usermodel_comparison.tsv"), sep = "\t")

