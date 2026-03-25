library(lavaan)
library(data.table)

output_dir <- "D:/文章/GS/GWAS/step3_cfa_results"
compare_dir <- file.path(output_dir, "final_model_comparison")
dir.create(compare_dir, showWarnings = FALSE, recursive = TRUE)

even_ldsc <- readRDS(file.path(output_dir, "ALPS5_EVEN_CFA_ldsc.rds"))
all_ldsc <- readRDS(file.path(output_dir, "ALPS5_ALL_CFA_ldsc.rds"))
search_summary <- fread(file.path(output_dir, "model_search", "model_search_summary.tsv"))
refined_summary <- fread(file.path(output_dir, "refined_final", "ALPS_refined_CFA_fit_summary.tsv"))

get_S <- function(covstruc) {
  S <- as.matrix(covstruc[[2]])
  rownames(S) <- colnames(S)
  S
}

bifactor_model <- "
G =~ aALPS + Left_ALPS + mALPS + pALPS + Right_ALPS
F1 =~ Left_ALPS + Right_ALPS
F2 =~ mALPS + pALPS
G ~~ 0*F1
G ~~ 0*F2
F1 ~~ 0*F2
"

fit_bifactor <- function(S, dataset_label) {
  fit <- tryCatch(
    sem(
      model = bifactor_model,
      sample.cov = S,
      estimator = "ML",
      sample.nobs = 200,
      std.lv = TRUE,
      sample.cov.rescale = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      dataset = dataset_label,
      model = "bifactor",
      converged = FALSE,
      chisq = NA_real_,
      df = NA_real_,
      p_chisq = NA_real_,
      cfi = NA_real_,
      srmr = NA_real_,
      aic = NA_real_,
      neg_resid_count = NA_integer_,
      min_resid = NA_real_,
      status = paste("error:", conditionMessage(fit)),
      stringsAsFactors = FALSE
    ))
  }

  pe <- parameterEstimates(fit, standardized = TRUE)
  fwrite(pe, file.path(compare_dir, paste0(dataset_label, "__bifactor__parameters.tsv")), sep = "\t")
  saveRDS(fit, file.path(compare_dir, paste0(dataset_label, "__bifactor__fit.rds")))

  resid_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs & !(pe$lhs %in% c("G", "F1", "F2")), ]

  data.frame(
    dataset = dataset_label,
    model = "bifactor",
    converged = lavInspect(fit, "converged"),
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    aic = fitMeasures(fit, "aic"),
    neg_resid_count = sum(resid_rows$est < 0, na.rm = TRUE),
    min_resid = min(resid_rows$est, na.rm = TRUE),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

bifactor_even <- fit_bifactor(get_S(even_ldsc), "EVEN")
bifactor_all <- fit_bifactor(get_S(all_ldsc), "ALL")
bifactor_summary <- rbindlist(list(bifactor_even, bifactor_all), fill = TRUE)
fwrite(bifactor_summary, file.path(compare_dir, "bifactor_summary.tsv"), sep = "\t")

selected_search <- search_summary[model %in% c("one_factor", "two_factor_cross")]
selected_search[, model_label := fifelse(model == "one_factor", "common_factor",
                                  fifelse(model == "two_factor_cross", "two_factor", model))]

selected_refined <- copy(refined_summary)
selected_refined[, dataset := sub("_refined$", "", analysis)]
selected_refined[, model_label := "refined_two_factor"]
selected_refined[, `:=`(
  chisq = as.numeric(chisq),
  df = as.numeric(df),
  p_chisq = as.numeric(p_chisq),
  cfi = as.numeric(cfi),
  srmr = as.numeric(srmr),
  aic = as.numeric(aic),
  converged = TRUE,
  neg_resid_count = NA_integer_,
  min_resid = NA_real_,
  status = "ok"
)]

selected_bifactor <- copy(bifactor_summary)
selected_bifactor[, model_label := "bifactor"]

comparison <- rbindlist(list(
  selected_search[, .(dataset, model_label, converged, chisq, df, p_chisq, cfi, srmr, aic, neg_resid_count, min_resid, status)],
  selected_refined[, .(dataset, model_label, converged, chisq, df, p_chisq, cfi, srmr, aic, neg_resid_count, min_resid, status)],
  selected_bifactor[, .(dataset, model_label, converged, chisq, df, p_chisq, cfi, srmr, aic, neg_resid_count, min_resid, status)]
), fill = TRUE, ignore.attr = TRUE)

setorder(comparison, dataset, model_label)
fwrite(comparison, file.path(compare_dir, "ALPS_CFA_final_model_comparison.tsv"), sep = "\t")

notes <- c(
  "Models compared:",
  "common_factor = one-factor model across all 5 ALPS indicators",
  "two_factor = cross-loading two-factor model from EFA",
  "refined_two_factor = final constrained model used for factor GWAS",
  "bifactor = general factor plus Left/Right and m/p group factors"
)
writeLines(notes, file.path(compare_dir, "ALPS_CFA_final_model_comparison_notes.txt"))
