library(lavaan)
library(data.table)

output_dir <- "D:/文章/GS/GWAS/step3_cfa_results"
search_dir <- file.path(output_dir, "model_search")
dir.create(search_dir, showWarnings = FALSE, recursive = TRUE)

even_ldsc <- readRDS(file.path(output_dir, "ALPS5_EVEN_CFA_ldsc.rds"))
all_ldsc <- readRDS(file.path(output_dir, "ALPS5_ALL_CFA_ldsc.rds"))

get_S <- function(covstruc) {
  S <- as.matrix(covstruc[[2]])
  rownames(S) <- colnames(S)
  S
}

fit_model <- function(model_name, model_syntax, S, dataset_label) {
  fit <- tryCatch(
    sem(
      model = model_syntax,
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
      model = model_name,
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

  pe <- parameterEstimates(fit)
  resid_rows <- pe[pe$op == "~~" & pe$lhs == pe$rhs & !grepl("^F[0-9]+$", pe$lhs), ]
  neg_resid_count <- sum(resid_rows$est < 0, na.rm = TRUE)
  min_resid <- min(resid_rows$est, na.rm = TRUE)

  params <- parameterEstimates(fit, standardized = TRUE)
  fwrite(
    params,
    file.path(search_dir, paste0(dataset_label, "__", model_name, "__parameters.tsv")),
    sep = "\t"
  )

  data.frame(
    dataset = dataset_label,
    model = model_name,
    converged = lavInspect(fit, "converged"),
    chisq = fitMeasures(fit, "chisq"),
    df = fitMeasures(fit, "df"),
    p_chisq = fitMeasures(fit, "pvalue"),
    cfi = fitMeasures(fit, "cfi"),
    srmr = fitMeasures(fit, "srmr"),
    aic = fitMeasures(fit, "aic"),
    neg_resid_count = neg_resid_count,
    min_resid = min_resid,
    status = "ok",
    stringsAsFactors = FALSE
  )
}

models <- list(
  one_factor = "
    F1 =~ aALPS + Left_ALPS + mALPS + pALPS + Right_ALPS
  ",
  two_factor_cross = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
  ",
  two_factor_a_on_f1 = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ mALPS + pALPS
    F1 ~~ F2
  ",
  two_factor_a_on_f2 = "
    F1 =~ Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
  ",
  two_factor_with_mp_resid = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
    mALPS ~~ pALPS
  ",
  two_factor_with_lr_resid = "
    F1 =~ aALPS + Left_ALPS + Right_ALPS
    F2 =~ aALPS + mALPS + pALPS
    F1 ~~ F2
    Left_ALPS ~~ Right_ALPS
  "
)

results <- list()
for (dataset_label in c("EVEN", "ALL")) {
  S <- if (dataset_label == "EVEN") get_S(even_ldsc) else get_S(all_ldsc)
  for (model_name in names(models)) {
    results[[paste(dataset_label, model_name, sep = "_")]] <-
      fit_model(model_name, models[[model_name]], S, dataset_label)
  }
}

summary_df <- rbindlist(results, fill = TRUE)
fwrite(summary_df, file.path(search_dir, "model_search_summary.tsv"), sep = "\t")

model_lines <- unlist(
  lapply(names(models), function(nm) c(paste0("### ", nm), models[[nm]], "")),
  use.names = FALSE
)
writeLines(model_lines, con = file.path(search_dir, "candidate_models.txt"))
