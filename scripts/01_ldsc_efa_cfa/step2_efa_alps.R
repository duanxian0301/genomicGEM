library(GenomicSEM)
library(Matrix)
library(data.table)

gwas_dir <- "D:/文章/GS/GWAS"
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
output_dir <- file.path(gwas_dir, "step2_efa_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(gwas_dir)

trait_names <- c("aALPS", "Left_ALPS", "mALPS", "pALPS", "Right_ALPS")
traits <- file.path(gwas_dir, paste0(trait_names, ".sumstats.gz"))

if (!all(file.exists(traits))) {
  stop("Missing .sumstats.gz files for one or more traits.")
}

write.model <- function(Loadings, S_LD, cutoff, fix_resid = TRUE,
                        bifactor = FALSE, mustload = FALSE, common = FALSE) {
  Model <- ""
  if (common == TRUE) {
    for (f in 1) {
      u <- 1
      Model1 <- ""
      for (i in 1:nrow(S_LD)) {
        if (u == 1) {
          linestart <- paste("F", f, "=~", colnames(S_LD)[i], sep = "")
          u <- u + 1
          linemid <- ""
        } else {
          linemid <- paste(linemid, " + ", colnames(S_LD)[i], sep = "")
        }
      }
    }
    Model <- paste(Model, linestart, linemid, " \n ", sep = "")
  } else {
    if (mustload == TRUE) {
      Mins <- apply(abs(Loadings), 1, max)
      for (i in 1:nrow(Loadings)) {
        for (f in 1:ncol(Loadings)) {
          if (Mins[i] == abs(Loadings[i, f])) {
            Loadings[i, f] <- sign(Loadings[i, f]) * (cutoff + .01)
          }
        }
      }
    }

    for (f in 1:ncol(Loadings)) {
      u <- 1
      linestart <- ""
      linemid <- ""
      for (i in 1:nrow(Loadings)) {
        if (abs(Loadings[i, f]) > cutoff) {
          if (u == 1) {
            linestart <- paste("F", f, "=~", colnames(S_LD)[i], sep = "")
            u <- u + 1
          } else {
            linemid <- paste(linemid, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      if (linestart != "") {
        Model <- paste(Model, linestart, linemid, " \n ", sep = "")
      }
    }

    if (bifactor == TRUE) {
      Model_bi <- ""
      u <- 1
      for (i in 1:ncol(S_LD)) {
        b <- grepl(colnames(S_LD)[i], Model, fixed = TRUE)
        if (b == TRUE) {
          if (u == 1) {
            linestart_bi <- paste("Common_F", "=~", colnames(S_LD)[i], sep = "")
            u <- u + 1
            linemid_bi <- ""
          } else {
            linemid_bi <- paste(linemid_bi, " + ", colnames(S_LD)[i], sep = "")
          }
        }
      }
      Model_bi <- paste(linestart_bi, linemid_bi, " \n ", sep = "")
      Factor_bi <- ""
      for (i in 1:ncol(Loadings)) {
        Factor_bi <- paste(
          Factor_bi,
          paste("Common_F~~0*F", i, " \n ", sep = ""),
          sep = ""
        )
      }

      Modelsat <- ""
      for (i in 1:(ncol(Loadings) - 1)) {
        Modelsat <- paste(Modelsat, paste("", "F", i, "~~0*F", i + 1, sep = ""), " \n ", sep = "")
        if (ncol(Loadings) - i >= 2) {
          for (j in (i + 2):ncol(Loadings)) {
            Modelsat <- paste(Modelsat, paste("", "F", i, "~~0*F", j, " \n ", sep = ""), sep = "")
          }
        }
      }
      Model <- paste(Model, Model_bi, Factor_bi, Modelsat)
    }
  }

  if (fix_resid == TRUE) {
    Model3 <- ""
    label_pool <- combn(letters, 4)
    if (ncol(S_LD) > ncol(label_pool)) {
      stop("Not enough residual labels available.")
    }
    for (i in 1:ncol(S_LD)) {
      if (grepl(colnames(S_LD)[i], Model, fixed = TRUE) == TRUE) {
        label <- paste(label_pool[, i], collapse = "")
        Model3 <- paste(
          Model3,
          paste(colnames(S_LD)[i], " ~~ ", label, "*", colnames(S_LD)[i], sep = ""),
          " \n ",
          paste(label, " > .0001", sep = ""),
          " \n ",
          sep = ""
        )
      }
    }
    Model <- paste(Model, Model3, sep = "")
  }

  Model
}

calc_factor_criteria <- function(S_LD) {
  eig <- eigen(cov2cor(S_LD), only.values = TRUE)$values
  nk <- length(eig)
  k <- 1:nk
  criteria <- mean(eig)
  nkaiser <- sum(eig >= rep(criteria, nk))
  aparallel <- rep(criteria, length(eig))
  af <- rep(NA_real_, nk)
  for (j in 2:(nk - 1)) {
    if (eig[j - 1] >= aparallel[j - 1]) {
      af[j] <- (eig[j + 1] - 2 * eig[j]) + eig[j - 1]
    }
  }
  naf <- which(af == max(af, na.rm = TRUE), TRUE) - 1

  pred.eig <- rep(NA_real_, nk)
  proportion <- eig / sum(eig)
  cond1 <- TRUE
  cond2 <- TRUE
  i <- 0
  noc <- NA_integer_
  while ((cond1 == TRUE) && (cond2 == TRUE) && (i < nk)) {
    i <- i + 1
    ind <- k[c(i + 1, nk)]
    vp.p <- lm(eig[c(i + 1, nk)] ~ ind)
    vp.prec <- pred.eig[i] <- sum(c(1, i) * coef(vp.p))
    cond1 <- (eig[i] >= vp.prec)
    cond2 <- (eig[i] >= aparallel[i])
    noc <- i - 1
    if (i == nk - 1) {
      break
    }
  }

  data.frame(
    criterion = c("eigenvalues", "kaiser_n", "acceleration_factor_n", "optimal_coordinates_n"),
    value = c(
      paste(round(eig, 6), collapse = "; "),
      nkaiser,
      ifelse(length(naf) == 0, NA, naf[1]),
      noc
    )
  )
}

message("Running ODD-chromosome LDSC for EFA...")
LDSCoutput_odd <- ldsc(
  traits = traits,
  sample.prev = rep(NA, length(traits)),
  population.prev = rep(NA, length(traits)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = trait_names,
  ldsc.log = file.path(output_dir, "ALPS5_ODD_EFA"),
  select = "ODD"
)

saveRDS(LDSCoutput_odd, file.path(output_dir, "ALPS5_ODD_EFA.rds"))

S_LD <- as.matrix((nearPD(LDSCoutput_odd$S, corr = FALSE))$mat)
write.csv(S_LD, file.path(output_dir, "ALPS5_ODD_S_matrix_smoothed.csv"), row.names = TRUE)
write.csv(cov2cor(S_LD), file.path(output_dir, "ALPS5_ODD_rg_matrix_smoothed.csv"), row.names = TRUE)

criteria_table <- calc_factor_criteria(S_LD)
fwrite(criteria_table, file.path(output_dir, "ALPS5_EFA_factor_criteria.tsv"), sep = "\t")

n_traits <- ncol(S_LD)
max_factors <- min(4, n_traits - 1)
fit_summary <- list()

for (nf in 1:max_factors) {
  message("Fitting EFA with ", nf, " factor(s)...")
  efa_fit <- tryCatch(
    factanal(
      factors = nf,
      covmat = S_LD,
      n.obs = 10000,
      rotation = if (nf == 1) "none" else "promax"
    ),
    error = function(e) e
  )

  if (inherits(efa_fit, "error")) {
    fit_summary[[nf]] <- data.frame(
      factors = nf,
      dof = NA,
      statistic = NA,
      pvalue = NA,
      status = paste("failed:", conditionMessage(efa_fit))
    )
    next
  }

  loadings_mat <- as.matrix(efa_fit$loadings[1:ncol(S_LD), 1:nf, drop = FALSE])
  rownames(loadings_mat) <- colnames(S_LD)
  colnames(loadings_mat) <- paste0("F", seq_len(nf))

  loadings_df <- data.frame(trait = rownames(loadings_mat), loadings_mat, check.names = FALSE)
  fwrite(loadings_df, file.path(output_dir, paste0("EFA_", nf, "factor_loadings.tsv")), sep = "\t")

  uniqueness_df <- data.frame(
    trait = names(efa_fit$uniquenesses),
    uniqueness = unname(efa_fit$uniquenesses)
  )
  fwrite(uniqueness_df, file.path(output_dir, paste0("EFA_", nf, "factor_uniqueness.tsv")), sep = "\t")

  model_text <- write.model(loadings_mat, S_LD, cutoff = 0.30, mustload = TRUE)
  writeLines(model_text, file.path(output_dir, paste0("EFA_", nf, "factor_model.txt")))

  fit_summary[[nf]] <- data.frame(
    factors = nf,
    dof = efa_fit$dof,
    statistic = efa_fit$STATISTIC,
    pvalue = efa_fit$PVAL,
    status = "ok"
  )
}

fit_summary_df <- rbindlist(fit_summary)
fwrite(fit_summary_df, file.path(output_dir, "EFA_fit_summary.tsv"), sep = "\t")

message("Step 2 EFA finished successfully.")
