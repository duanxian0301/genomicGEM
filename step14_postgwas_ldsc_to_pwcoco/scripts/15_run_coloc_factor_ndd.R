#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

root_out <- "D:/文章/GS/postgwas/07_coloc_factor_ndd"
dir.create(root_out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_out, "per_task"), recursive = TRUE, showWarnings = FALSE)

tasks_file <- "D:/文章/GS/postgwas/06_coloc_pwcoco_candidate_loci/conjfdr_factor_loci_region_pair_tasks.tsv"
f1_file <- "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F1_factorGWAS_native_standard.txt"
f2_file <- "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F2_factorGWAS_native_standard.txt"
ad_file <- "D:/文章/4NDD/NDDGWAS/AD.txt"
pd_file <- "D:/文章/4NDD/NDDGWAS/PD.txt"

trait_files <- list(F1 = f1_file, F2 = f2_file)
disease_files <- list(AD = ad_file, PD = pd_file)

# Disease case fractions follow the official published sample descriptions,
# using proxy-inclusive primary definitions where applicable.
disease_meta <- data.table(
  disease = c("AD", "PD"),
  disease_type = c("cc", "cc"),
  n_cases = c(39106 + 46828, 63555 + 17700),
  n_total = c(487511, 1835938)
)
disease_meta[, s_case_prop := n_cases / n_total]

required_cols <- c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N")
tasks <- fread(tasks_file)

comp_base <- c(A = "T", T = "A", C = "G", G = "C")
comp_vec <- function(x) {
  y <- comp_base[toupper(x)]
  y[is.na(y)] <- NA_character_
  unname(y)
}

is_palindromic <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

read_sumstats <- function(path) {
  dt <- fread(path, select = required_cols, showProgress = FALSE)
  dt[, `:=`(
    CHR = as.integer(CHR),
    BP = as.integer(BP),
    FREQ = as.numeric(FREQ),
    BETA = as.numeric(BETA),
    SE = as.numeric(SE),
    P = as.numeric(P),
    N = as.numeric(N),
    A1 = toupper(A1),
    A2 = toupper(A2)
  )]
  dt
}

harmonize_pair <- function(dt1, dt2) {
  merged <- merge(
    dt1, dt2,
    by = "SNP",
    suffixes = c(".1", ".2"),
    allow.cartesian = FALSE
  )
  if (!nrow(merged)) {
    return(list(ok = FALSE, data = NULL, reason = "No overlapping SNPs after merge"))
  }

  merged <- merged[CHR.1 == CHR.2 & BP.1 == BP.2]
  if (!nrow(merged)) {
    return(list(ok = FALSE, data = NULL, reason = "No overlapping SNPs with matched CHR/BP"))
  }

  merged[, `:=`(
    A1.2.comp = comp_vec(A1.2),
    A2.2.comp = comp_vec(A2.2)
  )]

  same <- merged$A1.1 == merged$A1.2 & merged$A2.1 == merged$A2.2
  swap <- merged$A1.1 == merged$A2.2 & merged$A2.1 == merged$A1.2
  comp_same <- merged$A1.1 == merged$A1.2.comp & merged$A2.1 == merged$A2.2.comp
  comp_swap <- merged$A1.1 == merged$A2.2.comp & merged$A2.1 == merged$A1.2.comp

  merged[, align_status := fifelse(
    same, "same",
    fifelse(swap, "swap",
      fifelse(comp_same, "comp_same",
        fifelse(comp_swap, "comp_swap", "drop")
      )
    )
  )]

  out <- merged[align_status != "drop"]
  if (!nrow(out)) {
    return(list(ok = FALSE, data = NULL, reason = "No alignable SNPs after allele harmonization"))
  }

  out[align_status %in% c("swap", "comp_swap"), `:=`(
    BETA.2 = -BETA.2,
    FREQ.2 = 1 - FREQ.2
  )]
  out[align_status %in% c("comp_same", "comp_swap"), `:=`(
    A1.2 = A1.1,
    A2.2 = A2.1
  )]

  out <- out[
    !is.na(FREQ.1) & !is.na(FREQ.2) &
      FREQ.1 > 0 & FREQ.1 < 1 &
      FREQ.2 > 0 & FREQ.2 < 1 &
      !is.na(BETA.1) & !is.na(BETA.2) &
      !is.na(SE.1) & !is.na(SE.2) &
      SE.1 > 0 & SE.2 > 0
  ]

  if (!nrow(out)) {
    return(list(ok = FALSE, data = NULL, reason = "No usable SNPs after numeric filtering"))
  }

  out[, `:=`(
    MAF.1 = pmin(FREQ.1, 1 - FREQ.1),
    MAF.2 = pmin(FREQ.2, 1 - FREQ.2),
    palindromic = is_palindromic(A1.1, A2.1)
  )]

  out <- out[MAF.1 > 0 & MAF.2 > 0]
  if (!nrow(out)) {
    return(list(ok = FALSE, data = NULL, reason = "No usable SNPs after MAF filtering"))
  }

  setorder(out, CHR.1, BP.1, SNP)
  list(ok = TRUE, data = out, reason = NA_character_)
}

run_abf <- function(task, trait_dt, disease_dt, disease_info) {
  region_trait <- trait_dt[
    CHR == task$chrnum &
      BP >= task$region_start &
      BP <= task$region_end
  ]
  region_disease <- disease_dt[
    CHR == task$chrnum &
      BP >= task$region_start &
      BP <= task$region_end
  ]

  pre_trait_n <- nrow(region_trait)
  pre_disease_n <- nrow(region_disease)
  if (!pre_trait_n || !pre_disease_n) {
    return(list(
      ok = FALSE,
      summary = data.table(
        region_id = task$region_id,
        pair = task$pair,
        trait = task$trait,
        disease = task$disease,
        chrnum = task$chrnum,
        region_start = task$region_start,
        region_end = task$region_end,
        priority = task$priority,
        n_trait_region = pre_trait_n,
        n_disease_region = pre_disease_n,
        n_overlap_raw = 0L,
        n_aligned = 0L,
        status = "failed",
        failure_reason = "No SNPs in one or both datasets for region"
      ),
      snp_results = NULL
    ))
  }

  harmonized <- harmonize_pair(region_trait, region_disease)
  if (!harmonized$ok) {
    return(list(
      ok = FALSE,
      summary = data.table(
        region_id = task$region_id,
        pair = task$pair,
        trait = task$trait,
        disease = task$disease,
        chrnum = task$chrnum,
        region_start = task$region_start,
        region_end = task$region_end,
        priority = task$priority,
        n_trait_region = pre_trait_n,
        n_disease_region = pre_disease_n,
        n_overlap_raw = length(intersect(region_trait$SNP, region_disease$SNP)),
        n_aligned = 0L,
        status = "failed",
        failure_reason = harmonized$reason
      ),
      snp_results = NULL
    ))
  }

  dat <- harmonized$data
  d1 <- list(
    snp = dat$SNP,
    beta = dat$BETA.1,
    varbeta = dat$SE.1^2,
    MAF = dat$MAF.1,
    N = unique(dat$N.1)[1],
    type = "quant"
  )
  d2 <- list(
    snp = dat$SNP,
    beta = dat$BETA.2,
    varbeta = dat$SE.2^2,
    MAF = dat$MAF.2,
    N = unique(dat$N.2)[1],
    type = disease_info$disease_type,
    s = disease_info$s_case_prop
  )

  check_dataset(d1)
  check_dataset(d2)
  res <- coloc.abf(d1, d2)
  sumres <- as.list(res$summary)

  top_idx <- which.max(res$results$SNP.PP.H4)
  top_snp <- if (length(top_idx) && !is.na(top_idx)) res$results$snp[top_idx] else NA_character_
  top_snp_pph4 <- if (length(top_idx) && !is.na(top_idx)) res$results$SNP.PP.H4[top_idx] else NA_real_

  summary_dt <- data.table(
    region_id = task$region_id,
    pair = task$pair,
    trait = task$trait,
    disease = task$disease,
    chrnum = task$chrnum,
    region_start = task$region_start,
    region_end = task$region_end,
    priority = task$priority,
    sentinel_snps = task$sentinel_snps,
    n_trait_region = pre_trait_n,
    n_disease_region = pre_disease_n,
    n_overlap_raw = length(intersect(region_trait$SNP, region_disease$SNP)),
    n_aligned = nrow(dat),
    n_palindromic_retained = sum(dat$palindromic, na.rm = TRUE),
    trait_N = unique(dat$N.1)[1],
    disease_N = unique(dat$N.2)[1],
    disease_case_prop = disease_info$s_case_prop,
    PP.H0.abf = unname(sumres$PP.H0.abf),
    PP.H1.abf = unname(sumres$PP.H1.abf),
    PP.H2.abf = unname(sumres$PP.H2.abf),
    PP.H3.abf = unname(sumres$PP.H3.abf),
    PP.H4.abf = unname(sumres$PP.H4.abf),
    nsnps = unname(sumres$nsnps),
    top_snp_h4 = top_snp,
    top_snp_pph4 = top_snp_pph4,
    coloc_interpretation = fifelse(
      unname(sumres$PP.H4.abf) >= 0.80, "Strong_H4",
      fifelse(
        unname(sumres$PP.H3.abf) >= 0.80, "Strong_H3",
        fifelse(
          unname(sumres$PP.H4.abf) >= 0.50, "Moderate_H4",
          fifelse(unname(sumres$PP.H3.abf) >= 0.50, "Moderate_H3", "Inconclusive")
        )
      )
    ),
    status = "ok",
    failure_reason = NA_character_
  )

  snp_dt <- as.data.table(res$results)
  snp_dt[, `:=`(
    region_id = task$region_id,
    pair = task$pair,
    trait = task$trait,
    disease = task$disease
  )]
  setcolorder(snp_dt, c("region_id", "pair", "trait", "disease", setdiff(names(snp_dt), c("region_id", "pair", "trait", "disease"))))

  list(ok = TRUE, summary = summary_dt, snp_results = snp_dt)
}

pair_order <- unique(tasks[, .(pair, trait, disease)])
all_summaries <- list()
manifest <- list()

for (i in seq_len(nrow(pair_order))) {
  pair_row <- pair_order[i]
  pair_tasks <- tasks[pair == pair_row$pair]
  trait_name <- pair_row$trait
  disease_name <- pair_row$disease
  disease_info <- disease_meta[disease == disease_name]

  cat(sprintf("\n[%d/%d] Processing %s\n", i, nrow(pair_order), pair_row$pair))
  trait_dt <- read_sumstats(trait_files[[trait_name]])
  disease_dt <- read_sumstats(disease_files[[disease_name]])

  pair_summary <- list()
  for (j in seq_len(nrow(pair_tasks))) {
    task <- pair_tasks[j]
    cat(sprintf("  - %s %s chr%s:%s-%s\n", task$region_id, task$pair, task$chrnum, task$region_start, task$region_end))
    out <- run_abf(task, trait_dt, disease_dt, disease_info)
    pair_summary[[j]] <- out$summary

    per_prefix <- file.path(root_out, "per_task", sprintf("%s_%s", task$pair, task$region_id))
    fwrite(out$summary, paste0(per_prefix, "_summary.tsv"), sep = "\t")
    manifest[[length(manifest) + 1L]] <- data.table(
      pair = task$pair,
      region_id = task$region_id,
      summary_file = paste0(per_prefix, "_summary.tsv"),
      snp_file = if (isTRUE(out$ok)) paste0(per_prefix, "_snp_results.tsv.gz") else NA_character_
    )

    if (isTRUE(out$ok)) {
      fwrite(out$snp_results, paste0(per_prefix, "_snp_results.tsv.gz"), sep = "\t")
    }
  }

  all_summaries[[length(all_summaries) + 1L]] <- rbindlist(pair_summary, fill = TRUE)
  rm(trait_dt, disease_dt)
  invisible(gc())
}

summary_dt <- rbindlist(all_summaries, fill = TRUE)
summary_dt[, pair := factor(pair, levels = c("F1_AD", "F2_AD", "F1_PD", "F2_PD"))]
setorder(summary_dt, pair, chrnum, region_start)
summary_dt[, pair := as.character(pair)]
manifest_dt <- rbindlist(manifest, fill = TRUE)

fwrite(summary_dt, file.path(root_out, "coloc_factor_ndd_summary.tsv"), sep = "\t")
fwrite(manifest_dt, file.path(root_out, "coloc_factor_ndd_manifest.tsv"), sep = "\t")

strong_h4 <- summary_dt[status == "ok" & PP.H4.abf >= 0.80]
moderate_h4 <- summary_dt[status == "ok" & PP.H4.abf >= 0.50 & PP.H4.abf < 0.80]
strong_h3 <- summary_dt[status == "ok" & PP.H3.abf >= 0.80]

lines <- c(
  "# coloc summary",
  "",
  "Workflow:",
  "- Official `coloc` R package (`coloc.abf`) with default priors.",
  "- Candidate regions were defined from factor-level conjFDR loci and merged into 500 kb region windows.",
  "- All analyses used harmonized summary statistics from the project-standardized GWAS files.",
  "- AD and PD were modeled as case-control traits using proxy-inclusive primary case fractions.",
  "",
  sprintf("Successful tasks: %d / %d", summary_dt[status == "ok", .N], nrow(summary_dt)),
  sprintf("Strong H4 (PP.H4 >= 0.80): %d", nrow(strong_h4)),
  sprintf("Moderate H4 (0.50 <= PP.H4 < 0.80): %d", nrow(moderate_h4)),
  sprintf("Strong H3 (PP.H3 >= 0.80): %d", nrow(strong_h3)),
  "",
  "Notes:",
  "- `coloc.abf` assumes at most one causal signal per trait per region.",
  "- The chr17 multi-sentinel region (R17) remains the top priority for multi-signal follow-up with `coloc.susie` and/or PWCoCo."
)
writeLines(lines, file.path(root_out, "README.md"))

cat("\nDone. Results written to:\n")
cat(root_out, "\n")
