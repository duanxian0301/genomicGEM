#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[hit[length(hit)] + 1L]
}

output_root <- get_arg("--output-root", "/home/shenjing/ctwas_paths/mctwas_work")
summary_dir <- file.path(output_root, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

success_traits <- c("F1", "F2", "AD", "PD")
all_traits <- c(success_traits, "LBD")
model_name <- "Broad_FUSION_GTEx_brain_hg19_9tissues"

parse_gene_symbol <- function(feature_id) {
  x <- sub("\\|.*$", "", feature_id)
  sub("^[^.]+\\.", "", x)
}

parse_region_coords <- function(region_id) {
  parts <- tstrsplit(region_id, "_", fixed = TRUE)
  data.table(
    chr = suppressWarnings(as.integer(parts[[1]])),
    region_start = suppressWarnings(as.integer(parts[[2]])),
    region_end = suppressWarnings(as.integer(parts[[3]]))
  )
}

feature_rows <- list()
variant_rows <- list()
gene_rows <- list()
tissue_rows <- list()
status_rows <- list()

for (trait in success_traits) {
  ctwas_file <- file.path(output_root, "single_trait", trait, "results", sprintf("%s_ctwas_full.rds", trait))
  summary_file <- file.path(output_root, "single_trait", trait, "results", sprintf("%s_run_summary.tsv", trait))
  x <- readRDS(ctwas_file)
  run_summary <- fread(summary_file)

  finemap <- as.data.table(x$finemap_res)
  finemap[, p := 2 * pnorm(-abs(z))]

  expr_dt <- copy(finemap[type == "expression"])
  expr_dt[, gene_symbol := parse_gene_symbol(id)]
  expr_dt[, tissue := context]
  expr_dt[, model := model_name]
  expr_dt[, FDR := p.adjust(p, method = "BH")]
  expr_dt[, prioritized_fdr := FDR < 0.05]
  expr_dt[, prioritized_pip := susie_pip >= 0.5]
  expr_dt[, prioritized_any := prioritized_fdr | prioritized_pip]
  expr_dt[, source_file := normalizePath(ctwas_file, winslash = "/", mustWork = TRUE)]
  expr_dt <- cbind(expr_dt, parse_region_coords(expr_dt$region_id))
  expr_dt[, trait := trait]

  feature_rows[[trait]] <- expr_dt[, .(
    trait,
    gene_symbol,
    molecular_feature = id,
    molecular_id,
    tissue,
    model,
    region_id,
    chr,
    region_start,
    region_end,
    z,
    p,
    FDR,
    pip = susie_pip,
    credible_set = cs,
    prioritized_fdr,
    prioritized_pip,
    prioritized_any,
    source_file
  )]

  variant_dt <- copy(finemap[type == "SNP"])
  variant_dt[, p := 2 * pnorm(-abs(z))]
  variant_dt[, FDR := p.adjust(p, method = "BH")]
  variant_dt[, prioritized_any := (FDR < 0.05) | (susie_pip >= 0.5)]
  variant_dt <- cbind(variant_dt, parse_region_coords(variant_dt$region_id))
  variant_dt[, trait := trait]
  variant_rows[[trait]] <- variant_dt[, .(
    trait,
    variant_id = id,
    region_id,
    chr,
    region_start,
    region_end,
    z,
    p,
    FDR,
    pip = susie_pip,
    credible_set = cs,
    prioritized_any
  )]

  expr_for_gene <- feature_rows[[trait]]
  expr_for_gene[, abs_z := abs(z)]
  setorder(expr_for_gene, -pip, FDR, -abs_z)
  gene_agg <- expr_for_gene[, .(
    best_feature = molecular_feature[1],
    best_tissue = tissue[1],
    best_region_id = region_id[1],
    chr = chr[1],
    region_start = region_start[1],
    region_end = region_end[1],
    best_z = z[1],
    best_p = p[1],
    min_FDR = min(FDR, na.rm = TRUE),
    max_pip = max(pip, na.rm = TRUE),
    n_tested_contexts = .N,
    n_prioritized_contexts = sum(prioritized_any, na.rm = TRUE),
    n_fdr_contexts = sum(prioritized_fdr, na.rm = TRUE),
    n_highpip_contexts = sum(prioritized_pip, na.rm = TRUE),
    prioritized_any = any(prioritized_any, na.rm = TRUE),
    prioritized_fdr = any(prioritized_fdr, na.rm = TRUE),
    prioritized_pip = any(prioritized_pip, na.rm = TRUE)
  ), by = .(trait, gene_symbol)]
  gene_agg[, priority_label := fifelse(
    prioritized_pip | min_FDR < 0.05, "high_priority",
    fifelse(max_pip >= 0.1 | min_FDR < 0.2, "moderate_priority", "exploratory")
  )]
  gene_rows[[trait]] <- gene_agg

  param <- x$param
  gp <- data.table(
    group = names(param$group_prior),
    group_prior = as.numeric(param$group_prior),
    group_prior_var = as.numeric(param$group_prior_var)
  )
  gp <- gp[group != "SNP"]
  gp[, c("tissue", "type") := tstrsplit(group, "|", fixed = TRUE)]
  tissue_support <- expr_for_gene[, .(
    n_features = .N,
    n_prioritized_features = sum(prioritized_any, na.rm = TRUE),
    n_fdr_features = sum(prioritized_fdr, na.rm = TRUE),
    n_highpip_features = sum(prioritized_pip, na.rm = TRUE),
    median_abs_z = median(abs(z), na.rm = TRUE),
    max_pip = max(pip, na.rm = TRUE)
  ), by = tissue]
  tissue_rows[[trait]] <- merge(
    gp[, .(tissue, type, group_prior, group_prior_var)],
    tissue_support,
    by = "tissue",
    all.x = TRUE
  )[, trait := trait][]

  status_rows[[trait]] <- data.table(
    trait = trait,
    status = "success",
    n_input = run_summary$n_input,
    n_harmonized = run_summary$n_harmonized,
    n_weights = run_summary$n_weights,
    result_dir = sub("^/home/shenjing/ctwas_paths/mctwas_work", "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd", run_summary$result_dir),
    note = NA_character_
  )
}

status_rows[["LBD"]] <- data.table(
  trait = "LBD",
  status = "failed",
  n_input = NA_integer_,
  n_harmonized = NA_integer_,
  n_weights = NA_integer_,
  result_dir = "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd/single_trait/LBD/results",
  note = "ctwas_sumstats failed twice during EM parameter estimation; group_prior became NA, so no stable gene- or tissue-level result was retained."
)

feature_dt <- rbindlist(feature_rows, use.names = TRUE, fill = TRUE)
variant_dt <- rbindlist(variant_rows, use.names = TRUE, fill = TRUE)
gene_dt <- rbindlist(gene_rows, use.names = TRUE, fill = TRUE)
tissue_dt <- rbindlist(tissue_rows, use.names = TRUE, fill = TRUE)
status_dt <- rbindlist(status_rows, use.names = TRUE, fill = TRUE)

fwrite(feature_dt, file.path(summary_dir, "ctwas_trait_level_feature_results.tsv"), sep = "\t")
fwrite(variant_dt, file.path(summary_dir, "ctwas_trait_level_variant_results.tsv"), sep = "\t")
fwrite(gene_dt, file.path(summary_dir, "ctwas_trait_level_gene_results.tsv"), sep = "\t")
fwrite(tissue_dt, file.path(summary_dir, "ctwas_trait_level_tissue_results.tsv"), sep = "\t")
fwrite(status_dt, file.path(summary_dir, "ctwas_trait_run_status.tsv"), sep = "\t")

message("Wrote rich cTWAS result tables to ", summary_dir)
