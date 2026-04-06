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

pair_defs <- list(
  c("F1", "AD"),
  c("F2", "AD"),
  c("F1", "PD"),
  c("F2", "PD"),
  c("F1", "LBD"),
  c("F2", "LBD")
)

trait_info <- list(
  F1 = list(gwas = "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F1_factorGWAS_native_standard.txt"),
  F2 = list(gwas = "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F2_factorGWAS_native_standard.txt"),
  AD = list(gwas = "D:/文章/4NDD/NDDGWAS/AD.txt"),
  PD = list(gwas = "D:/文章/4NDD/NDDGWAS/PD.txt"),
  LBD = list(gwas = "D:/文章/4NDD/NDDGWAS/LBD.txt")
)

success_traits <- c("F1", "F2", "AD", "PD")

model_name <- "Broad_FUSION_GTEx_brain_hg19_9tissues"

read_run_summary <- function(trait) {
  path <- file.path(output_root, "single_trait", trait, "results", sprintf("%s_run_summary.tsv", trait))
  if (!file.exists(path)) return(NULL)
  fread(path)
}

extract_trait_table <- function(trait) {
  rds_file <- file.path(output_root, "single_trait", trait, "results", sprintf("%s_ctwas_full.rds", trait))
  if (!file.exists(rds_file)) {
    stop("Missing ctwas result for trait: ", trait, call. = FALSE)
  }
  x <- readRDS(rds_file)
  finemap <- as.data.table(x$finemap_res)
  finemap <- finemap[type == "expression"]
  if (!nrow(finemap)) {
    return(data.table())
  }
  finemap[, p := 2 * pnorm(-abs(z))]
  finemap[, FDR := p.adjust(p, method = "BH")]
  finemap[, `whether significant` := FDR < 0.05]
  finemap[, trait := trait]
  finemap[, pair_source_file := normalizePath(rds_file, winslash = "/", mustWork = TRUE)]
  finemap[, `gene / feature` := molecular_id]
  finemap[, `tissue / model` := paste0(context, " / ", model_name)]
  finemap[, effect := z]
  finemap[, source_files := pair_source_file]
  finemap[, feature_key := paste(molecular_id, context, sep = "||")]
  finemap[, feature_id := id]
  finemap[, pip := susie_pip]
  finemap[, credible_set := cs]
  finemap[, .(
    trait,
    feature_key,
    `gene / feature`,
    tissue = context,
    model = model_name,
    `tissue / model`,
    z = z,
    effect = effect,
    p,
    FDR,
    `whether significant`,
    pip,
    credible_set,
    feature_id,
    source_files
  )]
}

trait_tables <- lapply(success_traits, extract_trait_table)
names(trait_tables) <- success_traits

single_trait_table <- rbindlist(trait_tables, use.names = TRUE, fill = TRUE)
single_trait_out <- file.path(summary_dir, "mctwas_single_trait_expression_results.tsv")
fwrite(single_trait_table, single_trait_out, sep = "\t")

run_status_rows <- list()
for (tr in c("F1", "F2", "AD", "PD", "LBD")) {
  rs <- read_run_summary(tr)
  if (!is.null(rs)) {
    run_status_rows[[tr]] <- data.table(
      trait = tr,
      status = "success",
      n_input = rs$n_input,
      n_harmonized = rs$n_harmonized,
      n_weights = rs$n_weights,
      input_gwas = trait_info[[tr]]$gwas,
      result_dir = sub("^/home/shenjing/ctwas_paths/mctwas_work", "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd", rs$result_dir),
      note = NA_character_
    )
  } else {
    note <- if (tr == "LBD") {
      "ctwas_sumstats failed twice at EM parameter estimation; group_prior became NA, so no stable M-cTWAS result was retained."
    } else {
      "No completed run summary found."
    }
    run_status_rows[[tr]] <- data.table(
      trait = tr,
      status = "failed",
      n_input = NA_integer_,
      n_harmonized = NA_integer_,
      n_weights = NA_integer_,
      input_gwas = trait_info[[tr]]$gwas,
      result_dir = sprintf("D:/文章/GS/postgwas/05_mctwas_f1f2_ndd/single_trait/%s/results", tr),
      note = note
    )
  }
}
run_status <- rbindlist(run_status_rows, use.names = TRUE, fill = TRUE)
run_status_out <- file.path(summary_dir, "mctwas_run_status.tsv")
fwrite(run_status, run_status_out, sep = "\t")

make_pair_table <- function(left, right) {
  pair_name <- sprintf("%s_vs_%s", left, right)
  left_ok <- left %in% names(trait_tables)
  right_ok <- right %in% names(trait_tables)

  if (!(left_ok && right_ok)) {
    fail_trait <- if (!left_ok) left else right
    fail_log <- sprintf("D:/文章/GS/postgwas/05_mctwas_f1f2_ndd/single_trait/%s/logs/%s_ctwas.log", fail_trait, fail_trait)
    return(data.table(
      pair = pair_name,
      `gene / feature` = NA_character_,
      `tissue / model` = NA_character_,
      `z / effect` = NA_character_,
      p = NA_character_,
      FDR = NA_character_,
      `whether significant` = NA_character_,
      source_files = fail_log,
      run_status = sprintf("%s failed; no stable pair-level comparison available", fail_trait)
    ))
  }

  x <- copy(trait_tables[[left]])
  y <- copy(trait_tables[[right]])
  setnames(x, old = c("z", "effect", "p", "FDR", "whether significant", "source_files"),
           new = c("z_left", "effect_left", "p_left", "FDR_left", "sig_left", "source_left"))
  setnames(y, old = c("z", "effect", "p", "FDR", "whether significant", "source_files"),
           new = c("z_right", "effect_right", "p_right", "FDR_right", "sig_right", "source_right"))

  merged <- merge(
    x[, .(feature_key, `gene / feature`, tissue, model, `tissue / model`, z_left, effect_left, p_left, FDR_left, sig_left, source_left)],
    y[, .(feature_key, `gene / feature`, tissue, model, `tissue / model`, z_right, effect_right, p_right, FDR_right, sig_right, source_right)],
    by = c("feature_key", "gene / feature", "tissue", "model", "tissue / model"),
    all = TRUE
  )

  merged[, pair := pair_name]
  merged[, `z / effect` := paste0(
    left, ":z=", fifelse(is.na(z_left), "NA", sprintf("%.6f", z_left)),
    ",effect=", fifelse(is.na(effect_left), "NA", sprintf("%.6f", effect_left)),
    "; ",
    right, ":z=", fifelse(is.na(z_right), "NA", sprintf("%.6f", z_right)),
    ",effect=", fifelse(is.na(effect_right), "NA", sprintf("%.6f", effect_right))
  )]
  merged[, p := paste0(
    left, ":", fifelse(is.na(p_left), "NA", formatC(p_left, format = "e", digits = 4)),
    "; ",
    right, ":", fifelse(is.na(p_right), "NA", formatC(p_right, format = "e", digits = 4))
  )]
  merged[, FDR := paste0(
    left, ":", fifelse(is.na(FDR_left), "NA", formatC(FDR_left, format = "e", digits = 4)),
    "; ",
    right, ":", fifelse(is.na(FDR_right), "NA", formatC(FDR_right, format = "e", digits = 4))
  )]
  merged[, `whether significant` := paste0(
    left, ":", fifelse(is.na(sig_left), "NA", as.character(sig_left)),
    "; ",
    right, ":", fifelse(is.na(sig_right), "NA", as.character(sig_right))
  )]
  merged[, source_files := paste0(
    left, ":", fifelse(is.na(source_left), "NA", source_left),
    "; ",
    right, ":", fifelse(is.na(source_right), "NA", source_right)
  )]
  merged[, run_status := "ok"]

  merged[, .(
    pair,
    `gene / feature`,
    `tissue / model`,
    `z / effect`,
    p,
    FDR,
    `whether significant`,
    source_files,
    run_status
  )]
}

pair_tables <- lapply(pair_defs, function(x) make_pair_table(x[1], x[2]))
pair_names <- vapply(pair_defs, function(x) sprintf("%s_vs_%s", x[1], x[2]), character(1))
names(pair_tables) <- pair_names

for (nm in names(pair_tables)) {
  fwrite(pair_tables[[nm]], file.path(summary_dir, sprintf("mctwas_%s.tsv", nm)), sep = "\t")
}

summary_table <- rbindlist(pair_tables, use.names = TRUE, fill = TRUE)
summary_out <- file.path(summary_dir, "mctwas_f1f2_ndd_summary.tsv")
fwrite(summary_table, summary_out, sep = "\t")

sig_by_trait <- single_trait_table[`whether significant` == TRUE, .N, by = trait][order(-N)]
if (!nrow(sig_by_trait)) {
  sig_by_trait <- data.table(trait = success_traits, N = 0L)
}

pair_sig_overlap <- rbindlist(lapply(pair_defs[1:4], function(x) {
  left <- x[1]
  right <- x[2]
  left_sig <- unique(trait_tables[[left]][`whether significant` == TRUE, feature_key])
  right_sig <- unique(trait_tables[[right]][`whether significant` == TRUE, feature_key])
  data.table(
    pair = sprintf("%s_vs_%s", left, right),
    shared_significant_features = length(intersect(left_sig, right_sig)),
    left_significant_features = length(left_sig),
    right_significant_features = length(right_sig)
  )
}), use.names = TRUE, fill = TRUE)

md_lines <- c(
  "# M-cTWAS Summary",
  "",
  "## Scheme",
  "Official single-trait cTWAS was run for F1, F2, AD and PD using hg19/GRCh37 summary statistics, Broad FUSION GTEx brain hg19 weights, and 9 brain tissues. Pair-level results were then generated by comparing matched expression features across traits.",
  "",
  "LBD was attempted twice, including a rerun with externally supplied EM initialization from F1/F2/AD/PD, but ctwas parameter estimation remained unstable and no stable LBD result was retained.",
  "",
  "## Reference Weights",
  paste0("- ", model_name),
  "- Brain_Caudate_basal_ganglia",
  "- Brain_Cerebellar_Hemisphere",
  "- Brain_Cerebellum",
  "- Brain_Cortex",
  "- Brain_Frontal_Cortex_BA9",
  "- Brain_Hippocampus",
  "- Brain_Hypothalamus",
  "- Brain_Nucleus_accumbens_basal_ganglia",
  "- Brain_Putamen_basal_ganglia",
  "",
  "## Key Parameters",
  "- ctwas 0.4.21",
  "- thin = 0.1",
  "- niter_prefit = 3",
  "- niter = 30",
  "- L = 5",
  "- group_prior_var_structure = shared_context",
  "- min_group_size = 100",
  "",
  "## Run Status",
  paste0("- ", run_status$trait, ": ", run_status$status, ifelse(is.na(run_status$note), "", paste0(" (", run_status$note, ")"))),
  "",
  "## Significant Feature Counts",
  paste0("- ", sig_by_trait$trait, ": ", sig_by_trait$N, " FDR<0.05 expression features"),
  "",
  "## Pair-Level Signal Snapshot",
  if (nrow(pair_sig_overlap)) paste0("- ", pair_sig_overlap$pair, ": shared=", pair_sig_overlap$shared_significant_features,
                                     ", left=", pair_sig_overlap$left_significant_features,
                                     ", right=", pair_sig_overlap$right_significant_features) else "- No stable pair overlap could be computed.",
  "",
  "## Consistency With Existing Results",
  "Within the successfully completed traits, AD and F2 provide the cleanest completed comparison set, while PD also yields stable cTWAS output for comparison with F1/F2. LBD remains unstable in this framework, which is directionally consistent with your prior note that LBD evidence is overall weaker.",
  "",
  "## Output Files",
  paste0("- Single-trait table: ", sub("^/home/shenjing/ctwas_paths/mctwas_work", "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd", single_trait_out)),
  paste0("- Pair summary table: ", sub("^/home/shenjing/ctwas_paths/mctwas_work", "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd", summary_out)),
  paste0("- Run status table: ", sub("^/home/shenjing/ctwas_paths/mctwas_work", "D:/文章/GS/postgwas/05_mctwas_f1f2_ndd", run_status_out))
)

md_out <- file.path(summary_dir, "mctwas_f1f2_ndd_summary.md")
writeLines(md_lines, md_out, useBytes = TRUE)

message("Wrote: ", summary_out)
