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

default_output_root <- "/home/shenjing/ctwas_paths/mctwas_work"
default_lib <- "/home/shenjing/R/ctwas_lib"

trait <- get_arg("--trait")
if (is.null(trait) || !nzchar(trait)) {
  stop("Use --trait with one of F1/F2/AD/PD/LBD", call. = FALSE)
}

gwas_file <- get_arg("--gwas-file")
if (is.null(gwas_file) || !file.exists(gwas_file)) {
  stop("Missing --gwas-file or file not found", call. = FALSE)
}

output_root <- get_arg("--output-root", default_output_root)
ctwas_lib <- get_arg("--ctwas-lib", default_lib)
ncore <- as.integer(get_arg("--ncore", "4"))
thin <- as.numeric(get_arg("--thin", "0.1"))
niter_prefit <- as.integer(get_arg("--niter-prefit", "3"))
niter <- as.integer(get_arg("--niter", "30"))
L <- as.integer(get_arg("--L", "5"))
seed <- as.integer(get_arg("--seed", "99"))
max_snp <- as.numeric(get_arg("--max-snp", "Inf"))
init_from_traits_arg <- get_arg("--init-from-traits", "")

.libPaths(c(ctwas_lib, .libPaths()))
suppressPackageStartupMessages({
  library(ctwas)
})

patch_ctwas_fusion_bestcv <- function() {
  fixed_fun <- function(wgt_rdata_file, wgt_ID, fusion_method = c("lasso",
                                                                  "enet",
                                                                  "top1",
                                                                  "blup",
                                                                  "bslmm",
                                                                  "best.cv"),
                        fusion_genome_version = NA) {
    fusion_method <- match.arg(fusion_method)
    snps <- wgt.matrix <- cv.performance <- NULL
    load(wgt_rdata_file)
    if (missing(wgt_ID)) {
      wgt_ID <- sub("\\.wgt\\.RDat$", "", basename(wgt_rdata_file))
    }
    snps <- as.data.frame(snps)
    colnames(snps) <- c("chrom", "rsid", "cm", "pos", "alt", "ref")
    snps$varID <- sprintf(
      "chr%s_%s_%s_%s_%s",
      snps$chrom, snps$pos, snps$ref, snps$alt, fusion_genome_version
    )
    snps[is.na(snps$rsid), "rsid"] <- snps[is.na(snps$rsid), "varID"]
    snps <- snps[, c("chrom", "rsid", "varID", "pos", "ref", "alt")]
    rownames(wgt.matrix) <- snps$rsid

    row.rsq <- grep("rsq", rownames(cv.performance))
    row.pval <- grep("pval", rownames(cv.performance))
    cv.rsq <- cv.performance[row.rsq, , drop = FALSE]

    if (fusion_method == "best.cv") {
      cv_cols <- colnames(cv.performance)
      if (is.null(cv_cols)) {
        cv_cols <- names(as.data.frame(cv.performance))
      }
      best.idx <- which.min(apply(cv.performance[row.pval, , drop = FALSE], 2, min, na.rm = TRUE))
      g.method <- cv_cols[best.idx]
    } else {
      g.method <- fusion_method
    }

    if (!length(g.method) || is.na(g.method) || !g.method %in% colnames(wgt.matrix)) {
      stop("Selected fusion method not found in wgt.matrix")
    }

    g.cv.rsq <- cv.rsq[, g.method, drop = TRUE]
    if (g.method == "top1") {
      wgt.matrix[, "top1"][-which.max(abs(wgt.matrix[, "top1"]))] <- 0
    }

    wgt.matrix <- wgt.matrix[wgt.matrix[, g.method] != 0, , drop = FALSE]
    wgt.matrix <- wgt.matrix[complete.cases(wgt.matrix), , drop = FALSE]

    if (nrow(wgt.matrix) > 0) {
      weight_table <- data.frame(
        gene = wgt_ID,
        rsid = rownames(wgt.matrix),
        weight = wgt.matrix[, g.method],
        stringsAsFactors = FALSE
      )
      weight_table <- merge(weight_table, snps, by = "rsid", all.x = TRUE, sort = FALSE)
      weight_table <- weight_table[, c("gene", "rsid", "varID", "ref", "alt", "weight")]
      colnames(weight_table) <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
    } else {
      weight_table <- data.frame(matrix(ncol = 6, nrow = 0))
      colnames(weight_table) <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
    }

    list(weight_table = weight_table, fusion_method = g.method, wgt_ID = wgt_ID, cv.rsq = g.cv.rsq)
  }

  assignInNamespace("load_fusion_wgt_data", fixed_fun, ns = "ctwas")
}

patch_ctwas_fusion_bestcv()

manifest_dir <- file.path(output_root, "manifests")
weights_manifest_file <- file.path(manifest_dir, "mctwas_hg19_brain9_weights.tsv")
ld_map_file <- file.path(manifest_dir, "mctwas_hg19_ld_map.tsv")
region_info_file <- file.path(manifest_dir, "mctwas_hg19_region_info.rds")

if (!file.exists(weights_manifest_file) || !file.exists(ld_map_file) || !file.exists(region_info_file)) {
  stop("Reference manifests are missing. Run 20_prepare_mctwas_hg19_refs.R first.", call. = FALSE)
}

weights_manifest <- fread(weights_manifest_file)
ld_map <- fread(ld_map_file)
if ("ld_file" %in% names(ld_map) && !"LD_file" %in% names(ld_map)) {
  ld_map[, LD_file := ld_file]
}
if ("snp_info_file" %in% names(ld_map) && !"SNP_file" %in% names(ld_map)) {
  ld_map[, SNP_file := snp_info_file]
}
region_info <- readRDS(region_info_file)
ref_snp_info <- ctwas::read_snp_info_files(ld_map$SNP_file)
snp_map_obj <- ctwas::create_snp_map(region_info, ref_snp_info, ncore = ncore)
region_info <- snp_map_obj$region_info
snp_map <- snp_map_obj$snp_map

trait_output_dir <- file.path(output_root, "single_trait", trait)
input_dir <- file.path(trait_output_dir, "input")
result_dir <- file.path(trait_output_dir, "results")
log_dir <- file.path(trait_output_dir, "logs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

gwas <- fread(gwas_file)
required_cols <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N")
missing_cols <- setdiff(required_cols, names(gwas))
if (length(missing_cols)) {
  stop("GWAS file missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

gwas <- gwas[!is.na(SNP) & !is.na(A1) & !is.na(A2) & !is.na(BETA) & !is.na(SE)]
gwas <- gwas[SE > 0]
gwas[, z := BETA / SE]
z_snp_raw <- gwas[, .(id = SNP, A1, A2, z)]
fwrite(z_snp_raw, file.path(input_dir, sprintf("%s_z_snp_raw.tsv.gz", trait)), sep = "\t")

z_snp <- ctwas::preprocess_z_snp(
  z_snp = z_snp_raw,
  snp_map = snp_map,
  drop_multiallelic = TRUE,
  drop_strand_ambig = TRUE,
  logfile = file.path(log_dir, sprintf("%s_preprocess_z_snp.log", trait))
)

fwrite(z_snp, file.path(input_dir, sprintf("%s_z_snp_harmonized.tsv.gz", trait)), sep = "\t")

weights_by_tissue <- vector("list", nrow(weights_manifest))
names(weights_by_tissue) <- weights_manifest$tissue

for (i in seq_len(nrow(weights_manifest))) {
  tissue <- weights_manifest$tissue[i]
  message("Preprocessing weights: ", tissue)
  weights_by_tissue[[i]] <- ctwas::preprocess_weights(
    weight_file = weights_manifest$extracted_dir[i],
    region_info = region_info,
    gwas_snp_ids = z_snp$id,
    snp_map = snp_map,
    LD_map = ld_map,
    type = "expression",
    context = tissue,
    weight_name = paste0(tissue, "_expression"),
    weight_format = "FUSION",
    drop_strand_ambig = TRUE,
    include_weight_LD = TRUE,
    fusion_method = "best.cv",
    fusion_genome_version = "hg19",
    LD_format = "rds",
    ncore = ncore,
    logfile = file.path(log_dir, sprintf("%s_weights_%s.log", trait, tissue)),
    verbose = FALSE
  )
}

weights_list <- do.call(c, weights_by_tissue)
saveRDS(weights_list, file.path(input_dir, sprintf("%s_weights_preprocessed.rds", trait)), compress = "xz")

init_group_prior <- NULL
init_group_prior_var <- NULL

if (nzchar(init_from_traits_arg)) {
  init_traits <- trimws(unlist(strsplit(init_from_traits_arg, ",", fixed = TRUE)))
  init_traits <- init_traits[nzchar(init_traits)]
  if (!length(init_traits)) {
    stop("No valid traits found in --init-from-traits", call. = FALSE)
  }

  param_list <- vector("list", length(init_traits))
  names(param_list) <- init_traits

  for (i in seq_along(init_traits)) {
    init_trait <- init_traits[i]
    init_file <- file.path(output_root, "single_trait", init_trait, "results",
                           sprintf("%s_ctwas_full.rds", init_trait))
    if (!file.exists(init_file)) {
      stop("Initialization file not found: ", init_file, call. = FALSE)
    }
    init_obj <- readRDS(init_file)
    if (is.null(init_obj$param$group_prior) || is.null(init_obj$param$group_prior_var)) {
      stop("Initialization file missing group_prior/group_prior_var: ", init_file, call. = FALSE)
    }
    param_list[[i]] <- init_obj$param
  }

  prior_names <- names(param_list[[1]]$group_prior)
  var_names <- names(param_list[[1]]$group_prior_var)

  if (any(vapply(param_list, function(x) !identical(names(x$group_prior), prior_names), logical(1)))) {
    stop("Initialization group_prior names do not match across traits", call. = FALSE)
  }
  if (any(vapply(param_list, function(x) !identical(names(x$group_prior_var), var_names), logical(1)))) {
    stop("Initialization group_prior_var names do not match across traits", call. = FALSE)
  }

  init_group_prior <- Reduce(`+`, lapply(param_list, function(x) x$group_prior)) / length(param_list)
  init_group_prior_var <- Reduce(`+`, lapply(param_list, function(x) x$group_prior_var)) / length(param_list)

  init_dt <- data.table(
    parameter = c(rep("group_prior", length(init_group_prior)),
                  rep("group_prior_var", length(init_group_prior_var))),
    group = c(names(init_group_prior), names(init_group_prior_var)),
    value = c(as.numeric(init_group_prior), as.numeric(init_group_prior_var)),
    source_traits = paste(init_traits, collapse = ",")
  )
  fwrite(
    init_dt,
    file.path(input_dir, sprintf("%s_init_params.tsv", trait)),
    sep = "\t"
  )
}

res <- ctwas::ctwas_sumstats(
  z_snp = z_snp,
  weights = weights_list,
  region_info = region_info,
  LD_map = ld_map,
  snp_map = snp_map,
  thin = thin,
  niter_prefit = niter_prefit,
  niter = niter,
  L = L,
  init_group_prior = init_group_prior,
  init_group_prior_var = init_group_prior_var,
  maxSNP = max_snp,
  group_prior_var_structure = "shared_context",
  min_group_size = 100,
  outputdir = result_dir,
  outname = trait,
  ncore = ncore,
  ncore_LD = max(1L, ncore - 1L),
  seed = seed,
  logfile = file.path(log_dir, sprintf("%s_ctwas.log", trait)),
  verbose = TRUE
)

saveRDS(res, file.path(result_dir, sprintf("%s_ctwas_full.rds", trait)), compress = "xz")

finemap_file <- file.path(result_dir, sprintf("%s.susieIrss.txt", trait))
if (file.exists(finemap_file)) {
  finemap <- fread(finemap_file)
  id_col <- intersect(c("variable", "id"), names(finemap))
  pip_col <- intersect(c("pip", "PIP"), names(finemap))
  cs_col <- intersect(c("cs", "CS"), names(finemap))
  out <- copy(finemap)
  if (length(id_col)) {
    out[, feature := get(id_col[1])]
    out[, tissue := fifelse(grepl("_expression$", feature), sub("_expression$", "", feature), NA_character_)]
    out[, model := fifelse(!is.na(tissue), "FUSION_GTEx_Broad_hg19", NA_character_)]
  }
  if (length(pip_col)) {
    setnames(out, pip_col[1], "pip")
  }
  if (length(cs_col)) {
    setnames(out, cs_col[1], "credible_set")
  }
  out[, trait := trait]
  fwrite(out, file.path(result_dir, sprintf("%s_finemap_annotated.tsv", trait)), sep = "\t")
}

summary_dt <- data.table(
  trait = trait,
  input_gwas = gwas_file,
  n_input = nrow(gwas),
  n_harmonized = nrow(z_snp),
  n_weights = length(weights_list),
  result_dir = result_dir,
  run_time = as.character(Sys.time())
)
fwrite(summary_dt, file.path(result_dir, sprintf("%s_run_summary.tsv", trait)), sep = "\t")

message("Completed cTWAS for ", trait)
