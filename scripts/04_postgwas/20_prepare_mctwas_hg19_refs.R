#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

default_output_root <- "/home/shenjing/ctwas_paths/mctwas_work"
default_lib <- "/home/shenjing/R/ctwas_lib"

get_arg <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) return(default)
  if (hit[length(hit)] == length(args)) {
    stop("Missing value for ", flag, call. = FALSE)
  }
  args[hit[length(hit)] + 1L]
}

flag_present <- function(flag) {
  flag %in% args
}

output_root <- get_arg("--output-root", default_output_root)
ctwas_lib <- get_arg("--ctwas-lib", default_lib)
plink_bin <- get_arg(
  "--plink-bin",
  file.path(output_root, "refs", "tools", "plink", "plink")
)
geno_prefix <- get_arg(
  "--geno-prefix",
  "/home/shenjing/ctwas_paths/gnova_eur1000g/EUR.QC"
)
ncore <- as.integer(get_arg("--ncore", "4"))
force <- flag_present("--force")

.libPaths(c(ctwas_lib, .libPaths()))
suppressPackageStartupMessages({
  library(ctwas)
})

brain_tissues <- c(
  "Brain_Caudate_basal_ganglia",
  "Brain_Cerebellar_Hemisphere",
  "Brain_Cerebellum",
  "Brain_Cortex",
  "Brain_Frontal_Cortex_BA9",
  "Brain_Hippocampus",
  "Brain_Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia",
  "Brain_Putamen_basal_ganglia"
)

refs_dir <- file.path(output_root, "refs")
fusion_dir <- file.path(refs_dir, "fusion_broad_gtex_old")
fusion_tar_dir <- file.path(fusion_dir, "tarballs")
fusion_extract_dir <- file.path(fusion_dir, "extracted")
ld_dir <- file.path(refs_dir, "ldref_1000g_eur_b37")
plink_chr_dir <- file.path(ld_dir, "plink_by_chr_ascii")
ctwas_ld_dir <- file.path(ld_dir, "ctwas_ld_all_ascii")
manifest_dir <- file.path(output_root, "manifests")

dir.create(refs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fusion_tar_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fusion_extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plink_chr_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ctwas_ld_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

region_info <- readRDS(
  system.file("extdata", "ldetect", "EUR.b37.ldetect.regions.RDS", package = "ctwas")
)

download_one_tissue <- function(tissue) {
  tar_name <- sprintf("GTEx.%s.tar.bz2", tissue)
  url <- sprintf(
    "https://storage.googleapis.com/broad-alkesgroup-public/FUSION/WGT/%s",
    tar_name
  )
  tar_path <- file.path(fusion_tar_dir, tar_name)
  tissue_dir <- file.path(fusion_extract_dir, sprintf("GTEx.%s", tissue))
  pos_file <- file.path(fusion_extract_dir, sprintf("GTEx.%s.pos", tissue))

  if (!file.exists(tar_path) || force) {
    message("Downloading ", tar_name)
    status <- system2("wget", c("-q", "-O", tar_path, url))
    if (!identical(status, 0L)) {
      stop("Failed downloading ", url, call. = FALSE)
    }
  }

  if (!dir.exists(tissue_dir) || force || !file.exists(pos_file)) {
    if (dir.exists(tissue_dir) && force) {
      unlink(tissue_dir, recursive = TRUE, force = TRUE)
    }
    dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
    message("Extracting ", tar_name)
    status <- system2("tar", c("-xjf", tar_path, "-C", fusion_extract_dir))
    if (!identical(status, 0L)) {
      stop("Failed extracting ", tar_path, call. = FALSE)
    }
  }

  if (!file.exists(pos_file)) {
    stop("Missing extracted .pos file for ", tissue, call. = FALSE)
  }

  data.table(
    tissue = tissue,
    tarball = tar_path,
    extracted_dir = tissue_dir,
    pos_file = pos_file
  )
}

split_chr_if_needed <- function(chr) {
  out_prefix <- file.path(plink_chr_dir, sprintf("1000G.EUR.QC.%s", chr))
  bed_path <- paste0(out_prefix, ".bed")
  if (file.exists(bed_path) && !force) {
    return(out_prefix)
  }
  message("Splitting LD reference chr", chr)
  status <- system2(
    plink_bin,
    c(
      "--bfile", geno_prefix,
      "--chr", chr,
      "--make-bed",
      "--out", out_prefix
    )
  )
  if (!identical(status, 0L)) {
    stop("PLINK failed on chr", chr, call. = FALSE)
  }
  out_prefix
}

tissue_manifest <- rbindlist(lapply(brain_tissues, download_one_tissue))
fwrite(tissue_manifest, file.path(manifest_dir, "mctwas_hg19_brain9_weights.tsv"), sep = "\t")

chroms <- sort(unique(region_info$chrom))
geno_files <- character(max(chroms))
varinfo_files <- character(max(chroms))

for (chr in chroms) {
  geno_chr_prefix <- split_chr_if_needed(chr)
  geno_files[chr] <- paste0(geno_chr_prefix, ".bed")
  varinfo_files[chr] <- paste0(geno_chr_prefix, ".bim")
}

message("Building ctwas LD blocks for all chromosomes")
ld_map <- as.data.table(
  ctwas::convert_geno_to_LD_matrix(
    region_info = region_info,
    genotype_files = geno_files,
    varinfo_files = varinfo_files,
    chrom = chroms,
    outputdir = ctwas_ld_dir,
    outname = "1000G_EUR_Phase3_b37",
    show_progress_bar = TRUE,
    verbose = FALSE
  )
)
setnames(ld_map, old = c("LD_file", "SNP_file"), new = c("ld_file", "snp_info_file"))

missing_ld <- ld_map[!file.exists(ld_file) | !file.exists(snp_info_file)]
if (nrow(missing_ld)) {
  stop("Missing LD outputs for ", nrow(missing_ld), " regions", call. = FALSE)
}

fwrite(ld_map, file.path(manifest_dir, "mctwas_hg19_ld_map.tsv"), sep = "\t")
saveRDS(region_info, file.path(manifest_dir, "mctwas_hg19_region_info.rds"))

message("Reference preparation complete")
message("Weights manifest: ", file.path(manifest_dir, "mctwas_hg19_brain9_weights.tsv"))
message("LD map manifest: ", file.path(manifest_dir, "mctwas_hg19_ld_map.tsv"))
