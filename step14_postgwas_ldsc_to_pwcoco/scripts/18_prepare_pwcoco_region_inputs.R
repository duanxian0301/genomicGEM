#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 18_prepare_pwcoco_region_inputs.R <region_id> <chr> <start> <end>")
}

region_id <- args[[1]]
region_chr <- as.integer(args[[2]])
region_start <- as.integer(args[[3]])
region_end <- as.integer(args[[4]])

out_dir <- file.path("D:/文章/GS/postgwas/08_pwcoco/inputs", region_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ref_prefix <- sprintf(
  "D:/Mixer/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s",
  region_chr
)
ref_bim <- paste0(ref_prefix, ".bim")

f1_file <- "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F1_factorGWAS_native_standard.txt"
f2_file <- "D:/文章/GS/GWAS/step11_factor_gwas_native_results/merged/standard_txt/ALPS_F2_factorGWAS_native_standard.txt"
ad_file <- "D:/文章/4NDD/NDDGWAS/AD.txt"
pd_file <- "D:/文章/4NDD/NDDGWAS/PD.txt"

trait_files <- list(F1 = f1_file, F2 = f2_file)
disease_files <- list(AD = ad_file, PD = pd_file)
disease_cases <- c(AD = 39106 + 46828, PD = 63555 + 17700)

comp_map <- c(A = "T", T = "A", C = "G", G = "C")
comp_vec <- function(x) {
  y <- comp_map[toupper(x)]
  y[is.na(y)] <- NA_character_
  unname(y)
}

read_region <- function(path, ref_snps) {
  dt <- fread(
    path,
    select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
    showProgress = FALSE
  )
  dt[
    CHR == region_chr &
      BP >= region_start &
      BP <= region_end &
      SNP %in% ref_snps
  ]
}

harmonize_to_ref <- function(dt, ref, freq_mode = c("maf", "a1freq")) {
  freq_mode <- match.arg(freq_mode)
  merged <- merge(dt, ref[, .(SNP, ref_A1 = A1, ref_A2 = A2)], by = "SNP")
  if (!nrow(merged)) {
    return(merged[0])
  }
  merged[, `:=`(
    A1 = toupper(A1),
    A2 = toupper(A2),
    ref_A1 = toupper(ref_A1),
    ref_A2 = toupper(ref_A2),
    A1_comp = comp_vec(A1),
    A2_comp = comp_vec(A2)
  )]

  same <- merged$A1 == merged$ref_A1 & merged$A2 == merged$ref_A2
  swap <- merged$A1 == merged$ref_A2 & merged$A2 == merged$ref_A1
  comp_same <- merged$A1_comp == merged$ref_A1 & merged$A2_comp == merged$ref_A2
  comp_swap <- merged$A1_comp == merged$ref_A2 & merged$A2_comp == merged$ref_A1

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
    return(out)
  }
  out[align_status %in% c("swap", "comp_swap"), BETA := -BETA]
  if (freq_mode == "a1freq") {
    out[align_status %in% c("swap", "comp_swap"), FREQ := 1 - FREQ]
  }
  out[, `:=`(A1 = ref_A1, A2 = ref_A2)]
  out <- out[!is.na(FREQ) & FREQ > 0 & FREQ < 1]
  setorder(out, BP, SNP)
  out
}

format_for_pwcoco <- function(dt, with_case = FALSE, n_case = NA_real_) {
  if (!nrow(dt)) {
    return(dt)
  }
  if (with_case) {
    dt[, .(
      SNP, A1, A2,
      A1_freq = FREQ,
      beta = BETA,
      se = SE,
      p = P,
      n = N,
      case = n_case
    )]
  } else {
    dt[, .(
      SNP, A1, A2,
      A1_freq = FREQ,
      beta = BETA,
      se = SE,
      p = P,
      n = N
    )]
  }
}

ref <- fread(ref_bim, header = FALSE, showProgress = FALSE)
setnames(ref, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
ref_region <- ref[
  BP >= region_start & BP <= region_end,
  .(SNP, A1, A2, BP)
]
ref_snps <- unique(ref_region$SNP)

manifest <- list()
for (trait in names(trait_files)) {
  raw <- read_region(trait_files[[trait]], ref_snps)
  out <- format_for_pwcoco(harmonize_to_ref(raw, ref_region, "maf"), with_case = FALSE)
  out_file <- file.path(out_dir, sprintf("%s_%s_pwcoco.tsv", trait, region_id))
  fwrite(out, out_file, sep = "\t", quote = FALSE, na = "NA")
  manifest[[length(manifest) + 1L]] <- data.table(
    dataset = trait, role = "trait", output_file = out_file, nsnps = nrow(out),
    n = if (nrow(out)) unique(out$n)[1] else NA_real_, case = NA_real_
  )
}

for (disease in names(disease_files)) {
  raw <- read_region(disease_files[[disease]], ref_snps)
  out <- format_for_pwcoco(harmonize_to_ref(raw, ref_region, "a1freq"), with_case = TRUE, n_case = disease_cases[[disease]])
  out_file <- file.path(out_dir, sprintf("%s_%s_pwcoco.tsv", disease, region_id))
  fwrite(out, out_file, sep = "\t", quote = FALSE, na = "NA")
  manifest[[length(manifest) + 1L]] <- data.table(
    dataset = disease, role = "disease", output_file = out_file, nsnps = nrow(out),
    n = if (nrow(out)) unique(out$n)[1] else NA_real_, case = if (nrow(out)) unique(out$case)[1] else NA_real_
  )
}

fwrite(rbindlist(manifest), file.path(out_dir, sprintf("%s_input_manifest.tsv", region_id)), sep = "\t")
cat("Prepared inputs for", region_id, "in", out_dir, "\n")
