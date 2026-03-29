#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- "D:/文章/GS/postgwas/08_pwcoco/inputs/R17"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

region_chr <- 17L
region_start <- 43070680L
region_end <- 44860349L
ref_prefix <- "D:/Mixer/mixer/reference/1000G_EUR_Phase3_plink/1000G.EUR.QC.17"
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

read_region <- function(path, with_case = FALSE, n_case = NA_real_, ref_snps) {
  dt <- fread(
    path,
    select = c("SNP", "CHR", "BP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"),
    showProgress = FALSE
  )
  dt <- dt[
    CHR == region_chr &
      BP >= region_start &
      BP <= region_end &
      SNP %in% ref_snps
  ]
  setorder(dt, BP, SNP)
  if (!nrow(dt)) {
    return(dt)
  }
  dt
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
    ref_A2 = toupper(ref_A2)
  )]
  merged[, `:=`(
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
  out[, `:=`(
    A1 = ref_A1,
    A2 = ref_A2
  )]
  out <- out[!is.na(FREQ) & FREQ > 0 & FREQ < 1]
  setorder(out, BP, SNP)
  out
}

format_for_pwcoco <- function(dt, with_case = FALSE, n_case = NA_real_) {
  if (!nrow(dt)) {
    return(dt)
  }
  if (with_case) {
    out <- dt[, .(
      SNP, A1, A2,
      A1_freq = FREQ,
      beta = BETA,
      se = SE,
      p = P,
      n = N,
      case = n_case
    )]
  } else {
    out <- dt[, .(
      SNP, A1, A2,
      A1_freq = FREQ,
      beta = BETA,
      se = SE,
      p = P,
      n = N
    )]
  }
  out
}

ref <- fread(ref_bim, header = FALSE, showProgress = FALSE)
setnames(ref, c("CHR", "SNP", "CM", "BP", "A1", "A2"))
ref_snps <- unique(ref[
  BP >= region_start & BP <= region_end,
  SNP
])
ref_region <- ref[
  BP >= region_start & BP <= region_end,
  .(SNP, A1, A2, BP)
]

manifest <- list()
for (trait in names(trait_files)) {
  trait_raw <- read_region(trait_files[[trait]], with_case = FALSE, ref_snps = ref_snps)
  trait_dt <- format_for_pwcoco(
    harmonize_to_ref(trait_raw, ref_region, freq_mode = "maf"),
    with_case = FALSE
  )
  trait_out <- file.path(out_dir, sprintf("%s_R17_pwcoco.tsv", trait))
  fwrite(trait_dt, trait_out, sep = "\t", quote = FALSE, na = "NA")
  manifest[[length(manifest) + 1L]] <- data.table(
    dataset = trait,
    role = "trait",
    output_file = trait_out,
    nsnps = nrow(trait_dt),
    n = if (nrow(trait_dt)) unique(trait_dt$n)[1] else NA_real_,
    case = NA_real_
  )
}

for (disease in names(disease_files)) {
  disease_raw <- read_region(
    disease_files[[disease]],
    with_case = TRUE,
    n_case = disease_cases[[disease]],
    ref_snps = ref_snps
  )
  disease_dt <- format_for_pwcoco(
    harmonize_to_ref(disease_raw, ref_region, freq_mode = "a1freq"),
    with_case = TRUE,
    n_case = disease_cases[[disease]]
  )
  disease_out <- file.path(out_dir, sprintf("%s_R17_pwcoco.tsv", disease))
  fwrite(disease_dt, disease_out, sep = "\t", quote = FALSE, na = "NA")
  manifest[[length(manifest) + 1L]] <- data.table(
    dataset = disease,
    role = "disease",
    output_file = disease_out,
    nsnps = nrow(disease_dt),
    n = if (nrow(disease_dt)) unique(disease_dt$n)[1] else NA_real_,
    case = if (nrow(disease_dt)) unique(disease_dt$case)[1] else NA_real_
  )
}

fwrite(rbindlist(manifest), file.path(out_dir, "R17_input_manifest.tsv"), sep = "\t")

pairs <- data.table(
  pair = c("F1_AD_R17", "F2_AD_R17", "F1_PD_R17", "F2_PD_R17"),
  sum_stats1 = file.path(out_dir, c("F1_R17_pwcoco.tsv", "F2_R17_pwcoco.tsv", "F1_R17_pwcoco.tsv", "F2_R17_pwcoco.tsv")),
  sum_stats2 = file.path(out_dir, c("AD_R17_pwcoco.tsv", "AD_R17_pwcoco.tsv", "PD_R17_pwcoco.tsv", "PD_R17_pwcoco.tsv")),
  bfile = ref_prefix,
  chr = region_chr
)
fwrite(pairs, file.path(out_dir, "R17_pwcoco_pairs.tsv"), sep = "\t")

cat("Prepared PWCoCo R17 inputs in:\n")
cat(out_dir, "\n")
