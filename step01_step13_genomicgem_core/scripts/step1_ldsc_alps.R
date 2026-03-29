library(GenomicSEM)
library(data.table)
library(stringr)

# -----------------------------
# ALPS GenomicSEM Step 1:
# 1) munge raw GWAS if needed
# 2) summarize univariate LDSC QC
# 3) run multivariate LDSC on retained traits
# -----------------------------

gwas_dir <- "D:/文章/GS/GWAS"
ld_ref_dir <- "D:/LDSC/ldsc-master/eur_w_ld_chr"
hm3 <- file.path(ld_ref_dir, "w_hm3.snplist")
output_dir <- file.path(gwas_dir, "step1_ldsc_results")

# If TRUE, re-run munge() from raw .txt files even when .sumstats.gz already exists.
force_munge <- FALSE
# If TRUE, re-run univariate LDSC even when the single-trait log already exists.
force_univariate_ldsc <- FALSE

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(gwas_dir)

txt_files <- sort(list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE))
if (length(txt_files) == 0) {
  stop("No .txt GWAS files found in: ", gwas_dir)
}

trait_names <- tools::file_path_sans_ext(basename(txt_files))
sumstats_files <- file.path(gwas_dir, paste0(trait_names, ".sumstats.gz"))

if (force_munge || !all(file.exists(sumstats_files))) {
  message("Running munge() on raw .txt files ...")
  munge(
    files = txt_files,
    hm3 = hm3,
    trait.names = trait_names,
    N = rep(NA, length(txt_files))
  )
} else {
  message("Reusing existing .sumstats.gz files.")
}

sumstats_files <- file.path(gwas_dir, paste0(trait_names, ".sumstats.gz"))
missing_sumstats <- sumstats_files[!file.exists(sumstats_files)]
if (length(missing_sumstats) > 0) {
  stop(
    "Missing .sumstats.gz files after munge/reuse:\n",
    paste(missing_sumstats, collapse = "\n")
  )
}

extract_ldsc_metric <- function(log_file, pattern) {
  lines <- readLines(log_file, warn = FALSE, encoding = "UTF-8")
  hit <- grep(pattern, lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  trimws(sub(".*?:\\s*", "", hit[1]))
}

run_single_trait_ldsc <- function(sumstats_file, trait_name) {
  log_name <- paste0(basename(sumstats_file), "_ldsc")
  log_path <- file.path(gwas_dir, paste0(log_name, ".log"))
  if (force_univariate_ldsc || !file.exists(log_path)) {
    ldsc(
      traits = sumstats_file,
      sample.prev = NA,
      population.prev = NA,
      ld = ld_ref_dir,
      wld = ld_ref_dir,
      trait.names = trait_name,
      ldsc.log = log_name
    )
  }
  data.frame(
    trait = trait_name,
    sumstats = basename(sumstats_file),
    intercept = extract_ldsc_metric(log_path, "^Intercept:"),
    ratio = extract_ldsc_metric(log_path, "^Ratio:"),
    h2 = extract_ldsc_metric(log_path, "^Total Observed Scale h2:"),
    h2_z = suppressWarnings(as.numeric(extract_ldsc_metric(log_path, "^h2 Z:"))),
    log_file = log_path,
    stringsAsFactors = FALSE
  )
}

message("Running univariate LDSC QC for each trait ...")
qc_list <- Map(run_single_trait_ldsc, sumstats_files, trait_names)
qc_table <- rbindlist(qc_list, fill = TRUE)
qc_table$pass_h2_z_ge_4 <- qc_table$h2_z >= 4
fwrite(qc_table, file.path(output_dir, "alps_univariate_ldsc_qc.tsv"), sep = "\t")

eligible_traits <- qc_table$trait[qc_table$pass_h2_z_ge_4 %in% TRUE]
eligible_sumstats <- file.path(gwas_dir, paste0(eligible_traits, ".sumstats.gz"))

if (length(eligible_traits) < 2) {
  stop(
    "Fewer than 2 traits passed h2 Z >= 4. Cannot run multivariate LDSC.\n",
    "See: ", file.path(output_dir, "alps_univariate_ldsc_qc.tsv")
  )
}

message("Running multivariate LDSC on retained traits: ",
        paste(eligible_traits, collapse = ", "))

multi_ldsc <- ldsc(
  traits = eligible_sumstats,
  sample.prev = rep(NA, length(eligible_sumstats)),
  population.prev = rep(NA, length(eligible_sumstats)),
  ld = ld_ref_dir,
  wld = ld_ref_dir,
  trait.names = eligible_traits,
  ldsc.log = file.path(output_dir, "ALPS_multivariate_ldsc")
)

saveRDS(multi_ldsc, file.path(output_dir, "ALPS_multivariate_ldsc.rds"))

S_matrix <- multi_ldsc$S
V_matrix <- multi_ldsc$V
rg_matrix <- cov2cor(S_matrix)

write.csv(S_matrix, file.path(output_dir, "ALPS_S_matrix.csv"), row.names = TRUE)
write.csv(V_matrix, file.path(output_dir, "ALPS_V_matrix.csv"), row.names = TRUE)
write.csv(rg_matrix, file.path(output_dir, "ALPS_rg_matrix.csv"), row.names = TRUE)

message("Step 1 finished successfully.")
message("QC table: ", file.path(output_dir, "alps_univariate_ldsc_qc.tsv"))
message("Multivariate LDSC R object: ", file.path(output_dir, "ALPS_multivariate_ldsc.rds"))
