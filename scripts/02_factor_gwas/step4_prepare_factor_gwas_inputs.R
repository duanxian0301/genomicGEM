library(GenomicSEM)
library(data.table)

gwas_dir <- "D:/文章/GS/GWAS"
ref_file <- "D:/文章/Mibiogen and AD/多变量LDSC/reference.1000G.maf.0.005.txt"
output_dir <- file.path(gwas_dir, "step4_factor_gwas_inputs")
input_dir <- file.path(output_dir, "input")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
setwd(gwas_dir)

files <- sort(list.files(gwas_dir, pattern = "\\.txt$", full.names = TRUE))
trait_names <- tools::file_path_sans_ext(basename(files))

if (length(files) == 0) {
  stop("No raw .txt GWAS files found in: ", gwas_dir)
}
if (!file.exists(ref_file)) {
  stop("Reference file not found: ", ref_file)
}

available_cores <- parallel::detectCores(logical = TRUE)
use_cores <- max(1, min(8, available_cores - 1))
chunk_size <- 50000
base_model <- "ALPS_refined2F"

message("Running GenomicSEM::sumstats() for factor GWAS input preparation...")
sumstats_obj <- sumstats(
  files = files,
  ref = ref_file,
  trait.names = trait_names,
  se.logit = rep(FALSE, length(files)),
  OLS = rep(FALSE, length(files)),
  linprob = rep(FALSE, length(files)),
  N = rep(NA, length(files)),
  betas = NULL,
  info.filter = 0.6,
  maf.filter = 0.01,
  keep.indel = FALSE,
  parallel = TRUE,
  cores = use_cores
)

saveRDS(sumstats_obj, file.path(output_dir, "ALPS_factorGWAS_sumstats.rds"))
fwrite(as.data.frame(sumstats_obj), file.path(output_dir, "ALPS_factorGWAS_sumstats.tsv.gz"), sep = "\t")

n_rows <- nrow(sumstats_obj)
chunk_index <- rep(seq_len(ceiling(n_rows / chunk_size)), each = chunk_size)[seq_len(n_rows)]
sumstats_list <- split(sumstats_obj, chunk_index)

chunk_manifest <- data.frame(
  chunk_id = seq_along(sumstats_list),
  file_name = paste0("GenSem_sub_", seq_along(sumstats_list), "_", base_model, ".tsv"),
  n_rows = vapply(sumstats_list, nrow, integer(1))
)

for (i in seq_along(sumstats_list)) {
  out_file <- file.path(input_dir, chunk_manifest$file_name[i])
  fwrite(sumstats_list[[i]], file = out_file, sep = "\t")
}

fwrite(chunk_manifest, file.path(output_dir, "chunk_manifest.tsv"), sep = "\t")

message("Step 4 input preparation finished successfully.")
