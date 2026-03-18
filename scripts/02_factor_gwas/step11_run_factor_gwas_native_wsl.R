.libPaths(c("/home/shenjing/R/genomicsem_fix_lib", .libPaths()))

library(GenomicSEM)
library(data.table)

root_dir <- "/home/shenjing/gs_paths/gwas"
input_dir <- file.path(root_dir, "step4_factor_gwas_inputs", "input")
manifest_file <- file.path(root_dir, "step4_factor_gwas_inputs", "chunk_manifest.tsv")
ldsc_file <- file.path(root_dir, "step3_cfa_results", "ALPS5_ALL_CFA_ldsc.rds")
output_dir <- file.path(root_dir, "step11_factor_gwas_native_results")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "F1"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "F2"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Q_SNP"), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(manifest_file)) {
  stop("Missing chunk manifest: ", manifest_file)
}
if (!file.exists(ldsc_file)) {
  stop("Missing LDSC object: ", ldsc_file)
}

manifest <- fread(manifest_file)
LDSCoutput <- readRDS(ldsc_file)

start_chunk <- as.integer(Sys.getenv("START_CHUNK", unset = "1"))
end_chunk <- as.integer(Sys.getenv("END_CHUNK", unset = as.character(nrow(manifest))))
use_cores <- as.integer(Sys.getenv("USE_CORES", unset = "4"))

if (is.na(start_chunk) || is.na(end_chunk) || is.na(use_cores)) {
  stop("START_CHUNK, END_CHUNK, and USE_CORES must be integers.")
}

start_chunk <- max(1L, start_chunk)
end_chunk <- min(nrow(manifest), end_chunk)
use_cores <- max(1L, use_cores)

model <- "
F1 =~ aALPS + Left_ALPS + Right_ALPS
F2 =~ aALPS + mALPS + pALPS
F1 ~~ F2
Left_ALPS ~~ lvar*Left_ALPS
mALPS ~~ mvar*mALPS
lvar > 0.0001
mvar > 0.0001
F1 ~ SNP
F2 ~ SNP
"

writeLines(model, file.path(output_dir, "ALPS_refined_factorGWAS_model_native.txt"))

for (i in seq.int(start_chunk, end_chunk)) {
  chunk_file <- file.path(input_dir, manifest$file_name[i])
  if (!file.exists(chunk_file)) {
    message("Skipping missing file: ", chunk_file)
    next
  }

  f1_out <- file.path(output_dir, "F1", paste0("F1_", i, ".tsv"))
  f2_out <- file.path(output_dir, "F2", paste0("F2_", i, ".tsv"))
  q_out <- file.path(output_dir, "Q_SNP", paste0("Q_SNP_", i, ".tsv"))
  # Treat a chunk as completed once both factor outputs exist.
  # Q_SNP is optional here because earlier successful runs did not always
  # persist a paired Q_SNP file for every completed chunk.
  if (file.exists(f1_out) && file.exists(f2_out)) {
    message("Skipping completed chunk ", i)
    next
  }

  message("Running native chunk ", i, " / ", nrow(manifest), ": ", manifest$file_name[i])
  snp_chunk <- fread(chunk_file)

  result <- tryCatch(
    userGWAS(
      covstruc = LDSCoutput,
      SNPs = snp_chunk,
      estimation = "DWLS",
      model = model,
      printwarn = TRUE,
      sub = c("F1~SNP", "F2~SNP"),
      cores = use_cores,
      parallel = use_cores > 1,
      GC = "none",
      MPI = FALSE,
      smooth_check = TRUE,
      fix_measurement = TRUE,
      Q_SNP = TRUE
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    message("Chunk failed: ", manifest$file_name[i], " | ", conditionMessage(result))
    next
  }

  fwrite(as.data.frame(result[[1]]), f1_out, sep = "\t")
  fwrite(as.data.frame(result[[2]]), f2_out, sep = "\t")

  if (length(result) >= 3 && !is.null(result[[3]])) {
    fwrite(as.data.frame(result[[3]]), q_out, sep = "\t")
  }
}

message("Native factor GWAS loop finished.")
