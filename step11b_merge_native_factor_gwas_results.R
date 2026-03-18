library(data.table)

root_dir <- "D:/文章/GS/GWAS"
native_dir <- file.path(root_dir, "step11_factor_gwas_native_results")
merged_dir <- file.path(native_dir, "merged")
dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)

merge_one <- function(subdir, outfile) {
  files <- list.files(file.path(native_dir, subdir), pattern = "\\.tsv$", full.names = TRUE)
  files <- sort(files)
  if (!length(files)) {
    message("No files found in ", file.path(native_dir, subdir), "; skipping.")
    return(invisible(NULL))
  }
  dt <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)
  fwrite(dt, outfile, sep = "\t")
}

merge_one("F1", file.path(merged_dir, "ALPS_F1_factorGWAS_native_merged.tsv.gz"))
merge_one("F2", file.path(merged_dir, "ALPS_F2_factorGWAS_native_merged.tsv.gz"))
merge_one("Q_SNP", file.path(merged_dir, "ALPS_QSNP_native_merged.tsv.gz"))
