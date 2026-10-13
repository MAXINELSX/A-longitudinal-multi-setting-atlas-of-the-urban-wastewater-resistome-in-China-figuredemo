read_nonempty_tsv <- function(path) {
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "", encoding = "UTF-8")
  keep <- if (nrow(x)) apply(x, 1, function(z) any(!is.na(z) & nzchar(trimws(as.character(z))))) else logical()
  attr(x, "empty_rows_removed") <- sum(!keep)
  out <- x[keep, , drop = FALSE]
  attr(out, "empty_rows_removed") <- sum(!keep)
  out
}

normalize_demo_setting <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("Wetmarket", "Wet Market", "Wet_market")] <- "Wet market"
  x
}

audit_fig3_inputs <- function(input_dir) {
  files <- c("species_types.tsv", "mag_burden.tsv",
             "contig_fractions.tsv", "species_abundance.tsv",
             "species_counts.tsv", "sample_metadata.tsv")
  paths <- file.path(input_dir, files)
  if (!all(file.exists(paths))) stop("Missing supplied inputs: ", paste(files[!file.exists(paths)], collapse = ", "))
  tables <- setNames(lapply(paths, read_nonempty_tsv), files)
  inventory <- data.frame(file = files, nonempty_rows = vapply(tables, nrow, integer(1)),
    columns = vapply(tables, ncol, integer(1)),
    empty_rows_removed = vapply(tables, function(x) attr(x, "empty_rows_removed"), integer(1)),
    row.names = NULL)
  mags <- tables[["mag_burden.tsv"]]
  mags$mag_uid <- paste(mags$Sample, mags$bin_id, sep = "__")
  duplicate_mag_records <- mags[duplicated(mags$mag_uid) | duplicated(mags$mag_uid, fromLast = TRUE), ]
  points <- tables[["contig_fractions.tsv"]]
  meta <- tables[["sample_metadata.tsv"]]
  j <- match(points$sample, meta$sample)
  metadata_comparison <- data.frame(sample = points$sample, city_in_fraction_table = points$City,
    city_in_heatmap_metadata = meta$City[j], setting_in_fraction_table = points$Sample_Type,
    setting_in_heatmap_metadata = meta$Sample_Type[j], city_equal = points$City == meta$City[j],
    setting_equal = normalize_demo_setting(points$Sample_Type) == normalize_demo_setting(meta$Sample_Type[j]))
  panel_readiness <- data.frame(panel = letters[1:5], status = "BLOCKED",
    reason_code = c("MISSING_SPECIES_COUNTS", "MISSING_CLASSIFIED_COUNTS", "CONTIG_UNIT_MISMATCH",
      "MISSING_SITE_DATE_DUPLICATE_IDS", "MISSING_MAG_REGISTRY"), stringsAsFactors = FALSE)
  list(inventory = inventory, duplicate_mag_records = duplicate_mag_records,
       metadata_comparison = metadata_comparison, panel_readiness = panel_readiness,
       input_hashes = data.frame(file = files, md5 = unname(tools::md5sum(paths))),
       tables = tables)
}

write_fig3_audit <- function(result, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in setdiff(names(result), "tables"))
    utils::write.csv(result[[name]], file.path(out_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) stop("Usage: Rscript audit_fig3_inputs.R INPUT_DIR OUTPUT_DIR")
  result <- audit_fig3_inputs(args[[1]])
  write_fig3_audit(result, args[[2]])
  print(result$panel_readiness)
  cat("PASS_INPUT_AUDIT\n")
}
