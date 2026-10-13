


prepare_fig3d_demo <- function(raw) {
  required <- c("Sample", "City", "Sample_Type", "bin_id", "species", "arg_burden")
  if (!all(required %in% names(raw))) stop("Missing MAG input columns.")
  if (any(!is.finite(raw$arg_burden)) || any(raw$arg_burden < 0)) stop("Invalid supplied burden values.")
  norm_setting <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c("Wetmarket", "Wet Market", "Wet_market")] <- "Wet market"
    x
  }
  meta <- unique(data.frame(original_sample = raw$Sample, City = raw$City,
                            Setting = norm_setting(raw$Sample_Type)))
  if (anyDuplicated(meta$original_sample)) stop("One source sample maps to multiple cities or settings.")
  if (anyNA(meta) || any(!nzchar(as.matrix(meta)))) stop("Incomplete embedded MAG metadata.")
  meta$Sample <- paste0("D_", meta$original_sample)
  meta$PhysicalSite <- paste("DEMO_SITE", meta$City, gsub(" ", "_", meta$Setting), sep = "_")
  visit <- ave(seq_len(nrow(meta)), meta$PhysicalSite, FUN = seq_along)
  dates <- seq(as.Date("2025-01-15"), by = "month", length.out = max(visit))
  meta$Sample_Date <- as.character(dates[visit])
  meta$SamplingMonth <- format(as.Date(meta$Sample_Date), "%Y-%m")
  meta$synthetic_metadata <- TRUE
  present <- !is.na(raw$species) & nzchar(trimws(raw$species))
  mapping <- data.frame(source_row = seq_len(nrow(raw)), original_sample = raw$Sample,
    original_bin_id = raw$bin_id,
    demo_mag_uid = paste0("D_", raw$Sample, "__DemoMAG", sprintf("%04d", seq_len(nrow(raw)))),
    included_in_species_analysis = present,
    handling = ifelse(present, "INCLUDED", "NO_SPECIES"),
    stringsAsFactors = FALSE)
  mags <- data.frame(site = norm_setting(raw$Sample_Type[present]), species = raw$species[present],
    mag_uid = mapping$demo_mag_uid[present], total_burden = raw$arg_burden[present],
    source_row = which(present), stringsAsFactors = FALSE)
  groups <- split(seq_len(nrow(mags)), interaction(mags$site, mags$species, drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(j) data.frame(
    site = mags$site[j[1]], species = mags$species[j[1]], n_mag = length(j),
    mean_burden = mean(mags$total_burden[j]), median_burden = median(mags$total_burden[j]),
    max_burden = max(mags$total_burden[j]), prop_pathogen_mags = NA_real_)))
  rownames(summary) <- NULL
  excluded <- cbind(source_row = which(!present), raw[!present, , drop = FALSE])
  list(metadata = meta, mags = mags, summary = summary, id_map = mapping, excluded = excluded)
}

run_fig3_demo <- function(root, input_dir = file.path(root, "input"), output_dir = file.path(root, "reference_output")) {
  started <- Sys.time()
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (identical(tolower(input_dir), tolower(output_dir))) stop("Input and output directories must differ.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  source(file.path(root, "audit_fig3_inputs.R"), local = TRUE)
  source(file.path(root, "upstream", "analysis_functions.R"), local = TRUE)
  packages <- c("data.table", "dplyr", "ggplot2", "multcompView", "readr", "scales")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

  suppressPackageStartupMessages(library(dplyr))
  audit <- audit_fig3_inputs(input_dir)
  audit$panel_readiness$legacy_input <- TRUE
  write_fig3_audit(audit, file.path(output_dir, "input_audit"))
  raw <- audit$tables[["mag_burden.tsv"]]
  d <- prepare_fig3d_demo(raw)
  adapted_dir <- file.path(output_dir, "adapted_fig3d_inputs_SYNTHETIC_metadata")
  dir.create(adapted_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- function(x, name) utils::write.table(x, file.path(adapted_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  tab(d$mags, "Fig3d_MAG_total_burden.tsv")
  tab(d$summary, "Fig3d_species_burden_summary.tsv")
  utils::write.csv(d$metadata, file.path(adapted_dir, "Fig3_sample_metadata.csv"), row.names = FALSE)
  utils::write.csv(d$id_map, file.path(adapted_dir, "demo_MAG_ID_map.csv"), row.names = FALSE)
  utils::write.csv(d$excluded, file.path(adapted_dir, "rows_without_species_assignment.csv"), row.names = FALSE)
  panel_d <- file.path(output_dir, "fig3d_SYNTHETIC_metadata")
  model <- run_fig3d(adapted_dir, panel_d, metadata_file = file.path(adapted_dir, "Fig3_sample_metadata.csv"), render = FALSE)
  p <- model$plot + ggplot2::theme(
    axis.text.y = ggplot2::element_text(face = 'italic', colour = '#c9332b', size = 8),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, size = 8),
    legend.title = ggplot2::element_text(size = 7),
    legend.text = ggplot2::element_text(size = 7),
    plot.margin = ggplot2::margin(2, 2, 2, 2))
  p$layers[[1]]$aes_params$size <- 4.4
  p$layers[[2]]$aes_params$size <- 2.4
  ggplot2::ggsave(file.path(panel_d, "Fig3d_SYNTHETIC_metadata_demo.pdf"), p, width = 5.1, height = 5.1)
  ggplot2::ggsave(file.path(panel_d, "Fig3d_SYNTHETIC_metadata_demo.png"), p, width = 5.1, height = 5.1, dpi = 180)

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  for (panel in c("a", "b", "c", "e")) {
    module <- paste0("panel_", panel, "_demo")
    result <- system2(rscript, c("--vanilla", shQuote(file.path(root, module, "run_demo.R")),
      "--output", shQuote(file.path(output_dir, module))))
    if (result != 0L) stop("Included processed-table demo failed: ", panel)
  }
  status <- data.frame(panel = letters[1:5],
    module = c("relative_abundance_display", "precomputed_fraction_display", "contig_sample_and_summary_displays",
      "original_Fig3d_functions", "precomputed_MAG_cell_display"),
    status = c(rep("PASS_DISPLAY_ONLY", 3), "PASS_SYNTHETIC_METADATA_DEMO", "PASS_DISPLAY_ONLY"),
    complete_manuscript_reproduction = FALSE, stringsAsFactors = FALSE)
  utils::write.csv(status, file.path(output_dir, "run_status.csv"), row.names = FALSE)
  utils::write.csv(data.frame(package = packages, version = vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1))),
    file.path(output_dir, "package_versions.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo())[!grepl("time zone:|tzcode source:", capture.output(sessionInfo()), ignore.case=TRUE)], file.path(output_dir, "sessionInfo.txt"))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  utils::write.csv(data.frame(elapsed_seconds = elapsed), file.path(output_dir, "run_time.csv"), row.names = FALSE)
  stopifnot(identical(audit$input_hashes$md5, unname(tools::md5sum(file.path(input_dir, audit$input_hashes$file)))))
  print(status)
  cat(sprintf("PASS_FIG3 elapsed_seconds=%.2f\n", elapsed))
  invisible(list(status = status, figure3d = model))
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(input = file.path(root, "input"), output = file.path(root, "reference_output"))
  if (length(args) %% 2L) stop("Usage: Rscript run_demo.R [--input DIR] [--output DIR]")
  if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
    if (!args[i] %in% c("--input", "--output")) stop("Unknown option: ", args[i])
    opts[[sub("^--", "", args[i])]] <- args[i + 1L]
  }
  run_fig3_demo(root, opts$input, opts$output)
}
