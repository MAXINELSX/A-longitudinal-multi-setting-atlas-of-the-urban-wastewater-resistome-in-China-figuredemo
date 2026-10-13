

read_demo_tsv <- function(path, required) {
  if (!file.exists(path)) stop("Missing input: ", basename(path), call. = FALSE)
  d <- read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                  colClasses = "character", quote = "", comment.char = "", fileEncoding = "UTF-8")
  absent <- setdiff(required, names(d))
  if (length(absent)) stop("Missing columns in ", basename(path), ": ", paste(absent, collapse = ", "))
  if (!nrow(d)) stop("Empty input: ", basename(path))
  d
}

write_demo_tsv <- function(d, path) {
  old <- options(digits = 17)
  on.exit(options(old), add = TRUE)
  write.table(d, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

save_demo_plot <- function(p, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height,
                  device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 160, bg = "white")
}

prepare_fig4d_demo <- function(path) {
  numeric_fields <- c("ANI_min", "fixed_snp_per_mbp_pair", "callable_fraction_pair",
                      "clinical_ARG_count", "mag_ARG_count", "shared_ARG_count")
  flag_fields <- c("direct_remap_supported_bool", "instrain_partial_support_bool",
                   "shared_carbapenemase_present_bool")
  required <- c("Pair_ID", "Clinical_ID", "MAG_ID", "setting", numeric_fields, flag_fields)
  d <- read_demo_tsv(path, required)[, required, drop = FALSE]
  if (anyNA(d$Pair_ID) || any(!nzchar(d$Pair_ID)) || anyDuplicated(d$Pair_ID)) stop("Unique Pair_ID values are required.")
  for (column in numeric_fields) {
    d[[column]] <- suppressWarnings(as.numeric(d[[column]]))
    if (any(!is.finite(d[[column]]))) stop("Non-finite values in ", column)
  }
  for (column in flag_fields) {
    if (anyNA(d[[column]]) || any(!d[[column]] %in% c("True", "False", ""))) {
      stop("Invalid supplied Boolean flag in ", column)
    }
  }
  if (any(d$fixed_snp_per_mbp_pair <= 0) || any(d$callable_fraction_pair < 0 | d$callable_fraction_pair > 1) ||
      any(d$ANI_min < 0 | d$ANI_min > 100)) stop("Invalid relatedness display values.")
  for (column in c("clinical_ARG_count", "mag_ARG_count", "shared_ARG_count")) {
    if (any(d[[column]] < 0 | d[[column]] != floor(d[[column]]))) stop("Invalid count in ", column)
  }
  if (any(d$shared_ARG_count > pmin(d$clinical_ARG_count, d$mag_ARG_count))) stop("Shared ARG count exceeds source counts.")
  if (any(!d$setting %in% c("Community", "Hospital", "Wet market"))) stop("Unexpected setting.")
  d
}

run_fig4d_demo <- function(path, output_dir, original) {
  d <- prepare_fig4d_demo(path)
  display <- data.frame(
    Pair_ID = d$Pair_ID,
    Clinical_display_ID = d$Clinical_ID,
    Clinical_specimen = rep(c("Sputum", "Oropharyngeal swab"), length.out = nrow(d)),
    Site = sprintf("Site%02d", seq_len(nrow(d))),
    Clinical_city = rep(c("DemoA", "DemoB"), length.out = nrow(d)),
    Wastewater_city = rep(c("DemoB", "DemoB"), length.out = nrow(d)),
    same_city = rep(c("No", "Yes"), length.out = nrow(d)),
    Wastewater_setting = d$setting,
    ANI_min = d$ANI_min,
    fixed_SNP_Mbp = d$fixed_snp_per_mbp_pair,
    callable_fraction = d$callable_fraction_pair,
    Direct_remap = as.integer(d$direct_remap_supported_bool == "True"),
    inStrain_partial = as.integer(d$instrain_partial_support_bool == "True"),
    Clinical_ARG_count = d$clinical_ARG_count,
    Wastewater_ARG_count = d$mag_ARG_count,
    Shared_ARG_count = d$shared_ARG_count,
    High_priority_ARG = rep("", nrow(d)),
    stringsAsFactors = FALSE
  )
  original_env <- environment(original$.fig4d$run)
  normalized <- get("normalize_pairs", envir = original_env)(display)
  draw <- get("draw_fig4d", envir = original_env)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(output_dir, "Fig4d_partial_demo")
  grDevices::png(paste0(stem, ".png"), width = 14.5, height = 9.3,
                 units = "in", res = 180, bg = "white", type = "cairo")
  tryCatch(draw(normalized, "bottom"), finally = grDevices::dev.off())
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 14.5, height = 9.3, bg = "white")
  tryCatch(draw(normalized, "bottom"), finally = grDevices::dev.off())
  write_demo_tsv(d, file.path(output_dir, "Fig4d_displayed_pairs.tsv"))
  list(data = d, display = normalized, plot = NULL)
}

prepare_fig4e_demo <- function(input_dir) {
  evidence <- read_demo_tsv(file.path(input_dir, "e_carbapenemase_overlap.tsv"), c("Evidence", "Allele", "Supported", "Support_score"))
  context <- read_demo_tsv(file.path(input_dir, "e_mobile_genetic_context.tsv"), c("Context", "Contigs", "Percent"))
  taxa <- read_demo_tsv(file.path(input_dir, "taxon_counts.tsv"), c("Taxon_group", "Contigs", "Percent_of_classified"))
  coverage <- read_demo_tsv(file.path(input_dir, "e_classification_coverage.tsv"), c("Category", "Contigs", "Percent_of_all_targets"))
  evidence$Support_score <- suppressWarnings(as.numeric(evidence$Support_score))
  keys <- paste(evidence$Evidence, evidence$Allele, sep = "|")
  expected <- as.vector(outer(c("Clinical dataset", "Wastewater dataset", "Pair-level shared", "Candidate KP MAG host"),
                             c("KPC-2", "NDM-5", "OXA-48"), paste, sep = "|"))
  if (anyDuplicated(keys) || !setequal(keys, expected) ||
      any(!evidence$Supported %in% c("Yes", "No")) || anyNA(evidence$Support_score) ||
      any(evidence$Support_score != as.integer(evidence$Supported == "Yes"))) stop("Invalid supplied evidence records.")
  validate_summary <- function(d, column) {
    d$Contigs <- suppressWarnings(as.numeric(d$Contigs))
    d[[column]] <- suppressWarnings(as.numeric(d[[column]]))
    if (any(!is.finite(d$Contigs)) || any(d$Contigs < 0 | d$Contigs != floor(d$Contigs)) || sum(d$Contigs) <= 0) stop("Invalid contig counts.")
    d$Percent_recalculated <- 100 * d$Contigs / sum(d$Contigs)
    if (any(!is.finite(d[[column]])) || any(abs(d[[column]] - d$Percent_recalculated) > 0.051)) stop("Rounded percentage disagrees with its summary denominator.")
    d
  }
  context <- validate_summary(context, "Percent")
  taxa <- validate_summary(taxa, "Percent_of_classified")
  coverage <- validate_summary(coverage, "Percent_of_all_targets")
  expected_coverage <- c("Classified / host-table matched", "No host-table match")
  if (nrow(coverage) != 2L || !setequal(coverage$Category, expected_coverage)) stop("Unexpected classification coverage categories.")
  if (sum(taxa$Contigs) != coverage$Contigs[coverage$Category == expected_coverage[1]]) stop("Taxonomy total differs from classified count.")
  if (sum(context$Contigs) != sum(coverage$Contigs)) stop("Context total differs from coverage total.")
  list(evidence = evidence, context = context, taxa = taxa, coverage = coverage,
       denominators = c(context = sum(context$Contigs), taxonomy = sum(taxa$Contigs), coverage = sum(coverage$Contigs)))
}

run_fig4e_demo <- function(input_dir, output_dir, original) {
  d <- prepare_fig4e_demo(input_dir)
  result <- original$.fig4e$run(input_dir, output_dir)
  save_demo_plot(result$plot, file.path(output_dir, "Fig4e"), 5.1, 8.2)
  write_demo_tsv(d$evidence, file.path(output_dir, "Fig4e_evidence_all.tsv"))
  write_demo_tsv(result$data$evidence, file.path(output_dir, "Fig4e_evidence_displayed.tsv"))
  write_demo_tsv(d$context, file.path(output_dir, "Fig4e_context_recalculated.tsv"))
  write_demo_tsv(d$taxa, file.path(output_dir, "Fig4e_taxa_recalculated.tsv"))
  write_demo_tsv(d$coverage, file.path(output_dir, "Fig4e_classification_coverage.tsv"))
  write_demo_tsv(data.frame(Summary = names(d$denominators), Contigs = unname(d$denominators)),
                 file.path(output_dir, "Fig4e_denominators.tsv"))
  list(data = d, plot = result$plot)
}

run_fig4_demo <- function(root, input_dir = file.path(root, "input"), output_dir = file.path(root, "reference_output")) {
  start <- proc.time()[["elapsed"]]
  optional_library <- Sys.getenv("DEMO_R_LIB")
  if (nzchar(optional_library)) .libPaths(c(optional_library, .libPaths()))
  dependencies <- c("ggplot2", "patchwork")
  missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))
  if (!capabilities("cairo")) stop("Cairo support is required for PDF output.")
  original <- new.env(parent = globalenv())
  sys.source(file.path(root, "source", "analysis_functions.R"), envir = original)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  c_result <- original$.fig4c$run(file.path(input_dir, "cluster_counts.tsv"), file.path(output_dir, "c"))
  c_result$plot <- c_result$plot +
    ggplot2::labs(x = "Genomes or MAGs  |  Wastewater members within ANI cluster") +
    ggplot2::theme(axis.title.x = ggplot2::element_text(size = 8, face = "bold",
                                                     margin = ggplot2::margin(t = 9)))
  save_demo_plot(c_result$plot, file.path(output_dir, "c", "Fig4c"), 5.3, 4.3)
  d_result <- run_fig4d_demo(file.path(input_dir, "pair_demo.tsv"), file.path(output_dir, "d_partial"), original)
  e_result <- run_fig4e_demo(input_dir, file.path(output_dir, "e"), original)
  session <- capture.output(sessionInfo())
  writeLines(session[!grepl("time zone:|tzcode source:", session, ignore.case = TRUE)], file.path(output_dir, "sessionInfo.txt"))
  elapsed <- proc.time()[["elapsed"]] - start
  manifest <- data.frame(module = c("Fig4c", "Fig4d_partial", "Fig4e_evidence", "Fig4e_context", "Fig4e_taxonomy", "Fig4e_coverage"),
    input_rows = c(nrow(c_result$data), nrow(d_result$data), 12, 3, 7, 2),
    plotted_rows = c(nrow(c_result$data), nrow(d_result$data), 8, 3, 7, 0),
    mode = c("cluster_summary", "pair_subset", "boolean_fixture",
             "simulated_context", "simulated_taxonomy", "simulated_coverage"))
  write_demo_tsv(manifest, file.path(output_dir, "run_manifest.tsv"))
  writeLines(sprintf("elapsed_seconds\t%.3f", elapsed), file.path(output_dir, "runtime.txt"))
  message(sprintf("PASS\tFig4\telapsed_seconds=%.3f", elapsed))
  invisible(list(c = c_result, d = d_result, e = e_result, elapsed = elapsed))
}

if (sys.nframe() == 0L) {
  full_args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", full_args[grepl("^--file=", full_args)][1])
  root <- dirname(normalizePath(script, winslash = "/"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 2L && identical(args[[1]], "--output")) args <- c(file.path(root, "input"), args[[2]])
  if (length(args) > 2L || any(args == "--output")) stop("Usage: Rscript run_demo.R [--output OUTPUT_DIRECTORY] or [INPUT_DIRECTORY OUTPUT_DIRECTORY]")
  input <- if (length(args) >= 1L) args[[1]] else file.path(root, "input")
  output <- if (length(args) >= 2L) args[[2]] else file.path(root, "reference_output")
  run_fig4_demo(root, input, output)
}
