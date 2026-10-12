

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

run_fig4d_demo <- function(path, output_dir) {
  d <- prepare_fig4d_demo(path)
  plot_data <- d
  plot_data$Pair_plot <- factor(d$Pair_ID, levels = rev(d$Pair_ID))
  plot_data$setting <- factor(d$setting, levels = c("Community", "Hospital", "Wet market"))
  palette <- c(Community = "#86BFC0", Hospital = "#C56D72", `Wet market` = "#E5BE4D")
  theme <- ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   plot.title = ggplot2::element_text(face = "bold"))
  p_snp <- ggplot2::ggplot(plot_data, ggplot2::aes(fixed_snp_per_mbp_pair, Pair_plot)) +
    ggplot2::geom_point(ggplot2::aes(fill = setting, size = callable_fraction_pair), shape = 21) +
    ggplot2::scale_x_log10() + ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::scale_size_continuous(limits = c(0, 1), breaks = c(0.5, 0.7, 0.9), range = c(1, 7), name = "Callable fraction") +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, order = 1), size = ggplot2::guide_legend(nrow = 1, order = 2)) +
    ggplot2::labs(x = "Fixed SNPs/Mbp (log scale)", y = "Pair", fill = "Setting", title = "Relatedness") + theme
  count_columns <- c("clinical_ARG_count", "mag_ARG_count", "shared_ARG_count")
  counts <- do.call(rbind, lapply(count_columns, function(column) {
    data.frame(Pair_ID = d$Pair_ID, Category = column, Count = d[[column]], stringsAsFactors = FALSE)
  }))
  counts$Pair_plot <- factor(counts$Pair_ID, levels = rev(d$Pair_ID))
  counts$Category <- factor(counts$Category, levels = count_columns, labels = c("Clinical", "MAG", "Shared"))
  p_counts <- ggplot2::ggplot(counts, ggplot2::aes(Count, Pair_plot, colour = Category)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 3) +
    ggplot2::scale_colour_manual(values = c(Clinical = "#B78387", MAG = "#6D82AC", Shared = "#4F8888")) +
    ggplot2::labs(x = "ARG count", y = NULL, title = "ARG counts") + theme
  flags <- c("direct_remap_supported_bool", "instrain_partial_support_bool", "shared_carbapenemase_present_bool")
  support <- do.call(rbind, lapply(flags, function(column) {
    data.frame(Pair_ID = d$Pair_ID, Field = column, Value = d[[column]], stringsAsFactors = FALSE)
  }))
  support$Pair_plot <- factor(support$Pair_ID, levels = rev(d$Pair_ID))
  support$Field <- factor(support$Field, levels = flags,
    labels = c("Direct\nremap", "inStrain\npartial", "Shared\ncarbapenemase"))
  support$Value[support$Value == ""] <- "Missing"
  p_flags <- ggplot2::ggplot(support, ggplot2::aes(Field, Pair_plot, shape = Value)) +
    ggplot2::geom_point(size = 3.5, colour = "#476878") +
    ggplot2::scale_shape_manual(values = c(True = 16, False = 1, Missing = 4), name = "Flag") +
    ggplot2::labs(x = NULL, y = NULL, title = "Flags") + theme
  p <- patchwork::wrap_plots(p_snp, p_counts, p_flags, ncol = 3, widths = c(1.2, 1, 1.15)) +
    patchwork::plot_annotation(title = "Figure 4d")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  save_demo_plot(p, file.path(output_dir, "Fig4d_partial_demo"), 13, 5.5)
  write_demo_tsv(d, file.path(output_dir, "Fig4d_displayed_pairs.tsv"))
  list(data = d, plot = p)
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
  result$plot <- result$plot + patchwork::plot_annotation(title = "Figure 4e")
  save_demo_plot(result$plot, file.path(output_dir, "Fig4e"), 7, 9.2)
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
    ggplot2::labs(title = "Figure 4c", x = "Genomes / MAGs | Wastewater (%)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"),
                   axis.title.x = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 10)),
                   plot.margin = ggplot2::margin(12, 24, 12, 12))
  save_demo_plot(c_result$plot, file.path(output_dir, "c", "Fig4c"), 11, 5.5)
  d_result <- run_fig4d_demo(file.path(input_dir, "pair_demo.tsv"), file.path(output_dir, "d_partial"))
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
