

extra_lib <- Sys.getenv("DEMO_R_LIB", "")
if (nzchar(extra_lib) && dir.exists(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
setting_order <- c("Hospital", "WWTP", "Community", "Wet market")
setting_palette <- c(Hospital = "#E15759", WWTP = "#4E79A7", Community = "#76B7B2", "Wet market" = "#E3BA22")

prepare_summary <- function(raw) {
  required <- c("Setting", "N_samples", "Mean", "Median", "Q1", "Q3", "Minimum", "Maximum", "SD")
  if (!identical(names(raw), required) || nrow(raw) != 4L || !setequal(raw$Setting, setting_order) || anyDuplicated(raw$Setting)) stop("Expected exactly four setting-summary rows.")
  if (!all(vapply(raw[-1], is.numeric, logical(1))) || any(!is.finite(as.matrix(raw[-1])))) stop("Summary values must be finite numeric values.")
  if (any(raw$N_samples <= 0 | raw$N_samples != floor(raw$N_samples))) stop("Invalid reported sample sizes.")
  fraction_cols <- c("Mean", "Median", "Q1", "Q3", "Minimum", "Maximum")
  if (any(as.matrix(raw[fraction_cols]) < 0 | as.matrix(raw[fraction_cols]) > 1) || any(raw$SD < 0)) stop("Invalid reported fractions or SD.")
  if (any(raw$Minimum > raw$Q1 | raw$Q1 > raw$Median | raw$Median > raw$Q3 | raw$Q3 > raw$Maximum | raw$Mean < raw$Minimum | raw$Mean > raw$Maximum)) stop("Inconsistent summary ordering.")
  raw
}

prepare_contigs <- function(raw) {
  counts <- c("total_arg_contigs", "mobile_arg_contigs", "plasmid_arg_contigs", "virus_arg_contigs", "chromosome_arg_contigs")
  fractions <- c("mobile_fraction", "plasmid_fraction", "virus_fraction", "chromosome_fraction")
  required <- c("Sample", "Setting", counts, fractions)
  if (!identical(names(raw), required) || !nrow(raw) || anyNA(raw$Sample) || any(!nzchar(raw$Sample)) || anyDuplicated(raw$Sample)) stop("Expected unique sample rows and contig columns.")
  if (anyNA(raw$Setting) || !setequal(raw$Setting, setting_order) || any(!raw$Setting %in% setting_order)) stop("Expected the four known settings.")
  if (!all(vapply(raw[c(counts, fractions)], is.numeric, logical(1))) || any(!is.finite(as.matrix(raw[c(counts, fractions)])))) stop("Invalid numeric contig data.")
  values <- as.matrix(raw[counts])
  if (any(values < 0 | values != floor(values)) || any(raw$total_arg_contigs <= 0)) stop("Contig counts must be nonnegative integers with positive totals.")
  if (any(raw$mobile_arg_contigs != raw$plasmid_arg_contigs + raw$virus_arg_contigs) || any(raw$total_arg_contigs != raw$mobile_arg_contigs + raw$chromosome_arg_contigs)) stop("Contig counts do not reconcile.")
  if (any(as.matrix(raw[fractions]) < 0 | as.matrix(raw[fractions]) > 1)) stop("Contig fractions must be in [0, 1].")
  expected <- as.matrix(raw[counts[-1]]) / raw$total_arg_contigs
  if (any(abs(as.matrix(raw[fractions]) - expected) > 1e-8)) stop("Supplied contig fractions disagree with the counts beyond rounding tolerance.")
  if (any(table(factor(raw$Setting, levels = setting_order)) < 2L)) stop("A violin illustration requires at least two supplied rows per setting.")
  raw
}

plot_setting_summary <- function(d) {
  d <- prepare_summary(d)
  d$Setting <- factor(d$Setting, levels = setting_order)
  ggplot2::ggplot(d, ggplot2::aes(x = Setting, y = Median, colour = Setting)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = Q1, ymax = Q3), linewidth = 3.2, alpha = .6) +
    ggplot2::geom_point(ggplot2::aes(shape = "Median"), size = 3, fill = "white", stroke = .85) +
    ggplot2::geom_point(ggplot2::aes(y = Mean, shape = "Mean"), size = 3, stroke = 1) +
    ggplot2::scale_colour_manual(values = setting_palette, guide = "none") +
    ggplot2::scale_shape_manual(values = c(Median = 21, Mean = 4), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::labs(title = "Figure 3c", x = NULL, y = "Mobile ARG-carrying contig fraction") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))
}

plot_contig_observations <- function(d) {
  d <- prepare_contigs(d)
  d$Setting <- factor(d$Setting, levels = setting_order)
  ggplot2::ggplot(d, ggplot2::aes(x = Setting, y = mobile_fraction, fill = Setting)) +
    ggplot2::geom_violin(trim = TRUE, scale = "width", alpha = .45, colour = NA, width = .8) +
    ggplot2::geom_boxplot(width = .15, outlier.shape = NA, fill = "white", linewidth = .5) +
    ggplot2::geom_point(ggplot2::aes(colour = Setting), position = ggplot2::position_jitter(width = .07, height = 0, seed = 31), size = 2.3) +
    ggplot2::scale_fill_manual(values = setting_palette, guide = "none") +
    ggplot2::scale_colour_manual(values = setting_palette, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::labs(title = "Figure 3c", x = NULL, y = "Mobile ARG-carrying contig fraction") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold"))
}

summarize_contigs <- function(contigs) {
  contigs <- prepare_contigs(contigs)
  do.call(rbind, lapply(setting_order, function(s) {
    values <- contigs$mobile_fraction[contigs$Setting == s]
    data.frame(Setting = s, N_samples = length(values), Mean = mean(values), Median = stats::median(values),
      Q1 = unname(stats::quantile(values, .25)), Q3 = unname(stats::quantile(values, .75)),
      Minimum = min(values), Maximum = max(values), SD = stats::sd(values))
  }))
}

run_panel_c_demo <- function(root, input_dir = file.path(root, "input"), output_dir = file.path(root, "reference_output")) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (identical(tolower(input_dir), tolower(output_dir))) stop("Input and output directories must differ.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing package: ggplot2")
  supplied_summary <- prepare_summary(utils::read.delim(file.path(input_dir, "setting_summary.tsv"), check.names = FALSE, stringsAsFactors = FALSE))
  contigs <- prepare_contigs(utils::read.delim(file.path(input_dir, "contig_observations.tsv"), check.names = FALSE, stringsAsFactors = FALSE))
  summary <- summarize_contigs(contigs)
  if (!isTRUE(all.equal(supplied_summary, summary, tolerance = 1e-12, check.attributes = FALSE))) stop("SUMMARY_MISMATCH")
  ps <- plot_setting_summary(summary)
  pc <- plot_contig_observations(contigs)
  descriptive <- summary
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- function(x, name) utils::write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
  tab(summary, "setting_summary_plot_input.tsv")
  tab(contigs, "contig_plot_input.tsv")
  tab(descriptive, "contig_descriptive_statistics.tsv")
  for (extension in c("png", "pdf")) {
    ggplot2::ggsave(file.path(output_dir, paste0("setting_summary.", extension)), ps, width = 8, height = 5.8, dpi = 180, bg = "white")
    ggplot2::ggsave(file.path(output_dir, paste0("contig_violin.", extension)), pc, width = 8, height = 5.8, dpi = 180, bg = "white")
  }
  status <- data.frame(module = c("contig_summary_display", "contig_sample_display", "annotation_LMM"), status = c("PASS_SUMMARY_DISPLAY_ONLY", "PASS_CONTIG_PLOT_ONLY", "NOT_RUN_ANNOTATION_INPUT_UNAVAILABLE"))
  utils::write.csv(status, file.path(output_dir, "run_status.csv"), row.names = FALSE)
  utils::write.csv(data.frame(package = "ggplot2", version = as.character(utils::packageVersion("ggplot2"))), file.path(output_dir, "package_versions.csv"), row.names = FALSE)
  print(status, row.names = FALSE)
  invisible(list(status = status, summary_plot = ps, contig_plot = pc, descriptive = descriptive))
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
  args <- commandArgs(TRUE)
  options <- list(input = file.path(root, "input"), output = file.path(root, "reference_output"))
  if (length(args) %% 2L) stop("Usage: Rscript run_demo.R [--input DIR] [--output DIR]")
  if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
    if (!args[i] %in% c("--input", "--output")) stop("Unknown option.")
    options[[sub("^--", "", args[i])]] <- args[i + 1L]
  }
  run_panel_c_demo(root, options$input, options$output)
}
