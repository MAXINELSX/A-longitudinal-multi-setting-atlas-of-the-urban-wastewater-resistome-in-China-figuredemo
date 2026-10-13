

extra_lib <- Sys.getenv("DEMO_R_LIB", "")
if (nzchar(extra_lib) && dir.exists(extra_lib)) .libPaths(c(extra_lib, .libPaths()))

prepare_fraction_paths <- function(raw) {
  required <- c("ARG_type", "Display_ARG_type", "Genetic_compartment", "Setting", "Fraction_of_all_classified_ORFs")
  if (!identical(names(raw), required) || !nrow(raw)) stop("Expected the five processed-fraction columns.")
  for (name in required[1:4]) if (anyNA(raw[[name]]) || any(!nzchar(trimws(raw[[name]])))) stop("Missing category labels.")
  values <- raw$Fraction_of_all_classified_ORFs
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0 | values > 1)) stop("Fractions must be finite values in [0, 1].")
  if (sum(values) > 1 + 1e-5) stop("Selected rows exceed their declared global denominator.")
  if (any(!raw$Genetic_compartment %in% c("Chromosome", "Phage", "Plasmid"))) stop("Unknown genetic compartment.")
  if (any(!raw$Setting %in% c("Hospital", "WWTP", "Community", "Wet market"))) stop("Unknown setting.")
  if (anyDuplicated(raw[c("ARG_type", "Genetic_compartment", "Setting")])) stop("Duplicate fraction path.")
  raw
}

plot_fraction_paths <- function(d) {
  d <- prepare_fraction_paths(d)
  labels <- c(aminoglycoside = "Aminoglycoside", antibacterial_fatty_acid = "Others", bacitracin = "Bacitracin", beta_lactam = "Beta-lactam")
  d$ARG_label <- factor(unname(labels[d$ARG_type]), levels = unname(labels))
  if (anyNA(d$ARG_label)) stop("This compact renderer supports the packaged four ARG types.")
  d$Genetic_compartment <- factor(d$Genetic_compartment, levels = c("Chromosome", "Plasmid", "Phage"))
  d$Setting <- factor(d$Setting, levels = c("Hospital", "WWTP", "Community", "Wet market"))
  palette <- c("Aminoglycoside" = "#658CAF", "Others" = "#C9C9C9", "Bacitracin" = "#CAB481", "Beta-lactam" = "#87A957",
    "Chromosome" = "#8EAD70", "Plasmid" = "#AA94CF", "Phage" = "#E4A17F",
    "Hospital" = "#E06C75", "WWTP" = "#658CAF", "Community" = "#6EC5C1", "Wet market" = "#E9C95B")
  ggplot2::ggplot(d, ggplot2::aes(y = Fraction_of_all_classified_ORFs, axis1 = ARG_label, axis2 = Genetic_compartment, axis3 = Setting)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = ARG_label), width = .065, alpha = .38, knot.pos = .38, colour = NA) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = ggplot2::after_stat(stratum)),
      width = .065, colour = "#777777", linewidth = .18) +
    ggalluvial::stat_stratum(geom = "text", ggplot2::aes(label = ggplot2::after_stat(as.character(stratum))),
      size = 3.2, hjust = -.06) +
    ggplot2::scale_x_discrete(limits = c("ARG types", "Genetic context", "Setting"), expand = c(.18, .18)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.01, .01))) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), plot.margin = ggplot2::margin(8, 28, 8, 8))
}

run_panel_b_demo <- function(root, input_dir = file.path(root, "input"), output_dir = file.path(root, "reference_output")) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (identical(tolower(input_dir), tolower(output_dir))) stop("Input and output directories must differ.")
  packages <- c("ggplot2", "ggalluvial")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
  raw <- utils::read.delim(file.path(input_dir, "processed_fractions.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  d <- prepare_fraction_paths(raw)
  p <- plot_fraction_paths(d)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tab <- function(x, name) utils::write.table(x, file.path(output_dir, name), sep = "\t", quote = FALSE, row.names = FALSE)
  tab(d, "fraction_plot_input.tsv")
  tab(data.frame(observed_paths = nrow(d), ARG_classes = length(unique(d$ARG_type)), selected_fraction_of_all_classified_ORFs = sum(d$Fraction_of_all_classified_ORFs)), "fraction_totals.tsv")
  ggplot2::ggsave(file.path(output_dir, "fraction_flow.png"), p, width = 10, height = 6.5, dpi = 180, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "fraction_flow.pdf"), p, width = 10, height = 6.5, bg = "white")
  status <- data.frame(module = c("adapted_processed_fraction_renderer", "original_integer_count_analysis"), status = c("PASS_PROCESSED_FRACTION_DEMO", "NOT_RUN_COUNT_INPUT_UNAVAILABLE"))
  utils::write.csv(status, file.path(output_dir, "run_status.csv"), row.names = FALSE)
  utils::write.csv(data.frame(package = packages, version = vapply(packages, function(x) as.character(utils::packageVersion(x)), character(1))), file.path(output_dir, "package_versions.csv"), row.names = FALSE)
  print(status, row.names = FALSE)
  invisible(list(status = status, plot = p, paths = d))
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
  run_panel_b_demo(root, options$input, options$output)
}
