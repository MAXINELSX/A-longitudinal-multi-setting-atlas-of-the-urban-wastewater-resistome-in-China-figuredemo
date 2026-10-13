

verify_saved_fig4 <- function(root, output_dir = file.path(root, "reference_output")) {
  read <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", quote = "", comment.char = "")
  near <- function(actual, expected) isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12))
  expected_files <- c("c/Fig4c.pdf", "c/Fig4c.png", "c/Fig4c_plot_data.tsv",
    "d_partial/Fig4d_partial_demo.pdf", "d_partial/Fig4d_partial_demo.png", "d_partial/Fig4d_displayed_pairs.tsv",
    "e/Fig4e.pdf", "e/Fig4e.png", "e/Fig4e_evidence_all.tsv", "e/Fig4e_evidence_displayed.tsv",
    "e/Fig4e_context_recalculated.tsv", "e/Fig4e_taxa_recalculated.tsv",
    "e/Fig4e_classification_coverage.tsv", "e/Fig4e_denominators.tsv",
    "run_manifest.tsv", "runtime.txt", "sessionInfo.txt")
  paths <- file.path(output_dir, expected_files)
  stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0))
  c_data <- read(file.path(output_dir, "c/Fig4c_plot_data.tsv"))
  stopifnot(nrow(c_data) == 3L, near(c_data$Clinical, c(1, 1, 3)),
            near(c_data$Wastewater, c(1, 1, 1)), near(c_data$WW_percent, c(50, 50, 25)))
  source_pairs <- read(file.path(root, "input/pair_demo.tsv"))
  d <- read(file.path(output_dir, "d_partial/Fig4d_displayed_pairs.tsv"))
  numeric_fields <- c("ANI_min", "fixed_snp_per_mbp_pair", "callable_fraction_pair",
                      "clinical_ARG_count", "mag_ARG_count", "shared_ARG_count")
  stopifnot(nrow(d) == 4L, identical(d$Pair_ID, sprintf("Pair%02d", 1:4)),
            near(d$ANI_min, c(99.6987, 99.5089, 99.4822, 99.4938)),
            near(d$clinical_ARG_count, c(11, 28, 26, 37)),
            near(d$mag_ARG_count, c(8, 11, 10, 8)), near(d$shared_ARG_count, c(8, 9, 9, 8)),
            setequal(names(d), names(source_pairs)),
            !any(c("Clinical_specimen", "High_priority_ARG", "displayed_in_Fig4d") %in% names(d)))
  for (column in names(source_pairs)) {
    if (column %in% numeric_fields) stopifnot(near(d[[column]], source_pairs[[column]]))
    else stopifnot(identical(d[[column]], source_pairs[[column]]))
  }
  evidence <- read(file.path(output_dir, "e/Fig4e_evidence_all.tsv"))
  source_evidence <- read(file.path(root, "input/e_carbapenemase_overlap.tsv"))
  displayed <- read(file.path(output_dir, "e/Fig4e_evidence_displayed.tsv"))
  stopifnot(identical(evidence, source_evidence), nrow(evidence) == 12L, nrow(displayed) == 8L,
            setequal(displayed$Allele, c("KPC-2", "OXA-48")),
            near(evidence$Support_score, as.integer(evidence$Supported == "Yes")))
  check_summary <- function(input, output, total, percentage) {
    source <- read(file.path(root, "input", input))
    saved <- read(file.path(output_dir, "e", output))
    stopifnot(nrow(source) == nrow(saved), identical(source[[1]], saved[[1]]),
              near(saved$Contigs, source$Contigs), sum(as.numeric(saved$Contigs)) == total,
              near(saved[[percentage]], source[[percentage]]),
              all(abs(as.numeric(saved[[percentage]]) - 100 * as.numeric(saved$Contigs) / total) <= 0.051),
              near(saved$Percent_recalculated, 100 * as.numeric(saved$Contigs) / total),
              near(sum(as.numeric(saved$Percent_recalculated)), 100))
  }
  check_summary("e_mobile_genetic_context.tsv", "Fig4e_context_recalculated.tsv", 120, "Percent")
  check_summary("taxon_counts.tsv", "Fig4e_taxa_recalculated.tsv", 80, "Percent_of_classified")
  check_summary("e_classification_coverage.tsv", "Fig4e_classification_coverage.tsv", 120, "Percent_of_all_targets")
  context <- read(file.path(output_dir, "e/Fig4e_context_recalculated.tsv"))
  taxa <- read(file.path(output_dir, "e/Fig4e_taxa_recalculated.tsv"))
  coverage <- read(file.path(output_dir, "e/Fig4e_classification_coverage.tsv"))
  stopifnot(near(context$Contigs, c(59, 58, 3)),
            near(taxa$Contigs, c(35, 15, 13, 7, 4, 4, 2)),
            near(coverage$Contigs, c(80, 40)),
            sum(as.numeric(context$Contigs)) == sum(as.numeric(coverage$Contigs)),
            sum(as.numeric(taxa$Contigs)) == as.numeric(coverage$Contigs[coverage$Category == "Classified / host-table matched"]))
  denominators <- read(file.path(output_dir, "e/Fig4e_denominators.tsv"))
  stopifnot(identical(denominators$Summary, c("context", "taxonomy", "coverage")),
            near(denominators$Contigs, c(120, 80, 120)))
  manifest <- read(file.path(output_dir, "run_manifest.tsv"))
  stopifnot(near(manifest$input_rows, c(3, 4, 12, 3, 7, 2)),
            near(manifest$plotted_rows, c(3, 4, 8, 3, 7, 0)))
  session <- readLines(file.path(output_dir, "sessionInfo.txt"), warn = FALSE)
  stopifnot(!any(grepl("time zone:|tzcode source:|[/\\\\]Users[/\\\\]", session, ignore.case = TRUE)),
            !file.exists(file.path(root, "input/pair_metrics.tsv")),
            !any(grepl("\\.rds$", list.files(output_dir, recursive = TRUE), ignore.case = TRUE)))
  cat("PASS\tFig4_saved_outputs\tfiles=17\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  full_args <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", full_args[grepl("^--file=", full_args)][1])
  root <- dirname(normalizePath(script, winslash = "/"))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0L) output <- file.path(root, "reference_output")
  else if (length(args) == 2L && identical(args[1], "--output")) output <- args[2]
  else if (length(args) == 1L) output <- args[1]
  else stop("Usage: Rscript verify_saved_outputs.R [--output OUTPUT_DIRECTORY]")
  verify_saved_fig4(root, output)
}
