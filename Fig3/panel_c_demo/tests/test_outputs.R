script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- if (exists("panel_demo_root", inherits = FALSE)) panel_demo_root else dirname(dirname(normalizePath(script, winslash = "/", mustWork = TRUE)))
args <- commandArgs(TRUE)
if (length(args) != 2L || args[1] != "--output") stop("Usage: Rscript tests/test_outputs.R --output DIR")
out <- args[2]
s <- read.delim(file.path(root, "input", "setting_summary.tsv"), check.names = FALSE)
c <- read.delim(file.path(root, "input", "contig_observations.tsv"), check.names = FALSE)
actual_s <- read.delim(file.path(out, "setting_summary_plot_input.tsv"), check.names = FALSE)
actual_c <- read.delim(file.path(out, "contig_plot_input.tsv"), check.names = FALSE)
stopifnot(nrow(actual_s) == 4L, isTRUE(all.equal(actual_s, s, check.attributes = FALSE)))
stopifnot(nrow(actual_c) == 12L, isTRUE(all.equal(actual_c, c, check.attributes = FALSE)))
stopifnot(sum(actual_s$N_samples) == 12L, all(actual_c$total_arg_contigs >= 40L & actual_c$total_arg_contigs <= 100L))
for (i in seq_len(nrow(actual_s))) {
  z <- actual_c$mobile_fraction[actual_c$Setting == actual_s$Setting[i]]
  expected <- c(length(z), mean(z), median(z), quantile(z, .25), quantile(z, .75), min(z), max(z), sd(z))
  stopifnot(max(abs(as.numeric(actual_s[i, -1]) - expected)) < 1e-12)
}
stopifnot(all(table(actual_c$Setting) == 3L), all(abs(actual_c$mobile_fraction - actual_c$mobile_arg_contigs / actual_c$total_arg_contigs) < 1e-8))
desc <- read.delim(file.path(out, "contig_descriptive_statistics.tsv"))
for (i in seq_len(nrow(desc))) {
  z <- actual_c$mobile_fraction[actual_c$Setting == desc$Setting[i]]
  stopifnot(desc$N_samples[i] == length(z), abs(desc$Mean[i] - mean(z)) < 1e-12, abs(desc$Median[i] - median(z)) < 1e-12)
}
status <- read.csv(file.path(out, "run_status.csv"))
stopifnot(identical(status$status, c("PASS_SUMMARY_DISPLAY_ONLY", "PASS_CONTIG_PLOT_ONLY", "NOT_RUN_ANNOTATION_INPUT_UNAVAILABLE")))
stopifnot(!any(grepl("p.value|p_value|Letters|estimate|contrast", names(desc), ignore.case = TRUE)))
files <- file.path(out, c("setting_summary.png", "setting_summary.pdf", "contig_violin.png", "contig_violin.pdf"))
stopifnot(all(file.exists(files)), all(file.info(files)$size > 1000L))
stopifnot(!any(file.exists(file.path(out, c("SCOPE.txt", "DEMO_SCOPE.txt", "warnings.txt")))))
cat("PASS_FIG3C_OUTPUTS\n")
