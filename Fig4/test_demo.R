
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)][1])
root <- dirname(normalizePath(script, winslash = "/"))
source(file.path(root, "run_demo.R"), local = TRUE)
if (!exists("prepare_fig4d_demo", mode = "function")) {
  stop("FAIL: Final-master partial Figure 4d adapter is not implemented.")
}
input <- file.path(root, "input")
d <- prepare_fig4d_demo(file.path(input, "pair_demo.tsv"))
stopifnot(nrow(d) == 4L, identical(d$Pair_ID, sprintf("Pair%02d", 1:4)),
          !any(c("Clinical_specimen", "High_priority_ARG", "displayed_in_Fig4d") %in% names(d)),
          identical(d$direct_remap_supported_bool, rep("True", 4)),
          identical(d$instrain_partial_support_bool, rep("True", 4)),
          identical(d$shared_carbapenemase_present_bool, rep("False", 4)))
stopifnot(isTRUE(all.equal(d$ANI_min, c(99.6987, 99.5089, 99.4822, 99.4938))),
          isTRUE(all.equal(d$fixed_snp_per_mbp_pair,
            c(22.903913151180273, 127.07604669745322, 164.52127023443416, 562.5572306296885))),
          isTRUE(all.equal(d$callable_fraction_pair,
            c(0.5314382566078127, 0.8237739060281446, 0.8613229689826616, 0.8069241446504305))),
          identical(d$clinical_ARG_count, c(11, 28, 26, 37)),
          identical(d$mag_ARG_count, c(8, 11, 10, 8)),
          identical(d$shared_ARG_count, c(8, 9, 9, 8)))
rejects <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
bad <- d
bad$direct_remap_supported_bool[1] <- "Unknown"
bad_file <- file.path(tempdir(), "fig4 invalid flag.tsv")
write_demo_tsv(bad, bad_file)
rejects(prepare_fig4d_demo(bad_file), "supplied Boolean")
bad <- d
bad$shared_ARG_count[1] <- 100
write_demo_tsv(bad, bad_file)
rejects(prepare_fig4d_demo(bad_file), "Shared ARG count")
e <- prepare_fig4e_demo(input)
stopifnot(nrow(e$evidence) == 12L, nrow(e$context) == 3L,
          nrow(e$taxa) == 7L, nrow(e$coverage) == 2L,
          identical(unname(e$denominators), c(120, 80, 120)),
          identical(e$context$Contigs, c(48, 60, 12)),
          identical(e$taxa$Contigs, c(20, 10, 10, 10, 10, 10, 10)),
          identical(e$coverage$Contigs, c(80, 40)),
          sum(e$coverage$Contigs) == sum(e$context$Contigs),
          sum(e$taxa$Contigs) == e$coverage$Contigs[1],
          isTRUE(all.equal(e$context$Percent_recalculated, c(40, 50, 10))),
          isTRUE(all.equal(e$taxa$Percent_recalculated, c(25, rep(12.5, 6)))),
          isTRUE(all.equal(e$coverage$Percent_recalculated, 100 * c(80, 40) / 120)))
bad_input <- file.path(tempdir(), "fig4 inconsistent denominator input")
dir.create(bad_input, showWarnings = FALSE)
invisible(file.copy(list.files(input, full.names = TRUE), bad_input, overwrite = TRUE))
bad_context <- e$context
bad_context$Contigs <- bad_context$Contigs / 2
write_demo_tsv(bad_context, file.path(bad_input, "e_mobile_genetic_context.tsv"))
rejects(prepare_fig4e_demo(bad_input), "Context total differs from coverage total")
out <- file.path(tempdir(), "fig4 demo tests with spaces")
result <- run_fig4_demo(root, input, out)
stopifnot(nrow(result$c$data) == 3L,
          sum(result$c$data$Clinical) == 5,
          sum(result$c$data$Wastewater) == 3,
          isTRUE(all.equal(result$c$data$WW_percent, c(50, 50, 25))))
stopifnot(identical(result$c$plot$labels$title, "Figure 4c"),
          is.null(result$c$plot$labels$subtitle), is.null(result$c$plot$labels$caption),
          identical(result$d$plot$patches$annotation$title, "Figure 4d"),
          is.null(result$d$plot$patches$annotation$subtitle),
          is.null(result$d$plot$patches$annotation$caption),
          identical(result$e$plot$patches$annotation$title, "Figure 4e"),
          is.null(result$e$plot$patches$annotation$subtitle),
          is.null(result$e$plot$patches$annotation$caption))
e_panels <- c(result$e$plot$patches$plots, list(result$e$plot))
stopifnot(all(vapply(e_panels, function(p) is.null(p$labels$subtitle) && is.null(p$labels$caption), logical(1))))
source(file.path(root, "verify_saved_outputs.R"), local = TRUE)
verify_saved_fig4(root, out)
cat("PASS\tFig4_tests\n")
