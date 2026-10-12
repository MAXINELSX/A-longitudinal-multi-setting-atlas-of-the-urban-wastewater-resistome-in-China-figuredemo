

script_root <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument)) return(dirname(normalizePath(sub("^--file=", "", argument[1]), mustWork = TRUE)))
  stop("Run this file with Rscript.")
}

run_demo <- function(root, output) {
  required <- c("readr", "dplyr", "ggplot2", "cowplot", "scales", "vegan")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({
    library(readr); library(dplyr); library(ggplot2); library(cowplot); library(scales); library(grid)
  })
  source(file.path(root, "source", "plot_functions.R"), local = environment())
  meta <- load_metadata(file.path(root, "inputs", "metadata.csv"))
  if (nrow(meta) != 96L || any(!grepl("^Sample[0-9]{3}$", meta$Sample)) ||
      any(!grepl("^City[0-9]{2}$", as.character(meta$City_Code)))) {
    stop("The bundled fixture must contain 96 pseudonymized samples and generic city codes.")
  }
  if (any(!is.finite(meta$Clean_Reads_GBp) | meta$Clean_Reads_GBp <= 0)) {
    stop("The bundled demonstration requires positive finite sequencing depth.")
  }
  if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory.")
  }
  result_paths <- c("fig1_bc_demo.pdf", "fig1_bc_demo.png", "precomputed_curves_demo.pdf", "precomputed_curves_demo.png", "monthly_counts.csv", "site_city_totals.csv",
                    "site_totals.csv", "summary.csv", "coverage.csv", "session_info.txt")
  if (any(file.exists(file.path(output, result_paths)))) {
    stop("Output files already exist; choose a new output directory to preserve existing results.")
  }

  monthly <- meta %>% count(Sample_Type, City_Code, YearMonth, name = "Sample_Count")
  grouped <- meta %>% group_by(Sample_Type, City_Code) %>%
    summarise(Sample_Count = n(), Total_GBp = sum(Clean_Reads_GBp),
              Total_TBp = sum(Clean_Reads_GBp) / 1000, .groups = "drop")
  by_site <- grouped %>% group_by(Sample_Type) %>%
    summarise(Sample_Count = sum(Sample_Count), Total_GBp = sum(Total_GBp),
              Total_TBp = sum(Total_TBp), .groups = "drop")
  summary <- data.frame(
    metric = c("n_samples", "n_cities", "n_settings", "n_monthly_groups", "total_clean_gbp", "total_clean_tbp"),
    value = c(nrow(meta), length(unique(meta$City_Code)), length(unique(meta$Sample_Type)),
              nrow(monthly), sum(meta$Clean_Reads_GBp), sum(meta$Clean_Reads_GBp) / 1000))

  reference <- read.csv(file.path(root, "expected", "summary.csv"), stringsAsFactors = FALSE)
  stopifnot(identical(summary$metric, reference$metric),
            isTRUE(all.equal(summary$value, reference$value, tolerance = 1e-12)))
  site_reference <- read.csv(file.path(root, "expected", "site_totals.csv"), stringsAsFactors = FALSE)
  observed <- as.data.frame(by_site)
  observed$Sample_Type <- as.character(observed$Sample_Type)
  stopifnot(isTRUE(all.equal(observed, site_reference, tolerance = 1e-12, check.attributes = FALSE)))

  curves <- list(setting = read.csv(file.path(root, "inputs", "precomputed_setting_curves.csv")),
                 city = read.csv(file.path(root, "inputs", "precomputed_city_curves.csv")))
  curve_summary <- do.call(rbind, lapply(names(curves), function(scope) {
    do.call(rbind, lapply(split(curves[[scope]], curves[[scope]]$Group), function(d) {
      stopifnot(all(is.finite(d$N)), all(is.finite(d$Subtypes)), all(d$N > 0),
                all(d$Subtypes >= 0), !anyDuplicated(d$N), all(diff(d$N) > 0))
      data.frame(Scope = scope, Group = d$Group[1], Rows = nrow(d), MaxN = max(d$N),
                 EndSubtypes = d$Subtypes[which.max(d$N)])
    }))
  }))
  rownames(curve_summary) <- NULL
  curve_reference <- read.csv(file.path(root, "expected", "precomputed_curves.csv"))
  stopifnot(isTRUE(all.equal(curve_summary, curve_reference, tolerance = 1e-12, check.attributes = FALSE)))

  bc <- plot_panels_bc(meta)
  title <- cowplot::ggdraw() + cowplot::draw_label("Figure 1b-c",
                                                x = 0.01, hjust = 0, size = 12)
  figure <- cowplot::plot_grid(title, bc$b, bc$c, bc$legends,
                               ncol = 1, rel_heights = c(0.14, 1, 1, 0.28))
  save_pdf_safely(figure, file.path(output, "fig1_bc_demo.pdf"), width = 14.2, height = 7.4)
  ggplot2::ggsave(file.path(output, "fig1_bc_demo.png"), figure, width = 14.2,
                  height = 7.4, units = "in", dpi = 150, bg = "white")
  curve_title <- cowplot::ggdraw() + cowplot::draw_label(
    "Figure 1d", x = 0.01, hjust = 0, size = 12)
  curve_figure <- cowplot::plot_grid(curve_title, plot_panel_d(curves),
                                     ncol = 1, rel_heights = c(0.12, 1))
  save_pdf_safely(curve_figure, file.path(output, "precomputed_curves_demo.pdf"), width = 10.8, height = 7)
  ggplot2::ggsave(file.path(output, "precomputed_curves_demo.png"), curve_figure,
                  width = 10.8, height = 7, units = "in", dpi = 150, bg = "white")
  write_csv(monthly, file.path(output, "monthly_counts.csv"))
  write_csv(grouped, file.path(output, "site_city_totals.csv"))
  write_csv(by_site, file.path(output, "site_totals.csv"))
  write_csv(summary, file.path(output, "summary.csv"))
  coverage <- data.frame(panel = c("b", "c", "d", "e"),
                         status = c("ran", "ran", "precomputed_plot_only", "separate_technical_branch"))
  write_csv(coverage, file.path(output, "coverage.csv"))

  versions <- c(paste("R", getRversion()),
                vapply(required, function(p) paste(p, as.character(packageVersion(p))), character(1)))
  writeLines(versions, file.path(output, "session_info.txt"))
  technical_env <- new.env(parent = environment())
  sys.source(file.path(root, "e_technical_demo", "run_demo.R"), envir = technical_env)
  technical_env$run_technical_demo(file.path(root, "e_technical_demo"), file.path(output, "e_technical_demo"))
  message("PASS: Fig1 b-e")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--help")) {
    cat("Usage: Rscript run_demo.R --output OUTPUT_DIRECTORY\n")
  } else {
    if (length(args) != 2L || args[1] != "--output" || !nzchar(args[2])) {
      stop("Usage: Rscript run_demo.R --output OUTPUT_DIRECTORY")
    }
    run_demo(script_root(), args[2])
  }
}
