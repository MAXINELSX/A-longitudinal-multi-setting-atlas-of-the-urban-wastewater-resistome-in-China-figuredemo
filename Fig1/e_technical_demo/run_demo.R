
load_technical_inputs <- function(root) {
  data <- read.csv(file.path(root, "inputs", "abundance.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  raw_meta <- read.csv(file.path(root, "inputs", "metadata.csv"), check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(names(data)[1] == "Subtype", !anyDuplicated(data$Subtype), nrow(data) == 2646L,
            ncol(data) == 73L, nrow(raw_meta) == 72L,
            all(raw_meta$Metadata_status == "SYNTHETIC_DEMO_ONLY"),
            all(grepl("^Sample[0-9]{3}$", names(data)[-1])),
            all(grepl("^City[0-9]{2}$", raw_meta$City)))
  matrix <- t(as.matrix(data[, -1, drop = FALSE]))
  storage.mode(matrix) <- "double"
  colnames(matrix) <- data$Subtype
  stopifnot(all(is.finite(matrix)), all(matrix >= 0), !anyDuplicated(rownames(matrix)))
  setting <- trimws(raw_meta$Setting)
  setting[tolower(setting) %in% c("wet market", "wetmarket")] <- "Wet Market"
  stopifnot(all(setting %in% c("Hospital", "WWTP", "Community", "Wet Market")))
  metadata <- data.frame(Sample = raw_meta$Sample, City_Code = raw_meta$City,
                         Sample_Type = setting, stringsAsFactors = FALSE)
  list(matrix = matrix, metadata = metadata)
}

run_technical_demo <- function(root, output) {
  required <- c("vegan", "dplyr", "ggplot2", "scales", "readr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({ library(ggplot2); library(grid) })
  source(file.path(root, "source", "analysis_functions.R"), local = environment())
  input <- load_technical_inputs(root)
  metadata <- align_metadata(input$matrix, input$metadata)
  if (!dir.exists(output) && !dir.create(output, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create technical output directory.")
  targets <- c("pcoa_technical_demo.pdf", "pcoa_technical_demo.png", "scores.csv", "permanova_technical.csv",
               "dispersion_technical.csv", "summary.csv")
  if (any(file.exists(file.path(output, targets)))) stop("Technical outputs already exist; choose a new output directory.")
  result <- compute_pcoa(input$matrix, metadata, test_permutations = 99L, seed = 123L)
  figure <- plot_panel_e(result) +
    labs(title = "Figure 1e") +
    theme(plot.title = element_text(size = 12), legend.position = "right")
  plotted_points <- ggplot2::ggplot_build(figure)$data[[1]]
  stopifnot(nrow(plotted_points) == 72L, !anyNA(plotted_points$shape), !anyNA(plotted_points$colour))
  ggplot2::ggsave(file.path(output, "pcoa_technical_demo.pdf"), figure, width = 10.8, height = 6.8,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(output, "pcoa_technical_demo.png"), figure, width = 10.8, height = 6.8,
                  units = "in", dpi = 150, bg = "white")
  result$scores$Setting_status <- "SYNTHETIC_DEMO_ONLY"
  summary <- data.frame(metric = c("n_samples", "n_subtypes", "n_nonzero_samples", "axis1_percent", "axis2_percent", "permutations"),
                         value = c(nrow(input$matrix), ncol(input$matrix), nrow(result$scores), result$var1, result$var2, 99))
  reference_path <- file.path(root, "expected", "summary.csv")
  reference <- read.csv(reference_path, stringsAsFactors = FALSE)
  stopifnot(identical(summary$metric, reference$metric),
            isTRUE(all.equal(summary$value, reference$value, tolerance = 1e-9)))
  readr::write_csv(result$scores, file.path(output, "scores.csv"))
  readr::write_csv(result$permanova, file.path(output, "permanova_technical.csv"))
  readr::write_csv(result$dispersion, file.path(output, "dispersion_technical.csv"))
  readr::write_csv(summary, file.path(output, "summary.csv"))
  message("PASS: Fig1 e")
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L || args[1] != "--output") stop("Usage: Rscript run_demo.R --output OUTPUT_DIRECTORY")
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
  root <- dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
  run_technical_demo(root, args[2])
}
