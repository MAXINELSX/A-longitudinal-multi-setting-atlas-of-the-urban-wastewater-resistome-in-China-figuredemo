args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
out <- if (length(args)) args[[1]] else file.path(root, "reference_output")
required <- c("run_status.csv", "input_audit/duplicate_mag_records.csv",
  "adapted_fig3d_inputs_SYNTHETIC_metadata/demo_MAG_ID_map.csv",
  "fig3d_SYNTHETIC_metadata/Fig3d_SYNTHETIC_metadata_demo.pdf",
  "fig3d_SYNTHETIC_metadata/Fig3d_SYNTHETIC_metadata_demo.png",
  "fig3d_SYNTHETIC_metadata/Fig3d_Dunn_BH_contrasts.csv",
  "fig3d_SYNTHETIC_metadata/Fig3d_Kruskal_Wallis_tests.csv")
paths <- file.path(out, required)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0))
status <- read.csv(file.path(out, "run_status.csv"))
stopifnot("Five included panel demos must run" =
  nrow(status) == 5L && all(grepl("^PASS", status$status)) &&
  identical(status$panel, letters[1:5]))
stopifnot("A demo must not claim complete manuscript reproduction" =
  all(!status$complete_manuscript_reproduction))
for (panel in c("a", "b", "c", "e")) {
  module <- paste0("panel_", panel, "_demo")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  code <- system2(rscript, c("--vanilla", shQuote(file.path(root, module, "test_outputs.R")),
    "--output", shQuote(file.path(out, module))))
  if (code != 0L) stop("Processed-table panel check failed: ", panel)
}
mapping <- read.csv(file.path(out, "adapted_fig3d_inputs_SYNTHETIC_metadata/demo_MAG_ID_map.csv"))
stopifnot(nrow(mapping) == 216L, !anyDuplicated(mapping$demo_mag_uid))
stopifnot(sum(mapping$included_in_species_analysis) == 171L)
contrast <- read.csv(file.path(out, "fig3d_SYNTHETIC_metadata/Fig3d_Dunn_BH_contrasts.csv"))
stopifnot(all(c("Z", "P.unadj", "P.adj", "n_sites_1", "n_sites_2") %in% names(contrast)))
stopifnot(all(is.finite(contrast$P.unadj)), all(contrast$P.unadj >= 0 & contrast$P.unadj <= 1))
stopifnot(!"Demo_scope" %in% names(contrast))
metadata <- read.csv(file.path(out, "adapted_fig3d_inputs_SYNTHETIC_metadata/Fig3_sample_metadata.csv"))
stopifnot(all(metadata$synthetic_metadata), !"metadata_provenance" %in% names(metadata))
products <- list.files(out, recursive = TRUE)
stopifnot(!any(basename(products) %in% c("SCOPE.txt", "DEMO_SCOPE.txt", "warnings.txt", "README.md")))
for (path in products[grepl("\\.(csv|tsv)$", products)]) {
  header <- readLines(file.path(out, path), n = 1L, warn = FALSE)
  stopifnot(!grepl("Demo_scope|metadata_provenance|measurement", header))
}
for (pdf in paths[grepl("\\.pdf$", paths)]) {
  con <- file(pdf, "rb")
  signature <- readChar(con, 5L, useBytes = TRUE)
  close(con)
  stopifnot(signature == "%PDF-")
}
cat("PASS_FIG3_OUTPUTS panels=5 burden_rows=216 species_assigned=171\n")
