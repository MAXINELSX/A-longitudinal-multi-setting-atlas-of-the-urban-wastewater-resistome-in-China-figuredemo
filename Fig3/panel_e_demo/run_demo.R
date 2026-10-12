run_panel_e_demo <- function(root, output_dir = file.path(root, 'reference_output')) {
  started <- proc.time()[['elapsed']]
  if (nzchar(Sys.getenv('DEMO_R_LIB'))) .libPaths(c(Sys.getenv('DEMO_R_LIB'), .libPaths()))
  packages <- c('ggplot2', 'scales')
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install required packages: ', paste(missing, collapse = ', '))
  source(file.path(root, 'adapter.R'), local = TRUE)
  output_dir <- normalizePath(output_dir, winslash = '/', mustWork = FALSE)
  input_dir <- normalizePath(file.path(root, 'input'), winslash = '/', mustWork = TRUE)
  if (identical(tolower(output_dir), tolower(input_dir))) stop('Output must differ from input.')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  warning_codes <- character()
  result <- withCallingHandlers({
    supplied <- read.delim(file.path(input_dir, 'precomputed_cells.tsv'), check.names = FALSE)
    cells <- prepare_panel_e(supplied)
    options(digits = 17)
    write.table(cells, file.path(output_dir, 'displayed_cells.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    plot <- plot_panel_e(cells)
    ggplot2::ggsave(file.path(output_dir, 'panel_e_demo.pdf'), plot, width = 12, height = 6.8, units = 'in', bg = 'white')
    ggplot2::ggsave(file.path(output_dir, 'panel_e_demo.png'), plot, width = 12, height = 6.8, units = 'in', dpi = 180, bg = 'white')
    cells
  }, warning = function(w) { warning_codes <<- c(warning_codes, class(w)[1]) })
  diagnostics <- data.frame(code = unique(warning_codes), count = as.integer(table(factor(warning_codes, levels = unique(warning_codes)))))
  write.table(diagnostics, file.path(output_dir, 'diagnostics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  versions <- data.frame(component = c('R', packages), version = c(as.character(getRversion()), vapply(packages, function(p) as.character(packageVersion(p)), character(1))))
  write.table(versions, file.path(output_dir, 'versions.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
  writeLines(c('status=complete', paste0('elapsed_seconds=', round(proc.time()[['elapsed']] - started, 3)),
    paste0('warnings=', length(warning_codes))), file.path(output_dir, 'run_summary.txt'))
  cat('PASS_FIG3E cells=150\n')
  invisible(result)
}

if (sys.nframe() == 0L) {
  cmd <- commandArgs(trailingOnly = FALSE)
  root <- dirname(normalizePath(sub('^--file=', '', cmd[grepl('^--file=', cmd)][1]), winslash = '/', mustWork = TRUE))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) && (length(args) != 2L || args[1] != '--output')) stop('Usage: Rscript run_demo.R [--output DIR]')
  run_panel_e_demo(root, if (length(args)) args[2] else file.path(root, 'reference_output'))
}
