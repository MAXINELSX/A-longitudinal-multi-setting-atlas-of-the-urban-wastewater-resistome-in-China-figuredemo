
prepare_panel_a <- function(abundance, metadata) {
  if (!all(c('sample', 'city', 'setting') %in% names(metadata)) || !'species' %in% names(abundance)) stop('Missing required columns.')
  if (anyNA(metadata) || any(!nzchar(as.matrix(metadata))) || anyDuplicated(metadata$sample) || anyDuplicated(abundance$species)) stop('Metadata and species identifiers must be complete and unique.')
  samples <- setdiff(names(abundance), 'species')
  if (!setequal(samples, metadata$sample)) stop('The sample matrix and metadata must match exactly.')
  values <- as.matrix(abundance[, samples, drop = FALSE])
  if (!is.numeric(values) || any(!is.finite(values)) || any(values < 0 | values > 1)) stop('Relative abundance must be finite and in [0, 1].')
  groups <- unique(metadata[, c('city', 'setting')])
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    ids <- metadata$sample[metadata$city == groups$city[i] & metadata$setting == groups$setting[i]]
    data.frame(species = abundance$species, city = groups$city[i], setting = groups$setting[i],
      n_samples = length(ids), mean_relative_abundance = rowMeans(as.matrix(abundance[, ids, drop = FALSE])),
      stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

plot_panel_a <- function(data) {
  data$species <- factor(data$species, levels = rev(unique(data$species)))
  data$setting <- factor(data$setting, levels = c('Hospital', 'WWTP', 'Community', 'Wet market'))
  ggplot2::ggplot(data, ggplot2::aes(city, species, size = mean_relative_abundance, colour = mean_relative_abundance)) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(. ~ setting) +
    ggplot2::scale_size_area(max_size = 8, name = 'Mean relative abundance', breaks = c(0.01, 0.02, 0.05), labels = scales::label_percent()) +
    ggplot2::scale_colour_gradientn(colours = c('#F7F7F7', '#F6D6B8', '#E98B4A', '#B30000'),
      values = c(0, 0.2, 0.4, 1), limits = c(0, 0.05), oob = scales::squish,
      breaks = c(0.01, 0.02, 0.05), labels = scales::label_percent(), name = 'Mean relative abundance') +
    ggplot2::labs(title = 'Figure 3a', x = 'City', y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(face = 'italic'),
      legend.position = 'bottom')
}
