
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

plot_panel_a <- function(data, composition = NULL, counts = NULL) {
  species <- unique(as.character(data$species))
  settings <- c('Hospital', 'WWTP', 'Community', 'Wet market')
  cities <- sort(unique(as.character(data$city)))
  columns <- expand.grid(city = cities, setting = settings, stringsAsFactors = FALSE)
  columns$x <- seq_len(nrow(columns))
  data$x <- columns$x[match(paste(data$city, data$setting), paste(columns$city, columns$setting))]
  data$y <- length(species) + 1 - match(data$species, species)
  setting_palette <- c(Hospital = '#E06C75', WWTP = '#4E79A7', Community = '#6EC5C1', 'Wet market' = '#E9C95B')
  city_palette <- setNames(c('#4E79A7', '#E06C75', '#6CBF84', '#6EC5C1')[seq_along(cities)], cities)
  columns$setting_fill <- unname(setting_palette[columns$setting])
  columns$city_fill <- unname(city_palette[columns$city])
  setting_blocks <- unique(columns[, c('setting', 'setting_fill')])
  setting_blocks$xmin <- vapply(settings, function(s) min(columns$x[columns$setting == s]) - .46, numeric(1))
  setting_blocks$xmax <- vapply(settings, function(s) max(columns$x[columns$setting == s]) + .46, numeric(1))
  palette <- c(multidrug = '#7E6BB5', 'macrolide-lincosamide-streptogramin' = '#74C2B3',
    beta_lactam = '#7DA34D', bacitracin = '#D6B97A', aminoglycoside = '#4E79A7',
    tetracycline = '#E8D98F', polymyxin = '#75C0C1')
  bars <- data.frame(xmin = numeric(), xmax = numeric(), y = numeric(), ARG = character())
  if (!is.null(composition)) {
    types <- setdiff(names(composition), 'species')
    rows <- lapply(species, function(s) {
      item <- composition[composition$species == s, types, drop = FALSE]
      if (!nrow(item)) return(NULL)
      values <- as.numeric(item[1, ])
      values <- values / sum(values)
      ends <- 8.78 + 1.18 * cumsum(values)
      data.frame(species = s, y = length(species) + 1 - match(s, species),
        ARG = types, xmin = c(8.78, head(ends, -1)), xmax = ends)
    })
    bars <- do.call(rbind, rows)
  }
  count_bars <- data.frame(xmin = numeric(), xmax = numeric(), y = numeric())
  if (!is.null(counts)) {
    count_bars <- counts[counts$species %in% species, , drop = FALSE]
    count_bars$y <- length(species) + 1 - match(count_bars$species, species)
    count_bars$xmin <- .58
    count_bars$xmax <- .58 + .28 * count_bars$total_count / max(count_bars$total_count)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = seq_along(species), colour = '#f3f3f3', linewidth = .22) +
    ggplot2::geom_vline(xintercept = cumsum(rep(length(cities), length(settings) - 1)) + .5,
      colour = '#666666', linewidth = .35) +
    ggplot2::geom_rect(data = count_bars, ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = y - .16, ymax = y + .16), fill = '#a7cfdf', inherit.aes = FALSE) +
    ggplot2::geom_point(data = data, ggplot2::aes(x, y, size = mean_relative_abundance,
      colour = mean_relative_abundance)) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = y - .18, ymax = y + .18, fill = ARG), inherit.aes = FALSE) +
    ggplot2::geom_rect(data = setting_blocks, ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = length(species) + .94, ymax = length(species) + 1.23),
      fill = setting_blocks$setting_fill, inherit.aes = FALSE) +
    ggplot2::geom_rect(data = columns, ggplot2::aes(xmin = x - .45, xmax = x + .45,
      ymin = length(species) + .65, ymax = length(species) + .89),
      fill = columns$city_fill, inherit.aes = FALSE) +
    ggplot2::scale_size_area(max_size = 5.4, limits = c(0, .05), oob = scales::squish,
      breaks = c(.01, .02, .05), labels = scales::label_percent(),
      name = 'Fraction of ARG-carrying contigs') +
    ggplot2::scale_colour_gradientn(colours = c('#fff8ee', '#f4b183', '#c12f0f'),
      limits = c(0, .05), oob = scales::squish, guide = 'none') +
    ggplot2::scale_fill_manual(values = palette, guide = 'none') +
    ggplot2::scale_x_continuous(limits = c(.45, nrow(columns) + 2.05), breaks = columns$x,
      labels = columns$city, position = 'top', expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(.5, length(species) + 1.25), breaks = seq_along(species),
      labels = rev(species), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = 'off') +
    ggplot2::theme_classic(base_size = 8.5) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(face = 'italic', colour = '#333333', size = 8),
      axis.text.x.top = ggplot2::element_text(size = 7, angle = 0, margin = ggplot2::margin(b = 1)),
      legend.position = 'bottom', legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7), plot.margin = ggplot2::margin(5, 4, 4, 8))
  p
}
