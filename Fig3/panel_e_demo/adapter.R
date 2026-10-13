
prepare_panel_e <- function(cells) {
  needed <- c('scope', 'species', 'arg_type', 'n_mag', 'n_arg_carrying_mags', 'pct_mags')
  if (!all(needed %in% names(cells)) || anyNA(cells[, needed])) stop('Missing required cell fields.')
  if (anyDuplicated(cells[, c('scope', 'species', 'arg_type')])) stop('Duplicate precomputed cells.')
  numeric_fields <- as.matrix(cells[, c('n_mag', 'n_arg_carrying_mags', 'pct_mags')])
  if (!is.numeric(numeric_fields) || any(!is.finite(numeric_fields))) stop('Cell numbers must be finite.')
  if (any(cells$n_mag <= 0) || any(cells$n_arg_carrying_mags < 0 | cells$n_arg_carrying_mags > cells$n_mag) || any(cells$pct_mags < 0 | cells$pct_mags > 100)) stop('Invalid supplied counts or percentages.')
  if (any(abs(cells$pct_mags - 100 * cells$n_arg_carrying_mags / cells$n_mag) > 1e-10)) stop('Supplied percentage and count fields disagree.')
  cells
}

plot_panel_e <- function(data) {
  settings <- c('Hospital', 'WWTP', 'Community', 'Wet market')
  types <- unique(data$arg_type)
  species <- rev(unique(data$species))
  data <- data[data$scope %in% settings, , drop = FALSE]
  data$scope <- factor(data$scope, levels = settings)
  data$x <- match(data$arg_type, types)
  data$y <- match(data$species, species)
  key <- split(seq_len(nrow(data)), interaction(data$scope, data$species, drop = TRUE))
  bars <- do.call(rbind, lapply(key, function(ix) {
    z <- data[ix, , drop = FALSE]
    z <- z[order(z$x), , drop = FALSE]
    fractions <- z$pct_mags / sum(z$pct_mags)
    ends <- 6.22 + 2.45 * cumsum(fractions)
    data.frame(scope = z$scope, species = z$species, y = z$y, arg_type = z$arg_type,
      xmin = c(6.22, head(ends, -1)), xmax = ends)
  }))
  palette <- c(aminoglycoside = '#4E79A7', antibacterial_fatty_acid = '#C9C9C9',
    bacitracin = '#D6B97A', beta_lactam = '#7DA34D', bicyclomycin = '#8BCB88')
  display_types <- c(aminoglycoside = 'Aminoglycoside', antibacterial_fatty_acid = 'Others',
    bacitracin = 'Bacitracin', beta_lactam = 'Beta-lactam', bicyclomycin = 'Bicyclomycin')
  bars$arg_label <- factor(unname(display_types[bars$arg_type]), levels = unname(display_types))
  ggplot2::ggplot(data, ggplot2::aes(x, y)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_mags, colour = pct_mags), alpha = .95) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax,
      ymin = y - .21, ymax = y + .21, fill = arg_label), inherit.aes = FALSE) +
    ggplot2::geom_vline(xintercept = 5.66, linewidth = .3, colour = '#777777') +
    ggplot2::facet_grid(scope ~ ., switch = 'y') +
    ggplot2::scale_size(range = c(.25, 2.7), limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100), name = '% MAGs') +
    ggplot2::scale_colour_gradientn(colours = c('#fff7d3', '#e4c45b', '#b88a06'),
      limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = '% MAGs') +
    ggplot2::scale_fill_manual(values = setNames(unname(palette[types]), unname(display_types[types])), name = 'ARG class') +
    ggplot2::scale_x_continuous(breaks = seq_along(types), labels = unname(display_types[types]),
      limits = c(.52, 8.92), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(species), labels = species,
      limits = c(.5, length(species) + .5), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = '#eeeeee', linewidth = .2),
      panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(linewidth = .35),
      strip.background = ggplot2::element_rect(fill = 'white', colour = '#666666', linewidth = .35),
      strip.placement = 'outside', strip.text.y.left = ggplot2::element_text(angle = 90, size = 8),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(face = 'italic', colour = '#bd2f26', size = 7),
      legend.position = 'right', legend.title = ggplot2::element_text(size = 7),
      legend.text = ggplot2::element_text(size = 7), panel.spacing.y = grid::unit(.05, 'lines'))
}
