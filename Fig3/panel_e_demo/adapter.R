
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
  data$species <- factor(data$species, levels = rev(unique(data$species)))
  data$arg_type <- factor(data$arg_type, levels = unique(data$arg_type))
  data$scope <- factor(data$scope, levels = c('Overall', 'Hospital', 'WWTP', 'Community', 'Wet market'))
  ggplot2::ggplot(data, ggplot2::aes(arg_type, species, size = pct_mags, colour = pct_mags)) +
    ggplot2::geom_point(alpha = 0.92) +
    ggplot2::facet_wrap(~ scope, ncol = 3, drop = FALSE) +
    ggplot2::scale_size(range = c(0.45, 3.7), limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = 'Supplied % MAGs') +
    ggplot2::scale_colour_gradientn(colours = c('#4E79A7', '#9ecae1', '#f7f7f7', '#f4a582', '#ca0020'),
      limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), name = 'Supplied % MAGs') +
    ggplot2::scale_x_discrete(labels = function(x) gsub('_', ' ', x), drop = FALSE) +
    ggplot2::labs(title = 'Figure 3e', x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(face = 'italic'), legend.position = 'bottom')
}
