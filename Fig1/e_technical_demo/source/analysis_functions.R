
city_palette <- c(City01 = "#4E79A7", City02 = "#E06C75", City03 = "#6CBF84")
setting_shapes <- c(Hospital = 17, WWTP = 18, Community = 16, `Wet Market` = 15)
site_order <- names(setting_shapes)
theme_fig1 <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35),
      axis.ticks.length = unit(0.10, "cm"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      panel.grid = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

align_metadata <- function(mat, meta) {
  if (is.null(rownames(mat)) || anyDuplicated(rownames(mat))) {
    stop("Abundance matrix needs unique sample row names.")
  }
  if (anyDuplicated(meta$Sample)) stop("Metadata has duplicate Sample IDs.")
  missing_meta <- setdiff(rownames(mat), meta$Sample)
  missing_mat <- setdiff(meta$Sample, rownames(mat))
  if (length(missing_meta) || length(missing_mat)) {
    stop("Sample IDs differ between abundance matrix and metadata (matrix-only: ",
         length(missing_meta), "; metadata-only: ", length(missing_mat), ").")
  }
  meta[match(rownames(mat), meta$Sample), , drop = FALSE]
}

compute_pcoa <- function(mat, meta, test_permutations = 999L, seed = 123L) {
  meta <- align_metadata(mat, meta)
  if (!is.numeric(mat) || any(!is.finite(mat)) || any(mat < 0)) {
    stop("PCoA requires a finite, nonnegative numeric matrix.")
  }
  nonzero <- rowSums(mat) > 0
  if (sum(nonzero) < 3L) stop("Fewer than three nonzero samples for PCoA.")
  meta <- meta[nonzero, , drop = FALSE]
  mat <- mat[nonzero, , drop = FALSE]
  dist_mat <- vegan::vegdist(mat, method = "bray")
  test_permutations <- as.integer(test_permutations)
  if (length(test_permutations) != 1L || is.na(test_permutations) || test_permutations < 1L)
    stop("test_permutations must be a positive integer.")
  meta$Sample_Type <- droplevels(factor(meta$Sample_Type))
  if (nlevels(meta$Sample_Type) < 2L || any(table(meta$Sample_Type) < 2L))
    stop("PERMANOVA and PERMDISP require at least two settings and two samples per setting.")
  set.seed(seed)
  perm <- vegan::adonis2(dist_mat ~ Sample_Type, data = meta, permutations = test_permutations)
  dispersion_model <- vegan::betadisper(dist_mat, meta$Sample_Type)
  set.seed(seed)
  dispersion_test <- vegan::permutest(dispersion_model, permutations = test_permutations)
  test_table <- function(x) {
    out <- data.frame(Term = rownames(x), as.data.frame(x), row.names = NULL, check.names = FALSE)
    names(out)[names(out) == "Pr(>F)"] <- "P"
    out
  }
  permanova <- test_table(perm)
  permanova$Term[permanova$Term == "Model"] <- "Sample_Type"
  pcoa <- stats::cmdscale(dist_mat, eig = TRUE, k = 2)
  eig_positive <- pcoa$eig[pcoa$eig > 0]
  if (length(eig_positive) < 2L) stop("PCoA has fewer than two positive axes.")
  scores <- data.frame(
    Sample = rownames(mat), City_Code = as.character(meta$City_Code),
    Sample_Type = as.character(meta$Sample_Type),
    PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )
  list(scores = scores,
       var1 = 100 * eig_positive[1] / sum(eig_positive),
       var2 = 100 * eig_positive[2] / sum(eig_positive),
       permanova = permanova, dispersion = test_table(dispersion_test$tab),
       test_permutations = test_permutations)
}

plot_panel_e <- function(pcoa_result) {
  points <- pcoa_result$scores
  points$City_Code <- factor(points$City_Code, levels = names(city_palette))
  points$Sample_Type <- factor(points$Sample_Type, levels = site_order)
  ggplot(points, aes(x = PCoA1, y = PCoA2, colour = City_Code,
                     shape = Sample_Type)) +
    geom_point(size = 1.45, alpha = 0.70, stroke = 0.18) +
    scale_colour_manual(values = city_palette, breaks = names(city_palette),
                        name = "City", drop = FALSE) +
    scale_shape_manual(values = setting_shapes, breaks = site_order,
                       name = "Setting", drop = FALSE) +
    coord_equal() +
    labs(x = sprintf("PCoA1(%.1f%%)", pcoa_result$var1),
         y = sprintf("PCoA2(%.1f%%)", pcoa_result$var2), tag = "e") +
    theme_fig1(11) +
    theme(legend.position = "right",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7.5),
          legend.key.height = unit(0.18, "cm"),
          legend.key.width = unit(0.25, "cm"),
          plot.tag = element_text(size = 20, face = "plain")) +
    guides(colour = guide_legend(ncol = 2, order = 1,
                                 override.aes = list(size = 2, alpha = 1)),
           shape = guide_legend(ncol = 1, order = 2,
                                override.aes = list(size = 2, alpha = 1)))
}

