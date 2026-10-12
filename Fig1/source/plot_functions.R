


city_palette <- c(City01 = "#4E79A7", City02 = "#E06C75", City03 = "#6CBF84")
city_to_code <- stats::setNames(names(city_palette), names(city_palette))
site_palette <- c(Hospital = "#E06C75", WWTP = "#4E79A7",
                  Community = "#6EC5C1", `Wet Market` = "#E9C95B")
site_order <- names(site_palette)

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

normalize_setting <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(raw))
  out[raw == "hospital"] <- "Hospital"
  out[raw == "wwtp"] <- "WWTP"
  out[raw == "community"] <- "Community"
  out[raw %in% c("wetmarket", "wet market")] <- "Wet Market"
  out
}

load_metadata <- function(csv_path) {
  if (!file.exists(csv_path)) stop("Metadata CSV not found: ", csv_path)
  raw <- readr::read_csv(csv_path, col_types = readr::cols(.default = "c"),
                         show_col_types = FALSE)
  required <- c("Sample", "Clean_Reads_GBp", "City", "Sample_Type", "Sample_Date")
  missing <- setdiff(required, names(raw))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse = ", "))

  out <- raw[, required]
  out$Sample <- trimws(out$Sample)
  if (anyNA(out$Sample) || any(out$Sample == "") || anyDuplicated(out$Sample)) {
    stop("Sample IDs must be nonblank and unique in metadata.")
  }

  city_raw <- trimws(out$City)
  out$City_Code <- unname(ifelse(city_raw %in% names(city_to_code),
                                 city_to_code[city_raw], city_raw))
  if (anyNA(out$City_Code) || any(!out$City_Code %in% names(city_palette))) {
    stop("Metadata contains an unknown city code/name.")
  }

  out$Sample_Type <- normalize_setting(out$Sample_Type)
  if (anyNA(out$Sample_Type)) stop("Metadata contains an unknown setting.")

  date_raw <- trimws(out$Sample_Date)
  dates <- as.Date(rep(NA_character_, length(date_raw)))
  compact <- grepl("^[0-9]{8}$", date_raw)
  iso <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", date_raw)
  dates[compact] <- as.Date(date_raw[compact], format = "%Y%m%d")
  dates[iso] <- as.Date(date_raw[iso], format = "%Y-%m-%d")
  if (anyNA(dates)) stop("Sample_Date contains invalid or unsupported dates.")
  out$Sample_Date <- dates
  out$YearMonth <- format(dates, "%Y-%m")
  out$Month_Plot_Date <- as.Date(paste0(out$YearMonth, "-15"))

  depth_raw <- trimws(out$Clean_Reads_GBp)
  depth <- suppressWarnings(as.numeric(depth_raw))
  bad_depth <- !is.na(depth_raw) & depth_raw != "" & !is.finite(depth)
  if (any(bad_depth)) stop("Clean_Reads_GBp contains nonnumeric values.")
  if (any(depth < 0, na.rm = TRUE)) stop("Clean_Reads_GBp contains negative values.")
  out$Clean_Reads_GBp <- depth
  out$Sample_Type <- factor(out$Sample_Type, levels = site_order)
  out$City_Code <- factor(out$City_Code, levels = names(city_palette))
  out
}


plot_panels_bc <- function(meta) {
  observed_cities <- names(sort(table(as.character(meta$City_Code)), decreasing = TRUE))
  meta$City_Plot <- factor(as.character(meta$City_Code), levels = rev(observed_cities))
  meta$Sample_Type <- factor(as.character(meta$Sample_Type), levels = site_order)

  monthly <- meta %>%
    count(Sample_Type, City_Plot, Month_Plot_Date, name = "Sample_Count")
  month_breaks <- sort(unique(monthly$Month_Plot_Date))
  month_labels <- function(x) {

    labels <- month.abb[as.integer(format(x, "%m"))]
    show_year <- format(x, "%m") == "01" | x == min(month_breaks)
    labels[show_year] <- paste0(labels[show_year], "\n", format(x[show_year], "%Y"))
    labels
  }
  p_month <- ggplot(monthly,
                    aes(x = Month_Plot_Date, y = City_Plot)) +
    geom_point(aes(size = Sample_Count, colour = Sample_Type), alpha = 0.85) +
    facet_wrap(~ Sample_Type, nrow = 1, drop = FALSE) +
    scale_colour_manual(values = site_palette, breaks = site_order,
                        name = "Setting", drop = FALSE) +
    scale_size_continuous(range = c(1.0, 3.4),
                          limits = c(1, max(15, monthly$Sample_Count)),
                          breaks = c(1, 5, 10, 15), name = "Sample count") +
    scale_x_date(breaks = month_breaks, labels = month_labels,
                 expand = expansion(mult = c(0.025, 0.025))) +
    labs(x = NULL, y = NULL, tag = "b") +
    theme_fig1(9.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5,
                                     size = 7.7),
          axis.text.y = element_text(size = 9.5),
          strip.text = element_text(size = 9.5),
          panel.spacing.x = unit(0.15, "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.key.width = unit(0.27, "cm"),
          plot.tag = element_text(size = 19, face = "plain")) +
    guides(colour = guide_legend(order = 1, nrow = 2,
                                 override.aes = list(size = 2.6, alpha = 1)),
           size = guide_legend(order = 2, nrow = 1,
                               override.aes = list(colour = "black", alpha = 1)))

  sample_stack <- meta %>%
    count(Sample_Type, City_Code, name = "Sample_Count") %>%
    mutate(City_Code = factor(as.character(City_Code), levels = names(city_palette)))
  sample_total <- sample_stack %>%
    group_by(Sample_Type) %>%
    summarise(Total = sum(Sample_Count), .groups = "drop")
  p_sample_bar <- ggplot(sample_stack,
                         aes(x = Sample_Type, y = Sample_Count, fill = City_Code)) +
    geom_col(width = 0.83) +
    geom_text(data = sample_total,
              aes(x = Sample_Type, y = Total, label = scales::comma(Total)),
              inherit.aes = FALSE, vjust = -0.30, size = 3.0) +
    scale_fill_manual(values = city_palette, breaks = names(city_palette),
                      name = "City", drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                       labels = scales::comma_format()) +
    labs(title = "Sample counts", x = NULL, y = NULL) +
    theme_fig1(9.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "bottom", legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.20, "cm")) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

  depth <- meta %>%
    filter(is.finite(Clean_Reads_GBp), Clean_Reads_GBp > 0) %>%
    mutate(Log_Depth = log10(Clean_Reads_GBp))
  if (!nrow(depth)) stop("No positive Clean_Reads_GBp values for panel c.")
  p_depth <- ggplot(depth,
                    aes(x = Log_Depth, y = City_Plot, colour = Sample_Type)) +
    geom_boxplot(width = 0.42, outlier.shape = NA, fill = NA, linewidth = 0.40) +
    geom_point(size = 0.75, alpha = 0.32,
               position = position_jitter(width = 0, height = 0.10, seed = 2614)) +
    facet_wrap(~ Sample_Type, nrow = 1, drop = FALSE) +
    scale_colour_manual(values = site_palette, breaks = site_order,
                        drop = FALSE) +
    scale_x_continuous(breaks = scales::breaks_width(0.2),
                       labels = function(x) sprintf("%.1f", x),
                       expand = expansion(mult = c(0.03, 0.05))) +
    scale_y_discrete(drop = FALSE) +
    labs(x = expression(log[10] * "(Clean bases [GBp])"), y = NULL, tag = "c") +
    theme_fig1(9.5) +
    theme(axis.text.y = element_text(size = 9.5),
          axis.text.x = element_text(size = 8.5),
          strip.text = element_text(size = 9.5),
          panel.spacing.x = unit(0.15, "cm"),
          legend.position = "none",
          plot.tag = element_text(size = 19, face = "plain"))

  depth_stack <- depth %>%
    group_by(Sample_Type, City_Code) %>%
    summarise(Total_TBp = sum(Clean_Reads_GBp) / 1000, .groups = "drop") %>%
    mutate(City_Code = factor(as.character(City_Code), levels = names(city_palette)))
  depth_total <- depth_stack %>%
    group_by(Sample_Type) %>%
    summarise(Total = sum(Total_TBp), .groups = "drop")
  p_depth_bar <- ggplot(depth_stack,
                        aes(x = Sample_Type, y = Total_TBp, fill = City_Code)) +
    geom_col(width = 0.83) +
    geom_text(data = depth_total,
              aes(x = Sample_Type, y = Total, label = sprintf("%.1f", Total)),
              inherit.aes = FALSE, vjust = -0.30, size = 3.0) +
    scale_fill_manual(values = city_palette, breaks = names(city_palette),
                      drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(title = "Total clean bases (TBp)", x = NULL, y = NULL) +
    theme_fig1(9.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
          legend.position = "none")

  shared_legend <- cowplot::get_legend(p_month)
  city_legend <- cowplot::get_legend(p_sample_bar)
  row_b <- cowplot::plot_grid(p_month + theme(legend.position = "none"),
                              p_sample_bar + theme(legend.position = "none"),
                              nrow = 1, rel_widths = c(3.55, 1.05),
                              align = "h", axis = "tb")
  row_c <- cowplot::plot_grid(p_depth, p_depth_bar, nrow = 1,
                              rel_widths = c(3.55, 1.05), align = "h", axis = "tb")
  legend_row <- cowplot::plot_grid(shared_legend, city_legend, nrow = 1,
                                   rel_widths = c(3.55, 1.05))
  list(b = row_b, c = row_c, legends = legend_row)
}


save_pdf_safely <- function(plot, output_path, width = 14.2, height = 10.8) {
  out_dir <- dirname(output_path)
  if (!dir.exists(out_dir) &&
      !dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create output directory: ", out_dir)
  }
  tmp <- tempfile(pattern = "fig1_be_", fileext = ".pdf")
  on.exit(unlink(tmp), add = TRUE)
  ggplot2::ggsave(tmp, plot = plot, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf, bg = "white")
  if (!file.copy(tmp, output_path, overwrite = TRUE)) {
    stop("Cannot write PDF: ", output_path)
  }
  invisible(output_path)
}


plot_panel_d <- function(curves) {
  setting <- curves$setting
  setting$Group <- factor(setting$Group, levels = site_order)
  p_main <- ggplot(setting, aes(x = N, y = Subtypes, colour = Group,
                                linetype = Group)) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = site_palette, breaks = site_order, name = "Setting") +
    scale_linetype_manual(values = c(Hospital = "solid", WWTP = "dashed",
                                    Community = "dashed", `Wet Market` = "solid"),
                          breaks = site_order, name = "Setting") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06)),
                       labels = scales::comma_format()) +
    labs(x = "Number of Samples", y = "Number of Subtypes", tag = "d") +
    theme_fig1(11) +
    theme(legend.position = c(0.17, 0.19),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          plot.tag = element_text(size = 20, face = "plain")) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(linewidth = 1)))

  city <- curves$city
  city$Group <- factor(city$Group, levels = names(city_palette))
  p_inset <- ggplot(city, aes(x = N, y = Subtypes, colour = Group)) +
    geom_line(linewidth = 0.47) +
    scale_colour_manual(values = city_palette, breaks = names(city_palette),
                        name = "City", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(labels = scales::comma_format(),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = "Number of Samples", y = "Number of Subtypes") +
    theme_fig1(8) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          axis.text = element_text(size = 6.5),
          axis.title = element_text(size = 7.3),
          legend.position = "right", legend.title = element_text(size = 7),
          legend.text = element_text(size = 6.5),
          legend.key.height = unit(0.13, "cm"),
          legend.key.width = unit(0.17, "cm"),
          plot.margin = margin(2, 2, 2, 2))
  cowplot::ggdraw(p_main) +
    cowplot::draw_plot(p_inset, x = 0.34, y = 0.21, width = 0.64, height = 0.59)
}

