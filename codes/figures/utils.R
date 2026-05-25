parse_clusters_to_ranges <- function(cluster_str) {
  if (is.na(cluster_str) || cluster_str == "" || trimws(cluster_str) == "NA") {
    return(data.frame(start = integer(), end = integer()))
  }

  values <- suppressWarnings(as.integer(strsplit(cluster_str, ",")[[1]]))
  values <- sort(values[!is.na(values)])
  if (length(values) == 0) {
    return(data.frame(start = integer(), end = integer()))
  }
  if (length(values) == 1) {
    return(data.frame(start = values[1], end = values[1]))
  }

  starts <- integer()
  ends <- integer()
  current_start <- values[1]
  current_end <- values[1]

  for (value in values[-1]) {
    if (value == current_end + 1) {
      current_end <- value
    } else {
      starts <- c(starts, current_start)
      ends <- c(ends, current_end)
      current_start <- value
      current_end <- value
    }
  }

  starts <- c(starts, current_start)
  ends <- c(ends, current_end)
  data.frame(start = starts, end = ends)
}

format_ic_plot <- function(scice_results, threshold = 1.005, y_limits = c(NA, 1.15)) {
  ic_plot <- plot_ic(
    scice_results,
    threshold = threshold,
    title = "Bootstrapped IC Dist",
    show_threshold = FALSE
  )

  ic_plot$layers <- Filter(
    function(layer) !inherits(layer$geom, c("GeomPoint", "GeomJitter", "GeomText", "GeomLabel")),
    ic_plot$layers
  )

  if (!is.null(ic_plot$scales$scales)) {
    ic_plot$scales$scales <- Filter(function(scale_obj) !"fill" %in% scale_obj$aesthetics, ic_plot$scales$scales)
  }

  for (idx in seq_along(ic_plot$layers)) {
    layer <- ic_plot$layers[[idx]]
    if (inherits(layer$geom, "GeomBoxplot")) {
      layer$mapping$fill <- NULL
      layer$aes_params$fill <- NA
      layer$aes_params$colour <- "#00468BFF"
      layer$aes_params$outlier.shape <- 16
      layer$aes_params$outlier.size <- 0.01
      layer$geom_params$outlier.size <- 1
      layer$geom_params$outlier.shape <- 16
      layer$geom_params$outlier.colour <- "#00468BFF"
      ic_plot$layers[[idx]] <- layer
    }
    if (inherits(layer$geom, "GeomHline")) {
      ic_plot$layers[[idx]] <- NULL
    }
  }

  ic_plot$layers <- ic_plot$layers[!vapply(ic_plot$layers, is.null, logical(1))]

  ic_plot +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(-0.1, "cm"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = y_limits)
}
