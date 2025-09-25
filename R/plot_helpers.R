

## File: R/plot_helpers.R
#' Scatter plot with Spearman correlation annotation
#'
#' @param df data.frame
#' @param x unquoted column name for x axis
#' @param y unquoted column name for y axis
#' @param xlab x axis label
#' @param ylab y axis label
#' @param title plot title
#' @import ggplot2
#' @importFrom stats cor.test
#' @import rlang
#' @export
scatter_plot <- function(df, x, y, xlab = NULL, ylab = NULL, title = NULL) {
  xq <- rlang::enquo(x)
  yq <- rlang::enquo(y)
  xval <- df[[rlang::as_name(xq)]]
  yval <- df[[rlang::as_name(yq)]]

  ct <- stats::cor.test(xval, yval, method = "spearman", use = "complete.obs")
  rho <- round(ct$estimate, 2)
  pval <- signif(ct$p.value, 2)
  label <- paste0("ρ = ", rho, "
P = ", pval)

  ggplot2::ggplot(df, ggplot2::aes(x = !!xq, y = !!yq)) +
    ggplot2::geom_point(alpha = 0.6, color = "steelblue") +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    ggplot2::annotate("text", x = Inf, y = Inf, label = label, hjust = 1.1, vjust = 1.5,
                      size = 4, color = "black") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = xlab %||% rlang::as_label(xq), y = ylab %||% rlang::as_label(yq), title = title)
}

#' Simple bar plot helper
#' @import ggplot2
#' @import rlang
#' @export
#' Simple bar plot helper (bars ordered by descending counts)
#' @import ggplot2
#' @import rlang
#' @export
bar_plot <- function(data, x, fill = "#4C72B0", title = NULL, rotate = FALSE) {
  xq <- rlang::enquo(x)
  var_name <- rlang::as_name(xq)
  
  # Reorder factor by descending counts
  data[[var_name]] <- factor(data[[var_name]], 
                             levels = names(sort(table(data[[var_name]]), decreasing = TRUE)))
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = !!xq)) +
    ggplot2::geom_bar(fill = fill) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = NULL, y = "Count")
  
  if (rotate) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1))
  }
  
  return(p)
}

#' Density plot helper
#' @import ggplot2
#' @import rlang
#' @export
density_plot <- function(data, x, fill = "#4C72B0", title = NULL) {
  xq <- rlang::enquo(x)
  ggplot2::ggplot(data, ggplot2::aes(x = !!xq)) +
    ggplot2::geom_density(fill = fill, alpha = 0.5) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = rlang::as_label(xq), y = "Density")
}

#' Histogram helper
#' @import ggplot2
#' @import rlang
#' @export
hist_plot <- function(data, x, fill = "#4C72B0", binwidth = NULL, title = NULL, xlab = NULL) {
  xq <- rlang::enquo(x)
  ggplot2::ggplot(data, ggplot2::aes(x = !!xq)) +
    ggplot2::geom_histogram(binwidth = binwidth, fill = fill, color = "black") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(title = title, x = xlab %||% rlang::as_label(xq), y = "Count")
}

#' Boxplot helper
#' @import ggplot2
#' @import rlang
#' @export
box_plot <- function(data, x, y, title = NULL, xlab = NULL, ylab = NULL, colors = NULL) {
  xq <- rlang::enquo(x)
  yq <- rlang::enquo(y)
  
  # Colors
  fillvals <- colors %||% grDevices::rainbow(length(unique(data[[rlang::as_name(xq)]])))
  
  # Compute counts per group
  counts <- data %>%
    dplyr::group_by(group = .data[[rlang::as_name(xq)]]) %>%
    dplyr::summarise(n = sum(!is.na(.data[[rlang::as_name(yq)]]))) %>%
    dplyr::mutate(label = paste0("n=", n))
  
  # Create boxplot with counts above boxes
  ggplot2::ggplot(data, ggplot2::aes(x = !!xq, y = !!yq, fill = as.factor(!!xq))) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_fill_manual(values = fillvals) +
    ggplot2::geom_text(data = counts, ggplot2::aes(x = group, y = max(data[[rlang::as_name(yq)]], na.rm = TRUE) * 1.05, 
                                                   label = label), inherit.aes = FALSE, size = 3) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggplot2::labs(title = title, x = xlab %||% rlang::as_label(xq), y = ylab %||% rlang::as_label(yq))
}

