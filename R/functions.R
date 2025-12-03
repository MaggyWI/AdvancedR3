#' Creating a table with descriptive stat
#'
#' @param data Give the data
#'
#' @returns A data.frame/tibble.
#'

create_table_descriptive_stats <- function(data) {
  tablestat <- data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd, median = median, iqr = IQR))) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), \(x) round(x, digits = 1))) |>
    dplyr::mutate(
      MeanSD = glue::glue("{value_mean} ({value_sd})"),
      MedianIQR = glue::glue("{value_median} ({value_iqr})")
    ) |>
    dplyr::select(Metabolite = metabolite, `Mean SD` = MeanSD, `Median IQR` = MedianIQR)
  return(tablestat)
}


#' Plot distribution metabolites value
#'
#' @param data table of data
#'
#' @returns a plot object
#'
create_plot_distributions <- function(data) {
  plotdistrib <- data |>
    ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free") +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "Lipidomics Values") +
    ggplot2::labs(
      x = "Metabolite Value",
      y = "Count"
    )
  return(plotdistrib)
}
