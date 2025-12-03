#' Creating a table with descriptive stat
#'
#' @param data Give the data
#'
#' @returns A data.frame/tibble.
#'

create_table_descriptive_stats <- function(data) {
  tablestat <- data |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), \(x) round(x, digits = 1))) |>
    dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
    dplyr::select(Metabolite = metabolite, `Mean SD` = MeanSD)
  return(tablestat)
}



#' Plot distribution metabolites value
#'
#' @param data table of data
#'
#' @returns histograms
#'
create_plot_distributions <- function(data) {
  plotdistrib <-
    ggplot2::ggplot(data, aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(vars(metabolite), scales = "free") +
    ggplot2::theme_minimal()
  return(plotdistrib)
}
