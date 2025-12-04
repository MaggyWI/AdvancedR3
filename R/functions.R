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


#' Do some cleaning to fix issues in the data.
#'
#' @param data The lipidomics data frame.
#'
#' @returns A data frame.
#'
clean <- function(data) {
  data |>
    dplyr::group_by(dplyr::pick(-value)) |>
    dplyr::summarise(value = mean(value), .groups = "keep") |>
    dplyr::ungroup()
}


#' Fix data to process it for model fitting.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
preprocess <- function(data) {
  data |>
    dplyr::mutate(
      class = as.factor(class),
      value = scale(value)
    )
}


#' Fit the model to the data and get the results.
#'
#' @param data The data to fit.
#' @param model The formula.
#'
#' @returns A data frame of the results.
#'
fit_model <- function(data, model) {
  glm(
    formula = model,
    data = data,
    family = binomial
  ) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(data$metabolite),
      model = format(model),
      .before = tidyselect::everything()
    )
}


#' Create model results for report.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame of model results.
#'
create_model_results <- function(data) {
  data |>
    dplyr::group_split(metabolite) |>
    purrr::map(preprocess) |>
    purrr::map(fit_all_models) |>
    purrr::list_rbind()
}


#' Fit all models to a given data frame.
#'
#' @param data The data frame to fit the models to.
#'
#' @returns A data frame with results from all models.
#'
fit_all_models <- function(data) {
  list(
    class ~ value,
    class ~ value + gender + age
  ) |>
    purrr::map(\(model) fit_model(data, model = model)) |>
    purrr::list_rbind()
}


#' Plot the estimates and standard errors of the model results.
#'
#' @param results The model results.
#'
#' @return A ggplot2 figure.
#'
create_plot_model_results <- function(results) {
  results |>
    dplyr::filter(term == "value", std.error <= 2, estimate <= 5) |>
    dplyr::select(metabolite, model, estimate, std.error) |>
    ggplot2::ggplot(ggplot2::aes(
      x = estimate,
      y = metabolite,
      xmin = estimate - std.error,
      xmax = estimate + std.error
    )) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed") +
    ggplot2::facet_grid(cols = ggplot2::vars(model))
}
