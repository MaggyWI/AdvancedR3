create_table_descriptive_stats <- function(data) {
  tablestat <- data |>
    group_by(metabolite) |>
    summarise(across(value, list(mean = mean, sd = sd))) |>
    mutate(across(where(is.numeric), \(x) round(x, digits = 1))) |>
    mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
    select(Metabolite = metabolite, `Mean SD` = MeanSD)
  return(tablestat)
}
