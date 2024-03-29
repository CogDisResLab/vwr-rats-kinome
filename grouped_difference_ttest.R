# Perfrom group-wise t-test on the difference between two groups


library(tidyverse)


datafiles <- list.files("results",
  pattern = "signal.*.csv", full.names = TRUE
) |>
  set_names(~ basename(.x) |>
    str_extract("signal_(.*-[SP]TK)\\.csv", 1))

significant_peptides <- list.files(
  "results",
  pattern = "significant_peptides",
  full.names = TRUE
) |>
  set_names(~ basename(.x) |>
    str_extract("significant_peptides_0.2_(.*-[SP]TK)\\.txt", 1)) |>
  map(read_lines)

computed_data <- datafiles |>
  imap(~ read_csv(.x) |>
    filter(Peptide %in% significant_peptides[[.y]]) |>
    select(Group, Peptide, slope) |>
    group_by(Group, Peptide) |>
    summarise(slope = mean(slope), .groups = "drop")) |>
  bind_rows(.id = "Comparison") |>
  nest(.by = Comparison) |>
  mutate(
    ttest = map(data, ~ t.test(.$slope ~ .$Group, data = .x)),
    ttest = map(ttest, broom::tidy),
    correlation = map(data, ~ .x |>
      select(Group, slope) |>
      pivot_wider(names_from = Group, values_from = slope) |>
      unnest() |>
      cor.test()),
    significant_peptides = case_when(
      str_detect(Comparison, "STK") ~ map_dbl(data, ~ nrow(.x) / 2),
      str_detect(Comparison, "PTK") ~ map_dbl(data, ~ nrow(.x) / 2),
      .default = NA
    )
  ) |>
  unnest_wider(ttest) |>
  select(-data)
