# Creedenzymatic for PTK - Striatum - Targeted

library(tidyverse)
library(creedenzymatic)

combined <- read_csv("results/STR-PTK-Combined_Tools.csv")

chosen_kinases <- c("CSK", "EGFR", "FAK", "RYK", "TRK", "VEGFR")

sig_kinases <- kinome_mp_file |>
  filter(krsa_id %in% chosen_kinases) |>
  pull(hgnc_symbol)

quartile_fig <- combined |>
  filter(hgnc_symbol %in% sig_kinases) |>
  quartile_figure() +
  guides(shape = "none") +
  scale_size_discrete(name = "Quartile")

quartile_fig |>
  ggsave("figures/str_ptk_targeted_quartile_figure.png", plot = _, width = 10, height = 2, units = "in", bg = "white")

quartile_fig |>
  ggsave("figures/str_ptk_targeted_quartile_figure.svg", plot = _, width = 10, height = 2, units = "in", bg = "white")
