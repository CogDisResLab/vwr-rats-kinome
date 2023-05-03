# Creedenzymatic for STK - Striatum - Targeted

library(tidyverse)
library(creedenzymatic)

combined <- read_csv("results/STR-STK-Combined_Tools.csv")

chosen_kinases <- c("PKD", "RIPK", "IKK", "COT", "CAMK2", "AKT2")

sig_kinases <- kinome_mp_file |>
  filter(krsa_id %in% chosen_kinases) |>
  pull(hgnc_symbol)

quartile_fig <- combined |>
  filter(hgnc_symbol %in% sig_kinases) |>
  quartile_figure() +
  guides(shape = "none") +
  scale_size_discrete(name = "Quartile")

quartile_fig |>
  ggsave("figures/str_stk_targeted_quartile_figure.png", plot = _, width = 8.5, height = 3, units = "in", bg = "white")

quartile_fig |>
  ggsave("figures/str_stk_targeted_quartile_figure.pdf", plot = _, width = 8.5, height = 3, units = "in", bg = "white")
