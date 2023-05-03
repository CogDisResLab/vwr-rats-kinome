# Creedenzymatic for STK - Hippocampus

library(tidyverse)
library(creedenzymatic)

krsa_df <-
  read_csv("results/HPC-krsa_table_full_Exer_HPC_CTL_HPC_STK.csv") |>
  select(Kinase, Score = AvgZ) |>
  unique() |>
  read_krsa(trns = "abs", sort = "desc")

uka_df <- read_tsv("kinome_data/UKA-STK/HPC/Summaryresults 20230404-1530.txt") |> select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
  unique() |>
  read_uka(trns = "abs", sort = "desc")

diff_peps <- read_csv("results/HPC-dpp_Exer_HPC_CTL_HPC-STK.csv") |>
  select(Peptide, Score = totalMeanLFC) |>
  unique()

kea_df <-
  read_kea(
    diff_peps,
    filter = T,
    cutoff = 0.4,
    cutoff_abs = T,
    sort = "asc",
    trns = "abs",
    rm_duplicates = T,
    method = "MeanRank",
    lib = "kinase-substrate"
  )

ptmsea_df <- diff_peps |>
  read_ptmsea()

combined <- combine_tools(krsa_df, uka_df, kea_df, ptmsea_df) |>
  write_csv("results/HPC-STK-Combined_Tools.csv")

top_kinases <- krsa_df |>
  filter(Score >= 1.5) |>
  pull(Kinase)

top_kinase_symbols <- kinome_mp_file |>
  filter(krsa_id %in% top_kinases) |>
  pull(hgnc_symbol)

sig_kinases <- combined |>
  filter(Perc >= 0.95) |>
  pull(hgnc_symbol) |>
  unique()

quartile_fig <- combined |>
  filter(hgnc_symbol %in% sig_kinases) |>
  quartile_figure() |>
  ggsave("figures/hpc_stk_quartile_figure.png", plot = _, width = 18, height = 6, units = "in", bg = "white")
