# Creedenzymatic for STK - Striatum

library(tidyverse)
library(creedenzymatic)

krsa_df <-
  read_csv("results/STR-krsa_table_full_Exer_STR_CTL_STR_STK.csv") |>
  select(Kinase, Score = AvgZ) |>
  read_krsa(trns = "abs", sort = "desc")

uka_df <- read_tsv("kinome_data/UKA-STK/STR/Summaryresults 20230404-1531.txt") |> select(Kinase = `Kinase Name`, Score = `Median Final score`) |>
  unique() |>
  read_uka(trns = "abs", sort = "desc")

diff_peps <- read_csv("results/STR-dpp_Exer_STR_CTL_STR-STK.csv") |>
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
  write_csv("results/STR-STK-Combined_Tools.csv")

sig_kinases <- combined |>
  filter(Perc >= 0.9) |>
  pull(hgnc_symbol) |>
  unique()

combined |>
  filter(hgnc_symbol %in% sig_kinases) |>
  quartile_figure()
