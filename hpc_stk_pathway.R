# Pathway Analysis STK - HPC Data

library(tidyverse)
library(enrichR)
library(writexl)

# Load the reference data
stk_id <- readRDS("reference_data/stk_id_map.Rds")
stk_hgnc <- readRDS("reference_data/stk_hgnc_map.Rds")

stk_map <- stk_id |>
    inner_join(stk_hgnc)

# Copy this section for each comparison
#################

# Change the name of the file you are loading
comparison <- read_csv("results/HPC-dpp_Exer_HPC_CTL_HPC-STK.csv") |>
    group_by(Peptide) |>
    filter(abs(LFC) == max(abs(LFC))) |>
    ungroup() |>
    select(Peptide, LFC) |>
    inner_join(stk_map) |>
    # Change the name of the file you are writing.
    write_csv("results/annotated_dpp_Exer_HPC_CTL_HPC-STK.csv")

genes <- comparison |>
    select(Gene, LFC) |>
    pull(Gene)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021",
         "GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")

enriched <- enrichr(genes, dbs) |>
    imap(~ write_csv(.x, str_glue("results/Exer_HPC_CTL_HPC-STK-{.y}-Pathways.csv")))

# Change the name of the file you are writing
write_xlsx(enriched, "results/HPC-STK-Pathways.xlsx")

#################
