# Generate targeted reverse KRSA figures for the paper
# Good morning Ali,
# I am currently working on the discussion portion of the Dr. Yuan paper. One of the paragraphs will be on
# Down regulation of STK vs. upregulation of PTK
# Would this be abundance plots? If so, would you be able to send me the reverse krsa plots chosen in the creedenzymatic. These are:
# CAMK2
# IKK
# PKD
# RIPK
# MAPK
# CSK
# EGFR
# FAK
# RYK
# TRK
# VEGFR
# If there are no interesting findings, is there a way to look at all
# the kinases and identify opposite directions in activity? Thank you.

# Load the necessary libraries
library(tidyverse)
library(KRSA)

kinases_of_interest <- c(
  "CAMK2", "IKK", "PKD", "RIPK", "JNK", "ERK", "P38",
  "CSK", "EGFR", "FAK", "RYK", "TRK", "VEGFR"
)

# Load the data
dpp_data <- list.files(path = "results", pattern = "dpp", full.names = TRUE) |>
  set_names(~ str_extract(.x, "(STR|HPC)_CTL(_(HPC|STR))?-[SP]TK") |>
    str_remove("_CTL(_(HPC|STR))?")) |>
  map(read_csv) |>
  bind_rows(.id = "Comparison") |>
  separate_wider_delim(Comparison, names = c("Tissue", "Chip"), delim = "-") |>
  select(Barcode, Tissue, Chip, Peptide, LFC)

stk_mapping <- KRSA_Mapping_STK_PamChip_87102_v1 |>
  rename(Peptide = Substrates, Kinase = Kinases) |>
  separate_longer_delim(Kinase, delim = " ") |>
  unique()

ptk_mapping <- KRSA_Mapping_PTK_PamChip_86402_v1 |>
  rename(Peptide = Substrates, Kinase = Kinases) |>
  separate_longer_delim(Kinase, delim = " ") |>
  unique()

chip_mapping <- stk_mapping |>
  bind_rows(ptk_mapping) |>
  filter(Kinase %in% kinases_of_interest)

dpp_kinase_data <- dpp_data |>
  inner_join(chip_mapping, by = "Peptide", relationship = "many-to-many") |>
  mutate(Significant = LFC < -0.2 | LFC > 0.2)

ordered_kinases <- dpp_kinase_data |>
  select(Chip, Kinase) |>
  unique() |>
  group_by(Chip) |>
  mutate(Kinase = sort(Kinase)) |>
  pull(Kinase)

hpc_g <- dpp_kinase_data |>
  filter(Tissue == "HPC") |>
  ggplot(aes(x = Kinase, y = LFC, color = Significant))

hpc_p <- hpc_g +
  # geom_boxplot() +
  geom_point(position = position_jitter(
    width = 0.1, seed = 1989L
  )) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  geom_hline(yintercept = c(0.2, -0.2), linetype = "dashed") +
  theme_linedraw() +
  xlab("Kinase Family") + ylab("Fold Change") +
  theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
  guides(color = "none") +
  facet_grid(cols = vars(Chip), scales = "free_x")

ggsave("figures/targeted_reverse_krsa_hpc.png", hpc_p, width = 12, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("figures/targeted_reverse_krsa_hpc.svg", hpc_p, width = 12, height = 8, units = "in", dpi = 300)

str_g <- dpp_kinase_data |>
  filter(Tissue == "STR") |>
  ggplot(aes(x = Kinase, y = LFC, color = Significant))

str_p <- str_g +
  # geom_boxplot() +
  geom_point(position = position_jitter(
    width = 0.1, seed = 1989L
  )) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  geom_hline(yintercept = c(0.2, -0.2), linetype = "dashed") +
  theme_linedraw() +
  xlab("Kinase Family") + ylab("Fold Change") +
  theme(axis.text.x = element_text(angle = 90L, hjust = 1L)) +
  guides(color = "none") +
  facet_grid(cols = vars(Chip), scales = "free_x")

ggsave("figures/targeted_reverse_krsa_str.png", hpc_p, width = 12, height = 8, units = "in", dpi = 300, bg = "white")
ggsave("figures/targeted_reverse_krsa_str.svg", hpc_p, width = 12, height = 8, units = "in", dpi = 300)

g <- dpp_kinase_data |>
mutate(Tissue = factor(Tissue, levels = c("HPC", "STR"), labels = c("Hippocampus", "Dorsal Striatum")),
Chip = factor(Chip, levels = c("PTK", "STK"), labels = c("Tyrosine Kinases", "Serine/Threonine Kinases"))) |>
  ggplot(aes(x = Kinase, y = LFC, color = Significant))

p <- g +
  # geom_boxplot() +
  geom_point(position = position_jitter(
    width = 0.1, seed = 1989L
  ), size = 3) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black")
  ) +
  scale_y_continuous(limits = c(-1, 1), oob = scales::squish) +
  geom_hline(yintercept = c(0.2, -0.2), linetype = "dashed") +
  theme_linedraw() +
  xlab("Kinase Family") + ylab("Fold Change") +
  theme(text = element_text(size = 24)) +
  guides(color = "none") +
  facet_grid(rows = vars(Tissue), cols = vars(Chip), scales = "free_x")

ggsave("figures/targeted_reverse_krsa_grid.png", p, width = 20, height = 16, units = "in", dpi = 1200, bg = "white")
ggsave("figures/targeted_reverse_krsa_grid.svg", p, width = 20, height = 16, units = "in", dpi = 1200)
