# Load packages and data
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidymodels)
  library(here)
})
source(here("src/funs.R"))

# Summarize all workflow metrics and write to file
metrics_files <- list.files(
  path = file.path(outputDir, "sequential", "merge_results", dataset),
  pattern = "metrics",
  full.names = TRUE
)
all_wflow_metrics <- metrics_files %>%
  set_names(gsub("wflow_(.*)_metrics.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow")

all_metrics_file <- file.path(
  outputDir,
  "sequential",
  "summarize_results",
  dataset,
  paste0("all_metrics_", dataset, ".rds")
)
saveRDS(all_wflow_metrics, all_metrics_file)

# Summarize all workflow variable importance ranks and write to file
vi_files <- list.files(
  path = file.path(outputDir, "sequential", "merge_results", dataset),
  pattern = "vi",
  full.names = TRUE
)
vi_ranked <- vi_files %>%
  set_names(gsub("(.*)_vi.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow") %>%
  summarize(Aggregated_Rank = mean(Mean_Importance), .by = Variable) %>%
  arrange(Aggregated_Rank)

# Only consider candidate genes not already in PrOTYPE and SPOT
candidates <- c("C10orf116", "GAD1", "TPX2", "KGFLP2", "EGFL6", "KLK7", "PBX1",
                "LIN28B", "TFF3", "MUC5B", "FUT3", "STC1", "BCL2", "PAX8", "GCNT3",
                "GPR64", "ADCYAP1R1", "IGKC", "BRCA1", "IGJ", "TFF1", "MET",
                "CYP2C18", "CYP4B1", "SLC3A1", "EPAS1", "HNF1B", "IL6", "ATP5G3",
                "DKK4", "SENP8", "CAPN2", "C1orf173", "CPNE8", "IGFBP1", "WT1",
                "TP53", "SEMA6A", "SERPINA5", "ZBED1", "TSPAN8", "SCGB1D2", "LGALS4",
                "MAP1LC3A")

vi_ranked_candidates <- vi_ranked %>%
  filter(Variable %in% candidates)

all_vi_file <- file.path(
  outputDir,
  "sequential",
  "summarize_results",
  dataset,
  paste0("all_vi_", dataset, ".rds")
)
saveRDS(vi_ranked_candidates, all_vi_file)
