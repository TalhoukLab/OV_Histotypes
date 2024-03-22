# Load packages and data
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidymodels)
  library(here)
})
source(here("src/funs.R"))

# Combine all best models and write to file
model_files <- list.files(
  path = file.path(outputDir, "sequential", "merge_results", dataset),
  pattern = "model",
  full.names = TRUE
)
model_order <- order(gsub(".*_s([0-9])_model.*", "\\1", model_files))
all_models <- model_files[model_order] %>%
  set_names(gsub("_model.*", "", basename(.))) %>%
  map(readRDS)

all_models_file <- file.path(
  outputDir,
  "sequential",
  "summarize_results",
  dataset,
  paste0("all_models_", dataset, ".rds")
)
saveRDS(all_models, all_models_file)

# Combine per-class metrics and write to file
metrics_files <- list.files(
  path = file.path(outputDir, "sequential", "merge_results", dataset),
  pattern = "per_class_metrics",
  full.names = TRUE
)
per_class_metrics <- metrics_files %>%
  set_names(gsub("_metrics.*", "", basename(.))) %>%
  map(readRDS) %>%
  list_rbind()

per_class_metrics_file <- file.path(
  outputDir,
  "sequential",
  "summarize_results",
  dataset,
  paste0("per_class_metrics_", dataset, ".rds")
)
saveRDS(per_class_metrics, per_class_metrics_file)

# Summarize all workflow variable importance ranks and write to file
vi_files <- list.files(
  path = file.path(outputDir, "sequential", "merge_results", dataset),
  pattern = "vi",
  full.names = TRUE
)
vi_ranked <- vi_files %>%
  set_names(gsub("_vi.*", "",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow") %>%
  separate(wflow, c("wflow", "Sequence"), sep = "_s(?=[^_]+$)", convert = TRUE) %>%
  arrange(Sequence) %>%
  nest(.by = Rank) %>%
  mutate(Gene_Lists = map(data, "Variable"),
         Gene_Order = accumulate(Gene_Lists, union)) %>%
  tail(1) %>%
  select(Gene_Order) %>%
  unnest(Gene_Order)

# Only consider candidate genes not already in PrOTYPE and SPOT
candidates <- c("C10orf116", "GAD1", "TPX2", "KGFLP2", "EGFL6", "KLK7", "PBX1",
                "LIN28B", "TFF3", "MUC5B", "FUT3", "STC1", "BCL2", "PAX8", "GCNT3",
                "GPR64", "ADCYAP1R1", "IGKC", "BRCA1", "IGJ", "TFF1", "MET",
                "CYP2C18", "CYP4B1", "SLC3A1", "EPAS1", "HNF1B", "IL6", "ATP5G3",
                "DKK4", "SENP8", "CAPN2", "C1orf173", "CPNE8", "IGFBP1", "WT1",
                "TP53", "SEMA6A", "SERPINA5", "ZBED1", "TSPAN8", "SCGB1D2", "LGALS4",
                "MAP1LC3A")

vi_ranked_candidates <- vi_ranked %>%
  filter(Gene_Order %in% candidates)

all_vi_file <- file.path(
  outputDir,
  "sequential",
  "summarize_results",
  dataset,
  paste0("all_vi_", dataset, ".rds")
)
saveRDS(vi_ranked_candidates, all_vi_file)
