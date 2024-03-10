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
  path = file.path(outputDir, "merge_results", dataset),
  pattern = "metrics",
  full.names = TRUE
)
all_wflow_metrics <- metrics_files %>%
  set_names(gsub("wflow_(.*)_metrics.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow")

all_metrics_file <- file.path(
  outputDir,
  "summarize_results",
  dataset,
  paste0("all_metrics_", dataset, ".rds")
)
saveRDS(all_wflow_metrics, all_metrics_file)

# Summarize all workflow variable importance ranks and write to file
vi_files <- list.files(
  path = file.path(outputDir, "merge_results", dataset),
  pattern = "vi",
  full.names = TRUE
)
all_wflow_vi <- vi_files %>%
  set_names(gsub("wflow_(.*)_vi.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow")

all_vi_file <- file.path(
  outputDir,
  "summarize_results",
  dataset,
  paste0("all_vi_", dataset, ".rds")
)
saveRDS(all_wflow_vi, all_vi_file)
