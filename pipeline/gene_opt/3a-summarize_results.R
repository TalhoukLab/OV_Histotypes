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
  path = file.path(outputDir, "gene_opt", "merge_results", dataset),
  full.names = TRUE
)
all_wflow_metrics <- metrics_files %>%
  set_names(gsub(".*add(.*)_metrics.*", "\\1", basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "Genes") %>%
  mutate(Genes = as.integer(Genes))

all_metrics_file <- file.path(
  outputDir,
  "gene_opt",
  "summarize_results",
  dataset,
  paste0("gene_opt_all_metrics_", dataset, ".rds")
)
saveRDS(all_wflow_metrics, all_metrics_file)
