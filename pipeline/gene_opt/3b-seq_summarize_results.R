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
  path = file.path(outputDir, "gene_opt", "sequential", "merge_results", dataset),
  full.names = TRUE
)

all_wflow_metrics <- metrics_files %>%
  set_names(gsub("_metrics.*", "\\1", basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow_genes") %>%
  separate(wflow_genes, c("wflow", "Genes"), sep = "_add", convert = TRUE) %>%
  separate(wflow, c("wflow", "Sequence"), sep = "_s(?=[^_]+$)", convert = TRUE) %>%
  arrange(Sequence, Genes)

all_metrics_file <- file.path(
  outputDir,
  "gene_opt",
  "sequential",
  "summarize_results",
  dataset,
  paste0("gene_opt_all_metrics_", dataset, ".rds")
)
saveRDS(all_wflow_metrics, all_metrics_file)
