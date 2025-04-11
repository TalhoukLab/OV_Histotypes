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
  path = file.path(outputDir, "retrain", "merge_results", dataset),
  pattern = "metrics",
  full.names = TRUE
)
all_wflow_metrics <- metrics_files %>%
  set_names(gsub("wflow_(.*)_metrics.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow")

all_metrics_file <- file.path(
  outputDir,
  "retrain",
  "summarize_results",
  dataset,
  paste0("all_metrics_", dataset, ".rds")
)
saveRDS(all_wflow_metrics, all_metrics_file)

# Summarize all workflow variable importance ranks and write to file
vi_files <- list.files(
  path = file.path(outputDir, "retrain", "merge_results", dataset),
  pattern = "vi",
  full.names = TRUE
)
all_wflow_vi <- vi_files %>%
  set_names(gsub("wflow_(.*)_vi.*", "\\1",  basename(.))) %>%
  map(readRDS) %>%
  list_rbind(names_to = "wflow")

all_vi_file <- file.path(
  outputDir,
  "retrain",
  "summarize_results",
  dataset,
  paste0("all_vi_", dataset, ".rds")
)
saveRDS(all_wflow_vi, all_vi_file)

# Find the top class of each sequence and generate a retraining subset for the
# next sequence.
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))
train_ref <- cbind(data, class = class)

# Number of classes used
n_class <- n_distinct(class)

# Top workflow by per-class F1-score out of n_class
seq_top <- all_wflow_metrics %>%
  distinct(pick(-fold_id, -.estimate)) %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = mean_estimate) %>%
  add_column(n_class = n_class)

saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))

# Create retrain data and class if more than 2 classes
if (n_class > 2) {
  retrain_ref <-
    train_ref %>% filter(!class %in% seq_top[["class_group"]])
  retrain_data <- retrain_ref %>% select(-class)
  retrain_class <- retrain_ref %>% pull(class)

  saveRDS(retrain_data, here("data", paste0("retrain_", n_class - 1, "_data.rds")))
  saveRDS(retrain_class, here("data", paste0("retrain_", n_class - 1, "_class.rds")))
}
