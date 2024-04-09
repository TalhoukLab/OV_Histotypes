# Load packages and data
suppressPackageStartupMessages({
  library(rlang)
  library(parallel)
  library(doParallel)
  library(future)
  library(themis)
  library(tidymodels)
  library(vip)
  library(here)
})
source(here("src/funs.R"))
source(here("pipeline/gene_opt/0a-setup_data.R"))
source(here("pipeline/0-model_specs.R"))

## Register parallel multicore
all_cores <- min(detectCores(logical = FALSE), 8L)
plan(multicore, workers = all_cores)

## Tune workflow across inner folds
inner_folds <- pluck(folds, "inner_resamples", as.numeric(fold_id))
tuned_set <- wflow_set %>%
  workflow_map(
    seed = 2024,
    resamples = inner_folds,
    grid = tuning_grid,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages() %>%
  suppressWarnings()

## Unregister parallel multicore
plan(sequential)

# Write tuned workflow to file
results_file <- file.path(
  outputDir,
  "gene_opt",
  "tune_wflows",
  dataset,
  paste0(top_wflow, "_add", ngene, "_", fold_id, "_", dataset, ".rds")
)
saveRDS(tuned_set, results_file)
