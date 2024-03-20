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
source(here("pipeline/gene_opt/0b-setup_seq_data.R"))

# Combined tuned sequential workflows across folds
seq_wflows <- readRDS(file.path(inputDir, paste0(dataset, "_wflows.rds")))
seq_wflow <- paste0(seq_wflows[nseq], "_s", nseq)

tune_wflows_files <- list.files(
  path = file.path(outputDir, "gene_opt", "sequential", "tune_wflows", dataset),
  pattern = paste0(seq_wflow, "_add", ngene, "_"),
  full.names = TRUE
)

tuned_set <- tune_wflows_files %>%
  map(readRDS) %>%
  do.call(rbind, .)
tuned_set$wflow_id <- paste0(tuned_set$wflow_id, seq_len(n_folds))

ranking <- tuned_set %>%
  rank_results(rank_metric = rank_metric, select_best = TRUE) %>%
  filter(.metric == rank_metric)

# Best workflow
best_wflow <- ranking[["wflow_id"]][1]

# Best tuning parameters of workflow with highest metric
best_params <- tuned_set %>%
  extract_workflow_set_result(best_wflow) %>%
  select_best(metric = rank_metric)

# Finalize model with best set of tuning parameters
best_model <- tuned_set %>%
  extract_workflow(best_wflow) %>%
  finalize_workflow(best_params)

# Train best model on outer training sets, evaluate on outer test sets
outer_folds <- folds$splits
test_results <- map(outer_folds, ~ {
  suppressMessages(last_fit(best_model, split = .x, metrics = mset))
})

# Test set predictions
test_preds <- test_results %>%
  map(collect_predictions) %>%
  list_rbind(names_to = "fold_id")

# Select predicted probability columns for multiclass/binary AUC
prob_cols <- test_preds %>%
  select(matches(".pred(?!_class)", perl = TRUE)) %>%
  names()
if (n_distinct(class) == 2) {
  prob_cols <- prob_cols[1]
}

# Calculate average of overall metrics across folds
overall_metrics <- test_preds %>%
  nest(.by = "fold_id") %>%
  mutate(
    metrics = data %>%
      map(~ {
        .x %>%
          mset(truth = class, prob_cols, estimate = .pred_class) %>%
          suppressWarnings()
      })
  ) %>%
  unnest(cols = metrics) %>%
  summarise(mean_estimate = mean(.estimate),
            .by = c(.metric, .estimator)) %>%
  add_column(class_group = "Overall") %>%
  arrange(.metric)

# Write all metrics to file
metrics_file <- file.path(
  outputDir,
  "gene_opt",
  "sequential",
  "merge_results",
  dataset,
  paste0(seq_wflow, "_add", ngene, "_metrics_", dataset, ".rds")
)
saveRDS(overall_metrics, metrics_file)
