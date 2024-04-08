# Load packages and data
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(rlang)
  library(themis)
  library(tidymodels)
  library(here)
})
source(here("src/funs.R"))
source(here("pipeline/sequential/0-setup_data.R"))

# Combined tuned sequential workflows across folds
seq_wflows <- readRDS(file.path(inputDir, paste0(dataset, "_wflows.rds")))
seq_wflow <- paste0(seq_wflows[nseq], "_s", nseq)

tune_wflows_files <- list.files(
  path = file.path(outputDir, "sequential", "tune_wflows", dataset),
  pattern = seq_wflow,
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

# Calculate average of per-class metrics across folds using one-vs-all predictions
per_class_metrics <- test_preds %>%
  nest(.by = fold_id) %>%
  mutate(metrics = map(
    data,
    ~ ova_metrics(
      x = .x,
      truth = class,
      estimate = .pred_class,
      metric_set = per_class_mset
    )
  ),
  .keep = "unused") %>%
  unnest(cols = metrics) %>%
  mutate(
    mean_estimate = mean(.estimate, na.rm = TRUE),
    .by = c(.metric, .estimator, class_group)
  ) %>%
  relocate(class_group, .after = mean_estimate) %>%
  filter(!grepl("non", class_group))

# Variable importance metrics
## Use model-specific metrics if available, otherwise calculate
## permutation-based variable importance
if (grepl("_(rf|xgb|mr)", best_wflow)) {
  vi_ranked <- test_results %>%
    map(~ {
      .x %>%
        extract_fit_parsnip() %>%
        vip::vi()
    }) %>%
    list_rbind(names_to = "fold_id") %>%
    summarize(Mean_Importance = mean(Importance), .by = "Variable") %>%
    mutate(Rank = dense_rank(-Mean_Importance)) %>%
    arrange(Rank)
} else if (grepl("_svm", best_wflow)) {
  svm_pfun <- function(object, newdata) {
    kernlab::predict(object, new_data = newdata)[[".pred_class"]]
  }
  set.seed(2024)
  vi_ranked <-
    map2(test_results, outer_folds, ~ {
      .x %>%
        extract_fit_parsnip() %>%
        vip::vi(
          train = analysis(.y),
          method = "permute",
          target = "class",
          metric = "accuracy",
          nsim = 5,
          pred_wrapper = svm_pfun
        )
    }) %>%
    list_rbind(names_to = "fold_id") %>%
    summarize(Mean_Importance = mean(Importance), .by = "Variable") %>%
    mutate(Rank = dense_rank(-Mean_Importance)) %>%
    arrange(Rank)
}

# Write best models to file
model_file <- file.path(
  outputDir,
  "sequential",
  "merge_results",
  dataset,
  paste0(seq_wflow, "_model_", dataset, ".rds")
)
saveRDS(best_model, model_file)

# Write per-class metrics to file
metrics_file <- file.path(
  outputDir,
  "sequential",
  "merge_results",
  dataset,
  paste0(seq_wflow, "_per_class_metrics_", dataset, ".rds")
)
saveRDS(per_class_metrics, metrics_file)

# Write variable importance to file
vi_file <- file.path(
  outputDir,
  "sequential",
  "merge_results",
  dataset,
  paste0(seq_wflow, "_vi_", dataset, ".rds")
)
saveRDS(vi_ranked, vi_file)
