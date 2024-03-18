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
source(here("pipeline/gene_opt/0a-setup_data.R"))

# Combine tuned workflows across folds
tune_wflows_files <- list.files(
  path = file.path(outputDir, "gene_opt", "tune_wflows", dataset),
  pattern = paste0("add", ngene, "_"),
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

# Calculate overall metrics
overall_metrics <- test_preds %>%
  mset(truth = class, prob_cols, estimate = .pred_class) %>%
  add_column(class_group = "Overall") %>%
  suppressWarnings()

# Calculate per-class metrics using one-vs-all predictions
per_class_metrics <- test_preds %>%
  mutate(
    pred_class_ova = map(.pred_class, ~ {
      ifelse(levels(.pred_class) %in% .x, as.character(.x), "class_0") %>%
        set_names(paste0(".pred_class_", levels(.pred_class)))
    }),
    class_ova = map(class, ~ {
      ifelse(levels(class) %in% .x, as.character(.x), "class_0") %>%
        set_names(paste0("class_", levels(class)))
    })
  ) %>%
  unnest_wider(col = c(pred_class_ova, class_ova)) %>%
  pivot_longer(
    matches("^.pred_class_.*"),
    names_to = ".pred_class_group",
    names_prefix = ".pred_class_",
    values_to = ".pred_class_value"
  ) %>%
  pivot_longer(
    matches("^class_.*"),
    names_to = "class_group",
    names_prefix = "class_",
    values_to = "class_value"
  ) %>%
  filter(class_group == .pred_class_group) %>%
  nest(.by = class_group) %>%
  mutate(
    data = data %>%
      map(~ mutate(.x, across(matches("class_value"),
                              ~ factor(.x, levels = unique(c(.pred_class_group, "class_0")))))) %>%
      map(per_class_mset, truth = class_value, estimate = .pred_class_value) %>%
      suppressWarnings()
  ) %>%
  unnest(cols = data)

# Combine all metrics
all_metrics <- bind_rows(overall_metrics, per_class_metrics) %>%
  arrange(.metric)

# Top class by selected metric
top_class <- all_metrics %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = .estimate)

# Write all metrics to file
metrics_file <- file.path(
  outputDir,
  "gene_opt",
  "merge_results",
  dataset,
  paste0(top_wflow, "_add", ngene, "_metrics_", dataset, ".rds")
)
saveRDS(all_metrics, metrics_file)
