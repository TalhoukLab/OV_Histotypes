# Load packages and data
suppressPackageStartupMessages({
  library(dplyr)
  library(tidymodels)
  library(here)
})
source(here("src/funs.R"))

# Outer fold used for assessment
id <- as.numeric(fold_id)
outer_fold <- folds$splits[[id]]

# Combine tuned workflows for each outer fold
tuned_set <- list.files(
  path = file.path(outputDir,
                   "tune_wflows",
                   dataset),
  pattern = paste0("_", fold_id, "_", dataset),
  full.names = TRUE
) %>%
  map(readRDS) %>%
  list_rbind()

# Ranking workflows for selected metric
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

# Train best model on training set, evaluate on test set
test_results <- best_model %>%
  last_fit(split = outer_fold, metrics = mset)

# Test set predictions
test_preds <- test_results %>%
  collect_predictions()

# Calculate overall metrics
overall_metrics <- test_results %>%
  collect_metrics() %>%
  select(-.config) %>%
  add_column(class_group = "Overall")

# Calculate per-class metrics using one-vs-all predictions
per_class_mset <- metric_set(accuracy, f_meas, kap, gmean)
per_class_metrics <- test_results %>%
  collect_predictions() %>%
  mutate(
    pred_max = pmap_dbl(select(., matches(".pred") & where(is.double)),
                        ~ pmap_dbl(list(...), max)),
    across(matches(".pred") & where(is.double),
           ~ factor(
             ifelse(.x == pred_max, as.character(.pred_class), "class_0")
           ),
           .names = "{gsub('.pred_', '.pred_class_', {col}, fixed = TRUE)}"),
    across(matches(".pred_class_.*"),
           ~ factor(
             ifelse(class %in% levels(.x), as.character(class), "class_0")
           ),
           .names = "{gsub('.pred_', '', {col}, fixed = TRUE)}")
  ) %>%
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
      map(droplevels) %>%
      map(per_class_mset, truth = class_value, estimate = .pred_class_value)
  ) %>%
  unnest(cols = data)

# Combine all metrics
all_metrics <- bind_rows(overall_metrics, per_class_metrics) %>%
  arrange(.metric)

# Top class by selected metric
top_class <- all_metrics %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = .estimate)

# Variable importance metrics
## Use model-specific metrics if available, otherwise calculate
## permutation-based variable importance
if (grepl("_(rf|xgb|mr)", best_wflow)) {
  vi_df <- test_results %>%
    extract_fit_parsnip() %>%
    vi()
} else if (grepl("_svm", best_wflow)) {
  svm_pfun <- function(object, newdata) {
    kernlab::predict(object, new_data = newdata)[[".pred_class"]]
  }
  set.seed(2024)
  vi_df <- test_results %>%
    extract_fit_parsnip() %>%
    vi(
      train = analysis(outer_fold),
      method = "permute",
      target = "class",
      metric = "accuracy",
      nsim = 5,
      pred_wrapper = svm_pfun
    )
}
