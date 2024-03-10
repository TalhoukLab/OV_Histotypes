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
source(here("pipeline/0-setup_data.R"))

# Combine tuned workflows across folds
wflow <- paste(samp, alg, sep = "_")
tune_wflows_files <- list.files(
  path = file.path(outputDir, "tune_wflows", dataset),
  pattern = wflow,
  full.names = TRUE
)

tuned_set <- tune_wflows_files %>%
  map(readRDS) %>%
  do.call(rbind, .)
tuned_set$wflow_id <- paste0(tuned_set$wflow_id, seq_len(n_folds))

# training_metrics <- tuned_set %>%
#   collect_metrics() %>%
#   summarize(
#     mean_estimate = mean(mean),
#     sd = sd(mean),
#     lower = mean_estimate - qnorm(0.975) * sd,
#     upper = mean_estimate + qnorm(0.975) * sd,
#     .by = .metric
#   )

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

# Don't interpret F-measure if some levels had no predicted events
overall_metrics <- tryCatch({
  test_preds %>%
    mset(truth = class, prob_cols, estimate = .pred_class) %>%
    add_column(class_group = "Overall")
}, warning = function(w) {
  overall_metrics %>%
    mutate(.estimate = ifelse(.metric == "f_meas", NA_real_, .estimate))
})

# Calculate per-class metrics using one-vs-all predictions
per_class_mset <- metric_set(accuracy, f_meas, kap, gmean)
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

# Variable importance metrics
## Use model-specific metrics if available, otherwise calculate
## permutation-based variable importance
if (grepl("_(rf|xgb|mr)", best_wflow)) {
  vi_ranked <- test_results %>%
    map(~ {
      .x %>%
        extract_fit_parsnip() %>%
        vip::vi(rank = TRUE)
    }) %>%
    list_rbind(names_to = "fold_id") %>%
    summarize(Mean_Importance = mean(Importance), .by = "Variable") %>%
    mutate(Mean_Importance = dense_rank(Mean_Importance)) %>%
    arrange(Mean_Importance)
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
          pred_wrapper = svm_pfun,
          rank = TRUE
        )
    }) %>%
    list_rbind(names_to = "fold_id") %>%
    summarize(Mean_Importance = mean(Importance), .by = "Variable") %>%
    mutate(Mean_Importance = dense_rank(Mean_Importance)) %>%
    arrange(Mean_Importance)
}

# Only consider candidate genes not already in PrOTYPE and SPOT
candidates <- c("C10orf116", "GAD1", "TPX2", "KGFLP2", "EGFL6", "KLK7", "PBX1",
                "LIN28B", "TFF3", "MUC5B", "FUT3", "STC1", "BCL2", "PAX8", "GCNT3",
                "GPR64", "ADCYAP1R1", "IGKC", "BRCA1", "IGJ", "TFF1", "MET",
                "CYP2C18", "CYP4B1", "SLC3A1", "EPAS1", "HNF1B", "IL6", "ATP5G3",
                "DKK4", "SENP8", "CAPN2", "C1orf173", "CPNE8", "IGFBP1", "WT1",
                "TP53", "SEMA6A", "SERPINA5", "ZBED1", "TSPAN8", "SCGB1D2", "LGALS4",
                "MAP1LC3A")

vi_ranked_candidates <- vi_ranked %>%
  filter(Variable %in% candidates)

# Write all metrics to file
metrics_file <- file.path(
  outputDir,
  "merge_results",
  dataset,
  paste0("wflow_", wflow, "_metrics_", dataset, ".rds")
)
saveRDS(all_metrics, metrics_file)

# Write variable importance to file
vi_file <- file.path(
  outputDir,
  "merge_results",
  dataset,
  paste0("wflow_", wflow, "_vi_", dataset, ".rds")
)
saveRDS(vi_ranked_candidates, vi_file)
