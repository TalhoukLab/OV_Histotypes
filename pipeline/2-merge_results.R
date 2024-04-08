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
  nest(.by = fold_id) %>%
  mutate(
    metrics = data %>%
      map(~ {
        .x %>%
          mset(truth = class, prob_cols, estimate = .pred_class) %>%
          suppressWarnings()
      }),
    .keep = "unused"
  ) %>%
  unnest(cols = metrics) %>%
  mutate(mean_estimate = mean(.estimate),
         .by = c(.metric, .estimator)) %>%
  add_column(class_group = "Overall")

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
  )

# Combine all metrics
all_metrics <-
  bind_rows(overall_metrics, per_class_metrics) %>%
  arrange(fold_id, .metric)

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
