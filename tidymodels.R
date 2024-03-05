# Server constants
inputDir <- "data"
dataset <- "train"
fold_id <- 1
rank_metric <- "accuracy"
alg <- "rf"
samp <- "hybrid"

# Load packages and data
suppressPackageStartupMessages({
  library(rlang)
  library(parallel)
  library(doParallel)
  library(themis)
  library(tidymodels)
  library(vip)
  library(here)
})
source(here("validation/cs_classifier.R"))

# Import training data and class labels
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))
train_ref <- cbind(data, class = factor(class))

# Nested resampling: 10-fold CV outside, 5-fold CV inside
set.seed(2024)
folds <- nested_cv(
  train_ref,
  outside = vfold_cv(v = 10, strata = class),
  inside = vfold_cv(v = 5, strata = class)
)

inner_folds <- folds$inner_resamples[[fold_id]]
outer_fold <- folds$splits[[fold_id]]

# Metrics
gmean <- new_class_metric(gmean, direction = "maximize")
mset <- metric_set(accuracy, roc_auc, f_meas, kap, gmean)

# Recipe
rec <- recipe(class ~ ., train_ref)

# Preprocessing
## Class imbalance subsampling methods
none_samp <- rec
down_samp <- step_downsample(rec, class, seed = 2024)
up_samp <- step_upsample(rec, class, seed = 2024)
smote_samp <- step_smote(rec, class, seed = 2024)
hybrid_samp <- rec %>%
  step_smote(class, over_ratio = 0.5, seed = 2024) %>%
  step_downsample(class, under_ratio = 1, seed = 2024)

preproc <- list(
  none = none_samp,
  down = down_samp,
  up = up_samp,
  smote = smote_samp,
  hybrid = hybrid_samp
)

# Models
## Random forest
rf_model <-
  rand_forest(
    mode = "classification",
    engine = "ranger",
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
  set_engine("ranger", importance = "impurity")

## XGBoost
xgb_model <-
  boost_tree(
    mode = "classification",
    engine = "xgboost",
    mtry = tune(),
    trees = 500,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    stop_iter = tune()
  )

## Support vector machine
svm_model <-
  svm_rbf(
    mode = "classification",
    engine = "kernlab",
    cost = tune(),
    rbf_sigma = tune()
  )

## Multinomial regression
mr_model <-
  multinom_reg(
    engine = "glmnet",
    penalty = tune(),
    mixture = tune()
  )

models <- list(rf = rf_model,
               xgb = xgb_model,
               svm = svm_model,
               mr = mr_model)

# Workflow sets
wflow_sets <- workflow_set(preproc, models)

# Hyperparameter tuning

## Algorithm-specific tuning setup
if (alg %in% c("rf", "xgb")) {
  wflow_set <- wflow_sets %>%
    filter(wflow_id == paste(samp, alg, sep = "_"))

  tuning_grid <- 10

} else if (alg %in% "svm") {
  svm_params <-
    parameters(cost(), rbf_sigma()) %>%
    update(
      cost = cost(c(0, 2)),
      rbf_sigma = rbf_sigma(c(-3, 0))
    )

  wflow_set <- wflow_sets %>%
    filter(wflow_id == paste(samp, alg, sep = "_")) %>%
    option_add(param_info = svm_params)

  tuning_grid <- 10

} else if (alg %in% "mr") {
  wflow_set <- wflow_sets %>%
    filter(wflow_id == paste(samp, alg, sep = "_"))

  set.seed(2024)
  tuning_grid <- crossing(
    grid_latin_hypercube(penalty(), size = 10),
    grid_regular(mixture(), levels = 10)
  )
}

## Register parallel multicore
all_cores <- min(detectCores(logical = FALSE), 32L)
cl <- makePSOCKcluster(all_cores)
registerDoParallel(cl)

## Tune workflow
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
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list = ls(name = env), pos = env)
}
unregister_dopar()

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
if (alg %in% c("rf", "xgb", "mr")) {
  vi_df <- test_results %>%
    extract_fit_parsnip() %>%
    vi()
} else if (alg %in% "svm") {
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
