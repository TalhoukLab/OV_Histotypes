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

# Nested resampling
set.seed(2024)
folds <- nested_cv(
  train_ref,
  outside = vfold_cv(v = folds, strata = class),
  inside = vfold_cv(v = 5, strata = class)
)

id <- as.numeric(fold_id)
inner_folds <- folds$inner_resamples[[id]]
outer_fold <- folds$splits[[id]]

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

# Write tuned workflow to file
results_file <- file.path(
  outputDir,
  "tune_wflows",
  dataset,
  paste0(samp, "_", alg, "_", fold_id, "_", dataset, ".rds")
)
saveRDS(tuned_set, results_file)
