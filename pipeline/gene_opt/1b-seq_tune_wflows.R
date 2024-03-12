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

# Inner folds used for tuning
id <- as.numeric(fold_id)
inner_folds <- folds$inner_resamples[[id]]

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

# Generate workflow set and select specified workflow
wflow_sets <- workflow_set(preproc, models)
sq_wflows <- readRDS(file.path(inputDir, paste0(dataset, "_wflows.rds")))
wflow <- sq_wflows[[nseq]]
wflow_set <- wflow_sets %>% filter(wflow_id == wflow)

# Hyperparameter tuning

## Algorithm-specific tuning setup
if (grepl("rf|xgb", wflow)) {
  tuning_grid <- 10

} else if (grepl("svm", wflow)) {
  svm_params <-
    parameters(cost(), rbf_sigma()) %>%
    update(
      cost = cost(c(0, 2)),
      rbf_sigma = rbf_sigma(c(-3, 0))
    )

  wflow_set <- wflow_set %>%
    option_add(param_info = svm_params)

  tuning_grid <- 10

} else if (grepl("mr", wflow)) {
  set.seed(2024)
  tuning_grid <- crossing(
    grid_latin_hypercube(penalty(), size = 10),
    grid_regular(mixture(), levels = 10)
  )
}

## Register parallel multicore
all_cores <- min(detectCores(logical = FALSE), 8L)
plan(multicore, workers = all_cores)

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
plan(sequential)

# Write tuned workflow to file
results_file <- file.path(
  outputDir,
  "gene_opt",
  "sequential",
  "tune_wflows",
  dataset,
  paste0(wflow, "_s", nseq, "_add", ngene, "_", fold_id, "_", dataset, ".rds")
)
saveRDS(tuned_set, results_file)
