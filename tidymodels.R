library(here)
library(rlang)
library(parallel)
library(themis)
library(rsample)
library(recipes)
library(parsnip)
library(yardstick)
library(workflowsets)
library(tune)
source(here("validation/cs_classifier.R"))

train_ref <- cbind(train_data, class = factor(train_class))

data_split <- bootstraps(data = train_ref,
                         times = 100,
                         strata = class)
train_data_list <- map(data_split[["splits"]], training)

set.seed(2024)
train_folds <- map(train_data_list, vfold_cv, v = 3, strata = class)
rank_metric <- "roc_auc"

gmean <- new_class_metric(gmean, direction = "maximize")
mset <- metric_set(accuracy, roc_auc, f_meas, kap, gmean)

# Recipe
id <- 1
rec <- recipe(class ~ ., train_data_list[[id]])

# Preprocessing for class imbalance
preproc_none <- rec
preproc_down <- step_downsample(rec, class, seed = 2024)
preproc_up <- step_upsample(rec, class, seed = 2024)
preproc_smote <- step_smote(rec, class, seed = 2024)
preproc <- list(none = preproc_none,
                down = preproc_down,
                up = preproc_up,
                smote = preproc_smote)

# Models
## Random forest
rf_model <- list(
  rf = rand_forest(
    mode = "classification",
    engine = "ranger",
    mtry = tune(),
    min_n = tune(),
    trees = 500
  ) %>%
    set_engine("ranger", importance = "impurity")
)

## XGBoost
xgb_model <- list(
  xgb = boost_tree(
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
)

# Support vector machine
svm_model <- list(
  svm = svm_rbf(
    mode = "classification",
    engine = "kernlab",
    cost = tune(),
    rbf_sigma = tune()
  )
)

# Multinomial logistic regression
mlr_model <- list(
  mlr = multinom_reg(
    engine = "glmnet",
    penalty = tune(),
    mixture = tune()
  )
)

# Workflow sets
rf_set <- workflow_set(preproc, rf_model)
xgb_set <- workflow_set(preproc, xgb_model)
svm_set <- workflow_set(preproc, svm_model)
mlr_set <- workflow_set(preproc, mlr_model)
wflow_sets <- bind_rows(rf_set, xgb_set, svm_set, mlr_set)

# Register parallel multicore
# all_cores <- min(detectCores(logical = FALSE), 8L)
# cl <- makePSOCKcluster(all_cores)
# registerDoParallel(cl)

# Hyperparameter tuning using 10 space-filling parameter grids run in parallel

## Random Forest
rf_tuned_set <- rf_set %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

rf_tuned_set2 <- wflow_sets %>%
  filter(grepl("rf", wflow_id)) %>%
  slice(2:4) %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

## XGBoost
xgb_tuned_set <- xgb_set %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()


## SVM
svm_params <-
  parameters(cost(), rbf_sigma()) %>%
  update(
    cost = cost(c(0, 2)),
    rbf_sigma = rbf_sigma(c(-3, 0))
  )

svm_tuned_set <- svm_set  %>%
  option_add(param_info = svm_params) %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

# wflow_tuned_set <- wflow_sets %>%
#   workflow_map(
#     seed = 2024,
#     resamples = train_folds[[as.numeric(id)]],
#     grid = 10,
#     metrics = mset,
#     control = control_grid(save_pred = TRUE, save_workflow = TRUE)
#   ) %>%
#   suppressMessages()

# Ranking workflows for selected metric
ranking <- wflow_tuned_set %>%
  rank_results(rank_metric = rank_metric) %>%
  filter(.metric == rank_metric) %>%
  select(model, wflow_id, f_meas = mean, rank)

# Best workflow
best_wflow <- head(ranking[["wflow_id"]], 1)

# Best tuning parameters of workflow with highest metric
best_params <- wflow_tuned_set %>%
  extract_workflow_set_result(best_wflow) %>%
  select_best(metric = rank_metric)

# Finalize model with best set of tuning parameters
best_model <- wflow_tuned_set %>%
  extract_workflow(best_wflow) %>%
  finalize_workflow(best_params)

# Train best model on training set, evaluate on test set
test_results <- best_model %>%
  last_fit(split = data_split$splits[[idx]], metrics = mset)
