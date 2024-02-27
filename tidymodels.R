library(rlang)
library(parallel)
library(themis)
library(tidymodels)
library(here)
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
rf_tuned_set <- wflow_sets %>%
  filter(grepl("rf", wflow_id)) %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

## XGBoost
xgb_tuned_set <- wflow_sets %>%
  filter(grepl("xgb", wflow_id)) %>%
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

svm_tuned_set <- wflow_sets %>%
  filter(grepl("svm", wflow_id)) %>%
  option_add(param_info = svm_params) %>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = 10,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

### MLR
set.seed(2024)
mlr_grid <- crossing(
  grid_latin_hypercube(penalty(), size = 10),
  grid_regular(mixture(), levels = 10)
)

mlr_tuned_set <- wflow_sets %>%
  filter(grepl("mlr", wflow_id))%>%
  workflow_map(
    seed = 2024,
    resamples = train_folds[[as.numeric(id)]],
    grid = mlr_grid,
    metrics = mset,
    control = control_grid(save_pred = TRUE, save_workflow = TRUE)
  ) %>%
  suppressMessages()

# Combined tuned sets
wflow_tuned_set <- bind_rows(rf_tuned_set,
                             xgb_tuned_set,
                             svm_tuned_set,
                             mlr_tuned_set)

tuned_set <- mlr_tuned_set

# Ranking workflows for selected metric
ranking <- tuned_set %>%
  rank_results(rank_metric = rank_metric) %>%
  filter(.metric == rank_metric) %>%
  select(model, wflow_id, f_meas = mean, rank)

# Best workflow
best_wflow <- head(ranking[["wflow_id"]], 1)

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
  last_fit(split = data_split$splits[[id]], metrics = mset)

# Calculate per-class metrics using one-vs-all predictions
mset_pc <- metric_set(accuracy, f_meas, kap, gmean)

metrics_per_class <- test_results %>%
  pluck(".predictions", 1) %>%
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
  mutate(data = data %>%
           map(droplevels) %>%
           map(mset_pc, truth = class_value, estimate = .pred_class_value)) %>%
  unnest(cols = data)

# Top class by selected metric
# metrics_per_class %>%
#   filter(.metric == "f_meas") %>%
#   slice_max(order_by = .estimate)
