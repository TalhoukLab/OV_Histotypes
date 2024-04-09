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
wflow <- paste(samp, alg, sep = "_")

# Hyperparameter tuning

## Algorithm-specific tuning setup
if (grepl("rf|xgb", wflow)) {
  wflow_set <- wflow_sets %>%
    filter(wflow_id == wflow)

  tuning_grid <- 10

} else if (grepl("svm", wflow)) {
  set.seed(2024)
  sigma_range <- as.vector(kernlab::sigest(as.matrix(data))[-2])

  svm_params <-
    parameters(cost(), rbf_sigma()) %>%
    update(
      cost = cost(c(0, 3)),
      rbf_sigma = rbf_sigma(sigma_range, trans = NULL)
    )

  wflow_set <- wflow_sets %>%
    filter(wflow_id == wflow) %>%
    option_add(param_info = svm_params)

  tuning_grid <- 10

} else if (grepl("mr", wflow)) {
  wflow_set <- wflow_sets %>%
    filter(wflow_id == wflow)

  set.seed(2024)
  tuning_grid <- crossing(
    grid_latin_hypercube(penalty(), size = 10),
    grid_regular(mixture(), levels = 10)
  )
}
