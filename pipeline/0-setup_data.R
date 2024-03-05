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
