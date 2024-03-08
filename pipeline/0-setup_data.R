# Import training data and class labels
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))
train_ref <- cbind(data, class = factor(class))

# Nested resampling
set.seed(2024)
folds <- nested_cv(
  train_ref,
  outside = vfold_cv(v = n_folds, strata = class),
  inside = vfold_cv(v = n_folds, repeats = 2, strata = class, pool = 0)
) %>%
  suppressWarnings()


# Metrics
gmean <- new_class_metric(gmean, direction = "maximize")
mset <- metric_set(accuracy, roc_auc, f_meas, kap, gmean)
