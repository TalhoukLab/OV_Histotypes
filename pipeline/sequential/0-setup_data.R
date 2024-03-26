# Import sequential training data and class labels
sq_data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
sq_class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))
data <- sq_data[[nseq]]
class <- sq_class[[nseq]]
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
mset <- metric_set(accuracy, f_meas, kap, gmean, roc_auc)
per_class_mset <- metric_set(accuracy, f_meas, kap, gmean)
