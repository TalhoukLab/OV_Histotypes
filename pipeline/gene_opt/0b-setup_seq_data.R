# Import sequential training data and class labels
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))[[nseq]]
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))[[nseq]]

# Genes from PrOTYPE and SPOT to keep
base_genes <- c("COL11A1", "CD74", "CD2", "TIMP3", "LUM", "CYTIP", "COL3A1",
                "THBS2", "TCF7L1", "HMGA2", "FN1", "POSTN", "COL1A2", "COL5A2",
                "PDZK1IP1", "FBN1", "HIF1A", "CXCL10", "DUSP4", "SOX17", "MITF",
                "CDKN3", "BRCA2", "CEACAM5", "ANXA4", "SERPINE1", "CRABP2", "DNAJC9")

# Import ranked variable importance for trained workflows
ranked_vi <- readRDS(file.path(
  outputDir,
  "sequential",
  "summarize_results",
  seq_data,
  paste0("all_vi_", seq_data, ".rds")
))

# Order of candidate genes for sequence
candidate_genes <- ranked_vi %>%
  slice_head(n = as.numeric(ngene)) %>%
  pull(Variable)

# Select genes
data <- data %>%
  select(all_of(c(base_genes, candidate_genes)))
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
mset <-
  metric_set(accuracy, sensitivity, specificity, f_meas, bal_accuracy, kap, roc_auc)
per_class_mset <-
  metric_set(accuracy, sensitivity, specificity, f_meas, bal_accuracy, kap)
