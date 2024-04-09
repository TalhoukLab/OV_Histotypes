# Find the top class of each sequence and generate a retraining subset for the
# next sequence. Set `dataset` to one of the following:
# c("train", "retrain_4", "retrain_3", "retrain_2")

# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(tibble)
library(here)
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))
train_ref <- cbind(data, class = class)

# Number of classes used and number of classes to retrain after removing top class
n_class <- n_distinct(class)

# Number of classes used and number of classes to retrain after removing top class
subDir <- ifelse(dataset == "train", "", "retrain")
metrics <- readRDS(file.path(
  outputDir,
  subDir,
  "summarize_results",
  dataset,
  paste0("all_metrics_", dataset, ".rds")
))
n_class <- n_distinct(class)

# Top workflow by per-class F1-score out of n_class
seq_top <- metrics %>%
  distinct(pick(-fold_id, -.estimate)) %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = mean_estimate) %>%
  add_column(n_class = n_class)

saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))

# Create retrain data and class if more than 2 classes
if (n_class > 2) {
  retrain_ref <-
    train_ref %>% filter(!class %in% seq_top[["class_group"]])
  retrain_data <- retrain_ref %>% select(-class)
  retrain_class <- retrain_ref %>% pull(class)

  saveRDS(retrain_data, here("data", paste0("retrain_", n_class - 1, "_data.rds")))
  saveRDS(retrain_class, here("data", paste0("retrain_", n_class - 1, "_class.rds")))
}
