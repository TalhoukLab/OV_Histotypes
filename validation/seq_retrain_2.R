# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(tibble)
library(here)
all_metrics_train <- readRDS(here("data/all_metrics_retrain_3.rds"))
train_data <- readRDS(here("data/retrain_3_data.rds"))
train_class <- readRDS(here("data/retrain_3_class.rds"))
train_ref <- cbind(train_data, class = train_class)

# Number of classes used and number of classes to retrain after removing top class
n_class <- n_distinct(train_class)
n_retrain <- n_class - 1

# Top workflow by per-class F1-score out of n_class
seq_top <- all_metrics_train %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = .estimate) %>%
  add_column(n_class = n_class)

# Create retrain data and class lists
retrain_data <- train_ref %>%
  filter(!class %in% seq_top[["class_group"]]) %>%
  select(-class)

retrain_class <- train_ref %>%
  filter(!class %in% seq_top[["class_group"]]) %>%
  pull(class)

# Save outputs
saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))
saveRDS(retrain_data, here("data", paste0("retrain_", n_retrain, "_data.rds")))
saveRDS(retrain_class, here("data", paste0("retrain_", n_retrain, "_class.rds")))
