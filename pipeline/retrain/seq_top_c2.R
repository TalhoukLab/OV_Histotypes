# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(tibble)
library(here)
all_metrics_train <- readRDS(here("data/all_metrics_retrain_2.rds"))
train_data <- readRDS(here("data/retrain_2_data.rds"))
train_class <- readRDS(here("data/retrain_2_class.rds"))
train_ref <- cbind(train_data, class = train_class)

# Number of classes used
n_class <- n_distinct(train_class)

# Top workflow by per-class F1-score out of n_class
seq_top <- all_metrics_train %>%
  filter(.metric == rank_metric, class_group != "Overall") %>%
  slice_max(order_by = .estimate) %>%
  add_column(n_class = n_class)

# Save outputs
saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))
