# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(here)
iv_summary_train <- readRDS(here("data/iv_summary_train.rds"))
train_data <- readRDS(here("data/train_data.rds"))
train_class <- readRDS(here("data/train_class.rds"))
train_df <- cbind(train_data, class = train_class)

# Top median F1-score per class and sampling
seq_summary <- iv_summary_train %>%
  filter(grepl("f1\\.", measure)) %>%
  transmute(
    class = gsub("f1\\.", "", measure),
    algorithm,
    sampling = factor(sampling, levels = c("none", "up", "down", "smote")),
    f1_median = percentile_50
  ) %>%
  slice_max(order_by = f1_median)

# Create retrain data and class lists
top_list <- seq_summary %>%
  select(class, algorithm, sampling) %>%
  as.list()

n_class <- n_distinct(train_class)
n_retrain <- n_class - 1

retrain_data <- train_df %>%
  filter(!class %in% top_list[["class"]]) %>%
  select(-class)

retrain_class <- train_df %>%
  filter(!class %in% top_list[["class"]]) %>%
  pull(class)

# # Save outputs
saveRDS(retrain_data, here("data", paste0("retrain_", retrain_n, "_data.rds")))
saveRDS(retrain_class, here("data", paste0("retrain_", retrain_n, "_class.rds")))
