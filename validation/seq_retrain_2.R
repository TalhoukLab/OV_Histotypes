# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(tibble)
library(here)
iv_summary_train <- readRDS(here("data/iv_summary_retrain_3.rds"))
train_data <- readRDS(here("data/retrain_3_data.rds"))
train_class <- readRDS(here("data/retrain_3_class.rds"))
train_df <- cbind(train_data, class = train_class)

# Number of classes used and number of classes to retrain after removing top class
n_class <- n_distinct(train_class)
n_retrain <- n_class - 1

# Top median F1-score
seq_summary <- iv_summary_train %>%
  filter(grepl("f1\\.", measure)) %>%
  transmute(
    class = gsub("f1\\.", "", measure),
    algorithm,
    sampling = factor(sampling, levels = c("none", "up", "down", "smote")),
    f1_median = percentile_50
  ) %>%
  slice_max(order_by = f1_median)

# Algorithm and subsampling method used for top class out of n=n_class
seq_top <- add_column(seq_summary, n_class = n_class)

# Create retrain data and class lists
retrain_data <- train_df %>%
  filter(!class %in% seq_top[["class"]]) %>%
  select(-class)

retrain_class <- train_df %>%
  filter(!class %in% seq_top[["class"]]) %>%
  pull(class)

# Save outputs
saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))
saveRDS(retrain_data, here("data", paste0("retrain_", n_retrain, "_data.rds")))
saveRDS(retrain_class, here("data", paste0("retrain_", n_retrain, "_class.rds")))
