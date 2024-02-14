# Read in multiclass data, class, and IV summary
library(dplyr)
library(purrr)
library(tibble)
library(here)
iv_summary_train <- readRDS(here("data/iv_summary_retrain_2.rds"))
train_data <- readRDS(here("data/retrain_2_data.rds"))
train_class <- readRDS(here("data/retrain_2_class.rds"))
train_df <- cbind(train_data, class = train_class)

# Number of classes used
n_class <- n_distinct(train_class)

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

# Save outputs
saveRDS(seq_top, here("data", paste0("seq_top_c", n_class, ".rds")))
