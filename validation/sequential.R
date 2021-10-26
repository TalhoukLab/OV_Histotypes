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
  split(.[["sampling"]]) %>%
  map(~ {
    group_by(.x, class) %>%
      filter(f1_median == max(f1_median)) %>%
      ungroup() %>%
      arrange(desc(f1_median))
  })

# Create sequential algorithm data and class lists
seq_ranks <- map(seq_summary, ~ {
  head(accumulate(.[["class"]], c), -1)
})

seq_data <- map_depth(seq_ranks, 2, ~ {
  train_df %>%
    filter(!class %in% head(.x, -1)) %>%
    select(-class)
})

seq_class <- map(seq_ranks, ~ {
  imap(.x, ~ {
    cl <- train_df %>%
      filter(!class %in% head(.x, -1)) %>%
      pull(class)
    if (.y != length(seq_ranks)) {
      ifelse(cl == tail(.x, 1), tail(.x, 1), "class_0")
    } else {
      cl
    }
  })
})

seq_algs <- map(seq_summary, ~ head(.[["algorithm"]], -1))

# Save outputs
saveRDS(seq_data, here("data/seq_data.rds"))
saveRDS(seq_class, here("data/seq_class.rds"))
saveRDS(seq_algs, here("data/seq_algs.rds"))

# Rank Aggregation summary table
