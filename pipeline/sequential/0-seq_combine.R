# Sequential algorithm setup
library(dplyr)
library(purrr)

# Combine workflows
seq_top <- map(5:2, ~ {
  readRDS(here("data", paste0("seq_top_c", .x, ".rds")))
}) %>%
  list_rbind()
seq_wflows <- seq_top[["wflow"]]
saveRDS(seq_wflows, here("data/seq_wflows.rds"))

# Combine training data
seq_data <- map(
  c(
    "train_data",
    "retrain_4_data",
    "retrain_3_data",
    "retrain_2_data"
  ),
  ~ readRDS(here("data", paste0(.x, ".rds")))
)
saveRDS(seq_data, here("data/seq_data.rds"))

# Combine training classes
seq_class <- map(
  c(
    "train_class",
    "retrain_4_class",
    "retrain_3_class",
    "retrain_2_class"
  ),
  ~ readRDS(here("data", paste0(.x, ".rds")))
) %>%
  map2(seq_top[["class_group"]], ~ {
    if (n_distinct(.x) > 2) {
      ifelse(.x == .y, .x, paste0("non_", .y))
    } else {
      .x
    }
  })
saveRDS(seq_class, here("data/seq_class.rds"))
