# 2-step algorithm setup
library(dplyr)
library(purrr)

# Combine workflows
seq_top_c5 <- readRDS(here("data/seq_top_c5.rds"))
step1_wflow <- seq_top_c5[["wflow"]]

all_metrics_retrain_4 <- readRDS(here("data/all_metrics_retrain_4.rds"))
rank_metric <- "f_meas"
step2_wflow <- all_metrics_retrain_4 %>%
  distinct(pick(-fold_id, -.estimate)) %>%
  filter(.metric == rank_metric, class_group == "Overall") %>%
  slice_max(order_by = mean_estimate) %>%
  pull(wflow)

two_step_wflows <- c(step1_wflow, step2_wflow)
saveRDS(two_step_wflows, here("data/two_step_wflows.rds"))
