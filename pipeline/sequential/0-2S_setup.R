# sequence
step2_seq <- iv_summary_retrain_4 %>%
  filter(measure == "macro_f1") %>%
  transmute(
    algorithm,
    sampling,
    f1_median = percentile_50
  ) %>%
  slice_max(order_by = f1_median) %>%
  pull(algorithm)

two_step_algs <- c(seq_top_c5[["algorithm"]], step2_seq)
saveRDS(two_step_algs, here("data/two_step_algs.rds"))

# data
two_step_data <- list(train_step1_data,
                      train_step2_data)
saveRDS(two_step_data, here("data/two_step_data.rds"))

# class
two_step_class <- list(train_step1_class,
                       train_step2_class)
saveRDS(two_step_class, here("data/two_step_class.rds"))
