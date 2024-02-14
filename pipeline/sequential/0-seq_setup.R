# sequence
seq_top <- bind_rows(seq_top_c5,
                     seq_top_c4,
                     seq_top_c3,
                     seq_top_c2)
seq_algs <- seq_top[["algorithm"]]
saveRDS(seq_algs, here("data/seq_algs.rds"))

# data
seq_data <- list(train_data,
                 retrain_4_data,
                 retrain_3_data,
                 retrain_2_data)
saveRDS(seq_data, here("data/seq_data.rds"))

# class
seq_class <- list(train_class,
                  retrain_4_class,
                  retrain_3_class,
                  retrain_2_class) %>%
  map2(seq_top[["class"]], ~ {
    ifelse(.x == .y, .x, "class_0")
  })
saveRDS(seq_class, here("data/seq_class.rds"))
