# Combine confirmation data for sequential algorithm with gene optimization

# Confirmdation data and class
data <- readRDS(file.path(inputDir, "conf_data.rds"))
class <- readRDS(file.path(inputDir, "conf_class.rds"))
conf_ref <- cbind(data, class = class)

# Sequential of workflows and classes from training data
seq_top_c5 <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
seq_top_c4 <- readRDS(file.path(inputDir, "seq_top_c4.rds"))
seq_top_c3 <- readRDS(file.path(inputDir, "seq_top_c3.rds"))
seq_top_c2 <- readRDS(file.path(inputDir, "seq_top_c2.rds"))

# Create confirmation data sequence
conf_data_s1 <- conf_ref %>%
  select(-class)
conf_data_s2 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  select(-class)
conf_data_s3 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c4, by = c("class" = "class_group")) %>%
  select(-class)
conf_data_s4 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c4, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c3, by = c("class" = "class_group")) %>%
  select(-class)
conf_seq_data <-
  list(conf_data_s1, conf_data_s2, conf_data_s3, conf_data_s4)
conf_two_step_data <-
  list(conf_data_s1, conf_data_s2)

# Create confirmation class sequence
conf_class_s1 <- conf_ref %>%
  mutate(class = ifelse(class == seq_top_c5[["class_group"]], class, paste0("non_", seq_top_c5[["class_group"]]))) %>%
  pull(class)
conf_class_s2 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  mutate(class = ifelse(class == seq_top_c4[["class_group"]], class, paste0("non_", seq_top_c4[["class_group"]]))) %>%
  pull(class)
conf_class_s3 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c4, by = c("class" = "class_group")) %>%
  mutate(class = ifelse(class == seq_top_c3[["class_group"]], class, paste0("non_", seq_top_c3[["class_group"]]))) %>%
  pull(class)
conf_class_s4 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c4, by = c("class" = "class_group")) %>%
  anti_join(seq_top_c3, by = c("class" = "class_group")) %>%
  pull(class)
conf_class_step2 <- conf_ref %>%
  anti_join(seq_top_c5, by = c("class" = "class_group")) %>%
  pull(class)
conf_seq_class <-
  list(conf_class_s1, conf_class_s2, conf_class_s3, conf_class_s4)
conf_two_step_class <-
  list(conf_class_s1, conf_class_step2)

# Save outputs
saveRDS(conf_seq_data, here("data/conf_seq_data.rds"))
saveRDS(conf_seq_class, here("data/conf_seq_class.rds"))

saveRDS(conf_two_step_data, here("data/conf_two_step_data.rds"))
saveRDS(conf_two_step_class, here("data/conf_two_step_class.rds"))
