# Process cohorts
source(here::here("validation/cs_process_cohorts.R"))
source(here::here("src/funs.R"))


# Remove Duplicates -------------------------------------------------------

# CS1: n=246
cs1_dedup <- cohorts %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE))

# CS2: n=795
cs2_dedup <- cohorts %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE))

# CS3: n=2011 (Vancouver site only and removed pools)
cs3_dedup <- cohorts %>%
  semi_join(hist_cs3_van, by = "col_name") %>%
  filter(col_name %in% cs3_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE))

# Reference Samples -------------------------------------------------------

# CS1: n=5
cs1_samples_R <- cs1_dedup %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  arrange(ottaID) %>%
  pull(col_name)

# CS2: n=5
cs2_samples_R <- cs2_dedup %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  arrange(ottaID) %>%
  pull(col_name)

# CS3: n=5
cs3_samples_R <- cs3_dedup %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  arrange(ottaID) %>%
  pull(col_name)


# Expression Samples ------------------------------------------------------

# CS1: n=241
cs1_samples_X <- cs1_dedup %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  pull(col_name)

# CS1 with duplicates: n=257
cs1_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  pull(col_name)

# CS2: n=790
cs2_samples_X <- cs2_dedup %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  pull(col_name)

# CS2 with duplicates: n=820
cs2_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  pull(col_name)

# CS3: n=2006
cs3_samples_X <- cs3_dedup %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  pull(col_name)

# CS3 with duplicates: n=2139
cs3_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  pull(col_name)


# Reference Datasets ------------------------------------------------------

# CS1: 5 samples by 256 genes
cs1_R <- select_samples(cs1_norm, cs1_samples_R)

# CS2: 5 samples by 365 genes
cs2_R <- select_samples(cs2_norm, cs2_samples_R)

# CS3: 5 samples by 513 genes
cs3_R <- select_samples(cs3_norm, cs3_samples_R)


# Expression Datasets -----------------------------------------------------

# CS1: 241 samples by 256 genes
cs1_X <- select_samples(cs1_norm, cs1_samples_X)

# CS1 with duplicates: 257 samples by 256 genes
cs1_all_X <- select_samples(cs1_norm, cs1_samples_all_X)

# CS2: 790 samples by 365 genes
cs2_X <- select_samples(cs2_norm, cs2_samples_X)

# CS2 with duplicates: 820 samples by 365 genes
cs2_all_X <- select_samples(cs2_norm, cs2_samples_all_X)

# CS3: 2006 samples by 513 genes
cs3_X <- select_samples(cs3_norm, cs3_samples_X)


# Normalization by Reference Method ---------------------------------------

# Normalizing CS1 to CS3 uses n=79 common genes
cs13_genes <- intersect(names(cs3_R), names(cs1_R))
cs1_train <- refMethod(Y = cs1_X[cs13_genes],
                       R1 = cs3_R[cs13_genes],
                       R2 = cs1_R[cs13_genes]) %>%
  as.data.frame() %>%
  select(all_of(common_genes123))

# Normalizing CS1 all to CS3 using n=79 common genes
cs1_all_train <- refMethod(Y = cs1_all_X[cs13_genes],
                           R1 = cs3_R[cs13_genes],
                           R2 = cs1_R[cs13_genes]) %>%
  as.data.frame()

# Normalizing CS2 to CS3 uses n=136 common genes
cs23_genes <- intersect(names(cs3_R), names(cs2_R))
cs2_train <- refMethod(Y = cs2_X[cs23_genes],
                       R1 = cs3_R[cs23_genes],
                       R2 = cs2_R[cs23_genes]) %>%
  as.data.frame() %>%
  select(all_of(common_genes123))

# Normalizing CS2 all to CS3 using n=136 common genes
cs2_all_train <- refMethod(Y = cs2_all_X[cs23_genes],
                           R1 = cs3_R[cs23_genes],
                           R2 = cs2_R[cs23_genes]) %>%
  as.data.frame()


# Normalization by Pools --------------------------------------------------

# CS3-USC normalized to CS3-VAN by pools
cs3_usc_norm <- normalize_pools(
  x = cs3_usc_X,
  x_pools = cs3_usc_R,
  ref_pools = ref_pools,
  p = 3,
  weigh = TRUE,
  same_codeset = TRUE
) |>
  mutate(FileName = gsub("^X", "", FileName)) |>
  column_to_rownames("FileName")

# CS3-AOC normalized to CS3-VAN by pools
cs3_aoc_norm <- normalize_pools(
  x = cs3_aoc_X,
  x_pools = cs3_aoc_R,
  ref_pools = ref_pools,
  p = 3,
  weigh = TRUE,
  same_codeset = TRUE
) |>
  mutate(FileName = gsub("^X", "", FileName)) |>
  column_to_rownames("FileName")

# Remove cross-site replicates, test sets, and pool samples from CS3
cs3_train <-
  list(VAN = cs3_X, USC = cs3_usc_norm, AOC = cs3_aoc_norm) |>
  map(~ {
    .x |>
      rownames_to_column("FileName") %>%
      mutate(col_name = paste0("X", FileName)) %>%
      inner_join(cohorts, by = "col_name") %>%
      filter(!cohort %in% c("TNCO", "DOVE4")) %>%
      column_to_rownames("FileName") %>%
      select(ottaID, all_of(common_genes123))
  }) |>
  list_rbind(names_to = "Site") |>
  filter(Site == "VAN", .by = ottaID) |>
  select(-c(Site, ottaID))


# Construct Training Set --------------------------------------------------

# Combined Training Set (CS1 + CS2 + CS3), n=241+790+470=1501
# Common genes n=72
train_ref_comb <-
  bind_rows(cs1_train, cs2_train, cs3_train) %>%
  rownames_to_column("FileName") %>%
  mutate(col_name = paste0("X", FileName)) %>%
  inner_join(cohorts, by = "col_name") %>%
  select(FileName, all_of(common_genes123), cohort) %>%
  inner_join(hist, by = "FileName") %>%
  column_to_rownames("FileName")

# Training set removed replicates (CS3 > CS2 > CS1), n=1501-244=1257
train_ref <- train_ref_comb %>%
  slice_tail(n = 1, by = ottaID)

train_data <- select(train_ref, where(is.double))
train_class <- train_ref[["hist_final"]]

saveRDS(train_data, here::here("data/train_data.rds"))
saveRDS(train_class, here::here("data/train_class.rds"))


# Two-Step Training Set ---------------------------------------------------

# Training sets for two-step classifier
train_step1_data <- train_data
train_step1_class <- train_ref |> pull(hist_gr)
train_step2_data <- train_ref |>
  filter(hist_gr != "HGSC") |>
  select(where(is.double))
train_step2_class <- train_ref |>
  filter(hist_gr != "HGSC") |>
  pull(hist_final)

two_step_data <- list(train_step1_data, train_step2_data)
saveRDS(two_step_data, here::here("data/two_step_data.rds"))

two_step_class <- list(train_step1_class, train_step2_class)
saveRDS(two_step_class, here::here("data/two_step_class.rds"))


# CS1 and CS2 with duplicates ---------------------------------------------

# CS1 all set, n=257, common genes with CS3 n=79
cs1_all_ref <- cs1_all_train %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  column_to_rownames("FileName")

cs1_all_data <- select(cs1_all_ref, where(is.double))
cs1_all_class <- cs1_all_ref[["hist_final"]]

saveRDS(cs1_all_data, here::here("data/cs1_all_data.rds"))
saveRDS(cs1_all_class, here::here("data/cs1_all_class.rds"))

# CS2 all set, n=820 common genes with CS3 n=136
cs2_all_ref <- cs2_all_train %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  column_to_rownames("FileName")

cs2_all_data <- select(cs2_all_ref, where(is.double))
cs2_all_class <- cs2_all_ref[["hist_final"]]

saveRDS(cs2_all_data, here::here("data/cs2_all_data.rds"))
saveRDS(cs2_all_class, here::here("data/cs2_all_class.rds"))


# Construct Test Sets -----------------------------------------------------

# Confirmation set, n=642 (TNCO)
conf_ref <- cs3_X %>%
  rownames_to_column("FileName") %>%
  mutate(col_name = paste0("X", FileName)) %>%
  inner_join(cohorts, by = "col_name") %>%
  select(FileName, all_of(common_genes123), cohort) %>%
  inner_join(hist, by = "FileName") %>%
  filter(cohort == "TNCO") %>%
  column_to_rownames("FileName")

conf_data <- select(conf_ref, where(is.double))
conf_class <- conf_ref[["hist_final"]]

saveRDS(conf_data, here::here("data/conf_data.rds"))
saveRDS(conf_class, here::here("data/conf_class.rds"))

# Validation set, n=894 (DOVE)
val_ref <- cs3_X %>%
  rownames_to_column("FileName") %>%
  mutate(col_name = paste0("X", FileName)) %>%
  inner_join(cohorts, by = "col_name") %>%
  select(FileName, all_of(common_genes123), cohort) %>%
  inner_join(hist, by = "FileName") %>%
  filter(cohort == "DOVE4") %>%
  column_to_rownames("FileName")

val_data <- select(val_ref, where(is.double))
val_class <- val_ref[["hist_final"]]

saveRDS(val_data, here::here("data/val_data.rds"))
saveRDS(val_class, here::here("data/val_class.rds"))
