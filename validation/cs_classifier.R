# Process cohorts
source(here::here("validation/cs_process_cohorts.R"))
source(here::here("src/funs.R"))

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand1 <- hist %>%
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Reference Samples (n=5)
cs1_samples_R <- cohorts %>%
  inner_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

cs2_samples_R <- cohorts %>%
  inner_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

cs3_samples_R <- cohorts %>%
  inner_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# Expression Samples
# CS1: n=263
cs1_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS1 with duplicates: n=279
cs1_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS2: n=832
cs2_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS2 with duplicates: n=882
cs2_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS3: n=2107
cs3_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# Reference Datasets
# CS1: 5 samples by 256 genes
cs1_R <- cs1_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs1_samples_R)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# CS2: 5 samples by 365 genes
cs2_R <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs2_samples_R)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# CS3: 5 samples by 513 genes
cs3_R <- cs3_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs3_samples_R)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# Expression Datasets

# CS1 transposed
cs1_norm_t <- cs1_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  spread(Name, value) %>%
  column_to_rownames("FileName")

# CS1: 263 samples by 256 genes
cs1_X <- cs1_norm_t[cs1_samples_X, ]

# CS1 with duplicates: 279 samples by 256 genes
cs1_all_X <- cs1_norm_t[cs1_samples_all_X, ]

# CS2 transposed
cs2_norm_t <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select(-c(Code.Class, Accession)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  spread(Name, value) %>%
  column_to_rownames("FileName")

# CS2: 833 samples by 365 genes
cs2_X <- cs2_norm_t[cs2_samples_X, ]

# CS2 with duplicates: 884 samples by 365 genes
cs2_all_X <- cs2_norm_t[cs2_samples_all_X, ]

# CS3: 2111 samples by 513 genes
cs3_X <- cs3_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs3_samples_X)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# Normalized Datasets
# Normalizing CS1 to CS3 uses n=79 common genes
cs13_genes <- intersect(names(cs3_R), names(cs1_R))
cs1_train <- as.data.frame(refMethod(cs1_X[cs13_genes],
                                     cs3_R[cs13_genes],
                                     cs1_R[cs13_genes])) %>%
  select(all_of(common_genes123))

# Normalizing CS1 all to CS3 using n=79 common genes
cs1_all_train <- as.data.frame(refMethod(cs1_all_X[cs13_genes],
                                         cs3_R[cs13_genes],
                                         cs1_R[cs13_genes]))

# Normalizing CS2 to CS3 uses n=136 common genes
cs23_genes <- intersect(names(cs3_R), names(cs2_R))
cs2_train <- as.data.frame(refMethod(cs2_X[cs23_genes],
                                     cs3_R[cs23_genes],
                                     cs2_R[cs23_genes])) %>%
  select(all_of(common_genes123))

# Normalizing CS2 all to CS3 using n=136 common genes
cs2_all_train <- as.data.frame(refMethod(cs2_all_X[cs23_genes],
                                         cs3_R[cs23_genes],
                                         cs2_R[cs23_genes]))

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_site_rand1 <- hist %>%
  filter(
    FileName %in%
      c(cs3_clean_van$FileName,
        cs3_clean_aoc$FileName,
        cs3_clean_usc$FileName)
  ) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
# Reference Datasets: each are 5 samples by 513 genes
cs3_usc_rand1 <- join_avg(cs3_clean_usc, hist_site_rand1, "ottaID", "keep")
cs3_aoc_rand1 <- join_avg(cs3_clean_aoc, hist_site_rand1, "ottaID", "keep")
cs3_van_rand1 <- join_avg(cs3_clean_van, hist_site_rand1, "ottaID", "keep")

# Remove common samples from CS1, preserving gene order
# Remove common samples and keep last run sample for duplicates
# Expression Datasets: each are 19 samples by 513 genes
cs3_usc_dist_counts1 <- cs3_clean_usc %>%
  anti_join(hist_site_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")
cs3_aoc_dist_counts1 <- cs3_clean_aoc %>%
  anti_join(hist_site_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")
cs3_van_dist_counts1 <- cs3_clean_van %>%
  anti_join(hist_site_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")

# Normalize by reference method using common samples, add histotypes from annot
cs3_norm_usc_rand1 <-
  refMethod(cs3_usc_dist_counts1, cs3_van_rand1, cs3_usc_rand1) %>%
  as.data.frame()

cs3_norm_aoc_rand1 <-
  refMethod(cs3_aoc_dist_counts1, cs3_van_rand1, cs3_aoc_rand1) %>%
  as.data.frame()

# Combine the two CS3-VAN used to normalize CS1/CS2 and normalize CS3-USC/CS3-AOC
# and remove duplicates
cs3_train <-
  list(cs3_X, cs3_van_dist_counts1, cs3_norm_usc_rand1, cs3_norm_aoc_rand1) %>%
  map_dfr(rownames_to_column, "FileName") %>%
  distinct() %>%
  filter(!duplicated(FileName, fromLast = TRUE)) %>%
  mutate(col_name = paste0("X", FileName)) %>%
  inner_join(cohorts, by = "col_name") %>%
  filter(!cohort %in% c("TNCO", "DOVE4"),
         !grepl("pool", col_name, ignore.case = TRUE)) %>%
  mutate(col_name = gsub("^X", "", col_name)) %>%
  column_to_rownames("col_name") %>%
  select(all_of(common_genes123))

# Training set, n=266+833+514-77=1536 (CS1 + CS2 + CS3 excluding TNCO and DOVE
# - other histotypes), common genes n=72
train_ref <-
  bind_rows(cs1_train, cs2_train, cs3_train) %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

train_data <- select(train_ref, where(is.double))
train_class <- train_ref[["revHist"]]

saveRDS(train_data, here::here("data/train_data.rds"))
saveRDS(train_class, here::here("data/train_class.rds"))

# Training set for 2-step process
train_step1_data <- train_data
train_step1_class <- ifelse(train_class == "HGSC", "HGSC", "non-HGSC")
train_step2_data <- train_data[train_step1_class != "HGSC", ]
train_step2_class <- train_class[train_step1_class != "HGSC"]

saveRDS(train_step1_data, here::here("data/train_step1_data.rds"))
saveRDS(train_step1_class, here::here("data/train_step1_class.rds"))
saveRDS(train_step2_data, here::here("data/train_step2_data.rds"))
saveRDS(train_step2_class, here::here("data/train_step2_class.rds"))

# CS1 all set, n=282-19=263 (CS1 - other histotypes), common genes with CS3 n=79
cs1_all_ref <- cs1_all_train %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs1_all_data <- select(cs1_all_ref, where(is.double))
cs1_all_class <- cs1_all_ref[["revHist"]]

saveRDS(cs1_all_data, here::here("data/cs1_all_data.rds"))
saveRDS(cs1_all_class, here::here("data/cs1_all_class.rds"))

# CS2 all set, n=884-69=815 (CS2 - other histotypes), common genes with CS3 n=136
cs2_all_ref <- cs2_all_train %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs2_all_data <- select(cs2_all_ref, where(is.double))
cs2_all_class <- cs2_all_ref[["revHist"]]

saveRDS(cs2_all_data, here::here("data/cs2_all_data.rds"))
saveRDS(cs2_all_class, here::here("data/cs2_all_class.rds"))

# Confirmation set, n=673-30=643 (TNCO - other histotypes)
conf_ref <- cs3_X %>%
  rownames_to_column("col_name") %>%
  mutate(col_name = paste0("X", col_name)) %>%
  inner_join(cohorts, by = "col_name") %>%
  filter(cohort == "TNCO") %>%
  mutate(col_name = gsub("^X", "", col_name)) %>%
  select(FileName = col_name, all_of(common_genes123)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

conf_data <- select(conf_ref, where(is.double))
conf_class <- conf_ref[["revHist"]]

saveRDS(conf_data, here::here("data/conf_data.rds"))
saveRDS(conf_class, here::here("data/conf_class.rds"))

# Validation set, n=940-31=909 (DOVE - other histotypes)
val_ref <- cs3_X %>%
  rownames_to_column("col_name") %>%
  mutate(col_name = paste0("X", col_name)) %>%
  inner_join(cohorts, by = "col_name") %>%
  filter(cohort == "DOVE4") %>%
  mutate(col_name = gsub("^X", "", col_name)) %>%
  select(FileName = col_name, all_of(common_genes123)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

val_data <- select(val_ref, where(is.double))
val_class <- val_ref[["revHist"]]

saveRDS(val_data, here::here("data/val_data.rds"))
saveRDS(val_class, here::here("data/val_class.rds"))
