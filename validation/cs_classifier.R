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
# CS1: n=270
cs1_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS2: n=840
cs2_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name) %>%
  gsub("^X", "", .)

# CS3: n=2267
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
# CS1: 270 samples by 256 genes
cs1_X <- cs1_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs1_samples_X)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# CS2: 840 samples by 365 genes
cs2_X <- cs2_norm %>%
  rename_all(~ gsub("^X", "", .)) %>%
  select_if(names(.) %in% c("Name", cs2_samples_X)) %>%
  mutate(Name = fct_inorder(Name)) %>%
  gather(FileName, value, -Name) %>%
  inner_join(hist, by = "FileName") %>%
  spread(Name, value) %>%
  select(-ottaID:-site) %>%
  column_to_rownames("FileName")

# CS3: 2267 samples by 513 genes
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
cs1_norm <- as.data.frame(refMethod(cs1_X, cs3_R, cs1_R))
# Normalizing CS2 to CS3 uses n=136 common genes
cs2_norm <- as.data.frame(refMethod(cs2_X, cs3_R, cs2_R))

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

# # Combine the two CS3-VAN used to normalize CS1/CS2 and normalize CS3-USC/CS3-AOC
# # and remove duplicates
# cs3_combined <- bind_rows(
#   rownames_to_column(cs3_dist_counts1, "FileName"),
#   rownames_to_column(cs3_van_dist_counts1, "FileName")
# ) %>%
#   distinct() %>%
#   column_to_rownames("FileName")
#
# # Combine all 5 datasets for classifier development
# cs_classifier <- bind_rows(
#   cs1_norm_rand1,
#   cs2_norm_rand1,
#   cs3_norm_usc_rand1,
#   cs3_norm_aoc_rand1,
#   cs3_combined
# )
