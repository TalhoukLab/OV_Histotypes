# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand1 <- hist %>%
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand1 <- join_avg(cs1_clean, hist_rand1, "ottaID", "keep")
cs2_rand1 <- join_avg(cs2_clean, hist_rand1, "ottaID", "keep")
cs3_rand1 <- join_avg(cs3_clean, hist_rand1, "ottaID", "keep")

# Remove common samples and keep last run sample for duplicates
cs1_dist_counts1 <- cs1_clean %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")
cs2_dist_counts1 <- cs2_clean %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")
cs3_dist_counts1 <- cs3_clean %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  select(-ottaID) %>%
  column_to_rownames("FileName")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand1 <-
  as.data.frame(refMethod(cs1_dist_counts1, cs3_rand1, cs1_rand1))
cs2_norm_rand1 <-
  as.data.frame(refMethod(cs2_dist_counts1, cs3_rand1, cs2_rand1))

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
cs3_usc_rand1 <- join_avg(cs3_clean_usc, hist_site_rand1, "ottaID", "keep")
cs3_aoc_rand1 <- join_avg(cs3_clean_aoc, hist_site_rand1, "ottaID", "keep")
cs3_van_rand1 <- join_avg(cs3_clean_van, hist_site_rand1, "ottaID", "keep")

# Remove common samples from CS1, preserving gene order
# Remove common samples and keep last run sample for duplicates
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
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Normalize by reference method using common samples, add histotypes from annot
cs3_norm_usc_rand1 <-
  refMethod(cs3_usc_dist_counts1, cs3_van_rand1, cs3_usc_rand1) %>%
  as.data.frame() %>%
  select(all_of(common_genes))

cs3_norm_aoc_rand1 <-
  refMethod(cs3_aoc_dist_counts1, cs3_van_rand1, cs3_aoc_rand1) %>%
  as.data.frame() %>%
  select(all_of(common_genes))

# Combine the two CS3-VAN used to normalize CS1/CS2 and normalize CS3-USC/CS3-AOC
# and remove duplicates
cs3_combined <- bind_rows(
  rownames_to_column(cs3_dist_counts1, "FileName"),
  rownames_to_column(cs3_van_dist_counts1, "FileName")
) %>%
  distinct() %>%
  column_to_rownames("FileName")

# Combine all 5 datasets for classifier development
cs_classifier <- bind_rows(
  cs1_norm_rand1,
  cs2_norm_rand1,
  cs3_norm_usc_rand1,
  cs3_norm_aoc_rand1,
  cs3_combined
)
