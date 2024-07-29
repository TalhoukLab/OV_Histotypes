# Process cohorts
source(here::here("validation/cs_process_cohorts.R"))
source(here::here("src/funs.R"))

# CS3 site mapping used for prioritizing Vancouver site and most recent date
hist_cs3 <- hist %>%
  mutate(
    col_name = paste0("X", FileName),
    site = factor(site, levels = c("Vancouver", "USC", "AOC")),
    .keep = "none"
  )

# Reference Samples
# CS1: n=5
cs1_samples_R <- cohorts %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  arrange(ottaID) %>%
  pull(col_name)

# CS2: n=5
cs2_samples_R <- cohorts %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  arrange(ottaID) %>%
  pull(col_name)

# CS3: n=5 (ensure Vancouver site)
cs3_samples_R <- cohorts %>%
  semi_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  inner_join(hist_cs3, by = "col_name") %>%
  group_by(ottaID) %>%
  arrange(site, desc(col_name)) %>%
  ungroup() %>%
  distinct(ottaID, .keep_all = TRUE) %>%
  arrange(ottaID) %>%
  pull(col_name)

# Expression Samples
# CS1: n=263
cs1_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name)

# CS1 with duplicates: n=279
cs1_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs1_samples) %>%
  pull(col_name)

# CS2: n=827
cs2_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name)

# CS2 with duplicates: n=876
cs2_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs2_samples) %>%
  pull(col_name)

# CS3: n=2094 (ensure Vancouver site)
cs3_samples_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  inner_join(hist_cs3, by = "col_name") %>%
  filter(site == "Vancouver") %>%
  filter(!duplicated(ottaID, fromLast = TRUE)) %>%
  pull(col_name)

# CS3 with duplicates: n=2264
cs3_samples_all_X <- cohorts %>%
  anti_join(hist_rand1, by = "ottaID") %>%
  filter(col_name %in% cs3_samples) %>%
  pull(col_name)

# Reference Datasets
# CS1: 5 samples by 256 genes
cs1_R <- cs1_norm %>%
  pivot_longer(
    cols = all_of(cs1_samples_R),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS2: 5 samples by 365 genes
cs2_R <- cs2_norm %>%
  pivot_longer(
    cols = all_of(cs2_samples_R),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS3: 5 samples by 513 genes
cs3_R <- cs3_norm %>%
  pivot_longer(
    cols = all_of(cs3_samples_R),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# Expression Datasets
# CS1: 263 samples by 256 genes
cs1_X <- cs1_norm %>%
  pivot_longer(
    cols = all_of(cs1_samples_X),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS1 with duplicates: 279 samples by 256 genes
cs1_all_X <- cs1_norm %>%
  pivot_longer(
    cols = all_of(cs1_samples_all_X),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS2: 827 samples by 365 genes
cs2_X <- cs2_norm %>%
  pivot_longer(
    cols = all_of(cs2_samples_X),
    names_to = "FileName",
    names_prefix = "X",
    values_to = "value"
  ) %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS2 with duplicates: 876 samples by 365 genes
cs2_all_X <- cs2_norm %>%
  pivot_longer(cols = all_of(cs2_samples_all_X),
               names_to = "FileName",
               names_prefix = "X",
               values_to = "value") %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# CS3: 2094 samples by 513 genes
cs3_X <- cs3_norm %>%
  pivot_longer(cols = all_of(cs3_samples_X),
               names_to = "FileName",
               names_prefix = "X",
               values_to = "value") %>%
  pivot_wider(id_cols = FileName,
              names_from = Name,
              values_from = value) %>%
  column_to_rownames("FileName")

# Normalized by Reference Method
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

# CS3-VAN reference pools setup
weights <- c("Pool1", "Pool2", "Pool3") %>%
  purrr::set_names() %>%
  purrr::map_dbl(~ ncol(dplyr::select(ref_pools, dplyr::matches(.))) /
                   ncol(ref_pools))

ref_mean_gx <-
  rowMeans(ref_pools) %>%
  tibble::enframe(name = "Name", value = "ref_exp")

# CS3-USC normalized to CS3-VAN by pools
cs3_usc_R_mean_gx <-
  weights %>%
  purrr::imap( ~ {
    df <- dplyr::select(cs3_usc_R, Name, dplyr::matches(.y)) %>%
      tibble::column_to_rownames("Name")
    tibble::enframe(.x * rowSums(df) / ncol(df), name = "Name", value = .y)
  })  %>%
  purrr::reduce(dplyr::inner_join, by = "Name") %>%
  dplyr::transmute(Name, norm_exp = rowSums(dplyr::select(., dplyr::contains("Pool"))))

merged_usc <- dplyr::inner_join(ref_mean_gx, cs3_usc_R_mean_gx, by = "Name") %>%
  dplyr::transmute(Name, be = ref_exp - norm_exp) %>%
  dplyr::inner_join(cs3_usc_X, by = "Name")

cs3_usc_norm <- merged_usc %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Name") %>%
  dplyr::select(-c(be, Code.Class, Accession)) %>%
  dplyr::rename_with(~ gsub("^X", "", .)) %>%
  apply(2, `+`, merged_usc[["be"]]) %>%
  t() %>%
  as.data.frame()

# CS3-USC normalized to CS3-VAN by pools
cs3_aoc_R_mean_gx <-
  weights %>%
  purrr::imap( ~ {
    df <- dplyr::select(cs3_aoc_R, Name, dplyr::matches(.y)) %>%
      tibble::column_to_rownames("Name")
    tibble::enframe(.x * rowSums(df) / ncol(df), name = "Name", value = .y)
  })  %>%
  purrr::reduce(dplyr::inner_join, by = "Name") %>%
  dplyr::transmute(Name, norm_exp = rowSums(dplyr::select(., dplyr::contains("Pool"))))

merged_aoc <- dplyr::inner_join(ref_mean_gx, cs3_aoc_R_mean_gx, by = "Name") %>%
  dplyr::transmute(Name, be = ref_exp - norm_exp) %>%
  dplyr::inner_join(cs3_aoc_X, by = "Name")

cs3_aoc_norm <- merged_aoc %>%
  as.data.frame() %>%
  tibble::column_to_rownames("Name") %>%
  dplyr::select(-c(be, Code.Class, Accession)) %>%
  dplyr::rename_with(~ gsub("^X", "", .)) %>%
  apply(2, `+`, merged_aoc[["be"]]) %>%
  t() %>%
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

cs3_train2 <-
  list(cs3_X, cs3_usc_norm, cs3_aoc_norm) %>%
  map_dfr(rownames_to_column, "FileName") %>%
  mutate(col_name = paste0("X", FileName)) %>%
  inner_join(cohorts, by = "col_name") %>%
  filter(!cohort %in% c("TNCO", "DOVE4", "POOL-1", "POOL-2", "POOL-3")) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes123))

# Training set, n=263+827+514-75-286=1243 (CS1 + CS2 + CS3 excluding TNCO and DOVE
# - other histotypes - duplicates), common genes n=72
train_ref_all <-
  bind_rows(cs1_train, cs2_train, cs3_train) %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  inner_join(transmute(cohorts, FileName = gsub("^X", "", col_name), cohort),
             by = "FileName") %>%
  column_to_rownames("FileName")

train_ref <- train_ref_all %>%
  mutate(
    CodeSet_Site = fct_cross(CodeSet, site, sep = "_") %>%
      fct_relevel(
        "CS3_Vancouver",
        "CS3_AOC",
        "CS3_USC",
        "CS2_Vancouver",
        "CS2_AOC",
        "CS1_Vancouver"
      )
  ) %>%
  arrange(CodeSet_Site) %>%
  filter(!duplicated(ottaID)) %>%
  select(-CodeSet_Site) %>%
  arrange(CodeSet)

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

# CS1 all set, n=279-19=260 (CS1 - other histotypes), common genes with CS3 n=79
cs1_all_ref <- cs1_all_train %>%
  rownames_to_column("FileName") %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs1_all_data <- select(cs1_all_ref, where(is.double))
cs1_all_class <- cs1_all_ref[["revHist"]]

saveRDS(cs1_all_data, here::here("data/cs1_all_data.rds"))
saveRDS(cs1_all_class, here::here("data/cs1_all_class.rds"))

# CS2 all set, n=876-67=807 (CS2 - other histotypes), common genes with CS3 n=136
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
  select(FileName = col_name, cohort, all_of(common_genes123)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

conf_data <- select(conf_ref, where(is.double))
conf_class <- conf_ref[["revHist"]]

saveRDS(conf_data, here::here("data/conf_data.rds"))
saveRDS(conf_class, here::here("data/conf_class.rds"))

# Validation set, n=924-29=895 (DOVE - other histotypes)
val_ref <- cs3_X %>%
  rownames_to_column("col_name") %>%
  mutate(col_name = paste0("X", col_name)) %>%
  inner_join(cohorts, by = "col_name") %>%
  filter(cohort == "DOVE4") %>%
  mutate(col_name = gsub("^X", "", col_name)) %>%
  select(FileName = col_name, cohort, all_of(common_genes123)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

val_data <- select(val_ref, where(is.double))
val_class <- val_ref[["revHist"]]

saveRDS(val_data, here::here("data/val_data.rds"))
saveRDS(val_class, here::here("data/val_class.rds"))
