# Setup -------------------------------------------------------------------

source(here::here("validation/cs_process_cohorts.R"))

# Pools Method ------------------------------------------------------------

# Pool set: pool samples from CS2
pool_samples <-
  gsub("-|\\.RCC", "", grep("POOL", pools[["CS2-FileName"]], value = TRUE))
cs2_pools <- dplyr::select(cs2_norm, Name, paste0("X", pool_samples))

# Locked down reference pools weights
weights <- c("Pool1", "Pool2", "Pool3") %>%
  purrr::set_names() %>%
  purrr::map_dbl(~ mean(grepl(., names(ref_pools), ignore.case = TRUE)))

# Weighted mean gene expression for CS2 pools (using reference pools weights)
cs2_pools_mgx <- weights %>%
  purrr::imap_dfc(~ .x * rowMeans(dplyr::select(cs2_pools, dplyr::matches(.y)))) %>%
  dplyr::transmute(Name = cs2_pools[["Name"]], cs2_exp = rowSums(.))

# Mean gene expression for CS3 reference pools
cs3_pools_mgx <-
  tibble::enframe(rowMeans(ref_pools), name = "Name", value = "cs3_exp")

# Extract cs2 norm counts (not pools)
cs2_norm_counts <- cs2_norm %>% dplyr::select(-dplyr::one_of(names(cs2_pools)[-1]))

# Normalize each gene by adding batch effect (diff in mean gx)
cs2_normalized_data_pools <-
  dplyr::inner_join(cs3_pools_mgx, cs2_pools_mgx, by = "Name") %>%
  dplyr::transmute(Name, be = cs3_exp - cs2_exp) %>%
  dplyr::inner_join(cs2_norm_counts, by = "Name") %>%
  tidyr::gather(FileName, exp, -1:-4) %>%
  dplyr::transmute(Name = forcats::fct_inorder(Name), FileName, exp = be + exp) %>%
  tidyr::spread(Name, exp)

# Extract cs3 norm counts (not pools)
cs3_norm_counts <- cs3_norm %>% dplyr::select(-dplyr::one_of(names(ref_pools)))

# Normalize each gene by adding batch effect (diff in mean gx)
cs3_normalized_data_pools <-
  dplyr::inner_join(cs2_pools_mgx, cs3_pools_mgx, by = "Name") %>%
  dplyr::transmute(Name, be = cs2_exp - cs3_exp) %>%
  dplyr::inner_join(cs3_norm_counts, by = "Name") %>%
  tidyr::gather(FileName, exp, -1:-4) %>%
  dplyr::transmute(Name = forcats::fct_inorder(Name), FileName, exp = be + exp) %>%
  tidyr::spread(Name, exp)

## Summary
# CS2 pools: 9 samples, 365 genes # dim(t(cs2_pools))
# CS3 pools: 22 samples, 513 genes # dim(t(ref_pools))
# CS2 normalized: 1214 samples (1223 original - 9 from pools), 136 common genes # dim(cs2_normalized_data_pools)
# CS3 normalized: 1295 samples (1317 original - 22 from pools), 136 common genes # dim(cs2_normalized_data_pools)
