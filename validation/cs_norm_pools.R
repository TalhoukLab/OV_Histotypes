# Setup -------------------------------------------------------------------

source(here::here("validation/cs_process_cohorts.R"))

# Pools Method ------------------------------------------------------------

# Pool set: pool samples from CS2
pool_samples <-
  gsub("-|\\.RCC", "", grep("POOL", pools[["CS2-FileName"]], value = TRUE))
cs2_pools <- select(cs2_norm, Name, paste0("X", pool_samples))

# Locked down reference pools weights
weights <- c("Pool1", "Pool2", "Pool3") %>%
  set_names() %>%
  map_dbl(~ mean(grepl(., names(ref_pools), ignore.case = TRUE)))

# Weighted mean gene expression for CS2 pools (using reference pools weights)
cs2_pools_mgx <- weights %>%
  imap_dfc(~ .x * rowMeans(select(cs2_pools, matches(.y)))) %>%
  transmute(Name = cs2_pools[["Name"]], cs2_exp = rowSums(.))

# Mean gene expression for CS3 reference pools
cs3_pools_mgx <-
  enframe(rowMeans(ref_pools), name = "Name", value = "cs3_exp")

# Extract cs2 norm counts (not pools)
cs2_norm_counts <- cs2_norm %>% select(-one_of(names(cs2_pools)[-1]))

# Normalize each gene by adding batch effect (diff in mean gx)
cs2_normalized_data_pools <-
  inner_join(cs3_pools_mgx, cs2_pools_mgx, by = "Name") %>%
  transmute(Name, be = cs3_exp - cs2_exp) %>%
  inner_join(cs2_norm_counts, by = "Name") %>%
  gather(FileName, exp, -1:-4) %>%
  transmute(Name = fct_inorder(Name), FileName, exp = be + exp) %>%
  spread(Name, exp)

# Extract cs3 norm counts (not pools)
cs3_norm_counts <- cs3_norm %>% select(-one_of(names(ref_pools)))

# Normalize each gene by adding batch effect (diff in mean gx)
cs3_normalized_data_pools <-
  inner_join(cs2_pools_mgx, cs3_pools_mgx, by = "Name") %>%
  transmute(Name, be = cs2_exp - cs3_exp) %>%
  inner_join(cs3_norm_counts, by = "Name") %>%
  gather(FileName, exp, -1:-4) %>%
  transmute(Name = fct_inorder(Name), FileName, exp = be + exp) %>%
  spread(Name, exp)

## Summary
# CS2 pools: 9 samples, 365 genes # dim(t(cs2_pools))
# CS3 pools: 22 samples, 513 genes # dim(t(ref_pools))
# CS2 normalized: 1214 samples (1223 original - 9 from pools), 136 common genes # dim(cs2_normalized_data_pools)
# CS3 normalized: 1295 samples (1317 original - 22 from pools), 136 common genes # dim(cs2_normalized_data_pools)
