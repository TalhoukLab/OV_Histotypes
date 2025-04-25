# Setup -------------------------------------------------------------------

source(here::here("validation/cs_process_cohorts.R"))


# Housekeeping Genes Normalization ----------------------------------------

cs2_norm_coh <- HKnorm(cs2_coh)
cs3_norm_coh <- HKnorm(cs3_coh)


# Pools Method ------------------------------------------------------------

# Pool set: pool samples from CS2
pool_samples <-
  gsub("-|\\.RCC", "", grep("POOL", pools[["CS2-FileName"]], value = TRUE))
cs2_pools <- select(cs2_norm_coh, Name, paste0("X", pool_samples))

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
cs2_norm_counts <- cs2_norm_coh %>% select(-one_of(names(cs2_pools)[-1]))

# Normalize each gene by adding batch effect (diff in mean gx)
cs2_normalized_data_pools <-
  inner_join(cs2_pools_mgx, cs3_pools_mgx, by = "Name") %>%
  transmute(Name, be = cs3_exp - cs2_exp) %>%
  inner_join(cs2_norm_counts, by = "Name") %>%
  gather(FileName, exp, -1:-4) %>%
  transmute(Name = fct_inorder(Name), FileName, exp = be + exp) %>%
  spread(Name, exp)

# Extract cs3 norm counts (not pools)
cs3_norm_counts <- cs3_norm_coh %>% select(-one_of(names(ref_pools)))

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


# Compare Common Samples --------------------------------------------------

# Keep only common samples between codesets, average out duplicates
tmp2 <- cs2_normalized_data_pools %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName")

tmp3 <- cs3_normalized_data_pools %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName")

cs2_common_pools <- tmp2 %>%
  semi_join(tmp3, by = "ottaID") %>%
  group_by(ottaID) %>%
  summarize_if(is.double, mean) %>%
  ungroup()

cs3_common_pools <- tmp3 %>%
  semi_join(tmp2, by = "ottaID") %>%
  group_by(ottaID) %>%
  summarize_if(is.double, mean) %>%
  ungroup()

# Keep common samples/genes for the non-normalized data
cs2_common_non <- cs2_norm_counts %>%
  gather(FileName, exp, -1:-3) %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  semi_join(cs2_common_pools, by = "ottaID") %>%
  filter(Name %in% names(cs2_common_pools)[-1]) %>%
  mutate(Name = factor(Name, levels = names(cs2_common_pools)[-1])) %>%
  select(Name, ottaID, exp) %>%
  group_by(Name, ottaID) %>%
  summarize(exp = mean(exp)) %>%
  ungroup() %>%
  spread(Name, exp)

cs3_common_non <- cs3_norm_counts %>%
  gather(FileName, exp, -1:-3) %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  semi_join(cs3_common_pools, by = "ottaID") %>%
  filter(Name %in% names(cs3_common_pools)[-1]) %>%
  mutate(Name = factor(Name, levels = names(cs3_common_pools)[-1])) %>%
  select(Name, ottaID, exp) %>%
  group_by(Name, ottaID) %>%
  summarize(exp = mean(exp)) %>%
  ungroup() %>%
  spread(Name, exp)

# Combined gene expression
all_comps <- tibble(
  `CS2Non_vs_CS3Non` = c("CS2Non", "CS3Non"),
  `CS2Pools_vs_CS3Pools` = c("CS2Pools", "CS3Pools")
)
common_gx <-
  set_names(
    list(
      cs2_common_non,
      cs3_common_non,
      cs2_common_pools,
      cs3_common_pools
    ),
    c("CS2Non", "CS3Non", "CS2Pools", "CS3Pools")
  ) %>%
  map(column_to_rownames, "ottaID")

# Concordance measures for all genes averaged across samples
metrics_pools <- all_comps %>%
  imap_dfr(~ {
    pmap_dfr(common_gx[.x], ~ {
      R2 <- cor(.x, .y) ^ 2
      ccc <- epiR::epi.ccc(.x, .y)
      Ca <- pluck(ccc, "C.b")
      Rc <- pluck(ccc, "rho.c", "est")
      lst(R2, Ca, Rc)
    }) %>%
      mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

# Plot concordance measure histograms
p_pools <- ggplot(metrics_pools, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 30, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Pools Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_pools)
