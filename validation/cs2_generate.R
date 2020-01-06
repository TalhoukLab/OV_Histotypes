# CodeSet Processing ------------------------------------------------------
source(here::here("validation/cs_process.R"))


# Pool Method -------------------------------------------------------------

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

## Summary
# CS3 pools: 22 samples, 513 genes # dim(t(ref_pools)
# CS2 pools: 9 samples, 365 genes # dim(t(cs2_pools))
# CS2 normalized: 1214 samples (1223 original - 9 from pools), 136 common genes # dim(cs2_normalized_data_pools)


# Reference Method: Common Samples ----------------------------------------

# Averaged gene expression within duplicate ottaID, ensure same gene order
cs2_avgd <- cs2_clean %>%
  dplyr::select(-FileName) %>%
  dplyr::group_by(ottaID) %>%
  dplyr::summarize_if(is.double, mean) %>%
  tibble::column_to_rownames("ottaID")

cs3_avgd <- cs3_clean %>%
  dplyr::select(-FileName) %>%
  dplyr::group_by(ottaID) %>%
  dplyr::summarize_if(is.double, mean) %>%
  tibble::column_to_rownames("ottaID") %>%
  magrittr::extract(match(names(cs2_avgd), names(.)))

# Common samples removed from CS2, preserving gene order
cs2_norm_counts <- cs2_norm %>%
  dplyr::select(-c(Code.Class, Accession, dplyr::matches(paste(
    cs2_clean[["FileName"]], collapse = "|"
  )))) %>%
  tidyr::gather(FileName, exp, -1) %>%
  tidyr::spread(Name, exp) %>%
  tibble::column_to_rownames("FileName") %>%
  dplyr::select(names(cs2_avgd))

# Normalize by reference method using common samples
cs2_normalized_data_common <-
  nanostringr::refMethod(cs2_norm_counts, cs2_avgd, cs3_avgd) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("FileName")

## Summary
# CS3 common samples: 78 aggregated samples (140 total), 72 genes # dim(cs3_avgd)
# CS2 common samples: 78 aggregated samples (87 total), 72 genes # dim(cs2_avgd)
# CS2 normalized: 1136 samples (1223 original - 87 total samples), 72 genes # dim(cs2_normalized_data_common)


# Concordance Histograms --------------------------------------------------

# Compare methods on common samples and genes (1127 by 72)
cs2_gx1 <- dplyr::semi_join(cs2_normalized_data_pools,
                            cs2_normalized_data_common,
                            by = "FileName") %>%
  dplyr::select(names(cs2_normalized_data_common))

cs2_gx2 <- dplyr::semi_join(cs2_normalized_data_common,
                            cs2_normalized_data_pools,
                            by = "FileName")

# No reference method normalization (only to HK genes)
cs2_gx0 <- cs2_norm %>%
  dplyr::select(Name, cs2_gx1[["FileName"]]) %>%
  dplyr::filter(Name %in% names(cs2_gx1)) %>%
  dplyr::mutate(Name = forcats::fct_inorder(Name)) %>%
  tidyr::gather(FileName, exp, -1) %>%
  tidyr::spread(Name, exp)

# Compare cs2_gx0, cs2_gx1, cs2_gx2 with each other
cs_combined <-
  tibble::lst(cs2_gx0, cs2_gx1, cs2_gx2) %>%
  dplyr::bind_rows(.id = "dataset") %>%
  tidyr::gather(gene, exp, -1:-2, factor_key = TRUE) %>%
  tidyr::spread(dataset, exp)

# Concordance measures between samples of same gene, for each method comparison
all_metrics <- cs_combined %>%
  tidyr::gather(Methods, cs2_val_norm, dplyr::matches("cs2")) %>%
  dplyr::left_join(., ., by = c("FileName", "gene")) %>%
  dplyr::filter(Methods.x < Methods.y) %>%
  dplyr::mutate_at(
    dplyr::vars(dplyr::matches("Methods")),
    dplyr::recode,
    cs2_gx0 = "None",
    cs2_gx1 = "Pools",
    cs2_gx2 = "Common"
  ) %>%
  tidyr::unite(Methods, dplyr::matches("Methods"), sep = "_vs_") %>%
  dplyr::group_by(Methods = forcats::fct_inorder(Methods), FileName) %>%
  dplyr::summarize(
    R2 = cor(cs2_val_norm.x, cs2_val_norm.y) ^ 2,
    Ca = epiR::epi.ccc(cs2_val_norm.x, cs2_val_norm.y)[["C.b"]],
    Rc = epiR::epi.ccc(cs2_val_norm.x, cs2_val_norm.y)[["rho.c"]][["est"]]
  ) %>%
  dplyr::ungroup() %>%
  tidyr::gather(Metric, Expression, c("R2", "Ca", "Rc"), factor_key = TRUE)

# Plot all combinations of cross-method concordance measure histograms
p <- ggplot(all_metrics, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  facet_grid(rows = vars(Methods), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "CS2 Cross-Method Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p)


# Add Histotype -----------------------------------------------------------

cs2_norm_hist <- cs2_normalized_data_pools %>%
  tibble::as_tibble() %>%
  dplyr::mutate(FileName = gsub("^X", "", FileName)) %>%
  dplyr::inner_join(hist, by = "FileName") %>%
  dplyr::filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  tibble::column_to_rownames("FileName")

# Counts
cs2_norm_hist %>% dplyr::count(revHist)
cs2_norm_hist %>% dplyr::count(hist_gr)


# Save Objects ------------------------------------------------------------

# Use pool method as it retains more samples and genes
# Split into training data and training labels, save objects
cs2_data <- dplyr::select_if(cs2_norm_hist, is.double)
cs2_class <- cs2_norm_hist[["revHist"]]

saveRDS(cs2_data, here::here("data/cs2_data.rds"))
saveRDS(cs2_class, here::here("data/cs2_class.rds"))
