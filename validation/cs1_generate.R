# Reference Method: Common Samples ----------------------------------------

# Averaged gene expression within duplicate ottaID, ensure same gene order
cs1_avgd <- cs1_clean %>%
  dplyr::select(-FileName) %>%
  dplyr::group_by(ottaID = forcats::fct_inorder(ottaID)) %>%
  dplyr::summarize_if(is.double, mean) %>%
  tibble::column_to_rownames("ottaID")

cs3_avgd <- cs3_clean %>%
  dplyr::select(-FileName) %>%
  dplyr::group_by(ottaID = forcats::fct_inorder(ottaID)) %>%
  dplyr::summarize_if(is.double, mean) %>%
  tibble::column_to_rownames("ottaID") %>%
  magrittr::extract(match(names(cs1_avgd), names(.)))

# Common samples removed from CS1, preserving gene order
cs1_norm_counts <- cs1_norm %>%
  dplyr::select(-c(Code.Class, Accession, dplyr::matches(paste(
    cs1_clean[["FileName"]], collapse = "|"
  )))) %>%
  tidyr::gather(FileName, exp, -1) %>%
  tidyr::spread(Name, exp) %>%
  tibble::column_to_rownames("FileName") %>%
  dplyr::select(names(cs3_avgd))

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_common <-
  nanostringr::refMethod(cs1_norm_counts, cs1_avgd, cs3_avgd) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("FileName") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(FileName = gsub("^X", "", FileName)) %>%
  dplyr::inner_join(hist, by = "FileName") %>%
  dplyr::filter(revHist %in% c("CCC", "ENOCa", "HGSC", "LGSC", "MUC")) %>%
  tibble::column_to_rownames("FileName")

# Counts
cs1_norm_common %>% dplyr::count(revHist)
cs1_norm_common %>% dplyr::count(hist_gr)

## Summary
# Number of genes: 72 # ncol(cs1_norm_common) - 3
# CS3 common samples: 78 aggregated samples (140 total) # nrow(cs3_avgd)
# CS1 common samples: 78 aggregated samples (93 total) # nrow(cs1_avgd)
# CS1 normalized: 240 samples (412 original - 93 total samples - 79 other hists) # nrow(cs1_norm_common)


# Reference Method: Reference Samples -------------------------------------

# CS1 gene order
cs1_genes <- names(cs1_clean)[-1:-2]

# Reference set: duplicate samples from common subset
cs1_ref <- cs1_clean %>%
  dplyr::semi_join(dplyr::filter(dplyr::count(., ottaID), n > 1), by = "ottaID")

# Validation set: reference set removed from all samples, preserving gene order
cs1_val <- cs1_norm %>%
  dplyr::select(-dplyr::one_of(paste0("X", cs1_ref[["FileName"]]), "Code.Class", "Accession")) %>%
  tidyr::gather(FileName, value, -1) %>%
  tidyr::spread(Name, value) %>%
  dplyr::transmute(FileName = gsub("^X", "", FileName), !!!rlang::syms(cs1_genes))

# Remove non-numeric column names
cs1_ref2 <- dplyr::select_if(cs1_ref, is.double)
cs1_val2 <- dplyr::select_if(cs1_val, is.double)

# Corresponding samples in CS3 found in CS1 reference set
cs3_ref <- cs3_clean %>%
  dplyr::semi_join(cs1_ref, by = "ottaID") %>%
  dplyr::select(cs1_genes)

# Normalize by reference method using reference samples, add histotypes from annot
cs1_norm_ref <-
  nanostringr::refMethod(cs1_val2, cs1_ref2, cs3_ref) %>%
  tibble::as_tibble() %>%
  tibble::add_column(FileName = cs1_val[["FileName"]], .before = 1) %>%
  dplyr::inner_join(hist, by = "FileName") %>%
  dplyr::filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  tibble::column_to_rownames("FileName")

# Counts
cs1_norm_ref %>% dplyr::count(revHist)
cs1_norm_ref %>% dplyr::count(hist_gr)

## Summary
# Number of genes: 72 # ncol(cs1_norm_ref) - 3
# CS3 reference samples: 20 samples (140 total) # nrow(cs3_ref)
# CS1 reference samples: 25 samples (93 total) # nrow(cs1_ref2)
# CS1 normalized: 304 samples (412 original - 25 reference samples - 83 other hists) # nrow(cs1_norm_ref)


# Save Objects ------------------------------------------------------------

# Use the reference samples as it retains more data
# Split into training data and training labels, save objects
cs1_data <- dplyr::select_if(cs1_norm_ref, is.double)
cs1_class <- cs1_norm_ref[["revHist"]]

# saveRDS(cs1_data, here::here("data/cs1_data.rds"))
# saveRDS(cs1_class, here::here("data/cs1_class.rds"))
