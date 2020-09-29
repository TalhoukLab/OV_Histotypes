# Reference Method: 3 Common Samples --------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand3 <- hist %>%
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) %>%
  group_by(CodeSet, revHist) %>%
  slice_sample(n = 3) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand <- cs1_clean %>%
  inner_join(hist_rand3, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_rand <- cs2_clean %>%
  inner_join(hist_rand3, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_rand <- cs3_clean %>%
  inner_join(hist_rand3, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Remove common samples from CS1, preserving gene order
cs1_norm_counts <- cs1_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs1_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_norm_counts <- cs2_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs2_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_norm_counts <- cs3_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs3_rand)))) %>%
  gather(FileName, exp,-1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand3 <-
  refMethod(cs1_norm_counts, cs1_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs2_norm_rand3 <-
  refMethod(cs2_norm_counts, cs2_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs3_norm_rand3 <-
  refMethod(cs3_norm_counts, cs3_rand, cs2_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")


# Reference Method: 2 Common Samples --------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand2 <- hist %>%
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) %>%
  group_by(CodeSet, revHist) %>%
  slice_sample(n = 2) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand <- cs1_clean %>%
  inner_join(hist_rand2, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_rand <- cs2_clean %>%
  inner_join(hist_rand2, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_rand <- cs3_clean %>%
  inner_join(hist_rand2, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Remove common samples from CS1, preserving gene order
cs1_norm_counts <- cs1_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs1_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_norm_counts <- cs2_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs2_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_norm_counts <- cs3_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs3_rand)))) %>%
  gather(FileName, exp,-1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand2 <-
  refMethod(cs1_norm_counts, cs1_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs2_norm_rand2 <-
  refMethod(cs2_norm_counts, cs2_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs3_norm_rand2 <-
  refMethod(cs3_norm_counts, cs3_rand, cs2_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")


# Reference Method: 1 Common Sample ---------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand1 <- hist %>%
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName)) %>%
  group_by(CodeSet, revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand <- cs1_clean %>%
  inner_join(hist_rand1, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_rand <- cs2_clean %>%
  inner_join(hist_rand1, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_rand <- cs3_clean %>%
  inner_join(hist_rand1, by = "FileName") %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Remove common samples from CS1, preserving gene order
cs1_norm_counts <- cs1_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs1_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs2_norm_counts <- cs2_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs2_rand)))) %>%
  gather(FileName, exp, -1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

cs3_norm_counts <- cs3_norm %>%
  select(-c(Code.Class, Accession, paste0("X", rownames(cs3_rand)))) %>%
  gather(FileName, exp,-1) %>%
  spread(Name, exp) %>%
  column_to_rownames("FileName") %>%
  select(all_of(common_genes))

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand1 <-
  refMethod(cs1_norm_counts, cs1_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs2_norm_rand1 <-
  refMethod(cs2_norm_counts, cs2_rand, cs3_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")

cs3_norm_rand1 <-
  refMethod(cs3_norm_counts, cs3_rand, cs2_rand) %>%
  as.data.frame() %>%
  rownames_to_column("FileName") %>%
  mutate(FileName = gsub("^X", "", FileName)) %>%
  inner_join(hist, by = "FileName") %>%
  filter(revHist %in% c("CCOC", "ENOC", "HGSC", "LGSC", "MUC")) %>%
  column_to_rownames("FileName")
