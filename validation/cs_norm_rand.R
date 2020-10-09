# Setup -------------------------------------------------------------------

source(here::here("src/funs.R"))

# Pairwise CodeSet comparisons
codesets <- c("CS1", "CS2", "CS3")
all_codesets <- combn(codesets, 2) %>%
  as_tibble() %>%
  set_names(map_chr(., paste, collapse = "_vs_"))


# Reference Method: 3 Common Samples --------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand3 <- hist %>%
  filter(ottaID %in% unique(c(
    cs1_clean$ottaID, cs2_clean$ottaID, cs3_clean$ottaID
  ))) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 3) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand3 <- join_avg(cs1_clean, hist_rand3, "ottaID", "keep")
cs2_rand3 <- join_avg(cs2_clean, hist_rand3, "ottaID", "keep")
cs3_rand3 <- join_avg(cs3_clean, hist_rand3, "ottaID", "keep")

# Remove common samples from CS1, preserving gene order
cs1_counts3 <- join_avg(cs1_clean, hist_rand3, "ottaID", "discard")
cs2_counts3 <- join_avg(cs2_clean, hist_rand3, "ottaID", "discard")
cs3_counts3 <- join_avg(cs3_clean, hist_rand3, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand3 <-
  refMethod(cs1_counts3, cs1_rand3, cs3_rand3) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs2_norm_rand3 <-
  refMethod(cs2_counts3, cs2_rand3, cs3_rand3) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs3_norm_rand3 <-
  refMethod(cs3_counts3, cs3_rand3, cs2_rand3) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

# Combined gene expression
counts3 <-
  set_names(list(cs1_counts3, cs2_counts3, cs3_counts3), codesets)
norm_rand3 <- list(cs1_norm_rand3, cs2_norm_rand3, cs3_norm_rand3) %>%
  set_names(codesets) %>%
  map(~ select(., -"revHist"))

# Concordance measures for all genes averaged across samples
metrics_non3 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts3[.x], ~ {
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

metrics_rand3 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand3[.x], ~ {
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

# Plot all combinations of cross-codeset concordance measure histograms
p_non3 <- ggplot(metrics_non3, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 30, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random3 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non3)

p_rand3 <- ggplot(metrics_rand3, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 15, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random3 Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_rand3)


# Reference Method: 2 Common Samples --------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand2 <- hist %>%
  filter(ottaID %in% unique(c(
    cs1_clean$ottaID, cs2_clean$ottaID, cs3_clean$ottaID
  ))) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 2) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand2 <- join_avg(cs1_clean, hist_rand2, "ottaID", "keep")
cs2_rand2 <- join_avg(cs2_clean, hist_rand2, "ottaID", "keep")
cs3_rand2 <- join_avg(cs3_clean, hist_rand2, "ottaID", "keep")

# Remove common samples from CS1, preserving gene order
cs1_counts2 <- join_avg(cs1_clean, hist_rand2, "ottaID", "discard")
cs2_counts2 <- join_avg(cs2_clean, hist_rand2, "ottaID", "discard")
cs3_counts2 <- join_avg(cs3_clean, hist_rand2, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand2 <-
  refMethod(cs1_counts2, cs1_rand2, cs3_rand2) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs2_norm_rand2 <-
  refMethod(cs2_counts2, cs2_rand2, cs3_rand2) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs3_norm_rand2 <-
  refMethod(cs3_counts2, cs3_rand2, cs2_rand2) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

# Combined gene expression
counts2 <-
  set_names(list(cs1_counts2, cs2_counts2, cs3_counts2), codesets)
norm_rand2 <- list(cs1_norm_rand2, cs2_norm_rand2, cs3_norm_rand2) %>%
  set_names(codesets) %>%
  map(~ select(., -"revHist"))

# Concordance measures for all genes averaged across samples
metrics_non2 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts2[.x], ~ {
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

metrics_rand2 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand2[.x], ~ {
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

# Plot all combinations of cross-codeset concordance measure histograms
p_non2 <- ggplot(metrics_non2, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 30, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random2 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non2)

p_rand2 <- ggplot(metrics_rand2, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 15, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random2 Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_rand2)


# Reference Method: 1 Common Sample ---------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand1 <- hist %>%
  filter(ottaID %in% unique(c(
    cs1_clean$ottaID, cs2_clean$ottaID, cs3_clean$ottaID
  ))) %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Gene expression from random common samples, preserving gene order
cs1_rand1 <- join_avg(cs1_clean, hist_rand1, "ottaID", "keep")
cs2_rand1 <- join_avg(cs2_clean, hist_rand1, "ottaID", "keep")
cs3_rand1 <- join_avg(cs3_clean, hist_rand1, "ottaID", "keep")

# Remove common samples from CS1, preserving gene order
cs1_counts1 <- join_avg(cs1_clean, hist_rand1, "ottaID", "discard")
cs2_counts1 <- join_avg(cs2_clean, hist_rand1, "ottaID", "discard")
cs3_counts1 <- join_avg(cs3_clean, hist_rand1, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand1 <-
  refMethod(cs1_counts1, cs1_rand1, cs3_rand1) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs2_norm_rand1 <-
  refMethod(cs2_counts1, cs2_rand1, cs3_rand1) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

cs3_norm_rand1 <-
  refMethod(cs3_counts1, cs3_rand1, cs2_rand1) %>%
  as.data.frame() %>%
  rownames_to_column("ottaID") %>%
  inner_join(hist %>% distinct(ottaID, revHist), by = "ottaID") %>%
  column_to_rownames("ottaID")

# Combined gene expression
counts1 <-
  set_names(list(cs1_counts1, cs2_counts1, cs3_counts1), codesets)
norm_rand1 <- list(cs1_norm_rand1, cs2_norm_rand1, cs3_norm_rand1) %>%
  set_names(codesets) %>%
  map(~ select(., -"revHist"))

# Concordance measures for all genes averaged across samples
metrics_non1 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts1[.x], ~ {
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

metrics_rand1 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand1[.x], ~ {
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

# Plot all combinations of cross-codeset concordance measure histograms
p_non1 <- ggplot(metrics_non1, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 30, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random1 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non1)

p_rand1 <- ggplot(metrics_rand1, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 15, label = Median),
            hjust = 0,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random1 Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_rand1)
