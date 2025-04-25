# Setup -------------------------------------------------------------------

source(here::here("src/funs.R"))

# Pairwise CodeSet comparisons
codesets <- c("CS1", "CS2", "CS3")
all_codesets <- combn(codesets, 2) %>%
  as_tibble(.name_repair = "unique") %>%
  set_names(map_chr(., paste, collapse = "_vs_"))

# Grouping dataset used for random selection of common samples
group_df <- hist |>
  filter(FileName %in% c(cs1_clean$FileName, cs2_clean$FileName, cs3_clean$FileName))


# Reference Method: 3 Common Samples --------------------------------------

# Random selection of common samples with equal number of histotypes
set.seed(2020)
hist_rand3 <- group_df %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 3) %>%
  ungroup()

# Remove common samples from CS1, preserving gene order
cs1_counts3 <- join_avg(cs1_clean, hist_rand3, "ottaID", "discard")
cs2_counts3 <- join_avg(cs2_clean, hist_rand3, "ottaID", "discard")
cs3_counts3 <- join_avg(cs3_clean, hist_rand3, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand3 <-
  normalize_random(
    x = cs1_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 3,
    strata = "revHist",
    seed = 2020
  )

cs2_norm_rand3 <-
  normalize_random(
    x = cs2_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 3,
    strata = "revHist",
    seed = 2020
  )

# Combined gene expression
counts3 <-
  set_names(list(cs1_counts3, cs2_counts3, cs3_counts3), codesets)
norm_rand3 <-
  set_names(list(cs1_norm_rand3, cs2_norm_rand3, cs3_counts3), codesets)

# Concordance measures for all genes averaged across samples
metrics_non3 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts3[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

metrics_rand3 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand3[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

# Plot all combinations of cross-codeset concordance measure histograms
p_non3 <- ggplot(metrics_non3, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 40, label = Median),
            hjust = 0,
            vjust = 1,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random3 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non3)

p_rand3 <- ggplot(metrics_rand3, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 40, label = Median),
            hjust = 0,
            vjust = 1,
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
hist_rand2 <- group_df %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 2) %>%
  ungroup()

# Remove common samples from CS1, preserving gene order
cs1_counts2 <- join_avg(cs1_clean, hist_rand2, "ottaID", "discard")
cs2_counts2 <- join_avg(cs2_clean, hist_rand2, "ottaID", "discard")
cs3_counts2 <- join_avg(cs3_clean, hist_rand2, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand2 <-
  normalize_random(
    x = cs1_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 2,
    strata = "revHist",
    seed = 2020
  )

cs2_norm_rand2 <-
  normalize_random(
    x = cs2_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 2,
    strata = "revHist",
    seed = 2020
  )

# Combined gene expression
counts2 <-
  set_names(list(cs1_counts2, cs2_counts2, cs3_counts2), codesets)
norm_rand2 <-
  set_names(list(cs1_norm_rand2, cs2_norm_rand2, cs3_counts2), codesets)

# Concordance measures for all genes averaged across samples
metrics_non2 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts2[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

metrics_rand2 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand2[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

# Plot all combinations of cross-codeset concordance measure histograms
p_non2 <- ggplot(metrics_non2, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 40, label = Median),
            hjust = 0,
            vjust = 1,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random2 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non2)

p_rand2 <- ggplot(metrics_rand2, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 50, label = Median),
            hjust = 0,
            vjust = 1,
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
hist_rand1 <- group_df %>%
  distinct(ottaID, revHist) %>%
  group_by(revHist) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Remove common samples from CS1, preserving gene order
cs1_counts1 <- join_avg(cs1_clean, hist_rand1, "ottaID", "discard")
cs2_counts1 <- join_avg(cs2_clean, hist_rand1, "ottaID", "discard")
cs3_counts1 <- join_avg(cs3_clean, hist_rand1, "ottaID", "discard")

# Normalize by reference method using common samples, add histotypes from annot
cs1_norm_rand1 <-
  normalize_random(
    x = cs1_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 1,
    strata = "revHist",
    seed = 2020
  )

cs2_norm_rand1 <-
  normalize_random(
    x = cs2_clean,
    ref = cs3_clean,
    group_df = group_df,
    n = 1,
    strata = "revHist",
    seed = 2020
  )

# Combined gene expression
counts1 <-
  set_names(list(cs1_counts1, cs2_counts1, cs3_counts1), codesets)
norm_rand1 <-
  set_names(list(cs1_norm_rand1, cs2_norm_rand1, cs3_counts1), codesets)

# Concordance measures for all genes averaged across samples
metrics_non1 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(counts1[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

metrics_rand1 <- all_codesets %>%
  imap_dfr(~ {
    pmap_dfr(norm_rand1[.x], ~ cor_stats(.x, .y)) %>% mutate(Sites = .y)
  }) %>%
  gather(key = "Metric", value = "Expression", -Sites) %>%
  group_by(Sites, Metric) %>%
  mutate(Median = paste0("Median = ", scales::number(median(Expression), accuracy = 0.01))) %>%
  ungroup()

# Plot all combinations of cross-codeset concordance measure histograms
p_non1 <- ggplot(metrics_non1, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 40, label = Median),
            hjust = 0,
            vjust = 1,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random1 Non-Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_non1)

p_rand1 <- ggplot(metrics_rand1, aes(Expression)) +
  geom_histogram(bins = 30, fill = "blue") +
  geom_text(aes(x = 0, y = 40, label = Median),
            hjust = 0,
            vjust = 1,
            check_overlap = TRUE) +
  facet_grid(rows = vars(Sites), cols = vars(Metric), scales = "free_x") +
  labs(y = "Count",
       title = "Random1 Normalized Concordance Measure Distributions") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
print(p_rand1)
