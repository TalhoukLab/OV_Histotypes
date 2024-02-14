# Best sampling method from classification of full training set
seq_top <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
samp <- as.character(seq_top[["sampling"]])

# All combinations
`%>%` <- magrittr::`%>%`
seqs <- seq_len(nseq) %>% purrr::set_names()
samps <- purrr::set_names(samp)

# All evaluation files
eval_files <-
  purrr::map(seqs, function(n) {
    list.files(
      path = file.path(outputDir, "sequential", "train_eval"),
      pattern = paste0("seq", n),
      full.names = TRUE
    )
  })

# Compute median + 95% CI of evaluations within subsampling, merge
# across sequences and resamples
eval_merged <- eval_files %>%
  purrr::map_depth(3, readRDS) %>%
  purrr::map_depth(2, ~ {
    purrr::transpose(.) %>%
      purrr::flatten() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("measure")
  }) %>%
  purrr::map(~ {
    purrr::reduce(.x, dplyr::full_join, by = "measure") %>%
      dplyr::filter(!grepl("class_0", measure)) %>%
      tibble::column_to_rownames("measure") %>%
      apply(1, quantile, probs = c(0.5, 0.05, 0.95), na.rm = TRUE)
  })

# Write all evaluations merged
saveRDS(
  eval_merged,
  file.path(outputDir, "sequential", "merge_eval", paste0("merge_eval_sequential.rds"))
)
