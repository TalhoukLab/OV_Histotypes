# All combinations
`%>%` <- magrittr::`%>%`
seqs <- seq_len(nseq) %>% purrr::set_names()

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
metric_df <- eval_files %>%
  purrr::map(~ purrr::map(., readRDS)) %>%
  purrr::map(purrr::transpose) %>%
  purrr::map_depth(2, ~ do.call(cbind, .x) %>%
                     set_names(seq_len(ncol(.)))) %>%
  purrr::list_flatten() %>%
  purrr::imap(~ {
    .x %>%
      tibble::rownames_to_column("measure") %>%
      set_names(c("measure", paste0("Seq", .y, "_B", names(.x))))
  }) %>%
  purrr::reduce(dplyr::full_join, by = "measure") %>%
  dplyr::filter(!grepl("class_0", measure)) %>%
  tidyr::pivot_longer(
    cols = where(is.numeric),
    names_to = c("Sequence", "Algorithm", "Bootstrap"),
    names_pattern = "(Seq.)_(.*)_(B.*)",
    values_to = "value"
  ) %>%
  dplyr::mutate(per_class = grepl("\\.", measure))

overall_metrics <- metric_df %>%
  dplyr::filter(!per_class) %>%
  dplyr::summarize(quants = list(quantile(
    value, probs = c(0.5, 0.05, 0.95), na.rm = TRUE
  )),
  .by = c(measure)) %>%
  tidyr::unnest_wider(col = quants)

per_class_metrics <- metric_df %>%
  dplyr::filter(per_class, !is.na(value)) %>%
  dplyr::summarize(quants = list(quantile(
    value, probs = c(0.5, 0.05, 0.95), na.rm = TRUE
  )),
  .by = c(measure)) %>%
  tidyr::unnest_wider(col = quants)

all_metrics <- dplyr::bind_rows(overall_metrics, per_class_metrics)

# Write all evaluations merged
saveRDS(
  all_metrics,
  file.path(outputDir, "sequential", "merge_eval", paste0("merge_eval_sequential.rds"))
)
