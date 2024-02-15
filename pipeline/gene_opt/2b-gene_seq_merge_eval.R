# All combinations
`%>%` <- magrittr::`%>%`
ngenes <- rlang::set_names(seq_len(ngenes))
nseq <- rlang::set_names(seq_len(nseq))

# All evaluation files
eval_files <-
  purrr::map(nseq, function(s) {
    purrr::map(ngenes, function(n) {
      list.files(
        path = file.path(outputDir, "gene_opt", "sequential", "train_eval"),
        pattern = paste0(sq, s, "_add", n, "_"),
        full.names = TRUE
      )
    })
  })

# Compute median + 95% CI of evaluations within algorithm & subsampling, merge
eval_merged <- eval_files %>%
  purrr::map(~ purrr::map(., ~ purrr::map(., ~ readRDS(.)))) %>%
  purrr::modify_depth(2, purrr::transpose) %>%
  purrr::map(purrr::transpose) %>%
  purrr::flatten()%>%
  purrr::map_depth(2, ~ do.call(cbind, .x) %>%
                     rlang::set_names(seq_len(ncol(.)))) %>%
  purrr::list_flatten() %>%
  purrr::imap(~ {
    .x %>%
      tibble::rownames_to_column("measure") %>%
      rlang::set_names(c("measure", paste0("Add_", .y, "_B", names(.x))))
  }) %>%
  purrr::reduce(dplyr::full_join, by = "measure") %>%
  dplyr::filter(!grepl("non-", measure)) %>%
  tidyr::pivot_longer(
    cols = where(is.numeric),
    names_to = c("Algorithm", "Genes", "Bootstrap"),
    names_pattern = "Add_(.*)_(.*)_(B.*)",
    values_to = "value"
  ) %>%
  dplyr::summarize(quants = list(quantile(
    value, probs = c(0.5, 0.05, 0.95), na.rm = TRUE
  )),
  .by = c(measure, Genes, Algorithm)) %>%
  tidyr::unnest_wider(col = quants) %>%
  dplyr::mutate(Sequence = dplyr::dense_rank(Algorithm),
                .before = Algorithm)

# Write all evaluations merged
saveRDS(
  eval_merged,
  file.path(outputDir, "gene_opt", "sequential", "merge_eval",
            paste0("merge_eval_", sq, ".rds"))
)
