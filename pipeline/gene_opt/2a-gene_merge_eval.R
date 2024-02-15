# All combinations
`%>%` <- magrittr::`%>%`
ngenes <- rlang::set_names(seq_len(ngenes))

# All evaluation files
eval_files <- purrr::map(ngenes, function(n) {
  list.files(
    path = file.path(outputDir, "gene_opt", "train_eval", dataset),
    pattern = paste0("add", n, "_"),
    full.names = TRUE
  )
})

# Compute median + 95% CI of evaluations within algorithm & subsampling, merge
eval_merged <- eval_files %>%
  purrr::map(~ purrr::map(., readRDS)) %>%
  purrr::map(purrr::transpose) %>%
  purrr::map_depth(2, ~ do.call(cbind, .x) %>%
                     rlang::set_names(seq_len(ncol(.)))) %>%
  purrr::list_flatten() %>%
  purrr::imap(~ {
    .x %>%
      tibble::rownames_to_column("measure") %>%
      rlang::set_names(c("measure", paste0("Add", .y, "_B", names(.x))))
  }) %>%
  purrr::reduce(dplyr::full_join, by = "measure") %>%
  tidyr::pivot_longer(
    cols = where(is.numeric),
    names_to = c("Genes", "Algorithm", "Bootstrap"),
    names_pattern = "(Add.*)_(.*)_(B.*)",
    values_to = "value"
  ) %>%
  dplyr::summarize(quants = list(quantile(
    value, probs = c(0.5, 0.05, 0.95), na.rm = TRUE
  )),
  .by = c(measure, Genes, Algorithm)) %>%
  tidyr::unnest_wider(col = quants)

# Write all evaluations merged
saveRDS(
  eval_merged,
  file.path(outputDir, "gene_opt", "merge_eval",
            paste0("merge_eval_", dataset, ".rds"))
)
