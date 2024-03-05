`%>%` <- magrittr::`%>%`

stats <- readRDS(list.files(
  path = file.path(outputDir, "merge_eval"),
  pattern = paste0(dataset, "\\.rds"),
  full.names = TRUE
))

df <- stats %>%
  purrr::imap_dfr(function(x, a) {
    purrr::imap_dfr(x, ~ {
      t(.x) %>%
        as.data.frame() %>%
        purrr::set_names(paste0("percentile_", gsub("%", "", names(.)))) %>%
        tibble::rownames_to_column("measure") %>%
        tibble::add_column(sampling = .y, .before = 1)
    }) %>%
      tibble::add_column(algorithm = a, .before = 1)
  }) %>%
  tibble::add_column(dataset = dataset, .before = 1) %>%
  tibble::as_tibble()

# write results to file
saveRDS(
  df,
  file.path(outputDir, "summary", dataset, paste0("iv_summary_", dataset, ".rds"))
)
