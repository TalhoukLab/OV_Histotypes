`%>%` <- magrittr::`%>%`

stats <- readRDS(list.files(
  path = file.path(outputDir, "sequential", "merge_eval"),
  full.names = TRUE
))

df <- stats %>%
  purrr::imap_dfr(~ {
    t(.x) %>%
      as.data.frame() %>%
      purrr::set_names(paste0("percentile_", gsub("%", "", names(.)))) %>%
      tibble::rownames_to_column("measure") %>%
      tibble::add_column(algorithm = "sequential", sampling = .y, .before = 1)
  }) %>%
  tibble::add_column(dataset = "train", .before = 1) %>%
  tibble::as_tibble()

# write results to file
saveRDS(
  df,
  file.path(outputDir, "sequential", "summary", "iv_summary_sequential.rds")
)
