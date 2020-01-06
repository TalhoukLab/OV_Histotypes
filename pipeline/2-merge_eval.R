# All combinations
`%>%` <- magrittr::`%>%`
algs <- purrr::set_names(algs)
samps <- purrr::set_names(samps)

# All evaluation files
eval_files <-
  purrr::map(algs, function(a) {
    purrr::map(samps, function(s) {
      list.files(
        path = file.path(outputDir, "train_eval", dataset),
        pattern = paste(a, s, sep = "_"),
        full.names = TRUE
      )
    })
  })

# Compute median + 95% CI of evaluations within algorithm & subsampling, merge
eval_merged <- eval_files %>%
  purrr::map(~ purrr::map(., ~ purrr::map(., readRDS))) %>%
  purrr::modify_depth(2, purrr::transpose) %>%
  purrr::map(purrr::transpose) %>%
  purrr::flatten() %>%
  purrr::map(~ purrr::map(., ~ apply(
    data.frame(.), 1, quantile, c(0.5, 0.05, 0.95), na.rm = TRUE
  )))

# Write all evaluations merged
saveRDS(
  eval_merged,
  file.path(outputDir, "merge_eval", paste0("merge_eval_", dataset, ".rds"))
)
