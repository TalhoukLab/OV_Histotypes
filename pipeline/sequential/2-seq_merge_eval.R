# All combinations
`%>%` <- magrittr::`%>%`
seqs <- seq_len(nseq) %>% purrr::set_names()
samps <- purrr::set_names(samps)

# All evaluation files
eval_files <-
  purrr::map(samps, function(s) {
    purrr::map(seqs, function(n) {
      list.files(
        path = file.path(outputDir, "train_eval"),
        pattern = paste0(s, "_seq", n),
        full.names = TRUE
      )
    })
  })

# Compute median + 95% CI of evaluations within subsampling, merge
# across sequences and resamples
eval_merged <- eval_files %>%
  purrr::map_depth(3, readRDS) %>%
  purrr::map_depth(2, ~ purrr::flatten(purrr::transpose(.))) %>%
  purrr::map(~ apply(
    data.frame(purrr::flatten(.)),
    1,
    quantile,
    probs = c(0.5, 0.05, 0.95),
    na.rm = TRUE
  ))

# Write all evaluations merged
saveRDS(
  eval_merged,
  file.path(outputDir, "merge_eval", paste0("merge_eval_sequential.rds"))
)
