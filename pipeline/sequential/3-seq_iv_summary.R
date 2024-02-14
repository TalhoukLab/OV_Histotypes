`%>%` <- magrittr::`%>%`

stats <- readRDS(list.files(
  path = file.path(outputDir, "sequential", "merge_eval"),
  pattern = seqData,
  full.names = TRUE
))

# Best sampling method from classification of full training set
seq_top <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
samp <- as.character(seq_top[["sampling"]])

df <- stats %>%
  dplyr::rename_with(~ paste0("percentile_", gsub("%", "", .)), where(is.numeric)) %>%
  tibble::add_column(
    dataset = "train",
    algorithm = seqData,
    sampling = samp,
    .before = 1
  ) %>%
  tibble::as_tibble()

# write results to file
saveRDS(
  df,
  file.path(outputDir, "sequential", "summary", paste0("iv_summary_", seqData, ".rds"))
)
