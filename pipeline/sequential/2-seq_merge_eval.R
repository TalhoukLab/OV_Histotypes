# All combinations
`%>%` <- magrittr::`%>%`
seqs <- seq_len(nseq) %>% rlang::set_names()

# All evaluation files
eval_files <-
  purrr::map(seqs, function(n) {
    list.files(
      path = file.path(outputDir, "sequential", "train_eval"),
      pattern = paste0(seqData, n),
      full.names = TRUE
    )
  })

# Compute median + 95% CI of evaluations within subsampling, merge
# across sequences and resamples
metric_df <- eval_files %>%
  purrr::map(~ purrr::map(., readRDS)) %>%
  purrr::map(purrr::transpose) %>%
  purrr::map_depth(2, ~ do.call(cbind, .x) %>%
                     rlang::set_names(seq_len(ncol(.)))) %>%
  purrr::list_flatten() %>%
  purrr::imap(~ {
    .x %>%
      tibble::rownames_to_column("measure") %>%
      rlang::set_names(c("measure", paste0("Seq", .y, "_B", names(.x))))
  }) %>%
  purrr::reduce(dplyr::full_join, by = "measure") %>%
  dplyr::filter(!grepl("class_0|non-", measure)) %>%
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
  file.path(outputDir, "sequential", "merge_eval", paste0(seqData, "_merge_eval.rds"))
)

# All variable importance files
vi_files <- purrr::map(seq_len(nseq), function(n) {
  list.files(
    path = file.path(outputDir, "sequential", "vi", seqData),
    pattern = paste0("vi_", seqData, n),
    full.names = TRUE
  )
})

# Merge variable importance results
vi_merged <- vi_files %>%
  purrr::map(~ purrr::map(., readRDS)) %>%
  purrr::modify_depth(2, ~ dplyr::bind_rows(.x, .id = "Algorithm")) %>%
  purrr::modify_depth(1, ~ dplyr::bind_rows(.x, .id = "Bootstrap")) %>%
  dplyr::bind_rows(.id = "Sequence")

# Rank variables using average importance scores per sampling method and algorithm
vi_ranked <- vi_merged %>%
  dplyr::group_by(Sequence, Algorithm, Variable) %>%
  dplyr::summarise(Mean_Importance = mean(Importance), .groups = "drop_last") %>%
  dplyr::arrange(Sequence, Algorithm, dplyr::desc(Mean_Importance)) %>%
  dplyr::mutate(Rank = dplyr::dense_rank(dplyr::desc(Mean_Importance))) %>%
  dplyr::ungroup()

# Rank aggregation across sequence
vi_rank_aggd <- vi_ranked %>%
  tidyr::pivot_wider(id_cols = Sequence,
                     names_from = "Rank",
                     values_from = "Variable") %>%
  tibble::column_to_rownames("Sequence") %>%
  as.matrix() %>%
  RankAggreg::RankAggreg(
    x = .,
    k = ncol(.),
    method = "GA",
    seed = 2024,
    verbose = FALSE
  )

# Only consider candidate genes not already in PrOTYPE and SPOT
candidates <- c("C10orf116", "GAD1", "TPX2", "KGFLP2", "EGFL6", "KLK7", "PBX1",
                "LIN28B", "TFF3", "MUC5B", "FUT3", "STC1", "BCL2", "PAX8", "GCNT3",
                "GPR64", "ADCYAP1R1", "IGKC", "BRCA1", "IGJ", "TFF1", "MET",
                "CYP2C18", "CYP4B1", "SLC3A1", "EPAS1", "HNF1B", "IL6", "ATP5G3",
                "DKK4", "SENP8", "CAPN2", "C1orf173", "CPNE8", "IGFBP1", "WT1",
                "TP53", "SEMA6A", "SERPINA5", "ZBED1", "TSPAN8", "SCGB1D2", "LGALS4",
                "MAP1LC3A")

vi_ranked_candidates <- vi_rank_aggd %>%
  purrr::pluck("top.list") %>%
  purrr::keep(~ . %in% candidates)

# Write ranked variable importance
saveRDS(
  vi_ranked_candidates,
  file.path(outputDir, "sequential", "ranked_vi", paste0("ranked_vi_", seqData, ".rds"))
)
