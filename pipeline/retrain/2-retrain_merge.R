# Best sampling method from classification of full training set
seq_top <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
samp <- seq_top[["sampling"]]

# Algorithms
`%>%` <- magrittr::`%>%`
algs <- purrr::set_names(algs)
samps <- purrr::set_names(samps)

# All evaluation files
eval_files <-
  purrr::map(algs, function(a) {
    purrr::map(samps, function(s) {
      list.files(
        path = file.path(outputDir, "retrain_eval", dataset),
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

# All variable importance files
vi_files <-
  purrr::map(algs, function(a) {
    purrr::map(samps, function(s) {
      list.files(
        path = file.path(outputDir, "retrain_vi", dataset),
        pattern = paste(a, s, sep = "_"),
        full.names = TRUE
      )
    })
  })

# Merge variable importance results
vi_merged <- vi_files %>%
  purrr::map(~ purrr::map(., ~ purrr::map(., readRDS))) %>%
  purrr::modify_depth(3, ~ dplyr::bind_rows(.x, .id = "Algorithm")) %>%
  purrr::modify_depth(2, ~ dplyr::bind_rows(.x, .id = "Bootstrap")) %>%
  purrr::map(~ .x %>%
               set_names(samps) %>%
               dplyr::bind_rows(.id = "Sampling")) %>%
  dplyr::bind_rows()

# Rank variables using average importance scores per sampling method and algorithm
vi_ranked <- vi_merged %>%
  dplyr::group_by(Sampling, Algorithm, Variable) %>%
  dplyr::summarise(Mean_Importance = mean(Importance), .groups = "drop_last") %>%
  dplyr::arrange(Sampling, Algorithm, dplyr::desc(Mean_Importance)) %>%
  dplyr::mutate(Rank = dplyr::dense_rank(dplyr::desc(Mean_Importance))) %>%
  dplyr::ungroup()

# Only consider candidate genes not already in PrOTYPE and SPOT
candidates <- c("C10orf116", "GAD1", "TPX2", "KGFLP2", "EGFL6", "KLK7", "PBX1",
                "LIN28B", "TFF3", "MUC5B", "FUT3", "STC1", "BCL2", "PAX8", "GCNT3",
                "GPR64", "ADCYAP1R1", "IGKC", "BRCA1", "IGJ", "TFF1", "MET",
                "CYP2C18", "CYP4B1", "SLC3A1", "EPAS1", "HNF1B", "IL6", "ATP5G3",
                "DKK4", "SENP8", "CAPN2", "C1orf173", "CPNE8", "IGFBP1", "WT1",
                "TP53", "SEMA6A", "SERPINA5", "ZBED1", "TSPAN8", "SCGB1D2", "LGALS4",
                "MAP1LC3A")

vi_ranked_candidates <- vi_ranked %>%
  dplyr::filter(Variable %in% candidates)

# Write ranked variable importance
saveRDS(
  vi_ranked_candidates,
  file.path(outputDir, "ranked_vi", paste0("ranked_vi_", dataset, ".rds"))
)
