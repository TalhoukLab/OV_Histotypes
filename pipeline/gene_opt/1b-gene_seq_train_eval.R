`%>%` <- magrittr::`%>%`
normalize <- function(data, norm_by = c("none", "genes", "samples"),
                      norm_type = "conventional", min_var = 0.5) {
  norm_by <- match.arg(norm_by)
  switch(
    norm_by,
    none = data %>%
      diceR::prepare_data(scale = FALSE,
                          min.var = min_var,
                          type = norm_type) %>%
      as.data.frame(),
    genes = data %>%
      diceR::prepare_data(scale = TRUE,
                          min.var = min_var,
                          type = norm_type) %>%
      as.data.frame(),
    samples = data %>%
      t() %>%
      diceR::prepare_data(scale = TRUE,
                          min.var = min_var,
                          type = norm_type) %>%
      t() %>%
      as.data.frame()
  )
}

# Import sequential training data, class, and algorithm by sampling method
sq_data <- readRDS(file.path(inputDir, paste0(sq, "_data.rds")))
sq_class <- readRDS(file.path(inputDir, paste0(sq, "_class.rds")))
sq_algs <- readRDS(file.path(inputDir, paste0(sq, "_algs.rds")))

data <- sq_data[[nseq]]
class <- sq_class[[nseq]]
alg <- sq_algs[[nseq]]

# Normalization
data <- normalize(data, norm_by, norm_type, min_var)

# Genes from PrOTYPE and SPOT to keep
base_genes <- c("COL11A1", "CD74", "CD2", "TIMP3", "LUM", "CYTIP", "COL3A1",
                "THBS2", "TCF7L1", "HMGA2", "FN1", "POSTN", "COL1A2", "COL5A2",
                "PDZK1IP1", "FBN1", "HIF1A", "CXCL10", "DUSP4", "SOX17", "MITF",
                "CDKN3", "BRCA2", "CEACAM5", "ANXA4", "SERPINE1", "CRABP2", "DNAJC9")

# Best sampling method from classification of full training set
seq_top <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
samp <- as.character(seq_top[["sampling"]])

# Import variable importance for selected algorithms and subsampling method
ranked_vi <- readRDS(file.path(outputDir, "sequential", "ranked_vi",
                               paste0("ranked_vi_", sq, ".rds")))

# Order of candidate genes
candidate_genes <- ranked_vi %>%
  dplyr::filter(Sequence == nseq) %>%
  dplyr::arrange(Rank) %>%
  dplyr::slice_min(order_by = Rank, n = ngene) %>%
  dplyr::pull(Variable)

# Select genes
data <- data %>%
  dplyr::select(all_of(c(base_genes, candidate_genes)))

# Supervised learning model output
suppressWarnings(
  sm <- splendid::splendid_model(
    data = data,
    class = class,
    algorithms = alg,
    n = 1,
    seed_boot = as.integer(reps),
    seed_samp = 2019,
    seed_alg = 2019,
    sampling = samp,
    stratify = TRUE,
    tune = TRUE
  )
)

# # Extract variable importance results
# vi_df <- sm[["models"]] %>%
#   purrr::imap(~ {
#     mod <- .x[[1]]
#     alg <- .y
#     if (alg %in% c("mlr_lasso", "mlr_ridge", "rf")) {
#       vip::vi(mod)
#     } else if (alg == "svm") {
#       pfun <- function(object, newdata) {
#         caret::predict.train(object, newdata = newdata, type = "prob")[, 1]
#       }
#       mod %>%
#         vip::vi_shap(pred_wrapper = pfun) %>%
#         dplyr::arrange(dplyr::desc(Importance))
#     } else if (alg == "adaboost") {
#       mod %>%
#         maboost::varplot.maboost(plot.it = FALSE,
#                                  type = "scores",
#                                  max.var.show = Inf) %>%
#         tibble::enframe(name = "Variable", value = "Importance")
#     }
#   })

# Write evaluations to file
outputFile <- file.path(
  outputDir, "gene_opt", "sequential", "train_eval",
  paste0(sq, nseq, "_add", ngene, "_", reps, ".rds")
)
saveRDS(sm[["evals"]], outputFile)

# # Write variable importance to file
# viFile <- file.path(outputDir, "sequential", "vi", sq,
#                     paste0("vi_", sq, nseq, "_", reps, ".rds"))
# saveRDS(vi_df, viFile)
