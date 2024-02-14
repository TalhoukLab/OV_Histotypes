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

# Import training data and class labels
data <- readRDS(file.path(inputDir, paste0(dataset, "_data.rds")))
class <- readRDS(file.path(inputDir, paste0(dataset, "_class.rds")))

# Normalization
data <- normalize(data, norm_by, norm_type, min_var)

# Genes from PrOTYPE and SPOT to keep
base_genes <- c("COL11A1", "CD74", "CD2", "TIMP3", "LUM", "CYTIP", "COL3A1",
                "THBS2", "TCF7L1", "HMGA2", "FN1", "POSTN", "COL1A2", "COL5A2",
                "PDZK1IP1", "FBN1", "HIF1A", "CXCL10", "DUSP4", "SOX17", "MITF",
                "CDKN3", "BRCA2", "CEACAM5", "ANXA4", "SERPINE1", "CRABP2", "DNAJC9")

# Import variable importance for selected algorithms and subsampling method
ranked_vi <- readRDS(file.path(outputDir, "ranked_vi", "ranked_vi_train.rds"))

# Order of candidate genes
candidate_genes <- ranked_vi %>%
  dplyr::filter(
    Sampling == bestSamp,
    Algorithm == bestAlg
  ) %>%
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
    algorithms = bestAlg,
    n = 1,
    seed_boot = as.integer(reps),
    seed_samp = 2019,
    seed_alg = 2019,
    sampling = bestSamp,
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
outputFile <- file.path(outputDir, "gene_train_eval", dataset,
                        paste0(bestAlg, "_", bestSamp, "_add", ngene, "_", reps, "_", dataset, ".rds"))
saveRDS(sm[["evals"]], outputFile)

# # Write variable importance to file
# viFile <- file.path(outputDir, "vi", dataset,
#                     paste0("vi_", alg, "_", samp, "_", reps, "_", dataset, ".rds"))
# saveRDS(vi_df, viFile)
