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

# Fast algorithms combined
if (alg == "combined") {
  a <- c("svm", "lda", "xgboost", "knn")
} else {
  a <- alg
}

# Supervised learning model output
sm <- splendid::splendid_model(
  data = data,
  class = class,
  algorithms = a,
  n = 1,
  seed_boot = as.integer(reps),
  seed_samp = 2019,
  seed_alg = 2019,
  sampling = samp,
  stratify = TRUE
)

# Write evaluations to file
outputFile <- file.path(outputDir, "train_eval", dataset,
                        paste0(alg, "_", samp, "_", reps, "_", dataset, ".rds"))
saveRDS(sm[["evals"]], outputFile)
