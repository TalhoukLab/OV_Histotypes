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
sq_class <- readRDS(file.path(inputDir, paste0(sq, "_classs.rds")))
sq_algs <- readRDS(file.path(inputDir, paste0(sq, "_algs.rds")))

data <- sq_data[[nseq]]
class <- sq_class[[nseq]]
alg <- sq_algs[[nseq]]

# Normalization
data <- normalize(data, norm_by, norm_type, min_var)

# Best sampling method from classification of full training set
seq_top <- readRDS(file.path(inputDir, "seq_top_c5.rds"))
samp <- as.character(seq_top[["sampling"]])

# Supervised learning model output
sm <- splendid::splendid_model(
  data = data,
  class = class,
  algorithms = alg,
  n = 1,
  seed_boot = as.integer(reps),
  seed_samp = 2019,
  seed_alg = 2019,
  sampling = samp,
  stratify = TRUE
)

# Write evaluations to file
outputFile <- file.path(outputDir, "sequential", "train_eval",
                        paste0(sq, nseq, "_", reps, ".rds"))
saveRDS(sm[["evals"]], outputFile)
