# Load packages, data, and common arguments
library(microbenchmark)
library(purrr)
library(splendid)
library(here)

cs1_data <- readRDS(here("data/cs1_all_data.rds"))
cs1_class <- readRDS(here("data/cs1_all_class.rds"))

cs2_data <- readRDS(here("data/cs2_all_data.rds"))
cs2_class <- readRDS(here("data/cs2_all_class.rds"))

sm_args <- list(
  n = 1,
  seed_boot = 1,
  seed_samp = 2019,
  seed_alg = 2019,
  stratify = TRUE
)

# CS1 ---------------------------------------------------------------------

# Test grouping structure
mb1 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "svm"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "rf"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "adaboost"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_lasso"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_ridge"),
  times = 2
)

# Test duration for each combination
mb_none1 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = c("svm", "rf", "adaboost"), sampling = "none"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_lasso", sampling = "none"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_ridge", sampling = "none"),
  times = 1
)

mb_up1 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = c("svm", "rf", "adaboost"), sampling = "up"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_lasso", sampling = "up"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_ridge", sampling = "up"),
  times = 1
)

mb_down1 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = c("svm", "rf", "adaboost"), sampling = "down"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_lasso", sampling = "down"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_ridge", sampling = "down"),
  times = 1
)

mb_smote1 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = c("svm", "rf", "adaboost"), sampling = "smote"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_lasso", sampling = "smote"),
  invoke(splendid_model, sm_args, data = cs1_data, class = cs1_class, alg = "mlr_ridge", sampling = "smote"),
  times = 1
)

# Time for one iteration of all algorithms and subsampling methods in CS1
mb_total1 <- list(mb_none1, mb_up1, mb_down1, mb_smote1) %>%
  purrr::map("time") %>%
  unlist() %>%
  magrittr::divide_by(1e9) %>%
  sum()

# 500 reps = ~ 16.5 hours in serial
mb_total1 * 500 / 60 / 60


# CS2 ---------------------------------------------------------------------

# Test grouping structure
mb2 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "svm"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "rf"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "adaboost"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_lasso"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_ridge"),
  times = 2
)

# Test duration for each combination
mb_none2 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = c("svm", "rf", "adaboost"), sampling = "none"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_lasso", sampling = "none"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_ridge", sampling = "none"),
  times = 1
)

mb_up2 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = c("svm", "rf", "adaboost"), sampling = "up"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_lasso", sampling = "up"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_ridge", sampling = "up"),
  times = 1
)

mb_down2 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = c("svm", "rf", "adaboost"), sampling = "down"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_lasso", sampling = "down"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_ridge", sampling = "down"),
  times = 1
)

mb_smote2 <- microbenchmark(
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = c("svm", "rf", "adaboost"), sampling = "smote"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_lasso", sampling = "smote"),
  invoke(splendid_model, sm_args, data = cs2_data, class = cs2_class, alg = "mlr_ridge", sampling = "smote"),
  times = 1
)

# Time for one iteration of all algorithms and subsampling methods in CS2
mb_total2 <- list(mb_none2, mb_up2, mb_down2, mb_smote2) %>%
  purrr::map("time") %>%
  unlist() %>%
  magrittr::divide_by(1e9) %>%
  sum()

# 500 reps = ~ 71 hours in serial
mb_total2 * 500 / 60 / 60
