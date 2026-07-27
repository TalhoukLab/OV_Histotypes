# Directories for I/O
inputDir <- "~/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"

# make all
dataset <- "train"
n_folds <- 5
fold_id <- "01"
alg <- "rf"
samp <- "smote"
rank_metric <- "f_meas"

# make tune_wflows
dataset <- "train"
n_folds <- 5
alg <- "rf"
samp <- "smote"
rank_metric <- "f_meas"

# make merge_results
dataset <- "train"
n_folds <- 5
alg <- "xgb"
samp <- "up"
rank_metric <- "f_meas"

# make retrain_merge
dataset <- "train"
n_folds <- 5
alg <- "rf"
samp <- "smote"
rank_metric <- "f_meas"

# make retrain_summarize
dataset <- "train"
n_folds <- 5
rank_metric <- "f_meas"

# make gene_opt
dataset <- "conf"
n_folds <- 5
fold_id <- "1"
ngene <- 5
gene_opt_wflow <- "smote_svm"
rank_metric <- "f_meas"

# make gene_opt_seq
dataset <- "conf_seq"
seq_data <- "seq"
n_folds <- 5
fold_id <- "1"
ngene <- 5
nseq <- 1
rank_metric <- "f_meas"
