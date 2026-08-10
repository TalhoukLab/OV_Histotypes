# Directories for I/O
inputDir <- "~/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"

## make train
# make tune_wflows
dataset <- "train"
n_folds <- 5
fold_id <- "01"
alg <- "rf"
samp <- "smote"

# make merge_results
dataset <- "train"
n_folds <- 5
alg <- "xgb"
samp <- "up"
rank_metric <- "f_meas"

# make summarize_results
dataset <- "train"
n_folds <- 5

## make retrain
# make retrain_tune
dataset <- "train"
n_folds <- 5
fold_id <- "01"
alg <- "rf"
samp <- "smote"

# make retrain_merge
dataset <- "train"
n_folds <- 5
alg <- "xgb"
samp <- "up"
rank_metric <- "f_meas"

# make retrain_summarize
dataset <- "train"
n_folds <- 5
rank_metric <- "f_meas"

## make sequential
# make seq_tune
dataset <- "seq"
nseq <- 4
n_folds <- 5
fold_id <- "01"

# make seq_merge
dataset <- "seq"
nseq <- 4
n_folds <- 5
rank_metric <- "f_meas"

# make seq_summarize
dataset <- "seq"
n_folds <- 5

## make gene_opt
# make gene_tune
dataset <- "conf"
n_folds <- 5
fold_id <- "01"
ngene <- 5
gene_opt_wflow <- "smote_svm"
gene_opt_rank <- "vi"

# make gene_merge
dataset <- "conf"
n_folds <- 5
ngene <- 5
rank_metric <- "f_meas"
gene_opt_wflow <- "smote_svm"
gene_opt_rank <- "vi"
ngenes <- 44

# make gene_summarize
dataset <- "conf"
nfolds <- 5

## make gene_opt_seq
# make gene_seq_tune
dataset <- "conf_seq"
n_folds <- 5
fold_id <- "1"
ngene <- 5
nseq <- 1
seq_data <- "seq"

# make gene_seq_merge
dataset <- "conf_seq"
n_folds <- 5
ngene <- 5
nseq <- 1
rank_metric <- "f_meas"
seq_data <- "seq"

# make gene_seq_summarize
dataset <- "conf_seq"
seq_data <- "seq"
nfolds <- 5
