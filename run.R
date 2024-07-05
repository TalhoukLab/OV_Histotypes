# make all
inputDir <- "~/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"
dataset <- "train"
n_folds <- 3
fold_id <- "01"
alg <- "rf"
samp <- "none"
rank_metric <- "f_meas"

# make retrain
inputDir <- "~/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"
dataset <- "retrain_4"
n_folds <- 5
fold_id <- "01"
alg <- "rf"

# make sequential gene_opt
inputDir <- "~/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"
dataset <- "seq"
n_folds <- 5
fold_id <- "1"
ngene <- 5
nseq <- 2
rank_metric <- "f_meas"

dataset <- "seq"
nseq <- 2
n_folds <- 5
fold_id <- "1"
inputDir <- "/home/dchiu/Projects/OV_Histotypes/data"
outputDir <- "/projects/dchiu_prj/results/OV_Histotypes/outputs"
