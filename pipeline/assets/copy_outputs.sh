#!/bin/bash

# Directory paths
user="$(whoami)"
projDir="/home"/$user"/Projects/OV_Histotypes"
inputDir=$projDir/"data"
resultsDir="/projects/"$user"_prj/results/OV_Histotypes"
outputDir=$resultsDir"/outputs"

# Training pipeline outputs
cp $outputDir/merge_results/train/wflow_up_rf_model_train.rds $inputDir/
cp $outputDir/merge_results/train/wflow_smote_svm_model_train.rds $inputDir/
cp $outputDir/merge_results/train/wflow_hybrid_xgb_model_train.rds $inputDir/
cp $outputDir/summarize_results/train/* $inputDir/

# Retraining pipeline outputs
cp $outputDir/retrain/summarize_results/retrain_4/*metrics* $inputDir/
cp $outputDir/retrain/summarize_results/retrain_3/*metrics* $inputDir/
cp $outputDir/retrain/summarize_results/retrain_2/*metrics* $inputDir/

# Sequential pipeline outputs
cp $outputDir/sequential/summarize_results/seq/* $inputDir/
cp $outputDir/sequential/summarize_results/two_step/* $inputDir/

# Gene optimization pipeline outputs
cp $outputDir/gene_opt/summarize_results/conf/* $inputDir/

# Sequential gene optimization pipeline outputs
cp $outputDir/gene_opt/sequential/summarize_results/conf_seq/* $inputDir/
cp $outputDir/gene_opt/sequential/summarize_results/conf_two_step/* $inputDir/
