#!/bin/bash

# Supervised parameters
dataSets=(train)
#reps=50
n_folds=10
algs=(rf xgb svm mr)
samps=(none up down smote hybrid)
rank_metric="accuracy"
#norm_by="none"
#norm_type="conventional"
#min_var=0
seqData=(seq)
nseq=4
bestAlg="adaboost"
bestSamp="up"
ngenes=44

# Directory parameters
user="$(whoami)"
rootDir="/home"/$user
projDir=$rootDir"/Projects/OV_Histotypes"
inputDir=$projDir"/data"
resultsDir=$rootDir"/results/OV_Histotypes"
scriptDir=$resultsDir"/scripts"
RDir=$scriptDir"/R"
shDir=$scriptDir"/sh"
outputDir=$resultsDir"/outputs"
logDir=$resultsDir"/logs"

# Bash parameters
RPath="$(which R | sort | tail -n 1)"
GREEN_TICK='\033[0;32m\xe2\x9c\x94\033[0m'
GREEN_BULLET='\033[0;32m\xe2\x80\xa2\033[0m'
BLUE_BULLET='\033[0;34m\xe2\x80\xa2\033[0m'
RED_CROSS='\033[0;31m\xe2\x9c\x96\033[0m'

# Queue parameters
maxQueueLength=8000
shouldWait=TRUE

. ./assets/check.sh
