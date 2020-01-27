#!/bin/bash

# Supervised parameters
dataSets=(cs1)
reps=500
algs=(rf adaboost combined)
samps=(none up down smote)
norm_by="none"
norm_type="conventional"
min_var=0

# Directory parameters
user="$(whoami)"
rootDir="/scratch/ovcare"/$user
inputDir=$rootDir"/Projects/histotype/data"
scriptDir=$rootDir"/results/histotype/scripts"
RDir=$scriptDir"/R"
shDir=$scriptDir"/sh"
outputDir=$rootDir"/results/histotype/outputs"
logDir=$rootDir"/results/histotype/logs"

# Bash parameters
RPath="$(which R | sort | tail -n 1)"
GREEN_TICK='\033[0;32m\xe2\x9c\x94\033[0m'
GREEN_BULLET='\033[0;32m\xe2\x80\xa2\033[0m'
BLUE_BULLET='\033[0;34m\xe2\x80\xa2\033[0m'
RED_CROSS='\033[0;31m\xe2\x9c\x96\033[0m'

# Queue parameters
maxQueueLength=8000
mem_free="4G"
mem_token="4G"
h_vmem="8G"

. ./assets/check.sh
