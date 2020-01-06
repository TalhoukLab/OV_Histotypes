#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=summary
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${dataSets[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    # Content of R file
    R_file=$RSubDir/$dataset/iv_summary.R
    echo 'outputDir <- "'$outputDir'"' > $R_file
    echo 'dataset <- "'$dataset'"' >> $R_file
    echo 'source("pipeline/3-iv_summary.R")' >> $R_file

    # Content of sh file
    sh_file=$shSubDir/$dataset/iv_summary.sh
    echo "Rscript $R_file" > $sh_file
    chmod +x $sh_file

    # Add to queue if qsub exists
    if command -v qsub &>/dev/null; then
       file_to_submit+=($sh_file)
       echo -e "$GREEN_TICK Added to queue: $sh_file"
    else
       bash $sh_file
    fi
done

# Submit to queue if qsub exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v qsub &>/dev/null; then
    . ./assets/submit_queue.sh
fi
