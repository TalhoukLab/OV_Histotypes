#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=merge_eval
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${dataSets[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir

    # Content of R file
    R_file=$RSubDir/$dataset/merge_eval.R
    echo 'outputDir <- "'$outputDir'"' > $R_file
    echo 'dataset <- "'$dataset'"' >> $R_file
    echo 'algs <- strsplit("'${algs[@]}'", " ")[[1]]' >> $R_file
    echo 'samps <- strsplit("'${samps[@]}'", " ")[[1]]' >> $R_file
    echo 'source("2-merge_eval.R")' >> $R_file

    # Content of sh file
    sh_file=$shSubDir/$dataset/merge_eval.sh
    echo "Rscript $R_file" > $sh_file
    chmod +x $sh_file

    # Add to queue if qsub exists
    if command -v qsub &>/dev/null; then
        file_to_submit+=($sh_file)
        echo -e "$GREEN_TICK Added to queue: $sh_file"
    else
        bash $shell_file
    fi
done

# Submit to queue if qsub exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v qsub &>/dev/null; then
    . ./assets/submit_queue.sh
fi
