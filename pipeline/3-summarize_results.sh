#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=summarize_results
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${dataSets[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    # Content of R file
    R_file=$RSubDir/$dataset/summarize_results.R
    echo 'outputDir <- "'$outputDir'"' > $R_file
    echo 'dataset <- "'$dataset'"' >> $R_file
    echo 'source("pipeline/3-summarize_results.R")' >> $R_file

    # Content of sh file
    job_file=$shSubDir/$dataset/summarize_results.sh
    cat ./assets/sbatch_params.sh > $job_file
    echo "cd $projDir" >> $job_file
    echo "Rscript $R_file" >> $job_file
    chmod +x $job_file

    # Add to queue if sbatch exists
    if command -v sbatch &>/dev/null; then
       file_to_submit+=($job_file)
       echo -e "$GREEN_TICK Added to queue: $job_file"
    else
       bash $job_file
    fi
done

# Submit to queue if sbatch exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v sbatch &>/dev/null; then
    . ./assets/submit_slurm.sh
fi
