#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=sequential/merge_results
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${seqData[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    for nseq in $(seq 1 $nseq); do
        # Content of R file
        R_file=$RSubDir/$dataset/"merge_"$dataset"_"$nseq.R
        echo 'dataset <- "'$dataset'"' > $R_file
        echo 'nseq <- '$nseq >> $R_file
        echo "n_folds <- $n_folds" >> $R_file
        echo 'rank_metric <- "'$rank_metric'"' >> $R_file
        echo 'inputDir <- "'$inputDir'"' >> $R_file
        echo 'outputDir <- "'$outputDir'"' >> $R_file
        echo 'source("pipeline/sequential/2-merge_results.R")' >> $R_file

        # Content of sh file
        job_file=$shSubDir/$dataset/"merge_"$dataset"_"$nseq.sh
        cat ./assets/sbatch_params.sh > $job_file
        echo "cd $projDir" >> $job_file
        echo "Rscript $R_file" >> $job_file
        chmod +x $job_file

        # Add to queue if sbatch exists
        if command -v sbatch &>/dev/null; then
            file_to_submit+=($job_file)
            echo -e "$GREEN_TICK Added to queue: $job_file"
        fi
    done
done

# Submit to queue if sbatch exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v sbatch &>/dev/null; then
    . ./assets/submit_slurm.sh
fi
