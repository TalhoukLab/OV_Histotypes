#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=retrain/tune_wflows
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${dataSets[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    for v in $(seq -f "%0${#n_folds}g" 1 $n_folds); do
        for alg in "${algs[@]}"; do
            for samp in "${samps[@]}"; do
                # Content of R file
                R_file=$RSubDir/$dataset/$samp"_"$alg"_"$v.R
                echo 'dataset <- "'$dataset'"' > $R_file
                echo "n_folds <- $n_folds" >> $R_file
                echo 'fold_id <- "'$v'"' >> $R_file
                echo 'alg <- "'$alg'"' >> $R_file
                echo 'samp <- "'$samp'"' >> $R_file
                echo 'inputDir <- "'$inputDir'"' >> $R_file
                echo 'outputDir <- "'$outputDir'"' >> $R_file
                echo 'source("pipeline/retrain/1-tune_wflows.R")' >> $R_file

                # Content of sh file
                job_file=$shSubDir/$dataset/$samp"_"$alg"_"$v.sh
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
    done
done

# Submit to queue if sbatch exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v sbatch &>/dev/null; then
    . ./assets/submit_slurm.sh
fi
