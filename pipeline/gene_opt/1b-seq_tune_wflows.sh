#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=gene_opt/sequential/tune_wflows
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${seqData[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    for nseq in $(seq 1 $nseq); do
        for v in $(seq -f "%0${#n_folds}g" 1 $n_folds); do
            for ng in $(seq 0 $ngenes); do
                # Content of R file
                R_file=$RSubDir/$dataset/$dataset"_s"$nseq"_"$v"_add"$ng.R
                echo 'dataset <- "'$dataset'"' > $R_file
                echo "n_folds <- $n_folds" >> $R_file
                echo 'fold_id <- "'$v'"' >> $R_file
                echo "ngene <- '$ng'" >> $R_file
                echo 'nseq <- '$nseq >> $R_file
                echo 'inputDir <- "'$inputDir'"' >> $R_file
                echo 'outputDir <- "'$outputDir'"' >> $R_file
                echo 'source("pipeline/gene_opt/1b-seq_tune_wflows.R")' >> $R_file

                # Content of sh file
                job_file=$shSubDir/$dataset/$dataset"_s"$nseq"_"$v"_add"$ng.sh
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
