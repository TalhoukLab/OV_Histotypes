#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=gene_opt/sequential/merge_eval
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for sq in "${seqData[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$sq
    mkdir -p $shSubDir/$sq
    mkdir -p $outputDir/$subDir

    # Content of R file
    R_file=$RSubDir/$sq/gene_opt/sequential/merge_eval.R
    echo 'outputDir <- "'$outputDir'"' > $R_file
    echo 'sq <- "'$sq'"' >> $R_file
    echo "nseq <- '$nseq'" >> $R_file
    echo "reps <- '$reps'" >> $R_file
    echo "ngenes <- '$ngenes'" >> $R_file
    echo 'source("pipeline/gene_opt/2-gene_merge_eval.R")' >> $R_file

    # Content of sh file
    job_file=$shSubDir/$sq/gene_opt/sequential/merge_eval.sh
    cat ./assets/sbatch_params.sh > $job_file
    echo "cd $projDir" >> $job_file
    echo "Rscript $R_file" >> $job_file
    chmod +x $job_file

    # Add to queue if sbatch exists
    if command -v sbatch &>/dev/null; then
        file_to_submit+=($job_file)
        echo -e "$GREEN_TICK Added to queue: $job_file"
    else
        bash $shell_file
    fi
done

# Submit to queue if sbatch exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v sbatch &>/dev/null; then
    . ./assets/submit_slurm.sh
fi
