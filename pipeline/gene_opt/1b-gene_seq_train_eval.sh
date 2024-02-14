#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=gene_opt/sequential/train_eval
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

# Make job and output directories
mkdir -p $RSubDir
mkdir -p $shSubDir
mkdir -p $outputDir/$subDir

for sq in "${seqData[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$sq
    mkdir -p $shSubDir/$sq
    mkdir -p $outputDir/$subDir/$sq
    mkdir -p $outputDir/sequential/vi/$sq

    for s in $(seq -f "%0${#reps}g" 1 $reps); do
        for nseq in $(seq 1 $nseq); do
            for ng in $(seq 1 $ngenes); do
                # Content of R file
                R_file=$RSubDir/$sq"_"$nseq"_add"$ng"_"$s.R
                echo 'reps <- "'$s'"' > $R_file
                echo "ngene <- '$ng'" >> $R_file
                echo 'sq <- "'$sq'"' >> $R_file
                echo 'nseq <- '$nseq >> $R_file
                echo 'inputDir <- "'$inputDir'"' >> $R_file
                echo 'outputDir <- "'$outputDir'"' >> $R_file
                echo 'norm_by <- "'$norm_by'"' >> $R_file
                echo 'norm_type <- "'$norm_type'"' >> $R_file
                echo "min_var <- '$min_var'" >> $R_file
                echo 'source("pipeline/sequential/1b-gene_seq_train_eval.R")' >> $R_file

                # Content of sh file
                job_file=$shSubDir/$sq"_"$nseq"_add"$ng"_"$s.sh
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
