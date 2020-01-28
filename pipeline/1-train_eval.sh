#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=train_eval
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

for dataset in "${dataSets[@]}"; do
    # Make job and output directories for dataset
    mkdir -p $RSubDir/$dataset
    mkdir -p $shSubDir/$dataset
    mkdir -p $outputDir/$subDir/$dataset

    for s in $(seq -f "%0${#reps}g" 1 $reps); do
        for alg in "${algs[@]}"; do
            for samp in "${samps[@]}"; do
                # Content of R file
                R_file=$RSubDir/$dataset/$alg"_"$samp"_"$s.R
                echo 'dataset <- "'$dataset'"' > $R_file
                echo 'reps <- "'$s'"' >> $R_file
                echo 'alg <- "'$alg'"' >> $R_file
                echo 'samp <- "'$samp'"' >> $R_file
                echo 'inputDir <- "'$inputDir'"' >> $R_file
                echo 'outputDir <- "'$outputDir'"' >> $R_file
                echo 'norm_by <- "'$norm_by'"' >> $R_file
                echo 'norm_type <- "'$norm_type'"' >> $R_file
                echo "min_var <- '$min_var'" >> $R_file
                echo 'source("1-train_eval.R")' >> $R_file

                # Content of sh file
                sh_file=$shSubDir/$dataset/$alg"_"$samp"_"$s.sh
                echo "Rscript $R_file" > $sh_file
                chmod +x $sh_file

                # Add to queue if qsub exists
                if command -v qsub &>/dev/null; then
                    file_to_submit+=($sh_file)
                    echo -e "$GREEN_TICK Added to queue: $sh_file"
                fi
            done
        done
    done
done

# Submit to queue if qsub exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v qsub &>/dev/null; then
    . ./assets/submit_queue.sh
fi
