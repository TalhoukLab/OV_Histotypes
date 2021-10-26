#!/bin/bash

. ./assets/params.sh

file_to_submit=()

# Make directories for R script, shell script
subDir=sequential/summary
RSubDir=$RDir/$subDir
shSubDir=$shDir/$subDir

# Make job and output directorie
mkdir -p $RSubDir
mkdir -p $shSubDir
mkdir -p $outputDir/$subDir

# Content of R file
R_file=$RSubDir/sequential_iv_summary.R
echo 'outputDir <- "'$outputDir'"' > $R_file
echo 'source("pipeline/sequential/3-seq_iv_summary.R")' >> $R_file

# Content of sh file
job_file=$shSubDir/sequential_iv_summary.sh
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

# Submit to queue if sbatch exists
logDir=$logDir/$subDir
outputDir=$outputDir/$subDir
if command -v sbatch &>/dev/null; then
    . ./assets/submit_slurm.sh
fi
