#!/bin/bash

#SBATCH --job-name=histotypes     # Job name
#SBATCH --partition=upgrade       # Partition requested (queue)
#SBATCH --nodes=1                 # Run all processes on a single node
#SBATCH --ntasks=1                # Run a single task
#SBATCH --cpus-per-task=8         # Number of CPU cores per task
#SBATCH --mem-per-cpu=8GB         # Memory per core
#SBATCH --output=o_%j.out         # Standard output
#SBATCH --error=e_%j.out          # Standard error
