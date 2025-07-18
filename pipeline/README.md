# Running Analyses

## Training

To run the ML pipeline on the **training set**, execute `make all` with parameters

``` bash
dataSets=(train)
```

## Retraining

To find the sequence of workflows to use in the sequential and two-step methods, we need to retrain the ML pipeline after each class is removed. This allows us to obtain an unbiased estimate of the internal validation score to use when assessing the next workflow to use.

Run `make retrain` with these parameters, one at a time:

``` bash
dataSets=(train)
dataSets=(retrain_4)
dataSets=(retrain_3)
dataSets=(retrain_2)
```

Do not supply an array to `dataSets` as the results of a later retraining step depends on the previous.

Run L14-16 of `copy_outputs.sh` to copy retraining outputs to this project's data directory.

## Sequential and Two-Step

First run these two scripts to generate the input data, classes, and workflows for the sequential and two-step methods:

``` r
source("pipeline/sequential/0-2S_combine.R")
source("pipeline/sequential/0-seq_combine.R")
```

Then run `make sequential` with these parameters:

``` bash
seqData=(seq)
nseq=4
```

``` bash
seqData=(two_step)
nseq=2
```

Run L19-20 of `copy_outputs.sh` to copy sequential outputs to this project's data directory.

## Gene Optimization

Gene optimization is performed on the confirmation set. First run this script to generate the input data, classes, and workflows for the confirmation set:

``` r
source("pipeline/gene_opt/0-conf_combine.R")
```

Then run `make gene_opt` with parameters:

``` bash
dataSets=(conf)
gene_opt_wflow="up_xgb"
```

to run gene optimization on the entire confirmation set.

Run `make gene_opt_seq` with parameters:

``` bash
dataSets=(conf_seq)
seqData="seq"
nseq=4
```

``` bash
dataSets=(conf_two_step)
seqData="two_step"
nseq=2
```

to run gene optimization for the sequential and two-step methods on the confirmation set.

# Obtaining Results

Results are stored at `resultsDir` specified in `pipeline/assets/params.sh`. It has three main directories:

-   `scriptDir`: contains the R and sh scripts used to submit jobs to SLURM
-   `logDir`: contains logs for standard output (extension `.o`) and standard error (extension `.e`) associated with the jobs, if any
-   `outputDir`: contains the outputs for the submitted jobs

The outputs in `outputDir` are structured first within the name of the task and within that subdirectories will use the data object name.
