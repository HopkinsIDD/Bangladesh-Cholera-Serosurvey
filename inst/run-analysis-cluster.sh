#!/bin/bash

#SBATCH --mem=30G
#SBATCH --time=100:00:00
#SBATCH --ntasks=5
#SBATCH --output=./generated_data/bgd-sero-est.o
#SBATCH --job-name=bgd-sero-est

## Steps 1-3 of inst/run-analysis.sh on a SLURM cluster. Then submit
## inst/run-boot-gridmap.sh (twice: main and season) and run steps 5-8 of
## inst/run-analysis.sh.

RSCRIPT=${RSCRIPT:-Rscript} # override with the path to your cluster's Rscript if needed

echo "Beginning of script"
date

RERUN_ANALYSIS=T
CHAINS=4

mkdir -p generated_data

echo "Make sure all packages are installed and data is downloaded"
$RSCRIPT --vanilla ./source/install-packages-dl-data.R > ./generated_data/install-packages-dl-data.Rout 2>&1

echo "Make analysis dataset with random forest predictions"
$RSCRIPT --vanilla ./source/generate-analyses-dataset.R $RERUN_ANALYSIS > ./generated_data/generate-analyses-dataset.Rout 2>&1

for DAYS in 365 200 100; do
    for DIAG in rf vib; do
        echo "Estimate cholera seroincidence in Bangladesh ($DAYS-day $DIAG)"
        $RSCRIPT --vanilla ./source/bgd-estimates-stan.R $DAYS $DIAG youden $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-$DIAG-$DAYS.Rout 2>&1
    done
done

echo "End of script"
date
