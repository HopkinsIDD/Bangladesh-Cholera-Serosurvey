#!/bin/bash

#SBATCH --mem=8G
#SBATCH --time=100:00:00
#SBATCH --array=1-1000%6
#SBATCH --ntasks=4
#SBATCH --output=./generated_data/sampling_uncert_gridmap_%A_%a.o
#SBATCH --job-name=sampling_uncert_gridmap

## usage: sbatch inst/run-boot-gridmap.sh main
##        sbatch inst/run-boot-gridmap.sh season

RSCRIPT=${RSCRIPT:-Rscript} # override with the path to your cluster's Rscript if needed
VARIANT=${1:-main}

echo "Beginning of script"
date
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

$RSCRIPT --vanilla ./source/sampling_uncert_gridmap.R $SLURM_ARRAY_TASK_ID $VARIANT > ./generated_data/sampling_uncert_gridmap_${VARIANT}_$SLURM_ARRAY_TASK_ID.Rout 2>&1

echo "End of script"
date
