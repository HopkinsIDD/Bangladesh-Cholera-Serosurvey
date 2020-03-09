#!/bin/bash

#SBATCH --mem=8G
#SBATCH --time=100:00:00
#SBATCH --array=1-1000%6
#SBATCH --ntasks=4
#SBATCH --output=./generated_data/sampling_uncert_gridmap_%A_%a.o
#SBATCH --job-name=sampling_uncert_gridmap

RSCRIPT=/opt/R/3.6.1/bin/Rscript

echo "Beginning of script"
date
# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# echo "Make sure all packages are installed and data is downloaded"
# $RSCRIPT --vanilla ./source/install-packages-dl-data.R > ./generated_data/install-packages-dl-data.Rout

$RSCRIPT --vanilla ./source/sampling_uncert_gridmap.R $SLURM_ARRAY_TASK_ID > ./generated_data/sampling_uncert_gridmap.Rout

echo "End of script"
date
