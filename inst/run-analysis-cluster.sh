#!/bin/bash

#SBATCH --mem=30G
#SBATCH --time=100:00:00
#SBATCH --ntasks=5
#SBATCH --output=./generated_data/bgd-sero-est.o
#SBATCH --job-name=bgd-sero-est

RSCRIPT=/opt/R/3.6.1/bin/Rscript

echo "Beginning of script"
date

RERUN_ANALYSIS=F
CHAINS=5

echo "Make sure all packages are installed and data is downloaded"
$RSCRIPT --vanilla ./source/install-packages-dl-data.R > ./generated_data/install-packages-dl-data.Rout

echo "Make analysis dataset with raster covariates, probability curves, and model predictions"
$RSCRIPT --vanilla ./source/generate-analyses-dataset.R $RERUN_ANALYSIS > ./generated_data/generate-analyses-dataset.Rout

echo "Estimate cholera incidence in Bangladesh (365-day random forest)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 365 "rf" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-rf-365.Rout

echo "Estimate cholera incidence in Bangladesh (365-day vibriocidal)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 365 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-365.Rout

echo "Estimate cholera incidence in Bangladesh (200-day random forest)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 200 "rf" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-rf-200.Rout

echo "Estimate cholera incidence in Bangladesh (200-day vibriocidal)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 200 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-200.Rout

echo "Estimate cholera incidence in Bangladesh (100-day random forest)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 100 "rf" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-rf-100.Rout

echo "Estimate cholera incidence in Bangladesh (100-day vibriocidal)"
$RSCRIPT --vanilla ./source/bgd-estimates-stan.R 100 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-100.Rout

echo "End of script"
date
