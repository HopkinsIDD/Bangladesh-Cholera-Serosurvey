#!/bin/bash

#SBATCH --mem=30G
#SBATCH --time=100:00:00
#SBATCH --ntasks=10
#SBATCH --output=./generated_data/bgd-sero-est.o
#SBATCH --job-name=bgd-sero-est

echo "Beginning of script"
date

RERUN_ANALYSIS=F
CHAINS=5

echo "Make sure all packages are installed and data is downloaded"
Rscript --vanilla ./source/install-packages-dl-data.R > ./generated_data/install-packages-dl-data.Rout

echo "Make analysis dataset with raster covariates, probability curves, and model predictions"
Rscript --vanilla ./source/generate-analyses-dataset.R $RERUN_ANALYSIS > ./generated_data/generate-analyses-dataset.Rout

echo "Estimate cholera incidence in Bangladesh (365-day random forest)"
Rscript --vanilla ./source/bgd-estimates-stan.R 365 "rf" "youden" $CHAINS 2 > ./generated_data/bgd-adjusted-estimate-rf-365.Rout

echo "Estimate cholera incidence in Bangladesh (365-day vibriocidal)"
Rscript --vanilla ./source/bgd-estimates-stan.R 365 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-365.Rout

echo "Estimate cholera incidence in Bangladesh (200-day random forest)"
Rscript --vanilla ./source/bgd-estimates-stan.R 200 "rf" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-rf-200.Rout

echo "Estimate cholera incidence in Bangladesh (200-day vibriocidal)"
Rscript --vanilla ./source/bgd-estimates-stan.R 200 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-200.Rout

echo "Estimate cholera incidence in Bangladesh (100-day random forest)"
Rscript --vanilla ./source/bgd-estimates-stan.R 100 "rf" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-rf-100.Rout

echo "Estimate cholera incidence in Bangladesh (100-day vibriocidal)"
Rscript --vanilla ./source/bgd-estimates-stan.R 100 "vib" "youden" $CHAINS 1 > ./generated_data/bgd-adjusted-estimate-vib-100.Rout

echo "Estimate risk factor effect sizes (this requires pandoc to run)"
Rscript --vanilla -e "rmarkdown::render('source/risk_factors_and_model_comp.Rmd')"

echo "End of script"
date
