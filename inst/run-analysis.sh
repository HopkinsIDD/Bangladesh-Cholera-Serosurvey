#!/bin/bash
## Run the full analysis pipeline on a single machine, from the repository root:
##   bash inst/run-analysis.sh
##
## Each step writes to generated_data/ (logs: generated_data/*.Rout) and skips
## work whose outputs already exist unless RERUN_ANALYSIS=T. See README.md for
## run times, memory needs and the order of dependencies.
## For a SLURM cluster see inst/run-analysis-cluster.sh and inst/run-boot-gridmap.sh.

set -e

echo "Beginning of script"
date

RERUN_ANALYSIS=T   # T = refit the random forests even if generated_data/dat_with_raster_covs.rds exists
CHAINS=4           # Stan chains (1,250 post-warmup draws each)
N_BOOT=1000        # number of INLA maps (posterior draws) per map type
BOOT_JOBS=6        # INLA maps to fit in parallel

mkdir -p generated_data/boot-maps generated_data/boot-maps-season figures tables

echo "1. Install packages and download the cohort data"
Rscript --vanilla ./source/install-packages-dl-data.R > ./generated_data/install-packages-dl-data.Rout 2>&1

echo "2. Build analysis dataset with random forest predictions"
Rscript --vanilla ./source/generate-analyses-dataset.R $RERUN_ANALYSIS > ./generated_data/generate-analyses-dataset.Rout 2>&1

echo "3. Estimate seroincidence with Stan (random forest and vibriocidal; 365, 200, 100 days)"
for DAYS in 365 200 100; do
    for DIAG in rf vib; do
        echo "   $DIAG, $DAYS days"
        Rscript --vanilla ./source/bgd-estimates-stan.R $DAYS $DIAG youden $CHAINS 1 \
            > ./generated_data/bgd-adjusted-estimate-$DIAG-$DAYS.Rout 2>&1
    done
done

echo "4. Fit $N_BOOT INLA grid-cell maps for the main and seasonality estimates ($BOOT_JOBS at a time; existing maps are skipped)"
for VARIANT in main season; do
    seq 1 $N_BOOT | xargs -P $BOOT_JOBS -I{} sh -c \
        "Rscript --vanilla ./source/sampling_uncert_gridmap.R {} $VARIANT >> ./generated_data/sampling_uncert_gridmap_${VARIANT}.Rout 2>&1"
done

echo "5. Leave-one-community-out cross validation of the spatial model"
Rscript --vanilla ./source/loocv-analysis.R > ./generated_data/loocv-analysis.Rout 2>&1

echo "6. Risk factor models (requires pandoc and LaTeX)"
Rscript --vanilla -e "rmarkdown::render('source/risk_factors_and_model_comp.Rmd')"

echo "7. Supplement; also summarises the maps into generated_data/grid-cell-data*.rds, INLA-seroincidence-est*.rds and bgd2-seroincidence*.rds (requires pandoc and LaTeX)"
Rscript --vanilla -e "rmarkdown::render('source/supplement.Rmd')"

echo "8. Main-text estimates, Figures 1-2 (requires pandoc and LaTeX)"
Rscript --vanilla -e "rmarkdown::render('source/manuscript-estimates.Rmd')"

echo "End of script"
date
