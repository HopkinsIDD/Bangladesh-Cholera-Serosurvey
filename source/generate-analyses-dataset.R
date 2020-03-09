## script to generate primary analysis dataset for Bangladesh serology paper
## includes bringing in lab data, the survey data and the raster data for mapping
## by Andrew Azman and Stephen A Lauer, August 2019

options(echo=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)
RERUN_CODE <- as.logical(args[1]) # re-run the analysis?

## load libraries and utilities
library(pacman)
p_load(here)
setwd(here())
source("source/utils.R")
reload_source()
set.seed(1)

data_missing <- (is.na(
    find_recent_file(name_start="serosurvey-with-rf-preds.rds",
                     path="generated_data/")))

if(data_missing | RERUN_CODE){
    ## read in serosurvey data
    serosurvey_dat <- readRDS("data/serosurvey-dat.rds")

    ## add randomForest predictions for serosurvey data (trained on cohort data)
    registerDoMC(detectCores()-2)
    cohort_dat <- readRDS("data/bgd-cohort.rds")
    ## fit randomForest at 100 days
    rf_100 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=serosurvey_dat,
                                   days_back=100,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 200 days
    rf_200 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=serosurvey_dat,
                                   days_back=200,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 365 days
    rf_365 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=serosurvey_dat,
                                   days_back=365,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## add rf predictions to analysis data
    serosurvey_dat$pred_100 <- rf_100$survey_preds
    serosurvey_dat$pred_200 <- rf_200$survey_preds
    serosurvey_dat$pred_365 <- rf_365$survey_preds

    ## save main dat
    saveRDS(serosurvey_dat,"generated_data/serosurvey-with-rf-preds.rds")
}
