## script to generate primary analysis dataset for Bangladesh serology paper
## brings together the (public) cleaned serosurvey data, community and grid-cell
## covariates, random forest predictions of recent infection, the daily
## probabilities of infection from the seasonality model and the age/sex
## distributions used for post-stratification
## by Andrew Azman and Stephen A Lauer, August 2019
##
## The raw lab, survey and sentinel-surveillance data are not public. The code
## that built the public inputs from them is in source/raw-data/.

options(echo=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)
RERUN_CODE <- as.logical(args[1]) # re-run the analysis?
if(is.na(RERUN_CODE)) RERUN_CODE <- FALSE
## use the random forest predictions from the published analysis? Refitting the
## random forests with current package versions reproduces ~99.5% of the
## individual predictions but not all of them, which shifts the estimates
## slightly. The random forests are refit either way (the 100- and 200-day cohort
## cross-validation predictions and the fitted forests are needed downstream).
USE_PUBLISHED_RF <- as.logical(args[2])
if(is.na(USE_PUBLISHED_RF)) USE_PUBLISHED_RF <- TRUE
age_cuts <- c(0,10,20,100)
my_labels <- paste("age", age_cuts[-length(age_cuts)], age_cuts[-1]-1, sep="_")

## load libraries and utilities
library(pacman)
p_load(here)
setwd(here())
source(here("source","utils.R"))
reload_source()
set.seed(1)

data_missing <- (is.na(find_recent_file(name_start="dat_with_raster_covs.rds",
                                        path="generated_data/")) |
                     is.na(find_recent_file(name_start="analysis_age_props.rds",
                                            path="generated_data/")))

if(data_missing | RERUN_CODE){

    ## cleaned serosurvey data (one row per participant); antibody columns are
    ## renamed to match the cohort data used to train the random forest
    labsurvey_data <- readRDS(here("data","serosurvey-dat.rds")) %>%
        rename(ptctxa_s=ctb_iga, pgctxb_s=ctb_igg,
               plpsa_s=lps_iga, pglsp_s=lps_igg) %>%
        ## vib_320 (max vibriocidal titre >= 320) is taken as provided: in the
        ## published analysis it was calculated from the originally reported titres
        mutate(age_cat=cut(age, age_cuts, my_labels, right=F),
               villid=as.factor(community_id)) %>%
        left_join(readRDS(here("data","community-data.rds")) %>%
                      select(-community_cell, -Division_name),
                  by="community_id")

    ## age and sex distribution of each grid cell (WorldPop 2015), aggregated
    ## to the age categories used in the model
    raw_age_props <- readRDS(here("data", "all_age_props.rds")) %>%
        ungroup()

    adj_age_props <- raw_age_props
    ## fill in NAs with nearby values
    age_prop_na_idx <- which(is.na(adj_age_props$age_prop))
    for(i in age_prop_na_idx){
        tmp <- adj_age_props %>%
            filter(sex == adj_age_props$sex[i], age == adj_age_props$age[i],
                   grid_id >= adj_age_props$grid_id[i]-1,
                   grid_id <= adj_age_props$grid_id[i]+1) %>%
            na.omit %>%
            summarize(subpop=mean(subpop))
        adj_age_props$subpop[i] <- tmp$subpop
    }

    age_props <- adj_age_props %>%
        group_by(grid_id) %>%
        mutate(age_prop=subpop/sum(subpop),
               age_grp=cut(as.numeric(age), age_cuts, my_labels, right=F)) %>%
        group_by(grid_id, sex, age_grp) %>%
        summarize(prop=sum(age_prop)) %>%
        ungroup()

    saveRDS(age_props, here("generated_data", "analysis_age_props.rds"))

    ## 5km x 5km grid with covariates
    my_grid <- readRDS(here("data", "grid_with_covs.rds"))

    ## add grid-cell covariates and grid centroids to the survey data
    analysis_data <- labsurvey_data %>%
        select(-grid_easting, -grid_northing) %>%
        left_join(my_grid %>% st_set_geometry(NULL), by="grid_id") %>%
        left_join(my_grid %>%
                      group_by(grid_id) %>%
                      st_centroid() %>%
                      st_coordinates %>%
                      data.frame %>%
                      mutate(grid_id = 1:n()) %>%
                      rename(grid_easting=X, grid_northing=Y),
                  by="grid_id")

    ## household coordinates in the Bangladesh Transverse Mercator projection
    ## (missing for households whose GPS location was imputed with the
    ## community centroid, imp_loc == 1)
    hhl_btm <- analysis_data %>%
        st_as_sf(coords=c("lon","lat"), crs="+proj=longlat +datum=WGS84 +no_defs") %>%
        transform_to_btm() %>%
        st_coordinates()
    analysis_data <- analysis_data %>%
        mutate(x_btm=ifelse(imp_loc == 1, NA, hhl_btm[,"X"]),
               y_btm=ifelse(imp_loc == 1, NA, hhl_btm[,"Y"]))

    ## add randomForest predictions for serosurvey data (trained on cohort data)
    ## the predictors are the columns shared by the two datasets, in this order
    rf_covs <- c("ptctxa_s", "pgctxb_s", "plpsa_s", "pglsp_s",
                 "vibinab", "vibogaw", "age", "sex")
    registerDoMC(max(1, detectCores()-2))
    cohort_dat <- readRDS("data/bgd-cohort.rds")
    ## fit randomForest at 100 days
    rf_100 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data[,rf_covs],
                                   days_back=100,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 200 days
    rf_200 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data[,rf_covs],
                                   days_back=200,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 365 days
    rf_365 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data[,rf_covs],
                                   days_back=365,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## add rf predictions to analysis data
    analysis_data$pred_100 <- rf_100$survey_preds
    analysis_data$pred_200 <- rf_200$survey_preds
    analysis_data$pred_365 <- rf_365$survey_preds

    if(USE_PUBLISHED_RF){
        published_preds <- readRDS(here("data", "published-rf-predictions.rds"))
        stopifnot(all(published_preds$sample_id == analysis_data$sample_id))
        cat("Agreement between refit and published predictions:\n")
        print(c(pred_100=mean(analysis_data$pred_100 == published_preds$pred_100, na.rm=T),
                pred_200=mean(analysis_data$pred_200 == published_preds$pred_200, na.rm=T),
                pred_365=mean(analysis_data$pred_365 == published_preds$pred_365, na.rm=T)))
        analysis_data$pred_100 <- published_preds$pred_100
        analysis_data$pred_200 <- published_preds$pred_200
        analysis_data$pred_365 <- published_preds$pred_365
        ## published 365-day cohort cross-validation predictions (used for the
        ## sensitivity and specificity); copied last so it is the most recent
        ## file that find_recent_file() picks up
        file.copy(here("data", "published-cohort-w-rf-preds-365-youden.csv"),
                  here("generated_data", "cohort-w-rf-preds-365-youden-published.csv"),
                  overwrite=TRUE)
        Sys.setFileTime(here("generated_data", "cohort-w-rf-preds-365-youden-published.csv"), Sys.time())
    }

    ## daily conditional probabilities of infection for the 100, 200 and 365
    ## days before sampling, estimated from sentinel cholera surveillance and
    ## acute watery diarrhoea data (see est_cholera_sentinel() and
    ## estimate_daily_probs() in source/raw-data/raw-data-utils.R)
    analysis_data <- analysis_data %>%
        left_join(readRDS(here("data", "daily-infection-probs.rds")),
                  by="sample_id")

    ## save main dat
    saveRDS(analysis_data, here("generated_data","dat_with_raster_covs.rds"))
}
