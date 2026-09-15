## script to generate primary analysis dataset for Bangladesh serology paper
## includes bringing in lab data, the survey data and the raster data for mapping
## by Andrew Azman and Stephen A Lauer, August 2019
## NOTE: this script cannot be run from this repository. It documents how the
## public inputs in data/ (serosurvey-dat.rds, daily-infection-probs.rds and
## community-data.rds) were built from the raw lab, questionnaire and sentinel
## surveillance data, which are not public (contact azman@jhu.edu).
## generate-raster-data.R similarly documents how data/grid_with_covs.rds and
## data/all_age_props.rds were built from WorldPop, Malaria Atlas Project and
## other public rasters.

options(echo=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)
RERUN_CODE <- as.logical(args[1]) # re-run the analysis?
age_cuts <- c(0,10,20,100)
my_labels <- paste("age", age_cuts[-length(age_cuts)], age_cuts[-1]-1, sep="_")

## load libraries and utilities
library(pacman)
p_load(here)
setwd(here())
source(here("source","utils.R"))
reload_source()
p_load(readxl, tidyxl, raster, exactextractr, fasterize, brms, scam)
source(here("source","raw-data","raw-data-utils.R"))
set.seed(1)

data_missing <- (is.na(find_recent_file(name_start="dat_with_raster_covs.rds",
                                 path="generated_data/")) |
                     is.na(find_recent_file(name_start="analysis_age_props.rds",
                                            path="generated_data/")))

if(data_missing | RERUN_CODE){

    ######################################################################################################
    ## brings in lab data and merges with survey data
    ##
    ## then adds pred model data and continous vibriocidals (experimental and not used in main analyses)
    ## does a little recoding
    ## note this takes a little while                                                     ##
    ## saves output to here("generated_data","full_data.rds")
    ##
    ######################################################################################################
    labsurvey_data <- clean_merge_and_save_labsurvey_data() %>%
        mutate(age_cat=cut(age, age_cuts, my_labels, right=F))

    ####################################################
    ## set up grid and bring in raster covariates etc ##
    ####################################################

    ## bring in the population raster from world pop
    bgd_pop <- here("data",
                    "Bangladesh 100m Population/BGD_ppp_2015_adj_v2.tif") %>%
        raster()

    ## read in age proportions
    raw_age_props <- readRDS(here("generated_data", "all_age_props.rds")) %>%
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

    ## read in raster data
    my_grid <- readRDS(here("generated_data", "grid_with_covs.rds"))

    ## load national polygon
    bgd2 <- load_bgd_shp(level=2) %>%
        transform_to_btm()

    ##############################################################################################
    ## around jagged water borders we have places that have population but no covariate or
    ## other places with no population but covariates. will only predict in places that both
    ## have people and covaraiates and need to make sure none of our sampled sites are in these
    ## get the grid cell ids for each community
    ##############################################################################################

    ## now figure out which villages go with which grid cells
    com_locs <- labsurvey_data %>%
        select(CommLong,CommLat,hhid,villid) %>%
        st_as_sf(coords = c("CommLong","CommLat"),
                 crs="+proj=longlat +datum=WGS84 +no_defs") %>%
        transform_to_btm()

    ## since we are missing some of the exact household GPS data
    ## we will impute the community centroid for these ~150 people
    ana_dat_gps <- labsurvey_data %>%
        mutate(imp_loc = ifelse(is.na(lon),1,0),
               lon = ifelse(is.na(lon),CommLong,lon),
               lat = ifelse(is.na(lat),CommLat,lat))

    hhl_locs <- ana_dat_gps %>%
        select(lon,lat,hhid,villid) %>%
        st_as_sf(coords = c("lon","lat"),
                 crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
        transform_to_btm()

    id_cells <- data.frame(
        community_cell = st_within(com_locs,my_grid) %>%
            unlist,
        grid_id = st_within(hhl_locs,my_grid) %>%
            unlist)

    ## quick look at instances where households are assigned to a different cell than
    ## their community centroid
    ## id_cells %>%
    ##    filter(community_cell != grid_id) %>%
    ##    distinct

    ## now adding grid id and population data to the original data.frame
    analysis_data <- bind_cols(ana_dat_gps,id_cells) %>%
        left_join(my_grid %>% st_set_geometry(NULL)) %>%
        ## adding grid cell centroids
        left_join(
            ## get grid centroids in btm
            my_grid %>%
                group_by(grid_id) %>%
                st_centroid() %>%
                st_coordinates %>%
                data.frame %>%
                mutate(grid_id = 1:n()) %>%
                rename(grid_easting=X, grid_northing=Y))

    ## add randomForest predictions for serosurvey data (trained on cohort data)
    registerDoMC(detectCores()-2)
    cohort_dat <- readRDS("data/bgd-cohort.rds")
    ## fit randomForest at 100 days
    rf_100 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data,
                                   days_back=100,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 200 days
    rf_200 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data,
                                   days_back=200,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## fit randomForest at 365 days
    rf_365 <- fit_rf_make_preds_cv(cohort_dat=cohort_dat,
                                   survey_dat=labsurvey_data,
                                   days_back=365,
                                   k_folds="id",
                                   cutoff_type="youden",
                                   ntree=1e4)
    ## add rf predictions to analysis data
    analysis_data$pred_100 <- rf_100$survey_preds
    analysis_data$pred_200 <- rf_200$survey_preds
    analysis_data$pred_365 <- rf_365$survey_preds

    ## estimate the amount of cholera given acute watery diarrhea and the
    ## number of cases based on AWD sentinel data
    ## read in sentinel surveillance data
    sentinel_data <- read_csv("data/df_subdistrict_cCh_2018.csv",
                              na=c("NULL","NA",""),
                              col_types=paste0(strrep("?",11),"c","c","?","?","c","?","c","c","?","c","?","c",strrep("?",22))) %>%
        rename(adm3name=subdistrict, date_sick=date) %>%
        mutate(year_sick=year(date_sick),
               month_sick=month(date_sick),
               under5=as.numeric(age<5)) %>%
        filter(year_sick<=2016, HOSP<=10) #sentinel surveillance halted from Jan 2016 - April 2016 because of funding
    ## check out how many surveillance cases there were in each month
    # cases_by_month <- sentinel_dat %>%
    #     group_by(year_sick, month_sick) %>%
    #     summarise(cases=n())
    ## read in map data
    # bgd2 <- readRDS("data/BGD_adm2.rds") %>%
    #     st_as_sf() %>%
    #     mutate(adm2name=NAME_2)
    ## read in AWD data
    awd_data <- read_csv("data/awd_dghs_2018.csv") %>%
        rename(adm2name=district, year_sick=year) %>%
        fix_adm2names() %>%
        filter(year_sick<=2016) %>%
        mutate(date_sick=myd(paste(date1,"01")),
               month_sick=month(date_sick))
    # ## check out the AWD by month in the dataset
    # awd_by_month <- awd %>%
    #     group_by(year_sick, month_sick, adm2name) %>%
    #     summarise(cases=sum(awd))
    options(mc.cores=4)
    season_dat <- est_cholera_sentinel(sentinel_data=sentinel_data,
                                       survey_data=analysis_data,
                                       my_grid=my_grid,
                                       bgd2=bgd2,
                                       bgd_pop=bgd_pop,
                                       awd_data=awd_data,
                                       iter=1500,
                                       warmup=250,
                                       control=list(adapt_delta=0.99,
                                                    max_treedepth=20))

    ## estimate daily conditional probabilities of infection (100 days)
    daily_probs_100 <- estimate_daily_probs(analysis_data,
                                            season_dat, 100) %>%
        t() %>% as_tibble() %>%
        mutate(id=seq(nrow(analysis_data))) %>%
        nest(-id, .key="daily_probs_100") %>%
        select(-id)
    ## estimate daily conditional probabilities of infection (200 days)
    daily_probs_200 <- estimate_daily_probs(analysis_data,
                                            season_dat, 200) %>%
        t() %>% as_tibble() %>%
        mutate(id=seq(nrow(analysis_data))) %>%
        nest(-id, .key="daily_probs_200") %>%
        select(-id)
    ## estimate daily conditional probabilities of infection (365 days)
    daily_probs_365 <- estimate_daily_probs(analysis_data,
                                            season_dat, 365) %>%
        t() %>% as_tibble() %>%
        mutate(id=seq(nrow(analysis_data))) %>%
        nest(-id, .key="daily_probs_365") %>%
        select(-id)

    analysis_data <- bind_cols(analysis_data,
                               daily_probs_100,
                               daily_probs_200,
                               daily_probs_365)

    ## save main dat
    saveRDS(analysis_data,here("generated_data","dat_with_raster_covs.rds"))

    ## just a little interactive visualization here
    #p_load(mapview)
    #mapview(my_grid_water %>% filter(is.na(poverty) & pop>100),zcol="pop")
    #mapview(my_grid_water %>% filter(is.na(travel_time) & pop>100),zcol="pop")
    #mapview(my_grid_water %>% filter(distance_to_water<.0001 & pop>100),zcol="pop")
    #mapview(my_grid_water %>% filter(grid_id %in% c(3438,3534)),zcol="pop")
}
