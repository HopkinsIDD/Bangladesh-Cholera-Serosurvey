## script to generate the raster data for mapping
## by Andrew Azman and Stephen A Lauer, June 2020

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

## load libraries and utilities
library(pacman)
p_load(here)
setwd(here())
source(here("source","utils.R"))
reload_source()
p_load(readxl, tidyxl, raster, exactextractr, fasterize, brms, scam)
set.seed(1)
age_cuts <- c(0,10,20,100)
my_labels <- paste("age", age_cuts[-length(age_cuts)], age_cuts[-1]-1, sep="_")

data_missing <- (is.na(find_recent_file(name_start="grid_with_covs.rds",
                                 path="generated_data/")) |
                     is.na(find_recent_file(name_start="all_age_props.rds",
                                            path="generated_data/")))

if(data_missing | RERUN_CODE){
    ## bring in the population raster from world pop
    bgd_pop <- here("data",
                    "Bangladesh 100m Population/BGD_ppp_2015_adj_v2.tif") %>%
        raster()

    ## load national polygon
    bgd0 <- load_bgd_shp(level=0) %>%
        transform_to_btm()

    ## making a 5 by 5 grid over bangladesh and then assigning ids
    my_grid <- st_make_grid(bgd0,cellsize=c(5000,5000)) %>%
        st_sf(grid_id = 1:length(.))

    ## using exact_extract() to extract pop data quickly
    pop_grid <-  bgd_pop %>%
        transform_to_btm() %>%
        exact_extract(my_grid, 'sum')

    ## adding population to grid cells
    my_grid_pop <- my_grid %>%
        mutate(pop = pop_grid,
               logpop = log(pop_grid)) %>%
        ## bring in distance to water body
        left_join(here("data","distance_to_major_water.shp") %>%
                      st_read() %>%
                      select(grid_id,distance) %>%
                      mutate(distance = ifelse(distance==0,0.0001,distance),
                             distance = distance/10000) %>%
                      rename(distance_to_water=distance) %>%
                      st_set_geometry(NULL)) %>%
        filter(pop > 0) %>%
        mutate(grid_id=as.numeric(as.factor(grid_id)))
    # plot(my_grid_pop["logpop"])
    # plot(my_grid_pop["distance_to_water"])

    ## add worldpop's poverty Wealth index from
    ## https://www.worldpop.org/visual/bgd_poverty
    ## Mapping poverty using mobile phone and satellite data, Steele et al, 2017?
    my_grid_pop$poverty <- here("data","world_pop_poverty/bgd2011wipov.tif") %>%
        raster() %>%
        transform_to_btm %>%
        exact_extract(my_grid_pop, 'mean')
    # plot(my_grid_pop["poverty"])

    ## let's also add the friction surface from Malaria atlas project
    my_grid_pop$friction <- find_recent_file("2015_friction_surface",
                                             path=here("data","getRaster/")) %>%
        raster() %>%
        transform_to_btm %>%
        exact_extract(my_grid_pop, 'max')
    # plot(my_grid_pop["friction"])

    travel_time_raster <- find_recent_file("2015_accessibility_to_cities",
                                           here("data","getRaster/")) %>%
        raster() %>%
        transform_to_btm
    travel_time_raster[travel_time_raster<0] <- NA #values in the ocean have negative values
    travel_time_raster[travel_time_raster==0] <- 0.0001  #for those in cities, setting this to a small number rather than having zeros
    my_grid_pop$travel_time <-  exact_extract(travel_time_raster,
                                              my_grid_pop,
                                              'mean') # mean and median about the same, but mean has built in function
    # plot(my_grid_pop["travel_time"])
    my_grid_pop$logtravel_time <- log(my_grid_pop$travel_time)

    ## get altitude
    my_grid_pop$altitude <- raster::getData("alt",country="BGD") %>%
        transform_to_btm %>%
        exact_extract(my_grid_pop, 'mean')
    # plot(my_grid_pop["altitude"])

    ## if missing, add age and sex demographic data
    if(is.na(find_recent_file(name_start="all_age_props.rds",
                              path=here("generated_data/")))){
        age_list <- list.files(path="data/BGD",
                               pattern="2015.tif")
        registerDoMC(4)
        time_in <- Sys.time()
        age_pop <- foreach(i=1:length(age_list), .combine=rbind, .inorder=F) %dopar% {
            tmp <- raster(here("data","BGD", age_list[i])) %>%
                transform_to_btm %>%
                exact_extract(my_grid_pop, 'sum')

            class_features <- strsplit(age_list[i], "_")[[1]]
            tibble(sex=class_features[2],
                   age=class_features[3],
                   subpop=tmp) %>%
                mutate(grid_id=1:n())
        }
        time_out <- Sys.time()
        time_out-time_in

        raw_age_props <- age_pop %>%
            group_by(grid_id) %>%
            mutate(age_prop=subpop/sum(subpop)) %>%
            ungroup()

        saveRDS(raw_age_props, here("generated_data", "all_age_props.rds"))
    }

    adj_age_props <- raw_age_props %>%
        ungroup()
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

    my_grid_pop$sexF <- (age_props %>%
                             group_by(grid_id, sex) %>%
                             summarize(sexF=sum(prop)) %>%
                             ungroup() %>%
                             filter(sex=="f"))$sexF

    my_grid_pop$age_0_9 <- (age_props %>%
                                group_by(grid_id, age_grp) %>%
                                summarize(age_prop=sum(prop)) %>%
                                ungroup() %>%
                                filter(age_grp == "age_0_9"))$age_prop

    my_grid_pop$age_10_19 <- (age_props %>%
                                group_by(grid_id, age_grp) %>%
                                summarize(age_prop=sum(prop)) %>%
                                ungroup() %>%
                                filter(age_grp == "age_10_19"))$age_prop

    my_grid_area <- my_grid_pop %>%
        st_intersection(bgd0) %>%
        select(grid_id)

    my_grid_area$area <- st_area(my_grid_area) %>% as.numeric()

    my_grid_pop$area <- NA
    for(i in 1:nrow(my_grid_pop)){
        if(!my_grid_pop$grid_id[i] %in% my_grid_area$grid_id){
            next
        }
        my_grid_pop$area[i] <- my_grid_area$area[which(my_grid_area$grid_id == my_grid_pop$grid_id[i])]
    }
    my_grid_pop$pop_density <- my_grid_pop$pop / my_grid_pop$area
    ## save grid with covariates
    saveRDS(my_grid_pop,here("generated_data", "grid_with_covs.rds"))
    ## my_grid_pop <- readRDS(here("generated_data", "grid_with_covs.rds"))
}
