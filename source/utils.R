## utlity functions for bangladesh serosurvey
## cholera-related analyses

## reloads pacakges and core functions (including this script)
reload_source <- function(){

    library(pacman)

    p_load(sf,
           conflicted,
           dplyr,
           stringr,
           purrr,
           tidyr,
           ggplot2,
           readr,
           randomForest,
           stringr,
           knitr,
           gridExtra,
           mgcv,
           RColorBrewer,
           ggsn,
           lubridate,
           here,
           doMC,
           ROCR,
           rstan)

    conflict_prefer("here", "here",quiet=T)
    conflict_prefer("filter", "dplyr",quiet=T)

    # here()

    # source(here("source","vibriocidal_extraction_funcs.R"))
    # source(here("source","utils.R"))

}

## aligns Zila names with gadm name standards
fix_adm2names <- function(dat){
    dat %>%
        mutate(adm2name = recode(adm2name,
                                 "Barguna" = "Borgona",
                                 "Munshiganj" = "Munshigonj",
                                 "Mymensingh" = "Nasirabad",
                                 "Netrokona" = "Netrakona",
                                 "Narsingdi" = "Narshingdi",
                                 "Kushtia" = "Kustia",
                                 "Chapai Nawabganj" = "Nawabganj",
                                 "Gaibandha" = "Gaibanda",
                                 "Rangpur" = "Rongpur",
                                 "Habiganj" = "Hobiganj",
                                 "Maulvibazar" = "Moulvibazar",
                                 "Cox's bazar"="Cox's Bazar",
                                 "Jhalokati"="Jhalakati",
                                 "Bandarban"="Bandarbon",
                                 "Khagrachhari"="Khagrachari",
                                 "Rangamati"="Parbattya Chattagram",
                                 "Gopalganj"="Gopalgonj",
                                 "Manikganj"="Manikgonj",
                                 "Narayanganj"="Naray Angonj",
                                 "Chuadanga"="Choua Danga",
                                 "Satkhira"="Shatkhira",
                                 "Joypurhat"="Jaipurhat",
                                 "Chapainawabganj"="Nawabganj",
                                 "Sirajganj"="Sirajgonj",
                                 "Sunamganj"="Sun Amgonj"))
}

##############################################
## some spatial helper files for this study ##
##############################################

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @return
##' @author Andrew A
load_bgd_shp <- function(level=2){

    if(level==2){
        rc <- readRDS(here("data","BGD_adm2.rds")) %>% sf::st_as_sf() %>% mutate(adm2name = NAME_2)
    } else if(level == 1){
        rc <- readRDS(here("data","BGD_adm1.rds")) %>% sf::st_as_sf() %>% mutate(adm1name = NAME_1)
    } else if(level == 0){
        rc <- readRDS(here("data","BGD_adm0.rds")) %>% sf::st_as_sf()
    }

    return(rc)
}

##' loads sampling locations as an sf object from cleaned study data
##' @title
##' @return
##' @author Andrew A
load_sampling_locs_sf <- function(){

    dat <- load_and_clean_data()

    return(
        sf::st_as_sf(dat %>%
                         filter(!is.na(lat)) %>%
                         dplyr::select(lon,lat),coords = c("lon","lat"),crs="+proj=longlat +datum=WGS84 +no_defs")
    )
}

## labels for variable names
variable_labs <- ggplot2::as_labeller(
    c(`ageSumm2` = "11-20 years",
      `ageSumm3` = "21-30 years",
      `ageSumm4` = "31-40 years",
      `ageSumm5` = "41-50 years",
      `ageSumm6` = "51-60 years",
      `ageSumm7` = ">60 years",
      `age_group2`="5-14 years",
      `age_group3`="\u2265 15 years",
      `sex1` = "male",
      `travel1` = "trav. last week",
      `travel2` = "trav. last month",
      `travel3` = "trav. last 6 months",
      `elecTRUE` = "electricity in house",
      `ownerHHTRUE` = "owns home",
      `landTRUE` = "owns land",
      `urbanTRUE` = "urban",
      `educationSumm2` = "secondary school",
      `educationSumm3` = "primary school",
      `educationSumm4` = "no school",
      `incomeSumm2` = "7,000-9,999TK",
      `incomeSumm3` = "10,000-20,000TK",
      `incomeSumm4` = ">20,000TK",
      `logpop` = "population density (log)",
      `travel_time` = "travel time to nearest city (log)",
      `poverty` = "poverty index",
      `distance_to_water` = "distance to major water body (per 10km)"
    ))


recode_risk_factors <- function(x) {

    factor_levels <-
        c(
            'age_group1',
            'age_group2',
            'age_group3',
            'ageSumm1',
            'ageSumm2',
            'ageSumm3',
            'ageSumm4',
            'ageSumm5',
            'ageSumm6',
            'ageSumm7',
            'sexMale',
            'travel1',
            'travel2',
            'travel3',
            'elecTRUE',
            'ownerHHTRUE',
            'landTRUE',
            'urbanTRUE',
            'educationSumm1',
            'educationSumm2',
            'educationSumm3',
            'educationSumm4',
            'incomeSumm1',
            'incomeSumm2',
            'incomeSumm3',
            'incomeSumm4',
            'HHSize',
            'logpop',
            'travel_time',
            'poverty',
            'altitude',
            'distance_to_water')

    factor(x,levels=factor_levels) %>%
        recode(.,!!!risk_factor_label_lookup)

}


## corrects num positive by sens and spec
correct_sero_misclass <- function(num_pos,num_neg,sens=.806,spec=.83){
    pmax((num_pos + (num_pos+num_neg)*(spec -1)) / (sens + spec - 1),0)
}

## takes known sensitivity and specificity of test
## and returns proportion of sample that 'true' positive
##' @param p_A - proportion of positives by imperfect test
##' @param N - number tested
##' @param sens - senstivitity
##' @param spec - specifcity
##' @return numeric vector
correct_sero_misclass_p <- function(p_A,sens=.891,spec=.792){
    pmin(pmax((p_A + (spec-1))/(sens+spec -1),0),1)
}


## label lookup function for plotting and tables
risk_factor_label_lookup <-c(`ageSumm2` = "11-20 years",
                             `ageSumm3` = "21-30 years",
                             `ageSumm4` = "31-40 years",
                             `ageSumm5` = "41-50 years",
                             `ageSumm6` = "51-60 years",
                             `ageSumm7` = ">60 years",
                             `age_group1` = "0-4 years",
                             `age_group2` = "5-14 years",
                             `age_group3` = ">14 years",
                             `sexMale` = "male",
                             `travel1` = "travel last week",
                             `travel2` = "travel last month",
                             `travel3` = "travel last 6 months",
                             `elecTRUE` = "electricity in house",
                             `ownerHHTRUE` = "owns home",
                             `landTRUE` = "owns land",
                             `urbanTRUE` = "urban",
                             `anyMosqContTRUE` = "mosquito control",
                             `educationSumm1` = "post-secondary education (head hhl)",
                             `educationSumm2` = "secondary school (head hhl)",
                             `educationSumm3` = "primary school (head hhl)",
                             `educationSumm4` = "no school (head hhl)",
                             `incomeSumm1` = "<7,000TK",
                             `incomeSumm2` = "7,000-9,999TK",
                             `incomeSumm3` = "10,000-20,000TK",
                             `incomeSumm4` = ">20,000TK",
                             `logpop` = "population density (log)",
                             `logtravel_time` = "travel time to nearest city (log)",
                             `travel_time` = "travel time to nearest city (min)",
                             `aedesAeg` = "Aedes aegypti mosquitos captured",
                             `aedesAlbo` = "Aedes albopictus mosquitos captured",
                             `poverty` = "poverty index",
                             `altitude` = "altitude",
                             `distance_to_water` = "distance to major water body (per 10km)")


#' @param name_start character string, first letters in file name
#' @param path character string, path to folder of interest, end with "/"
#' @param exclude character string, patterns to exclude from the file names of interest
#'
#' @return character string, path to most recent file
#' @export
#'
#' @examples
find_recent_file <- function(name_start, path, exclude=NULL){
    if(substring(path, nchar(path))!="/")
        warning("Path does not end with a '/', problems may ensue.")
    ## view all files of that name at that path
    file_list <- list.files(path=path,
                            pattern=paste0(name_start, "*"))
    ## remove files with unwanted patterns
    if(!is.null(exclude)){
        for(i in 1:length(exclude))
            file_list <- file_list[!grepl(pattern = exclude[i], file_list)]
    }
    if(length(file_list)==0){
        warning('File not found')
        return(NA)
    }
    ## view file info
    file_info <- file.info(paste0(path, file_list))
    ## find most recent file
    most_recent_file <- paste0(path,
                               file_list[which.max(file_info$mtime)])

    cat(sprintf("Loaded file: \n %s last updated on \n %s \n",most_recent_file,file_info$mtime[which.max(file_info$mtime)]))

    return(most_recent_file)
}

#' Makes age-sex pyramid for data
#' @param dat
#' @return
make_age_pyramid <- function(dat){
    cdat <- read_csv(here("data","census_data_sex_age.csv"),skip=1) %>% filter(Age != "100+" & Age != "95-99") %>%
        mutate(prop_male = 100*`Male Population`/sum(`Both Sexes Population`),
               prop_female = 100*`Female Population`/sum(`Both Sexes Population`)
        )

    dat <- dat %>% mutate(age_cat = cut(age,seq(0,100,by=5),right = FALSE) %>% ordered %>% droplevels)

    ind_summary <- dat %>% summarize(mean_age = mean(age,na.rm=T),
                                     prop_u5 = mean(age<5,na.rm=T))

    age_sex_dat <- dat %>% filter(!is.na(age_cat)) %>% group_by(sex,age_cat) %>%
        summarize(count = n()) %>%
        ungroup %>%
        mutate(prop = 100* count / sum(count))


    cdat$age_cat <-  ordered(levels(age_sex_dat$age_cat),levels = levels(age_sex_dat$age_cat))

    age_sex_dat %>% ggplot(aes(x=age_cat)) +
        geom_histogram(data = age_sex_dat %>% filter(sex=="Male"),aes(x=age_cat,y=prop),
                       binwidth = 1,fill="steelblue",stat="identity") +
        geom_histogram(data = age_sex_dat %>% filter(sex=="Female") %>% filter(!is.na(age_cat)),aes(y=-prop),
                       binwidth = 1,fill="pink",stat="identity") +
        geom_point(data = cdat,aes(x=age_cat,y=prop_male)) +
        geom_point(data = cdat,aes(x=age_cat,y=-prop_female)) +
        coord_flip() +
        scale_y_continuous(breaks = -10:10,labels = abs(-10:10),limits=c(-7,7)) +
        theme_minimal() +
        xlab("age category") +
        ylab("percent of total population")

}

## transforms sf file to Bangladesh Transverse Mercator projection
transform_to_btm <- function(my_obj){

    if ("sf" %in% class(my_obj)){
        rc <- st_transform(my_obj,
                           crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs")
    } else if(class(my_obj) == "RasterLayer"){
        rc <- raster::projectRaster(my_obj,
                                    crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs",
                                    method = "bilinear" )
    } else {
        stop("can only transform rasterlayers and sf objects at the moment")
    }

}

#' Find ROC curve cutoff
#' Finds the optimal cutoff for probabilities using either the Youden statistic or the point with the maximum sensitivity given the maximum specificity
#'
#' @param preds numeric vector of predictions, each between 0 and 1
#' @param labs numeric or logical vector of outcomes, each equaling 1/0 or T/F
#' @param cutoff_type character, either "youden" or "spec"
#'
#' @return numeric cutoff
#' @export
#'
#' @examples
find_roc_cutoff <- function(preds,
                            labs,
                            cutoff_type=c("youden","spec")){
    cutoff_type <- match.arg(cutoff_type)
    roc_sens_spec <- prediction(preds, labs) %>%
        performance("sens","spec")
    if(cutoff_type=="youden"){
        roc_dat <- data.frame(spec=roc_sens_spec@x.values[[1]] %>% unlist,
                              sens=roc_sens_spec@y.values[[1]] %>% unlist,
                              cutoff=roc_sens_spec@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+spec-1)
        cutoff <- roc_dat$cutoff[which.max(roc_dat$youden)]
    }
    if(cutoff_type=="spec"){
        roc_dat <- data.frame(spec=roc_sens_spec@x.values[[1]] %>% unlist,
                              sens=roc_sens_spec@y.values[[1]] %>% unlist,
                              cutoff=roc_sens_spec@alpha.values[[1]] %>% unlist) %>%
            filter(cutoff<=1) %>%
            filter(spec==max(spec))
        cutoff <- roc_dat$cutoff[which.max(roc_dat$sens)]
    }
    return(cutoff)
}

#' Fit randomForest model to cohort data and predict survey data, using cross validation to find sensitivity and specificity
#'
#' @param cohort_dat data.frame, cohort study data
#' @param survey_dat data.frame, cross-sectional serosurvey data
#' @param days_back numeric, time frame of interest for seropositivity
#' @param k_folds character or numeric, how should the cross validation be carried out? If "loo" then leave-one-out cross validation is performed. If a character string that matches a column name, then folds are determined by that covariate (e.g. "id" will create a new fold for each id). If a number is given, then that is the number of CV folds, e.g. 10-fold CV.
#' @param cutoff_type character, use the Youden statistic cutoff ("youden") or maximize the specificity ("spec")
#' @param save_results logical, save the individual results within the function?
#' @param ... randomForest parameters
#'
#' @return list with randomForest objects and diagnostics as well as predictions for the survey data
#' @export
#'
#' @examples
fit_rf_make_preds_cv <- function(cohort_dat,
                                 survey_dat,
                                 days_back,
                                 k_folds=c("id", "loo"),
                                 cutoff_type=c("youden","spec"),
                                 save_results=T,
                                 ...){
    cutoff_type <- match.arg(cutoff_type)
    ## assign observations to folds
    if(k_folds=="loo"){
        fold <- seq(1:nrow(cohort_dat))
    } else
        if(k_folds %in% colnames(cohort_dat)){
            ids <- unique(cohort_dat[,k_folds,T])
            fold <- rep(NA, nrow(cohort_dat))
            for(i in 1:length(ids)){
                fold[cohort_dat[,k_folds,T]==ids[i]] <- i
            }
        } else
            if(as.numeric(k_folds)){
                fold <- sample(k_folds, nrow(cohort_dat), replace=T)
            } else
                stop("k_fold should be 'loo', the name of a column in cohort_dat, or a number.")

    ## list the covariates to be used in the RF model and make into a formula
    shared_covs <- colnames(survey_dat)[which(colnames(survey_dat)%in%colnames(cohort_dat))]
    full_set <- paste(shared_covs, collapse="+")
    rf_formula <- paste0("factor(inf_",days_back,") ~",full_set) %>% as.formula

    ## we run k-fold cross validation with random forests leaving out
    ## one fold, fitting on the remaining observations, then predicting
    ## the left-out fold
    k_preds <- foreach(i=1:max(fold),
                       .combine=rbind) %dopar%{
        print(paste("Rep", i, "of", length(unique(fold))))
        set.seed(i)
        k_idx <- fold==i
        tmp_forest <- randomForest(rf_formula,data=cohort_dat[!k_idx,],...)
        tmp_cutoff <- find_roc_cutoff(preds=tmp_forest$votes[,2],
                                      labs=tmp_forest$y)
        votes <- predict(tmp_forest, cohort_dat[k_idx,], type="prob")[,2]

        return(cohort_dat[k_idx,] %>%
                   mutate(rf_cutoff=tmp_cutoff,
                          rf_votes=votes,
                          rf_preds=votes>tmp_cutoff))
    }
    if(save_results)
        write_csv(k_preds, path=paste0("generated_data/cohort-w-rf-preds-",days_back,"-",cutoff_type,"-",Sys.Date(),".csv"))

    ## fit and save full random forest model
    set.seed(1)
    rf_model <- randomForest(rf_formula, data=cohort_dat,...)
    if(save_results)
        saveRDS(rf_model, paste0("generated_data/forest_no_ogrp_",days_back,".rds"))
    ## find the optimum cutoff and predict serosurvey
    full_cutoff <- find_roc_cutoff(preds=rf_model$votes[,2],
                                   labs=rf_model$y)
    survey_preds <- predict(rf_model, survey_dat, type="prob")[,2] > full_cutoff
    return(list(rf_formula=rf_formula,
                cohort_preds=k_preds,
                rf_model=rf_model,
                survey_preds=survey_preds,
                days_back=days_back,
                k_folds=fold,
                cutoff_type=cutoff_type))
}

#' Make a bootstrapped map from sampled data
#'
#' @param my_rep numeric, simulation number
#' @param full_dat data.frame, containing all data
#' @param mesh inla.mesh object,
#' @param my_grid
#' @param s_index
#' @param A.pred
#'
#' @return
#' @export
#'
#' @examples
samp_map <- function(my_rep,
                     full_dat,
                     mesh,
                     my_grid,
                     s_index,
                     A.pred){
    set.seed(my_rep)

    tdat <- make_tmp_dat(my_rep,full_dat)

    out <- make_boot_map(tdat,mesh,my_grid,s_index,A.pred)

    return(out)
}

## helper lookup
assay_names_labels <- assay_names <- c(`lps_iga` = "Anti-LPS IgA",
                                       `ctb_iga` = "Anti-CTB IgA",
                                       `vibinab` = "Vibriocidal\nInaba",
                                       `vibogaw` = "Vibriocidal\nOgawa",
                                       `lps_igg` = "Anti-LPS IgG",
                                       `ctb_igg` = "Anti-CTB IgG",
                                       `lps_igm` = "Anti-LPS IgM",
                                       `ctb_igm` = "Anti-CTB IgM",
                                       `vibo_to_lps`= "Vibriocidal to anti-CTB IgG",
                                       `vibo_to_ctx` = "Vibriocidal to anti-LPS IgG",
                                       `lps_iga_inaba`="Anti-LPS IgA\nInaba",
                                       `lps_iga_ogawa`="Anti-LPS IgA\nOgawa",
                                       `lps_igg_inaba`="Anti-LPS IgG\nInaba",
                                       `lps_igg_ogawa`="Anti-LPS IgG\nOgawa",
                                       `maxvib` = "max(Inaba,Ogawa) Vibriocidal",
                                       `titer_interp_simp_inaba` = "Interp. Vibriocidal Inaba (simp)",
                                       `titer_interp_simp_ogawa` = "Interp. Vibriocidal Ogawa (simp)",
                                       `titer_top_interp_inaba` = "Interp. Vibriocidal Inaba",
                                       `titer_top_interp_ogawa` = "Interp. Vibriocidal Ogawa",
                                       `titer_interp_simp_max` = "Interp. Vibriocidal (max,simp)",
                                       `titer_top_interp_max` = "Interp. Vibriocidal (max)",
                                       `2` = "Day 2",
                                       `7` = "Day 7",
                                       `30` = "Day 30",
                                       `90` = "Day 90",
                                       `180` = "Day 180",
                                       `360` = "Day 365",
                                       `365` = "Day 365",
                                       `540` = "Day 540",
                                       `720` = "Day 720",
                                       `900` = "Day 900",
                                       `10` = "Day 10",
                                       `28` = "Day 28",
                                       `170` = "Day 170",
                                       `0` = "Day 0"
)


## function to make all the data-specific SDPE bits and
## to estimate the map
make_boot_map <- function(tmp_dat,
                          mesh,
                          my_grid,
                          s_index,
                          A.pred){

    ## observation/prediction matrix
    A_est_multi <- inla.spde.make.A(mesh=mesh,
                                    loc=tmp_dat %>%
                                        select(grid_easting,
                                               grid_northing) %>%
                                        as.matrix)

    ## model matrix using only spatial correlation
    model_matrix <- model.matrix(~logpop+distance_to_water,
                                 tmp_dat)[, -1]


    ## creating a stack called "est"
    stack_est_cor <- inla.stack(
        data=list(y=tmp_dat[,'n_pos_rf_cor'] %>% unlist,
                  n=tmp_dat[,'n_tested'] %>% unlist), ## response
        A=list(A_est_multi,1), ## projector matrix
        effects=list(s_index,
                     data.frame(Intercept=rep(1,nrow(tmp_dat)),
                                model_matrix)),
        tag="est")

    ## mean imputation for missing covaraiates (due to clipping issues and imperfect alignment of rasters)
    ## the idea is to predict everywhere but then to exclude cells with missing covariate data from analyses
    stack_pred <- inla.stack(
        data=list(y=NA,n=1),
        A=list(A.pred,1,1,1),
        effects=list(s_index,
                     data.frame(Intercept=rep(1,nrow(my_grid))),
                     list(logpop=ifelse(is.na(my_grid$logpop),
                                        mean(my_grid$logpop),
                                        my_grid$logpop)),
                     list(distance_to_water=
                              ifelse(is.na(my_grid$distance_to_water),
                                     mean(my_grid$distance_to_water),
                                     my_grid$distance_to_water))),
        tag="pred")

    est_and_pred_stack_cor <- inla.stack(stack_pred,stack_est_cor)

    ## formula for model
    all_covs_f <- y ~ -1 + Intercept + logpop + distance_to_water +
        f(spatial.field, model=spde)

    ## run the model
    inla_fit <- inla(all_covs_f,
                     data=inla.stack.data(est_and_pred_stack_cor),
                     family="binomial",
                     quantiles=c(0.025, 0.5, .975),
                     verbose=FALSE,
                     control.predictor=list(
                         A=inla.stack.A(est_and_pred_stack_cor),link=1,
                         compute=FALSE),
                     control.compute=list(dic=FALSE,
                                          config=TRUE,
                                          hyperpar=FALSE,
                                          mlik=FALSE),
                     Ntrials=n
    )

    inla_preds <- inla.posterior.sample(result=inla_fit)
    ## get ids for each cell
    id.prd <- inla.stack.index(est_and_pred_stack_cor, "pred")$data
    posterior_means <- inla_fit$summary.fitted.values$mean[id.prd]
    posterior_sds <- inla_fit$summary.fitted.values$sd[id.prd]
    posterior_lb <- inla_fit$summary.fitted.values$`0.025quant`[id.prd]
    posterior_ub <- inla_fit$summary.fitted.values$`0.975quant`[id.prd]
    posterior_median <- inla_fit$summary.fitted.values$`0.5quant`[id.prd]
    posterior_sample <- plogis(inla_preds[[1]]$latent[id.prd,1])

    rc <- data.frame(
        grid_id=id.prd,
        mean=posterior_means,
        sd=posterior_sds,
        rr=posterior_means/
            weighted.mean(posterior_means,w=my_grid$pop),
        median=posterior_median,
        lb=posterior_lb,
        ub=posterior_ub,
        sample_rate=posterior_sample,
        sample_total=posterior_sample*my_grid$pop,
        sample_rr=posterior_sample/
            weighted.mean(posterior_sample,w=my_grid$pop)
    )

    return(rc)
}
