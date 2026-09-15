## utlity functions for bangladesh serosurvey
## cholera-related analyses

## reloads packages and core functions
## (packages needed only for building the dataset from raw data are loaded in
## source/raw-data/generate-analyses-dataset-from-raw.R)
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
           knitr,
           gridExtra,
           mgcv,
           RColorBrewer,
           ggspatial,
           lubridate,
           here,
           doMC,
           ROCR,
           rstan)

    conflict_prefer("here", "here",quiet=T)
    conflict_prefer("filter", "dplyr",quiet=T)
    ## resolve all remaining masked functions up front; newer versions of
    ## conflicted error on *any* unresolved conflict when rstan::stan_model()
    ## scans the search path
    conflict_prefer("accumulate", "purrr",quiet=T)
    conflict_prefer("when", "purrr",quiet=T)
    conflict_prefer("collapse", "dplyr",quiet=T)
    conflict_prefer("combine", "dplyr",quiet=T)
    conflict_prefer("lag", "dplyr",quiet=T)
    conflict_prefer("extract", "rstan",quiet=T)
    conflict_prefer("margin", "ggplot2",quiet=T)

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





## helper lookup
assay_names_labels <- assay_names <- c(`plpsa_s` = "Anti-LPS IgA",
                                       `ptctxa_s` = "Anti-CTB IgA",
                                       `vibinab` = "Vibriocidal\nInaba",
                                       `vibogaw` = "Vibriocidal\nOgawa",
                                       `pglsp_s` = "Anti-LPS IgG",
                                       `pgctxb_s` = "Anti-CTB IgG",
                                       `plpsm_s` = "Anti-LPS IgM",
                                       `pmctxb_s` = "Anti-CTB IgM",
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



#----------------------------------
# TAKEN from Ben Arnold: https://github.com/ben-arnold/enterics-seroepi/blob/master/R/SI-File6-Fig3-Fig4-agecurve-analyses.Rmd
# simulataneous CIs for GAMs
# estimated by resampling the
# Baysian posterior estimates of
# the variance-covariance matrix
# assuming that it is multivariate normal
# the function below also estimates
# the unconditional variance-covariance
# matrix, Vb=vcov(x,unconditional=TRUE),
# which allows for undertainty in the actual
# estimated mean as well
# (Marra & Wood 2012 Scandinavian Journal of Statistics,
#  Vol. 39: 53–74, 2012, doi: 10.1111/j.1467-9469.2011.00760.x )
# simultaneous CIs provide much better coverage than pointwise CIs
# see: http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#----------------------------------
gamCI <- function(m,newdata,nreps=10000) {
    require(mgcv)
    require(dplyr)
    Vb <- vcov(m,unconditional = TRUE)
    pred <- predict(m, newdata, se.fit = TRUE)
    fit <- pred$fit
    se.fit <- pred$se.fit
    BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
    Cg <- predict(m, newdata, type = "lpmatrix")
    simDev <- Cg %*% t(BUdiff)
    absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
    masd <- apply(absDev, 2L, max)
    crit <- quantile(masd, prob = 0.95, type = 8)
    pred <- data.frame(newdata,fit=pred$fit,se.fit=pred$se.fit)
    pred <- mutate(pred,
                   uprP = fit + (2 * se.fit),
                   lwrP = fit - (2 * se.fit),
                   uprS = fit + (crit * se.fit),
                   lwrS = fit - (crit * se.fit)
    )
    return(pred)
}


## takes in dat and outcome column name
#' @param vgm - optional varigram object
#' @param my_dat - sp SpatialPointsDataFrame
#' @param outcome
ordinary_krig <- function(my_dat,outcome,vgm){
    if(is.null(vgm)){
        print("fitting variogram")
        my_formula <- paste0(outcome," ~ 1") %>% as.formula
        vgm <- fit.variogram(variogram(my_formula,data=my_dat),vgm(c("Exp","Sph","Mat"),fit.kappa=FALSE))
        print(vgm)

    }


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
        ## note that this bit gets broken with newer versions of the tibble package you may need to downgrade!
        rc <- readRDS(here("data","BGD_adm2.rds")) %>% sf::st_as_sf() %>% mutate(adm2name = NAME_2)
    } else if(level == 1){
        ## note that this bit gets broken with newer versions of the tibble package you may need to downgrade!
        rc <- readRDS(here("data","BGD_adm1.rds")) %>% sf::st_as_sf() %>% mutate(adm1name = NAME_1)
    } else if(level == 0){
        rc <- readRDS(here("data","BGD_adm0.rds")) %>% sf::st_as_sf()
    }

    return(rc)
}


## not needed anymore
load_bgd_national_shp <- function(){
    return(
        readRDS(here("data","BGD_adm0.rds")) %>% sf::st_as_sf()
    )

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
      `age_group2`="10-19 years",
      `age_group3`="\u2265 19 years",
      `sex1` = "male",
      `travel1` = "trav. last week",
      `travel2` = "trav. last month",
      `travel3` = "trav. last 6 months",
      `elecTRUE` = "electricity in house",
      `ownerHHTRUE` = "owns home",
      `landTRUE` = "owns land",
      `urbanTRUE` = "urban",
      `anyMosqContTRUE` = "mosquito control",
      `educationSumm2` = "secondary school",
      `educationSumm3` = "primary school",
      `educationSumm4` = "no school",
      `incomeSumm2` = "7,000-9,999TK",
      `incomeSumm3` = "10,000-20,000TK",
      `incomeSumm4` = ">20,000TK",
      `logpop` = "population density (log)",
      `travel_time` = "travel time to nearest city (log)",
      `aedesAeg` = "Aedes aegypti mosquitos captured",
      `aedesAlbo` = "Aedes albopictus mosquitos captured",
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
            'pop_density',
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
                             `age_group1` = "0-9 years",
                             `age_group2` = "10-19 years",
                             `age_group3` = ">19 years",
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
                             `logpop` = "population (log)",
                             `logtravel_time` = "travel time to nearest city (log)",
                             `travel_time` = "travel time to nearest city (min)",
                             `aedesAeg` = "Aedes aegypti mosquitos captured",
                             `aedesAlbo` = "Aedes albopictus mosquitos captured",
                             `poverty` = "poverty index",
                             `altitude` = "altitude",
                             `distance_to_water` = "distance to major water body (per 10km)",
                             `pop_density` = "population density (per sq. km)")


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

#' Adjust sensitivity
#'
#' @param sensitivity a dataframe of estimates of sensitivity with columns as days since infection and rows as iterations (in this case, output from Stan file)
#' @param prob_mtx a matrix of conditional probabilities of infection (given infection took place in the last year) where rows are days since infection and columns are individuals
#'
#' @return a vector where each entry represents one simulated estimate of average sensitivity across individuals and days since infection
#' @export
#'
#' @examples
adjust_sensitivity <- function(sensitivity,
                               prob_mtx){
    sensitivity_mtx <- as.matrix(sensitivity) %*% prob_mtx
    apply(sensitivity_mtx, 1, function(x) mean(x,na.rm=T))
}

#' Adjust incidence estimate
#'
#' @param sensitivity a dataframe of estimates of sensitivity with columns as days since infection and rows as iterations
#' @param prob_mtx a matrix of conditional probabilities of infection (given infection took place recently) where rows are days since infection and columns are individuals
#' @param specificity a vector of specificity estimates
#' @param p_A numeric, original predicted incidence estimate to be adjusted; otherwise a vector of the predicted serostatus for each individual
#'
#' @return a vector of adjusted incidence estimates
#' @export
#'
#' @examples
adjust_incidence <- function(sensitivity,
                             prob_mtx,
                             specificity,
                             p_A,
                             B,
                             id1=analysis_data$villid,
                             id2=analysis_data$hhid,
                             id3=analysis_data$sample_id,
                             name1="villid",
                             name2="hhid",
                             name3="sample_id",
                             ...){
    ## if p_A is a vector of greater than length 1 then make bootstrap reps
    ## (called p_reps)
    if(length(p_A)>1){
        analysis_data <- tibble(id1=id1, id2=id2, id3=id3, p_A=p_A)
        colnames(analysis_data) <- c(name1, name2, name3, "p_A")
        tmp_boot <- foreach(i=seq(B),.combine=rbind) %dopar%{
            set.seed(i)
            # sample(p_A, length(p_A), replace=T)
            tmp <- random_clusters(id1, id2, id3, name1, name2, name3, ...)
            tmp$boot <- i
            suppressWarnings(
                suppressMessages(
                    output <- left_join(tmp, analysis_data)
                    ))
            return(output)
        }
        # p_reps <- apply(tmp_boot,2,function(x) mean(x,na.rm=T))
        boot_dat <- tmp_boot %>%
            group_by(boot) %>%
            summarise(boot_p_A=mean(p_A, na.rm=T)) %>%
            ungroup()
        p_reps <- boot_dat$boot_p_A
    } else{
        ## if p_A is a single number, then take draws from a binomial
        ## with size=ncol(prob_mtx)
        p_reps <- rbinom(n=B, size=ncol(prob_mtx), p=p_A)/ncol(prob_mtx)
    }
    ## if sensitivity is a matrix, then it needs to be adjusted by the
    ## daily conditional probabilities of infection
    if(class(sensitivity)=="matrix")
        sensitivity <- adjust_sensitivity(sensitivity=sensitivity,
                                          prob_mtx=prob_mtx)
    ## make the adjusted estimate by taking the original p_A value and the
    ## median sensitivity and specificity values
    adjusted_incidence <- correct_sero_misclass_p(p_A=mean(p_A,na.rm=T),
                                                  sens=median(sensitivity,na.rm=T),
                                                  spec=median(specificity,na.rm=T))
    ## p_reps, sensitivity and specificity all must be same length
    ## if not, then resample to match longest vector
    max_length <- max(c(length(p_reps),length(sensitivity),length(specificity)))
    if(length(p_reps)<max_length){
        p_reps <- sample(p_reps, max_length, T)
    }
    if(length(sensitivity)<max_length){
        sensitivity <- sample(sensitivity, max_length, T)
    }
    if(length(specificity)<max_length){
        specificity <- sample(specificity, max_length, T)
    }
    return(
        list(original_incidence=mean(p_A,na.rm=T),
             unadjusted_reps=p_reps,
             adjusted_incidence=adjusted_incidence,
             adjusted_reps=correct_sero_misclass_p(p_A=p_reps,
                                                   spec=specificity,
                                                   sens=sensitivity),
             adjusted_sensitivity=sensitivity,
             specificity=specificity)
    )
}

#' Find the adjusted incidence by subset with variance estimated by bootstrap
#'
#' @param p_A numeric, original predicted incidence estimate to be adjusted; otherwise a vector of the predicted serostatus for each individual
#' @param prob_mtx matrix, with daily conditional probabilities
#' @param sensitivity matrix, where columns are times since infection and rows are estimates of sensitivity for each time
#' @param specificity vector, of specificity estimates, length should equal the number of rows in sensitivity
#' @param subset_var vector, of the values to subset p_A and prob_mtx by
#' @param subset_vals vector or matrix with 2 columns. If a vector, then the values of subset_var to subset by. If a matrix, then the columns represent the upper and lower bound for the subsets that occur in each row. If NULL, then the default is to subset over all unique values of subset_var.
#' @param B numeric, number of bootstrap simulations
#' @param var_name character, name of subset variable
#'
#' @return data frame with columns for the value that defines the subsetted group, the number of observations, the unadjusted incidence, the adjusted incidence, and the adjusted sensitivity of that group.
#' @export
#'
#' @examples
adjust_incidence_subset <- function(p_A,
                                    prob_mtx,
                                    sensitivity,
                                    specificity,
                                    subset_var,
                                    subset_vals,
                                    B,
                                    var_name=NULL,
                                    id1,
                                    id2,
                                    id3,
                                    name1="villid",
                                    name2="hhid",
                                    name3="sample_id"){
    if(missing(subset_vals))
        ## set subset values based on unique levels of the chosen variable
        subset_vals <- unique(subset_var)
    if(class(subset_vals)=="factor"){
        subset_vals <- as.character(subset_vals)
    }
    if(class(subset_vals)=="matrix"){
        if(ncol(subset_vals)>2)
            stop("subset_vals should be a vector or a matrix with 2 columns")
        subset_var <- cut(subset_var,
                          breaks=c(subset_vals[,1],subset_vals[nrow(subset_vals),2]),
                          right=F)
        subset_vals <- cut(subset_vals[,1],
                           breaks=c(subset_vals[,1],subset_vals[nrow(subset_vals),2]),
                           right=F)
    }
    num_subsets <- length(subset_vals)

    foreach(i=seq(num_subsets),.combine=rbind) %dopar%{
        print(paste("subset", i, "of", num_subsets))
        idx <- which(subset_var==subset_vals[i])
        tmp_output <- adjust_incidence(sensitivity=sensitivity,
                                       prob_mtx=prob_mtx[,idx],
                                       specificity=specificity,
                                       p_A=p_A[idx],
                                       B=B,
                                       id1=id1[idx],
                                       id2=id2[idx],
                                       id3=id3[idx])
        output <- tibble(subset_var=subset_vals[i],
                         obs=length(idx),
                         original_incidence=tmp_output$original_incidence,
                         adjusted_incidence=tmp_output$adjusted_incidence,
                         unadjusted_reps=tmp_output$unadjusted_reps,
                         adjusted_reps=tmp_output$adjusted_reps,
                         adjusted_sensitivity=tmp_output$adjusted_sensitivity,
                         specificity=tmp_output$specificity)
        if(!is.null(var_name))
            colnames(output)[1] <- var_name
        return(output)
    }
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

#' Combine objects within multiple lists in foreach
#' Within foreach .combine=combine_custom will allow the user to bind multiple objects within a list.
#'
#' @param LL1 LL2, two lists that are to be combined within foreach
#'
#' @return A list that combines each of the elements of LL1 and LL2
#' @export
#'
#' @examples
combine_custom <- function(LL1, LL2) {
    combo_spec <- c(LL1$spec, LL2$spec)
    combo_mat <- cbind(LL1$mat, LL2$mat)
    combo_sens <- c(LL1$sens, LL2$sens)
    combo_cutoff <- c(LL1$cutoff, LL2$cutoff)
    return(list(spec = combo_spec, sens=combo_sens, cutoff=combo_cutoff, mat = combo_mat))
}

#' Make bootstrap predictions with randomForest fit
#'
#' @param rf_fit randomForest fit object
#' @param pred_data data.frame to predict
#' @param B numeric, number of bootstraps to
#'
#' @return matrix, rows are observations, columns are bootstraps
#' @export
#'
#' @examples
rf_boot_pred <- function(rf_fit,
                         pred_data,
                         B,
                         cutoff){
    if(missing(cutoff)){
        cutoff <- 0.5
    }
    if(!is.numeric(cutoff))
        cutoff <- find_rf_cutoff(rf_fit=rf_fit,
                                 cutoff_type=cutoff)
    if(!is.numeric(cutoff))
        stop('cutoff needs to either be specified as numeric or estimated using "youden" or "spec".')
    rf_preds <- predict(rf_fit, pred_data, type="prob")[,2]>cutoff
    ## make a matrix of predictions for each observation and bootstrap
    pred_mtx <- foreach(i=seq(B),.combine=cbind) %dopar%{
        set.seed(i)
        sample(rf_preds, length(rf_preds), replace=T)
    }
    return(pred_mtx)
}

#' Find randomForest cutoff
#' Finds the optimal cutoff for randomForest votes using either the Youden statistic or the point with the maximum sensitivity given the maximum specificity
#'
#' @param rf_fit randomForest object
#' @param cutoff_type character, either "youden" or "spec"
#'
#' @return numeric cutoff
#' @export
#'
#' @examples
find_rf_cutoff <- function(rf_fit,
                           cutoff_type=c("youden","spec")){
    cutoff_type <- match.arg(cutoff_type)
    rf_roc <- prediction(rf_fit$votes[,2], rf_fit$y) %>%
        performance("sens","spec")
    if(cutoff_type=="youden"){
        roc_dat <- data.frame(spec=rf_roc@x.values[[1]] %>% unlist,
                              sens=rf_roc@y.values[[1]] %>% unlist,
                              cutoff=rf_roc@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+spec-1)
        cutoff <- roc_dat$cutoff[which.max(roc_dat$youden)]
    }
    if(cutoff_type=="spec"){
        roc_dat <- data.frame(spec=rf_roc@x.values[[1]] %>% unlist,
                              sens=rf_roc@y.values[[1]] %>% unlist,
                              cutoff=rf_roc@alpha.values[[1]] %>% unlist) %>%
            filter(cutoff<=1) %>%
            filter(spec==max(spec))
        cutoff <- roc_dat$cutoff[which.max(roc_dat$sens)]
    }
    return(cutoff)
}

#' Fit randomForest model to cohort data and predict survey data
#'
#' @param cohort_dat data.frame, cohort study data
#' @param survey_dat data.frame, cross-sectional serosurvey data
#' @param days_back numeric, time frame of interest for seropositivity
#' @param n_iter numeric, number of repetitions to estimate the sensitivity and specificity of the random forest model
#' @param cutoff_type character, use the Youden statistic cutoff ("youden") or maximize the specificity ("spec")
#' @param save_results logical, save the individual results within the function?
#' @param ... randomForest parameters
#'
#' @return list with randomForest objects and diagnostics as well as predictions for the survey data
#' @export
#'
#' @examples
fit_rf_make_preds <- function(cohort_dat,
                              survey_dat,
                              days_back,
                              n_iter,
                              cutoff_type=c("youden","spec"),
                              save_results=T,
                              ...){
    cutoff_type <- match.arg(cutoff_type)
    ## list the covariates to be used in the RF model and make into a formula
    shared_covs <- colnames(survey_dat)[which(colnames(survey_dat)%in%colnames(cohort_dat))]
    full_set <- paste(shared_covs, collapse="+")
    rf_formula <- paste0("factor(inf_",days_back,") ~",full_set) %>% as.formula

    ## to estimate the specificity of the RF model
    ## we re-fit the RF model n_iter times, find the cutoff and specificity for each sim
    ## we also save the predictions for use in calculating time-varying sensitivity
    cohort_preds <- foreach(i=seq(n_iter),.combine=combine_custom) %dopar%{
        print(paste("Rep", i, "of", n_iter))
        set.seed(i)
        tmp_forest <- randomForest(rf_formula,data=cohort_dat,...)
        tmp_roc <- prediction(tmp_forest$votes[,2], tmp_forest$y) %>%
            performance("sens","spec")
        if(cutoff_type=="youden"){
            tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                     sens=tmp_roc@y.values[[1]] %>% unlist,
                                     cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                mutate(youden=sens+spec-1)
            tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
        }
        if(cutoff_type=="spec"){
            tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                     sens=tmp_roc@y.values[[1]] %>% unlist,
                                     cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                filter(cutoff<=1) %>%
                filter(spec==max(spec))
            tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$sens)]
        }
        return(list(spec=tmp_roc_df$spec[tmp_roc_df$cutoff==tmp_cutoff],
                    sens=tmp_roc_df$sens[tmp_roc_df$cutoff==tmp_cutoff],
                    cutoff=tmp_cutoff,
                    mat=tmp_forest$votes[,2]>tmp_cutoff))
    }
    ## save specificity
    if(save_results)
        saveRDS(cohort_preds$spec, paste0("generated_data/specificity-rf-",days_back,"-",cutoff_type,"-",Sys.Date(),".rds"))

    ## add the predictions with the median sensitivity to the cohort data
    median_idx <- which.min(abs(cohort_preds$sens-median(cohort_preds$sens)))
    rf_dat <- cohort_dat %>%
        mutate(rf_pred=cohort_preds$mat[,median_idx])

    if(save_results)
        write_csv(rf_dat, path=paste0("generated_data/cohort-w-rf-preds-",days_back,"-",cutoff_type,"-",Sys.Date(),".csv"))

    ## fit and save random forest model
    set.seed(median_idx)
    rf_model <- randomForest(rf_formula, data=cohort_dat,...)

    if(save_results)
        saveRDS(rf_model, paste0("generated_data/forest_no_ogrp_",days_back,".rds"))

    survey_preds <- predict(rf_model, survey_dat, type="prob")[,2] > cohort_preds$cutoff[median_idx]
    return(list(rf_formula=rf_formula,
                cohort_preds=cohort_preds,
                rf_model=rf_model,
                survey_preds=survey_preds,
                days_back=days_back,
                n_iter=n_iter,
                cutoff_type=cutoff_type))
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
                                 cutoff_type=c("youden", "youden_plus","spec"),
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
        tmp_roc <- prediction(tmp_forest$votes[,2], tmp_forest$y) %>%
            performance("sens","spec")
        if(cutoff_type=="youden"){
            tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                     sens=tmp_roc@y.values[[1]] %>% unlist,
                                     cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                mutate(youden=sens+spec-1)
            tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
        }
        if(cutoff_type=="youden_plus"){
            tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                     sens=tmp_roc@y.values[[1]] %>% unlist,
                                     cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                mutate(youden=sens+1.1*spec-1)
            tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
        }
        if(cutoff_type=="spec"){
            tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                     sens=tmp_roc@y.values[[1]] %>% unlist,
                                     cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                filter(cutoff<=1) %>%
                filter(spec==max(spec))
            tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$sens)]
        }
        votes <- predict(tmp_forest, cohort_dat[k_idx,], type="prob")[,2]

        return(cohort_dat[k_idx,] %>%
                   mutate(rf_cutoff=tmp_cutoff,
                          rf_votes=votes,
                          rf_preds=votes>tmp_cutoff))
    }
    if(save_results)
        write_csv(k_preds, file=paste0("generated_data/cohort-w-rf-preds-",days_back,"-",cutoff_type,"-",Sys.Date(),".csv"))

    ## fit and save full random forest model
    set.seed(1)
    rf_model <- randomForest(rf_formula, data=cohort_dat,...)
    if(save_results)
        saveRDS(rf_model, paste0("generated_data/forest_no_ogrp_",days_back,".rds"))
    ## find the optimum cutoff and predict serosurvey
    full_roc <- prediction(rf_model$votes[,2], rf_model$y) %>%
        performance("sens","spec")
    if(cutoff_type=="youden"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+spec-1)
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$youden)]
    }
    if(cutoff_type=="youden_plus"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+1.1*spec-1)
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$youden)]
    }
    if(cutoff_type=="spec"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            filter(cutoff<=1) %>%
            filter(spec==max(spec))
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$sens)]
    }
    survey_preds <- predict(rf_model, survey_dat, type="prob")[,2] > full_cutoff
    return(list(rf_formula=rf_formula,
                cohort_preds=k_preds,
                rf_model=rf_model,
                survey_preds=survey_preds,
                days_back=days_back,
                k_folds=fold,
                cutoff_type=cutoff_type))
}

#' Fit ranger model to cohort data and predict survey data, using cross validation to find sensitivity and specificity
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
fit_ranger_make_preds_cv <- function(cohort_dat,
                                     survey_dat,
                                     days_back,
                                     k_folds=c("id", "loo"),
                                     cutoff_type=c("youden", "youden_plus","spec"),
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
                           tmp_forest <- ranger(rf_formula,
                                                data=cohort_dat[!k_idx,],
                                                probability=T,
                                                ...)
                           tmp_roc <- prediction(tmp_forest$predictions[,2],
                                                 cohort_dat[!k_idx,paste0("inf_",
                                                                          days_back)]) %>%
                               performance("sens","spec")

                           if(cutoff_type=="youden"){
                               tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                                        sens=tmp_roc@y.values[[1]] %>% unlist,
                                                        cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                                   mutate(youden=sens+spec-1)
                               tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
                           }
                           if(cutoff_type=="youden_plus"){
                               tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                                        sens=tmp_roc@y.values[[1]] %>% unlist,
                                                        cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                                   mutate(youden=sens+1.1*spec-1)
                               tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$youden)]
                           }
                           if(cutoff_type=="spec"){
                               tmp_roc_df <- data.frame(spec=tmp_roc@x.values[[1]] %>% unlist,
                                                        sens=tmp_roc@y.values[[1]] %>% unlist,
                                                        cutoff=tmp_roc@alpha.values[[1]] %>% unlist) %>%
                                   filter(cutoff<=1) %>%
                                   filter(spec==max(spec))
                               tmp_cutoff <- tmp_roc_df$cutoff[which.max(tmp_roc_df$sens)]
                           }
                           votes <- predict(tmp_forest,
                                            cohort_dat[k_idx,],
                                            type="response")$predictions[,2]

                           return(cohort_dat[k_idx,] %>%
                                      mutate(rf_cutoff=tmp_cutoff,
                                             rf_votes=votes,
                                             rf_preds=votes>tmp_cutoff))
                       }
    if(save_results)
        write_csv(k_preds, file=paste0("generated_data/cohort-w-rf-preds-",days_back,"-",cutoff_type,"-",Sys.Date(),".csv"))

    ## fit and save full random forest model
    set.seed(1)
    rf_model <- ranger(rf_formula,
                       data=cohort_dat,
                       probability=T,
                       ...)
    ## find the optimum cutoff and predict serosurvey
    full_roc <- prediction(rf_model$predictions[,2],
                           cohort_dat[paste0("inf_", days_back)]) %>%
        performance("sens","spec")

    if(save_results)
        saveRDS(rf_model, paste0("generated_data/forest_no_ogrp_",days_back,".rds"))
    if(cutoff_type=="youden"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+spec-1)
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$youden)]
    }
    if(cutoff_type=="youden_plus"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            mutate(youden=sens+1.1*spec-1)
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$youden)]
    }
    if(cutoff_type=="spec"){
        full_roc_df <- data.frame(spec=full_roc@x.values[[1]] %>% unlist,
                                  sens=full_roc@y.values[[1]] %>% unlist,
                                  cutoff=full_roc@alpha.values[[1]] %>% unlist) %>%
            filter(cutoff<=1) %>%
            filter(spec==max(spec))
        full_cutoff <- full_roc_df$cutoff[which.max(full_roc_df$sens)]
    }
    ## predict the serosurvey data
    survey_preds <- tibble(votes=predict(rf_model, survey_dat,
                                         type="response")$predictions[,2],
                           preds=votes>full_cutoff)
    return(list(rf_formula=rf_formula,
                cohort_preds=k_preds,
                rf_model=rf_model,
                survey_preds=survey_preds,
                days_back=days_back,
                k_folds=fold,
                cutoff_type=cutoff_type))
}



prep_inla_est <- function(eqn="logpop + sexF + age_0_9 + age_10_19 + distance_to_water + poverty + friction + travel_time + logtravel_time + altitude",
                          dat,
                          loc_btm,
                          ran_eff){
    ## MESH
    ## Build mesh around household coordinates
    bnd <- inla.nonconvex.hull(loc_btm,convex=-0.1) ## inla mesh segment

    #meshbuilder() ## shiny app
    mesh <- inla.mesh.2d(boundary = bnd,
                         loc=loc_btm,
                         offset=c(-0.05, -0.05),
                         cutoff = 2000, ## if households are less than 2km apart, builds only a single vertex
                         max.edge=c(30000,50000))

    ## let's have a quick look at the mesh
    # plot(mesh,asp=1,main="com")
    # points(loc_btm[,"X"],loc_btm[,"Y"],pch=19,cex=0.5,col="orange")

    ## now make the SPDE model on the mesh
    spde <- inla.spde2.matern(mesh=mesh, alpha=2) # alpha is Fractional operator order
    s_index <- inla.spde.make.index(name="spatial.field",n.spde=spde$n.spde)

    ## observation/prediction matrix
    A_est_multi <- inla.spde.make.A(mesh=mesh,
                                    loc=loc_btm)

    ## create model matrix of data
    model_mtx <- model.matrix(as.formula(paste("y~",eqn)), dat)[,-1, drop=F]

    ## create stack for estimates
    effect_list <- list(s_index,
                        data.frame(Intercept=rep(1,nrow(model_mtx)),
                                   model_mtx))
    if(!missing(ran_eff)){
        for(i in 1:length(ran_eff)){
            effect_list[[length(effect_list)+1]] <-
                list(as.vector(unlist(dat[,ran_eff[i]])) %>% as.factor())
            names(effect_list[[length(effect_list)]]) <- ran_eff[[i]]
        }
    }
    A_list <- list(A_est_multi,1)
    if(!missing(ran_eff)){
        for(i in 1:length(ran_eff)){
            A_list[length(A_list)+1] <- 1
        }
    }
    est_stack = inla.stack(
        data=list(y=dat[,'y'] %>% unlist,
                  Ntrials=dat[,'obs'] %>% unlist), ## response
        A=A_list, ## projector matrix
        effects=effect_list,
        tag="est")
    ## get stack data
    est_stack_data <- inla.stack.data(est_stack)

    return(list(bnd=bnd,
                mesh=mesh,
                spde=spde,
                s_index=s_index,
                A_est_multi=A_est_multi,
                model_mtx=model_mtx,
                stack=est_stack,
                stack_data=est_stack_data))
}


prep_inla_pred <- function(eqn="logpop + sexF + age_0_9 + age_10_19 + distance_to_water + poverty + friction + travel_time + logtravel_time + altitude",
                           est_dat,
                           loc_btm,
                           ran_eff,
                           pred_dat=my_grid){
    ## prepare the estimation data
    inla_prep_est <- prep_inla_est(eqn=eqn,
                                   dat=est_dat,
                                   loc_btm=loc_btm,
                                   ran_eff=ran_eff)
    ## find centroids for prediction data
    pred_cents <- pred_dat %>%
        st_centroid() %>%
        st_coordinates

    ## calculate a lattice projection from the estimated mesh to the
    ## predicted centroids
    pred_projector <- inla.mesh.projector(inla_prep_est$mesh,
                                          loc=pred_cents)
    ## get A matrix
    A_pred <- pred_projector$proj$A

    ## create prediction stack
    stack_pred <- inla.stack(
        data=list(y=NA,Ntrials=1),
        A=list(A_pred,1),
        effects=list(inla_prep_est$s_index,
                     data.frame(Intercept=rep(1,nrow(pred_dat)))),
        tag="pred")

    ## combine stacks
    est_and_pred_stack <- inla.stack(stack_pred,inla_prep_est$stack)

    ## get stack data
    est_and_pred_data <- inla.stack.data(est_and_pred_stack)

    return(list(est_prep=inla_prep_est,
                pred_projector=pred_projector,
                A_pred=A_pred,
                stack_pred=stack_pred,
                est_and_pred_stack=est_and_pred_stack,
                est_and_pred_data=est_and_pred_data))
}
