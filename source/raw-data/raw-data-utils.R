## functions used only to build the analysis dataset from the raw (non-public) lab,
## survey and sentinel-surveillance data; see source/raw-data/generate-analyses-dataset-from-raw.R
## Original study sample ids (used to resolve two duplicated samples) are redacted
## as SAMPLE_ID_1 and SAMPLE_ID_2.

##' loads lab data, merges with survey data and
##' does a little recoding etc.
clean_merge_and_save_labsurvey_data <- function(){

    ## anti-CTB IgG
    ctb_igg <- read_xlsx(here("data","ELISA GM1CT Oct 2018.xlsx"),
                         sheet="GM1-CT IgG",
                         range = "B2:S2934",
                         col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`)) %>%
        filter(sample_id != "SAMPLE_ID_1-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_ctb_igg=str_split(string=Date,pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_ctb_igg = str_split(string=Date,pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>% rename(pgctxb_s = titer) %>%
        select(sample_id,pgctxb_s,date_of_test_ctb_igg,test_plate_ctb_igg)

    ## anti-CTB IgA
    ctb_iga <- read_xlsx(here("data","ELISA GM1CT Oct 2018.xlsx"),
                         sheet="GM1-CT IgA",
                         range = "B2:S2934",
                         col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`)) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_ctb_iga=str_split(string=Date,pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_ctb_iga = str_split(string=Date,pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>% rename(ptctxa_s = titer) %>%
        dplyr::select(sample_id,ptctxa_s,
                      date_of_test_ctb_iga,test_plate_ctb_iga)

    ## anti-LPS IgG Ogawa
    lps_igg_ogawa <- read_xlsx(here("data","ELISA LPS Ogawa Oct 2018.xlsx"),
                               sheet="LPS og IgG",
                               range="A2:S2934",
                               col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`)) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_lps_igg_ogawa = str_split(string=Date,
                                                      pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_lps_igg_ogawa = str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>%
        rename(lps_igg_ogawa = titer) %>%
        select(sample_id,lps_igg_ogawa,date_of_test_lps_igg_ogawa,
               test_plate_lps_igg_ogawa)

    ## anti-LPS IgG Inaba
    lps_igg_inaba <- read_xlsx(here("data","ELISA LPS Inaba Oct 2018.xlsx"),
                               sheet="LPS-IgG Inaba",
                               range="B2:T2934",
                               col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`)) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_lps_igg_inaba = str_split(string=Date,
                                                      pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_lps_igg_inaba = str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>%
        rename(lps_igg_inaba = titer) %>%
        select(sample_id,lps_igg_inaba,date_of_test_lps_igg_inaba,
               test_plate_lps_igg_inaba)

    ## anti-LPS IgA Ogawa
    lps_iga_ogawa <- read_xlsx(here("data","ELISA LPS Ogawa Oct 2018.xlsx"),
                               sheet="LPS og IgA",
                               range="A2:T2934",
                               col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`))%>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_lps_iga_ogawa=str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_lps_iga_ogawa = str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>%
        rename(lps_iga_ogawa = titer) %>%
        select(sample_id,lps_iga_ogawa,date_of_test_lps_iga_ogawa,
               test_plate_lps_iga_ogawa)

    ## anti-LPS IgA Inaba
    lps_iga_inaba <- read_xlsx(here("data","ELISA LPS Inaba Oct 2018.xlsx"),
                               sheet="LPS-IgA Inaba",
                               range="B2:T2934",
                               col_types=c('text')) %>%
        rename(sample_id = `Sample ID`) %>%
        mutate(titer = 100*as.numeric(`Mean Sample`)/as.numeric(`Mean PS`)) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1")) %>%
        ## while we sort out duplicates we will simply take the first of each ID number
        group_by(sample_id) %>% top_n(1) %>% ungroup() %>%
        mutate(date_of_test_lps_iga_inaba=str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(1) %>%
                   unlist %>%
                   as.Date(format="%d%m%y"),
               test_plate_lps_iga_inaba = str_split(string=Date,
                                                    pattern="-",2) %>%
                   map(2) %>%
                   unlist) %>%
        rename(lps_iga_inaba = titer) %>%
        select(sample_id,lps_iga_inaba,date_of_test_lps_iga_inaba,
               test_plate_lps_iga_inaba)

    ## vibriocidal Inaba
    vibinab <- read_xlsx(here("data","Natoinal serosurvey Vibriocidal assay Inaba 17 Oct 2018.xlsx"),
                         range="A2:H2934",
                         col_types = c("text","numeric",'text','text','text','numeric','numeric','text'),
                         na = c("N/S",NA)) %>%
        select(sample_id = `Sample ID`,
               vibinab =`Vibriocidal titer`,
               date_of_test_vibinab=Date,
               test_plate_vibinab="...3") %>%
        tidyr::fill(test_plate_vibinab) %>%
        mutate(date_of_test_vibinab = as.Date(date_of_test_vibinab,
                                              format="%d.%m.%y")) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1"),
               sample_id = recode(sample_id,`SAMPLE_ID_2-B`="SAMPLE_ID_2"))

    ## vibriocidal Ogawa
    vibogaw <- read_xlsx(here("data","Natoinal serosurvey Vibriocidal assay Ogawa 16 Oct 2018.xlsx"),
                         range="A2:H2934",
                         col_types = c("text","numeric",'text','text','text','numeric','numeric','text'),
                         na=c("N/S",NA)) %>%
        select(sample_id = `Sample ID`,
               vibogaw =`Vibriocidal Titer`,
               date_of_test_vibogaw=Date,
               test_plate_vibogaw=`# plate`) %>%
        tidyr::fill(test_plate_vibogaw) %>%
        mutate(date_of_test_vibogaw = as.Date(date_of_test_vibogaw,
                                              format="%d.%m.%y")) %>%
        filter(sample_id != "SAMPLE_ID_1-C",sample_id != "SAMPLE_ID_2-C") %>%
        mutate(sample_id = recode(sample_id,`SAMPLE_ID_1-A`="SAMPLE_ID_1"),
               sample_id = recode(sample_id,`SAMPLE_ID_2-B`="SAMPLE_ID_2"))

    ## Bringing it all together
    study_dat <- full_join(vibinab,vibogaw) %>%
        full_join(lps_iga_inaba) %>%
        full_join(lps_iga_ogawa) %>%
        full_join(lps_igg_inaba) %>%
        full_join(lps_igg_ogawa) %>%
        full_join(ctb_iga) %>%
        full_join(ctb_igg) %>%
        mutate(plpsa_s = pmax(lps_iga_ogawa,lps_iga_inaba),
               pglsp_s = pmax(lps_igg_ogawa,lps_igg_inaba),
               maxvib = pmax(vibinab,vibogaw))

    ## now loading the raw OD data from vibriocidal experiments
    ## we will use these data to re-derive titers for all
    vibs <- load_vibriocidal_data(study_dat=study_dat,
                                  re_read_og_files = FALSE)
    # vibs_interp <- add_nplr_vibs(vibs)
    # vibs_interp <- vibs_interp %>%
    #   mutate(titer_top_interp = ifelse(titer_og == 5,5,titer_top_interp)) %>%
    #   mutate(titer_top_interp = ifelse(titer_top_interp > 20000,titer_interp_simp,titer_top_interp))

    vibs_tmp <- vibs %>%
      select(sample_id,serotype,titer_c1,titer_og) %>%
        distinct()

    ## read in meta-data
    col_dates <- readr::read_csv(here("data","Y2Dates.csv")) %>%
        mutate(sample_id = as.character(sample_id),
               date_of_collection  = date_of_collection %>%
                   as.Date(.,format="%d/%m/%Y"))

    sample_meta_dat <- readr::read_csv(here("data","CombinedDataY1Y2.csv")) %>%
        rename(sample_id = indid, community_id = CommunityID) %>%
        mutate(sample_id = as.character(sample_id),
               year = substr(indidyr,nchar(indidyr),nchar(indidyr))) %>%
        filter(year ==2) %>%
        mutate(age = ifelse(age == 999,NA,age)) %>%
        left_join(col_dates)

    ## adding titers to meta-data
    rc <- left_join(study_dat,
                     vibs_tmp %>%
                        filter(serotype == "inaba") %>%
                        select(-serotype) %>%
                        set_names(c("sample_id",
                                    "titer_c1_inaba",
                                    "titer_og_inaba"#,
                                    # "titer_interp_simp_inaba",
                                    # "titer_top_interp_inaba"
                        ))) %>%
        left_join(vibs_tmp %>%
                      filter(serotype == "ogawa") %>%
                      select(-serotype) %>%
                      set_names(c("sample_id",
                                  "titer_c1_ogawa",
                                  "titer_og_ogawa"#,
                                  # "titer_interp_simp_ogawa",
                                  # "titer_top_interp_ogawa"
                      ))) %>%
        #rc <- rc %>% mutate(titer_interp_simp_max = pmax(titer_interp_simp_inaba,titer_interp_simp_ogawa),
        #                      titer_top_interp_max = pmax(titer_top_interp_inaba,titer_top_interp_ogawa))

        ## we will go with the cleaned vibriocidal data
        rename(vibogaw = titer_c1_ogawa, vibogaw_og = vibogaw,
               vibinab = titer_c1_inaba, vibinab_og = vibinab) %>%
    ## just a check
        ## rc$sample_id[!rc$sample_id %in% sample_meta_dat$sample_id]
        left_join(sample_meta_dat) %>%
        mutate(adm2name = str_replace_all(Zila_name," Zila","")) %>%
        fix_adm2names() %>%
        ## a little recoding
        mutate(villid = villid %>% as.factor,
               hhid = hhid %>% as.factor,
               ageSumm = ageSumm %>% factor,
               incomeSumm = incomeSumm %>% factor(levels=1:4),
               educationSumm = educationSumm %>% factor(levels=1:4),
               sex = ifelse(sex==2,"Male","Female") %>% as.factor,
               travel = travel %>% factor(levels=c(4,1,2,3)),
               vib_320 = ifelse(maxvib>=320,1,0),
               vib_640 = ifelse(maxvib>=640,1,0))

    ## adding BTM coordinate data
    sampling_locs <- sf::st_as_sf(rc %>%
                                      filter(!is.na(lat)) %>%
                                      dplyr::select(lon,lat),
                                  coords = c("lon","lat"),
                                  crs="+proj=longlat +datum=WGS84 +no_defs")
    sampling_locs_btm <- sampling_locs %>%
        transform_to_btm %>%
        st_coordinates

    rc$x_btm <- NA
    rc$y_btm <- NA
    rc$x_btm[which(!is.na(rc$lat))] <- sampling_locs_btm[,"X"]
    rc$y_btm[which(!is.na(rc$lat))] <- sampling_locs_btm[,"Y"]


    print("saving generated_data/full_data.rds")
    saveRDS(rc,here("generated_data","full_data.rds"))

    return(rc)
}

load_bangladesh_app <- function(){
    library(lubridate)
    library(sf)
    library(stringr)
    library(purrr)

    ## bring in serolocations
    sero_locs <- readr::read_csv(here("data","sampled_communities.csv")) %>%
        mutate(NAME_2=str_replace_all(Zila_name,pattern = "\\s+Zila",""),                                                           NAME_1=str_replace_all(Division_name,pattern = "\\s+Division",""))

    ## convert into spatial object
    map2(sero_locs$Lat, sero_locs$Long, ~st_point(c(.x, .y))) %>%
        st_sfc(crs=4326) %>%
        st_sf(sero_locs[,-c(13:14)], .) -> sero_locs_sf


    ## locations of Qadri surveillance sites
    chol_surv <- read_xls(here("data","Distribution of cholera - FBIS sites- Dec2014-Dec2015.xls"),skip=2) %>%
        filter(!`Surveillance Sites and months` %in% c("Total","Grand Total"))

    loc_names <- chol_surv %>%  filter(str_detect(`Surveillance Sites and months`,"\\(")) %>%  rename(site=`Surveillance Sites and months`) %>%
        mutate(site = str_remove_all(site,"\\([0-9]*\\)") %>% str_trim()) %>%
        dplyr::select(site) %>%
        unlist %>% rep(.,each=13)


    chol_surv <- chol_surv %>% filter(!str_detect(`Surveillance Sites and months`,"\\(")) %>% rename(month=`Surveillance Sites and months`) %>%
        mutate(site = loc_names)

    chol_surv <- chol_surv %>%
        mutate(month = str_replace_all(month,"April","Apr")) %>%
        mutate(month = str_replace_all(month,"June","Jun")) %>%
        mutate(month = str_replace_all(month,"July","Jul")) %>%
        mutate(month = parse_date_time(month,"%b-%y") %>% as.Date) %>%
        rename(diarrhea = `Line listing (Diarrhea)`,
               culture_done = `Culture (Vibrio cholerae) Done`,
               culture_pos_n = `Culture (Vibrio cholerae) (+)`,
               culture_pos_pct = `(%)`)

    chol_surv <- chol_surv %>% mutate(adm2name = recode(site,
                                                        "Narshindi" = "Narshingdi",
                                                        "Habiganj" = "Hobiganj",
                                                        "DMCH, Dhaka" = "Dhaka",
                                                        "UAMC&H, Dhaka" = "Dhaka",
                                                        "BITID, Chittagang" = "Chittagong"),
                                      adm3name = site,
                                      adm3name = recode(adm3name,
                                                        "Habiganj" = "Habiganj S.",
                                                        "Narshindi" = "Narsingdi S.",
                                                        "Cox's Bazar" = "Cox Bazar S",
                                                        "Naogaon" = "Naogaon S.",
                                                        "Patuakhali" = "Patuakhali S.",
                                                        "Thakurgaon" = "Thakurgaon S." ,
                                                        "Shatkhira" = "Satkhira S.",
                                                        "BITID, Chittagang" = "Sitakunda", ## this is an issue as there is not Chattagong thanas as adm3 so just putting this as placeholdser
                                                        "UAMC&H, Dhaka" = "Savar", ## need to verify
                                                        "DMCH, Dhaka" = "Tejgaon Circle Upazila"))


    # bit of a hack, once we get full codes we can modify
    zl_codes <- bind_rows(sero_locs %>% filter(NAME_2 %in% c(key_adm2,"Habiganj","Narsingdi")) %>% group_by(NAME_2) %>% summarize(zl[1]) %>% rename(adm2name = NAME_2,zl=`zl[1]`),
                          data.frame(adm2name = c("Patuakhali","Shatkhira","Narshingdi","Hobiganj"),zl=c(78,87,68,36)))

    chol_surv <- left_join(chol_surv,zl_codes)

    ## bring in census data so we can get the population data per zilla
    census <- readr::read_csv(here("data","BDcensus2011.csv"))
    pop_by_zilla <- census %>% group_by(zl) %>% summarise(pop=sum(pop_both))

    ## adding population
    chol_surv <- left_join(chol_surv,pop_by_zilla)

    chol_surv_agg <- chol_surv %>% group_by(adm2name) %>% summarize(dia = sum(diarrhea),
                                                                    culture_done = sum(culture_done),
                                                                    culture_pos = sum(culture_pos_n),
                                                                    implied_cholera = dia * culture_pos/culture_done,
                                                                    implied_inc_per10k = 10000*implied_cholera/pop[1] ) %>%
        arrange(desc(implied_inc_per10k))

}

#' Estimate daily conditional probabilities of infection
#'
#' @param survey_dat data.frame, cross-sectional serology data
#' @param grid_dat list, each element is a monthly estimate of the probability of cholera given AWD for grid cells across Bangladesh
#' @param days_back numeric, time frame of interest for seropositivity
#'
#' @return matrix, rows represent the days since the survey was conducted and columns represent individual samples
#' @export
#'
#' @examples
estimate_daily_probs <- function(survey_dat,
                                 grid_dat,
                                 days_back){
    ## find the seasonality for each grid cell with a community
    grid_seasonality_dat <- c()
    for(i in seq(length(grid_dat))){
        tmp <- grid_dat[[i]] %>%
            as_tibble() %>%
            select(grid_id,year_sick,month_sick,chol_est_case) %>%
            mutate(date_sick=ymd(paste(year_sick,month_sick,"01",sep="-")))
        grid_seasonality_dat <- bind_rows(grid_seasonality_dat, tmp)
    }

    ## AWD and cholera data by date and cell
    grid_awd <- grid_seasonality_dat %>%
        group_by(grid_id) %>%
        transmute(date_sick=date_sick,
                  last_date=ymd(paste(year_sick, month_sick,
                                      days_in_month(month_sick),sep="-")),
                  chol_est_case=chol_est_case,
                  chol_cum=cumsum(chol_est_case)) %>%
        bind_rows(tibble(grid_id=unique(grid_seasonality_dat$grid_id),
                         date_sick=min(grid_seasonality_dat$date_sick)-1,
                         last_date=min(grid_seasonality_dat$date_sick)-1,
                         chol_est_case=0,
                         chol_cum=0)) %>%
        arrange(grid_id,date_sick) %>%
        mutate(days_sick=as.numeric(last_date-min(date_sick))) %>%
        ungroup()

    ## fit a monotonically increasing model for cumulative proportion of cholera for each cell
    ## then make predictions for each cell and bind into single dataframe
    all_grid_dat <- foreach(i=seq(length(unique(grid_awd$grid_id))),.combine="rbind") %dopar%{
        tmp <- filter(grid_awd, grid_id==unique(grid_awd$grid_id)[i])
        scam_fit <- scam(chol_cum~s(days_sick, bs="mpi", k=nrow(tmp), m=5), data=tmp)
        pred_dat <- data.frame(grid_id=unique(grid_awd$grid_id)[i],
                               days_sick=seq(0,max(grid_awd$days_sick))) %>%
            mutate(last_date=min(grid_awd$date_sick)+days_sick)
        pred_dat$chol_cum_pred <- predict(scam_fit, pred_dat)
        pred_dat$chol_pred <- 0
        pred_dat$chol_pred[2:nrow(pred_dat)] <- diff(pred_dat$chol_cum_pred)*days_in_month(month(pred_dat$last_date[2:nrow(pred_dat)]))
        pred_dat[-1,]
    }

    ## based on the community for each individual find the probability of being infected
    ## on each day in that cell
    prob_grid_mtx <- matrix(NA, nrow=days_back, ncol=nrow(survey_dat))
    for(i in seq(nrow(survey_dat))){
        grid_tmp <- filter(all_grid_dat, grid_id==survey_dat$community_cell[i])
        date_idx <- which(grid_tmp$last_date==survey_dat$date_of_collection[i])
        tmp <- grid_tmp$chol_pred[(date_idx-1):(date_idx-days_back)]
        tmp_sum <- sum(tmp)
        prob_grid_mtx[,i] <- tmp/tmp_sum
    }
    return(prob_grid_mtx)
}

#' Estimate cholera incidence at the grid level using sentinel data
#'
#' @param sentinel_data data.frame, sentinel surveillance data of cholera and acute watery diarrhea (AWD)
#' @param survey_data data.frame, cross-sectional serosurvey data used for identifying grid cells of interest and collection dates
#' @param my_grid sf data.frame, gridded cell map of Bangladesh
#' @param bgd2 sf data.frame, map of Bangladesh with adm2 regions
#' @param bgd_pop raster layer, WorldPop data for Bangladesh
#' @param awd_data data.frame, AWD cases by adm2 region in Bangladesh
#' @param save_results logical, do you want to save the output?
#' @param ... other options for brm()
#'
#' @return
#' @export
#'
#' @examples
est_cholera_sentinel <- function(sentinel_data,
                                 survey_data,
                                 my_grid,
                                 bgd2,
                                 bgd_pop,
                                 awd_data,
                                 save_results=T,
                                 ...){
    ## conflicted package breaks code during brms unless all conflicts
    ## resolved first
    # conflict_scout()
    conflict_prefer("accumulate", "foreach",quiet=T)
    conflict_prefer("ar", "brms",quiet=T)
    conflict_prefer("collapse", "nlme",quiet=T)
    conflict_prefer("combine", "gridExtra",quiet=T)
    conflict_prefer("extract", "rstan",quiet=T)
    conflict_prefer("intersect", "base",quiet=T)
    conflict_prefer("lag", "stats",quiet=T)
    conflict_prefer("lowess", "stats",quiet=T)
    conflict_prefer("margin", "ggplot2",quiet=T)
    conflict_prefer("Position", "base",quiet=T)
    conflict_prefer("s", "brms",quiet=T)
    conflict_prefer("setdiff", "base",quiet=T)
    conflict_prefer("t2", "brms",quiet=T)
    conflict_prefer("union", "base",quiet=T)
    conflict_prefer("when", "foreach",quiet=T)

    ## model monthly likelihood of cholera given awd
    chol_month_dat <- sentinel_data %>%
        mutate(year_sick=as.factor(year_sick)) %>%
        filter(!is.na(count)) %>%
        group_by(adm3name,year_sick,month_sick,X,Y) %>%
        summarize(chol=sum(count),
                  awd=sum(total)) %>%
        mutate(prop_cholera=chol/awd,
               logit_prop=log(prop_cholera/(1-prop_cholera)))

    chol_locs <- st_as_sf(chol_month_dat %>%
                              ungroup() %>%
                              dplyr::select(X,Y),
                          coords = c("X","Y"),crs="+proj=longlat +datum=WGS84 +no_defs") %>%
        transform_to_btm() %>%
        st_coordinates

    chol_month_dat$X <- chol_locs[,"X"]
    chol_month_dat$Y <- chol_locs[,"Y"]

    ## model spatial GAM using Bayesian GLM
    cholera_brm <- brm(chol|trials(awd)~s(X,Y,bs='tp')+s(month_sick,bs='cc')+
                           s(year_sick,bs='re'),
                       family="binomial", data=chol_month_dat, ...)

    ## predict monthly likelihood of cholera given awd
    ## create data frame to predict estimates
    my_grid <- filter(my_grid, grid_id %in% survey_data$community_cell)
    my_grid_xy <- my_grid %>%
        st_centroid() %>%
        st_coordinates %>%
        data.frame
    grid_coords <- as.data.frame(my_grid) %>%
        mutate(X=my_grid_xy[,"X"],
               Y=my_grid_xy[,"Y"]) %>%
        select(grid_id,X,Y)

    pred_dat <- crossing(year_sick=chol_month_dat$year_sick,
                         month_sick=chol_month_dat$month_sick,
                         grid_coords,
                         awd=1) %>%
        mutate(date_sick=ymd(paste(year_sick,month_sick,"01",sep="-")))

    #draw params from posterior and then make predictions
    brm_preds <- predict(cholera_brm,newdata=pred_dat,type="response")
    pred_dat$est <- brm_preds[,1]

    ## extract population by adm2 region
    pop_dist <- exact_extract(bgd_pop, bgd2, 'sum')
    # bgd_pop_v <- velox(bgd_pop)
    # pop_dist <- bgd_pop_v$extract(sp=bgd2, fun=function(x) sum(x, na.rm=T))
    pop_dist_df <- data.frame(OBJECTID=seq(length(pop_dist)),
                              pop=pop_dist)

    ## make maps for each month with estimated awd and cholera in each grid cell
    ## construct one map for each month
    yr_mth_dat <- pred_dat %>%
        select(year_sick, month_sick) %>%
        distinct()
    map_month <- foreach(i=seq(nrow(yr_mth_dat))) %dopar% {
        print(paste("start: map",i,"of",nrow(yr_mth_dat)))
        ## look only at observations for one month
        tmp_pred <- filter(pred_dat,
                           year_sick==yr_mth_dat$year_sick[i],
                           month_sick==yr_mth_dat$month_sick[i])
        ## find the awd_incidence for that month
        awd_map <- left_join(bgd2,
                             awd_data %>%
                                 filter(month_sick==yr_mth_dat$month_sick[i],
                                        year_sick==yr_mth_dat$year_sick[i]) %>%
                                 group_by(adm2name) %>%
                                 summarize(awds = sum(awd,na.rm=T))) %>%
            left_join(pop_dist_df, by="OBJECTID") %>%
            mutate(awd_incidence=awds/pop)
        ## make a raster grid of awd incidence
        awd_raster <- raster(awd_map, res = 5000)
        awd_fasterized <- fasterize(awd_map, awd_raster,
                                    field = "awd_incidence", fun="sum") %>%
            transform_to_btm()
        ## extract awd incidence across my_grid
        awd_month <- exact_extract(awd_fasterized, my_grid, 'mean')
        # bgd_awd <-  velox(awd_fasterized)
        # awd_month <- bgd_awd$extract(sp=my_grid,
        #                              fun=function(x)
        #                                  ifelse(all(is.na(x)),0,mean(x,na.rm=T)))
        ## make cholera incidence estimates for my_grid
        tmp <- my_grid %>%
            mutate(year_sick=tmp_pred$year_sick,
                   month_sick=tmp_pred$month_sick,
                   chol_surv_pred=tmp_pred$est,
                   chol_est_case=awd_month*chol_surv_pred*pop,
                   awd_incidence_rate=awd_month*1000,
                   chol_incidence_rate=awd_incidence_rate*chol_surv_pred,
                   rr=chol_incidence_rate/mean(chol_incidence_rate))
        return(tmp)
    }

    ## save rds
    if(save_results)
        saveRDS(map_month, paste0("generated_data/gam-cholera-bgd-month-",Sys.Date(),".rds"))
    return(map_month)
}

