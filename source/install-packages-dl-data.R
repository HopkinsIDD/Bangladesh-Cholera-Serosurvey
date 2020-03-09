## install packages and download data
## by Stephen Lauer, July 2019

options(echo=TRUE)
Sys.time()

if(!require(pacman)){
    install.packages("pacman", repos='https://cloud.r-project.org')
    library(pacman)
}

source("source/utils.R")
reload_source()

if(!require(INLA)){
    install.packages("INLA", repos=c('https://cloud.r-project.org', INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
    library(INLA)
}

## download data and adjust age by days since infection
bgd_cohort <- read_csv("https://raw.githubusercontent.com/scottyaz/cross-sectional-cholera-sero/master/data/public_bangladesh_data.csv") %>%
    mutate(age=round(age+lbday2/365.25),
           sex=ifelse(sex==1,"Male","Female") %>% as.factor) %>%
    rename(lps_iga=plpsa_s, ctb_iga=ptctxa_s, lps_igg=pglsp_s, ctb_igg=pgctxb_s,
           lps_igm=plpsm_s, ctb_igm=pmctxb_s)
## save into data folder
saveRDS(bgd_cohort, "data/bgd-cohort.rds")

## install other packages for generating manuscript and supplement
p_load(bookdown)
