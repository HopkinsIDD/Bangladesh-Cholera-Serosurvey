## install packages and download data
## by Stephen Lauer, July 2019

options(echo=TRUE, timeout=1000)
Sys.time()

repos <- c(CRAN="https://cloud.r-project.org")

if(!require(pacman)){
    install.packages("pacman", repos=repos)
    library(pacman)
}

if(!require(INLA)){
    install.packages("INLA",
                     repos=c(repos, INLA="https://inla.r-inla-download.org/R/stable"),
                     dep=TRUE)
    library(INLA)
}

## INLAutils (backwards model selection for the LOOCV analysis) is only on GitHub
if(!require(INLAutils)){
    p_load(remotes)
    remotes::install_github("timcdlucas/INLAutils", upgrade="never", dependencies=FALSE)
}

source("source/utils.R")
reload_source()

## download the cohort data used to train the random forest models
## (Azman et al. 2019, Science Translational Medicine) and adjust age by days
## since infection
bgd_cohort <- read_csv("https://raw.githubusercontent.com/scottyaz/cross-sectional-cholera-sero/master/data/public_bangladesh_data.csv") %>%
    mutate(age=round(age+lbday2/365.25),
           sex=ifelse(sex==1,"Male","Female") %>% as.factor)
## save into data folder
saveRDS(bgd_cohort, "data/bgd-cohort.rds")

## other packages for the INLA models, manuscript and supplement
p_load(foreign, broom, ggpubr, ggnewscale, kableExtra, flextable, officer,
       bookdown, rticles, sn, fasterize, exactextractr)
## captioner is no longer on CRAN
if(!require(captioner)){
    remotes::install_github("adletaw/captioner", upgrade="never")
}
