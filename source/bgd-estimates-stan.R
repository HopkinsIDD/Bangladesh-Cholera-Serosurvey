## make adjusted seropositivity estimates for Bangladesh, using Stan
## by Stephen Lauer, November 2019

options(echo=TRUE, warnings=TRUE, message=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)

days_back <- as.numeric(args[1]) # how recently was person infected to be positive
diagnostic <- as.character(args[2]) # using random forest diagnostic or vibriocidal standard? ("rf","vib")
cutoff_type <- as.character(args[3]) # use cutoff related to Youden statistic or maximum specificity ("youden", "spec")
n_chains <- as.numeric(args[4]) # number of chains to use for analysis
stan_seed <- as.numeric(args[5])

## load libraries and utilities
library(pacman)
p_load(here, foreign)
setwd(here())
source("source/utils.R")
reload_source()
set.seed(1)
# conflict_prefer("accumulate", "foreach",quiet=T)
# conflict_prefer("collapse", "nlme",quiet=T)
conflict_prefer("combine", "gridExtra",quiet=T)
conflict_prefer("extract", "rstan",quiet=T)
conflict_prefer("intersect", "base",quiet=T)
conflict_prefer("lag", "stats",quiet=T)
conflict_prefer("lowess", "stats",quiet=T)
conflict_prefer("margin", "ggplot2",quiet=T)
conflict_prefer("Position", "base",quiet=T)
# conflict_prefer("s", "brms",quiet=T)
conflict_prefer("setdiff", "base",quiet=T)
# conflict_prefer("t2", "brms",quiet=T)
conflict_prefer("union", "base",quiet=T)
# conflict_prefer("when", "foreach",quiet=T)
options(mc.cores=n_chains)

n_iter <- 1e2
n_warmup <- n_iter-1e1
p_adapt_delta <- 0.8
n_max_treedepth <- 10

## read in division codes
division_codes <- read.csv("data/division-zila-codes.csv")

## read in serodata
survey_dat <- readRDS("generated_data/serosurvey-with-rf-preds.rds") %>%
    left_join(readRDS("data/community-data.rds")) %>%
    filter(!is.na(pred_365)) %>%
    mutate(hh_id_min=dense_rank(hhid),
           comm_id_min=dense_rank(community_id),
           cell_id_min=dense_rank(community_cell),
           div_name=gsub(" Division", "", Division_name),
           vib_320=pmax(vibinab, vibogaw)>=320) %>%
    left_join(division_codes %>%
                  select(-zl) %>%
                  distinct,
              by=c("div_name"="Division"))

## load cohort_dat object
cohort_dat <- find_recent_file(paste("cohort-w-rf-preds",days_back,cutoff_type,sep="-"),
                               "generated_data/") %>%
    read_csv() %>%
    ungroup() %>%
    mutate(vib_320=pmax(vibinab, vibogaw)>=320)

## select the appropriate predictions from the cohort data
if(diagnostic=="rf"){
    cohort_dat$cohort_preds <- cohort_dat$rf_preds
}
if(diagnostic=="vib"){
    cohort_dat$cohort_preds <- cohort_dat$vib_320
}

N_cohort_pos <- sum(cohort_dat[,paste0("inf_",days_back)]==1)
t <- cohort_dat$lbday2[cohort_dat[,paste0("inf_",days_back)]==1]
z_cohort_pos <- cohort_dat$cohort_preds[cohort_dat[,paste0("inf_",days_back)]==1]
N_cohort_neg <- sum(cohort_dat[,paste0("inf_",days_back)]==0)
z_cohort_neg <- cohort_dat$cohort_preds[cohort_dat[,paste0("inf_",days_back)]==0]

## select the appropriate predictions from the serosurvey
if(diagnostic=="rf"){
    survey_dat$sero_preds <- survey_dat[,paste0("pred_",days_back),drop=T]
}
if(diagnostic=="vib"){
    survey_dat$sero_preds <- survey_dat$vib_320
}

## aggregate predictions to household level
hh_dat <- survey_dat %>%
    group_by(dv, community_cell, community_id, cell_id_min, pop, comm_id_min, hhid, hh_id_min) %>%
    summarise(hh_obs=n(),
              hh_inc=sum(sero_preds)) %>%
    ungroup() %>%
    arrange(hh_id_min)

## read in census data to find size of unsampled population
census_dat <- read.dta("data/census_2011_3.dta") %>%
    mutate(zl=as.numeric(zl),
           pop_both=ifelse(is.na(pop_both),1,pop_both),
           hh_total=as.numeric(hh_total)) %>%
    left_join(division_codes) %>%
    select(dv, zl, uz, un_wa, mz_mh, vill, rmo, community, pop_both, hh_total)

## for the sampled communities, find the size of the unsampled population
sampled_comms <- read_csv("data/sampled_communities.csv") %>%
    left_join(census_dat) %>%
    select(community_id=CommunityID, dv, zl, uz, un_wa, mz_mh, vill, rmo, community, pop_both, hh_total) %>%
    left_join(hh_dat %>%
                  group_by(community_id) %>%
                  summarize(sampled_hh=n(),
                            sampled_pop=sum(hh_obs),
                            sampled_inc=sum(hh_inc))
    ) %>%
    mutate(unsampled_hh=hh_total-sampled_hh,
           unsampled_pop=pop_both-sampled_pop,
           comm_id_min=dense_rank(community_id))
## separate out the unsampled communities
unsampled_comms <- anti_join(census_dat, sampled_comms) %>%
    mutate(community_id=100+row_number(),
           sampled_hh=0,
           sampled_pop=0,
           sampled_inc=0,
           dv=ifelse(is.na(dv), 1, dv))

rm(census_dat)

## fit a model to find the overall seroincidence across all observations
overall_model <- stan_model("source/overall-analysis.stan")

overall_est <- sampling(overall_model,
                        data=list(N_survey=nrow(hh_dat),
                                  N_comm=max(hh_dat$comm_id_min),
                                  comm=hh_dat$comm_id_min,
                                  comm_pop=sampled_comms$pop_both,
                                  hh_obs=hh_dat$hh_obs,
                                  z_survey=hh_dat$hh_inc,
                                  N_cohort_pos=N_cohort_pos,
                                  T_max=days_back,
                                  t=t,
                                  z_cohort_pos=z_cohort_pos,
                                  N_cohort_neg=N_cohort_neg,
                                  z_cohort_neg=z_cohort_neg,
                                  pop_unobs=sampled_comms$unsampled_pop,
                                  N_unsamp=nrow(unsampled_comms),
                                  pop_unsamp=unsampled_comms$pop_both
                        ),
                        pars=c("sens", "sens_tv", "spec", "y_survey", "y_unobs", "y_unsamp"),
                        chains=n_chains,
                        iter=n_iter,
                        warmup=n_warmup,
                        control=list(adapt_delta=p_adapt_delta,
                                     max_treedepth=n_max_treedepth),
                        seed=stan_seed,
                        save_warmup=F)

## file too large to save...
# saveRDS(overall_est, file=paste0("generated_data/overall-est-stan-", Sys.Date(), ".rds"))
#
# ## get likelihood information
# ll_w_div <- loo::extract_log_lik(overall_est)
#
# loo::loo_compare(list(loo(ll_w_div), loo(ll_no_div)))

## extract important bits of information
y_survey <- extract(overall_est, pars="y_survey")[[1]] %>%
    data.frame() %>%
    gather("hh_id_min", "estimate") %>%
    mutate(hh_id_min=gsub("X", "", hh_id_min) %>%
               as.character %>% as.numeric) %>%
    group_by(hh_id_min) %>%
    mutate(sim=1:n()) %>%
    left_join(hh_dat) %>%
    mutate(p_est=estimate/hh_obs)
y_unobs <- extract(overall_est, pars="y_unobs")[[1]] %>%
    data.frame() %>%
    gather("comm_id_min", "unobs_est") %>%
    mutate(comm_id_min=gsub("X", "", comm_id_min) %>%
               as.character %>% as.numeric) %>%
    group_by(comm_id_min) %>%
    mutate(sim=1:n()) %>%
    left_join(sampled_comms) %>%
    ungroup()

loc_incidence <- y_survey %>%
    group_by(sim, community_id, comm_id_min) %>%
    summarize(samp_inc=sum(hh_inc),
              obs=sum(hh_obs),
              samp_est=sum(estimate)) %>%
    # left_join(y_unobs) %>%
    group_by(rep_ind=sim, community_id, obs) %>%
    transmute(original_incidence=samp_inc/obs,
              adjusted_reps=samp_est/obs) %>%
    group_by(community_id) %>%
    mutate(adjusted_incidence=median(adjusted_reps)) %>%
    ungroup()

seroincidence_dat <- data.frame(unsampled_incidence=extract(overall_est, pars="y_unsamp")[[1]] %>%
                                    apply(1, sum),
                                unsampled_pop=sum(unsampled_comms$pop_both),
                                unobs_incidence=extract(overall_est, pars="y_unobs")[[1]] %>%
                                    apply(1, sum),
                                unobs_pop=sum(sampled_comms$unsampled_pop),
                                sampled_incidence=extract(overall_est, pars="y_survey")[[1]] %>%
                                    apply(1, sum),
                                sampled_pop=sum(sampled_comms$sampled_pop)) %>%
    mutate(p_survey=sampled_incidence/sampled_pop,
           total_incidence=unsampled_incidence+unobs_incidence+sampled_incidence,
           total_pop=unsampled_pop+unobs_pop+sampled_pop,
           p_overall=total_incidence/total_pop,
           sim=1:n())

## add adjusted results to output list
output_list <- list(original_incidence=mean(survey_dat$sero_preds),
                    adjusted_incidence=median(seroincidence_dat$p_overall),
                    adjusted_reps=seroincidence_dat$p_overall,
                    adjusted_sensitivity=extract(overall_est, pars="sens")[[1]],
                    sens_tv=extract(overall_est, pars="sens_tv")[[1]],
                    specificity=extract(overall_est, pars="spec")[[1]],
                    loc_adjusted_incidence=loc_incidence,
                    seroincidence_dat=seroincidence_dat)
## save output
if(diagnostic=="rf")
    saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,"-",cutoff_type,".rds"))
if(diagnostic=="vib")
    saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,".rds"))

if(diagnostic=="rf"){
    ## remove objects to free memory
    rm(overall_est)
    rm(y_survey)
    rm(y_unobs)
    rm(overall_model)
    rm(loc_incidence)
    rm(seroincidence_dat)

    ## fit a model to find the unadjusted seroincidence across all observations
    unadj_model <- stan_model("source/unadj-analysis.stan")

    unadj_est <- sampling(unadj_model,
                          data=list(N_survey=nrow(hh_dat),
                                    N_comm=max(hh_dat$comm_id_min),
                                    comm=hh_dat$comm_id_min,
                                    comm_pop=sampled_comms$pop_both,
                                    hh_obs=hh_dat$hh_obs,
                                    z_survey=hh_dat$hh_inc,
                                    pop_unobs=sampled_comms$unsampled_pop,
                                    N_unsamp=nrow(unsampled_comms),
                                    pop_unsamp=unsampled_comms$pop_both),
                          pars=c("y_survey", "pi_overall"),
                          chains=n_chains,
                          iter=n_iter,
                          warmup=n_warmup,
                          control=list(adapt_delta=p_adapt_delta,
                                       max_treedepth=n_max_treedepth),
                          seed=stan_seed,
                          save_warmup=F)

    ## extract important bits of information
    unadjusted_reps <- extract(unadj_est, pars="pi_overall")[[1]]

    loc_unadj_incidence <- extract(unadj_est, pars="y_survey")[[1]] %>%
        data.frame() %>%
        gather("hh_id_min", "estimate") %>%
        mutate(hh_id_min=gsub("X", "", hh_id_min) %>%
                   as.character %>% as.numeric) %>%
        group_by(hh_id_min) %>%
        mutate(rep_ind=1:n()) %>%
        left_join(hh_dat) %>%
        group_by(rep_ind, community_id, comm_id_min) %>%
        summarize(samp_inc=sum(hh_inc),
                  obs=sum(hh_obs),
                  samp_est=sum(estimate)) %>%
        group_by(rep_ind, community_id, obs) %>%
        transmute(original_incidence=samp_inc/obs,
                  unadjusted_reps=samp_est/obs)

    ## remove objects to free up memory
    rm(unadj_model)
    rm(unadj_est)

    ## add unadj results to output
    output_list[["unadjusted_reps"]] <- unadjusted_reps
    output_list[["loc_unadj_incidence"]] <- loc_unadj_incidence

    ## save output
    if(diagnostic=="rf")
        saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,"-",cutoff_type,".rds"))
    if(diagnostic=="vib")
        saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,".rds"))

    if(days_back==365){
        ## remove objects to free up memory
        rm(unadjusted_reps)
        rm(loc_unadj_incidence)

        ## find the seroincidence for each grid cell
        cell_dat <- hh_dat %>%
            mutate(pop=ceiling(pop)) %>%
            select(dv, community_cell, cell_id_min, pop) %>%
            distinct() %>%
            left_join(hh_dat %>%
                          group_by(community_cell) %>%
                          summarise(sampled_hh=n(),
                                    sampled_pop=sum(hh_obs),
                                    sampled_inc=sum(hh_inc))) %>%
            mutate(unsampled_pop=pop-sampled_pop) %>%
            arrange(cell_id_min)

        ## fit a model to find the seroincidence for each grid cell
        cell_model <- stan_model("source/cell-analysis.stan")
        cell_est <- sampling(cell_model,
                             data=list(N_survey=nrow(hh_dat),
                                       N_cell=max(hh_dat$cell_id_min),
                                       cell=hh_dat$cell_id_min,
                                       cell_pop=cell_dat$sampled_pop,
                                       hh_obs=hh_dat$hh_obs,
                                       z_survey=hh_dat$hh_inc,
                                       N_cohort_pos=N_cohort_pos,
                                       T_max=days_back,
                                       t=t,
                                       z_cohort_pos=z_cohort_pos,
                                       N_cohort_neg=N_cohort_neg,
                                       z_cohort_neg=z_cohort_neg),
                             pars=c("pi_cell"),
                             chains=n_chains,
                             iter=n_iter,
                             warmup=n_warmup,
                             control=list(adapt_delta=p_adapt_delta,
                                          max_treedepth=n_max_treedepth),
                             seed=stan_seed,
                             save_warmup=F)

        cell_preds <- extract(cell_est, pars="pi_cell")[[1]] %>%
            data.frame() %>%
            gather("cell_id_min", "adjusted_reps") %>%
            mutate(cell_id_min=gsub("X", "", cell_id_min) %>%
                       as.character %>% as.numeric) %>%
            group_by(cell_id_min) %>%
            mutate(rep_ind=1:n()) %>%
            left_join(cell_dat) %>%
            ungroup() %>%
            select(rep_ind, community_cell, obs=sampled_pop, adjusted_reps)

        ## remove objects to free up memory
        rm(cell_model)
        rm(cell_est)

        ## add cell results to output
        output_list[["cell_adjusted_incidence"]] <- cell_preds

        ## save output
        if(diagnostic=="rf")
            saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,"-",cutoff_type,".rds"))
        if(diagnostic=="vib")
            saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,".rds"))
    }
}
