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
options(mc.cores=n_chains)

## 1,250 post-warmup draws per chain (4 chains = 5,000 posterior draws)
n_iter <- 1500
n_warmup <- 250
p_adapt_delta <- 0.8
n_max_treedepth <- 20

## read in division codes
division_codes <- read.csv("data/division-zila-codes.csv", fileEncoding="UTF-8-BOM")

## read in serodata
survey_dat <- readRDS("generated_data/dat_with_raster_covs.rds") %>%
    filter(!is.na(pred_365)) %>%
    mutate(hh_id_min=dense_rank(hhid),
           comm_id_min=dense_rank(community_id),
           cell_id_min=dense_rank(community_cell),
           sample_id_min=dense_rank(sample_id),
           div_name=gsub(" Division", "", Division_name)) %>%
    left_join(division_codes %>%
                  select(-zl) %>%
                  distinct,
              by=c("div_name"="Division")) %>%
    arrange(sample_id_min)

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
t_ort <- poly(log(t), 3)
t_new_ort <- predict(t_ort, log(1:days_back))
z_cohort_pos <- cohort_dat$cohort_preds[cohort_dat[,paste0("inf_",days_back)]==1]
id_cohort_pos <- cohort_dat$id[cohort_dat[,paste0("inf_",days_back)]==1] %>% as.factor() %>% as.numeric()
N_cohort_neg <- sum(cohort_dat[,paste0("inf_",days_back)]==0)
z_cohort_neg <- cohort_dat$cohort_preds[cohort_dat[,paste0("inf_",days_back)]==0]
id_cohort_neg <- cohort_dat$id[cohort_dat[,paste0("inf_",days_back)]==0] %>% as.factor() %>% as.numeric()

## select the appropriate predictions from the serosurvey
if(diagnostic=="rf"){
    survey_dat$sero_preds <- survey_dat[,paste0("pred_",days_back),drop=T]
}
if(diagnostic=="vib"){
    survey_dat$sero_preds <- survey_dat$vib_320
}

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
    left_join(survey_dat %>%
                  select(community_id, community_cell, comm_id_min) %>%
                  distinct())

rm(census_dat)

## fit a model to find the overall seroincidence across all observations
sens_model <- stan_model("source/sens-spec.stan")
# beepr::beep(2)

time_in <- Sys.time()
sens_sample <- sampling(sens_model,
                        data=list(N_cohort_pos=N_cohort_pos,
                                  N_cohort_pos_ids=max(id_cohort_pos),
                                  id_cohort_pos=id_cohort_pos,
                                  z_cohort_pos=z_cohort_pos,
                                  T_max=days_back,
                                  K=ncol(t_ort),
                                  t_ort=as.matrix(t_ort),
                                  t_new_ort=t_new_ort,
                                  N_cohort_neg=N_cohort_neg,
                                  N_cohort_neg_ids=max(id_cohort_neg),
                                  id_cohort_neg=id_cohort_neg,
                                  z_cohort_neg=z_cohort_neg
                        ),
                        chains=n_chains,
                        iter=n_iter,
                        warmup=n_warmup,
                        control=list(adapt_delta=p_adapt_delta,
                                     max_treedepth=n_max_treedepth),
                        seed=stan_seed,
                        save_warmup=F)
Sys.time() - time_in
# beepr::beep(2)
# loo::extract_log_lik(sens_sample) %>% loo::loo()
# loo::extract_log_lik(sens_sample) %>% loo::loo() %>% plot
# print(sens_sample, pars=c("sens", "spec"))
# print(sens_sample, pars=c("alpha_spec", "alpha_sens", "beta", "sigma_sens_id", "sigma_spec_id"))

## sensitivity and specificity are relatively logit-normal, get params
sens_logit <- qlogis(extract(sens_sample, "sens")[[1]])
sens_logit_mean <- mean(sens_logit)
sens_logit_sd <- sd(sens_logit)

spec_logit <- qlogis(extract(sens_sample, "spec")[[1]])
spec_logit_mean <- mean(spec_logit)
spec_logit_sd <- sd(spec_logit)

X <- model.matrix(sero_preds ~ age_cat + sex,
                  data=survey_dat %>%
                      mutate(age_cat=relevel(age_cat, "age_20_99")))[,-1]

## fit a model to find the overall seroincidence across all observations
overall_model <- stan_model("source/random-effects.stan")
# beepr::beep(2)

time_in <- Sys.time()
overall_est <- sampling(overall_model,
                        data=list(N_survey=nrow(survey_dat),
                                  N_comm=max(survey_dat$comm_id_min),
                                  N_hh=max(survey_dat$hh_id_min),
                                  comm=survey_dat$comm_id_min,
                                  hh=survey_dat$hh_id_min,
                                  comm_pop=sampled_comms$pop_both,
                                  z_survey=survey_dat$sero_preds,
                                  K=ncol(X),
                                  X=X,
                                  sens_logit_mean=sens_logit_mean,
                                  sens_logit_sd=sens_logit_sd,
                                  spec_logit_mean=spec_logit_mean,
                                  spec_logit_sd=spec_logit_sd
                        ),
                        pars=c("beta", "gamma_0", "alpha_0", "sigma_comm", "sigma_hh", "alpha_comm", "sens", "spec", "y_survey", "log_lik"),
                        chains=n_chains,
                        iter=n_iter,
                        warmup=n_warmup,
                        control=list(adapt_delta=p_adapt_delta,
                                     max_treedepth=n_max_treedepth),
                        seed=stan_seed,
                        save_warmup=F)
Sys.time() - time_in
# beepr::beep(2)
# loo::extract_log_lik(overall_est) %>% loo::loo()
# loo::extract_log_lik(overall_est) %>% loo::loo() %>% plot()
# print(overall_est, pars=c("sens", "spec"))
# print(overall_est, pars=c("beta","gamma_0", "alpha_0", "sigma_comm", "sigma_hh"))

stan_pars <- extract(overall_est)
rm(overall_est)
rm(overall_model)

loc_incidence <- stan_pars$y_survey %>%
    data.frame() %>%
    gather("sample_id_min", "estimate") %>%
    mutate(sample_id_min=gsub("X", "", sample_id_min) %>%
               as.character %>% as.numeric) %>%
    group_by(sample_id_min) %>%
    mutate(sim=1:n()) %>%
    left_join(survey_dat %>%
                  select(sample_id_min, community_id, comm_id_min, community_cell, sero_preds)) %>%
    mutate(p_est=estimate) %>%
    group_by(sim, community_id, comm_id_min, comm_id_min, grid_id=community_cell) %>%
    summarize(samp_inc=sum(sero_preds),
              obs=n(),
              samp_est=sum(estimate)) %>%
    group_by(rep_ind=sim, community_id, comm_id_min, grid_id, obs) %>%
    transmute(original_incidence=samp_inc/obs,
              adjusted_reps=samp_est/obs) %>%
    group_by(community_id, comm_id_min, grid_id, obs) %>%
    mutate(adjusted_incidence=median(adjusted_reps)) %>%
    ungroup()

## post-stratify community results
analysis_age_props <- readRDS("generated_data/analysis_age_props.rds") %>%
    filter(grid_id %in% loc_incidence$grid_id) %>%
    full_join(sampled_comms %>%
                  select(community_id, comm_id_min, grid_id=community_cell, pop_both))

p_comm_mat <- model.matrix(prop ~ age_cat + sex,
                           data=analysis_age_props %>%
                               mutate(age_cat=relevel(as.factor(age_grp),
                                                      "age_20_99"),
                                      sex=relevel(as.factor(sex), "f")))[,-1]

registerDoMC(n_chains)
time_in <- Sys.time()
poststrat_reps <- foreach(k=1:nrow(loc_incidence), .combine=c, .inorder=T)  %dopar% {
    i <- loc_incidence$comm_id_min[k]
    j <- loc_incidence$rep_ind[k]
    tmp_dat <- analysis_age_props[analysis_age_props$community_id == loc_incidence$community_id[k],]
    tmp_mat <- p_comm_mat[analysis_age_props$community_id == loc_incidence$community_id[k],]
    tmp_dat$p_est <- NA
    for(m in 1:nrow(tmp_dat)){
        tmp_dat$p_est[m] <- integrate(function(x)
            plogis(qnorm(x, stan_pars$alpha_comm[j,i] + stan_pars$beta[j,,drop=F] %*% t(tmp_mat[m,,drop=F]),
                         stan_pars$sigma_hh[j])), 0,1,
            stop.on.error=F)[[1]]
    }
    sum(tmp_dat$prop*tmp_dat$p_est)
}
Sys.time() - time_in

# beepr::beep(2)
loc_incidence$poststrat_reps <- poststrat_reps
rm(poststrat_reps)
loc_incidence <- loc_incidence %>%
    group_by(community_id) %>%
    mutate(poststrat_incidence=median(poststrat_reps)) %>%
    ungroup

loc_incidence %>%
    select(community_id, original_incidence, adjusted_incidence, poststrat_incidence) %>%
    distinct %>%
    pivot_longer(cols=-community_id,
                 names_to="type",
                 values_to="est") %>%
    ggplot(aes(x=est, y=community_id, color=type)) +
    geom_point(alpha=0.6) +
    scale_color_discrete("",
                         labels=c("adjusted", "original", "poststrat")) +
    theme_bw() +
    theme(legend.position="bottom")

seroincidence_dat <- loc_incidence %>%
    group_by(rep_ind) %>%
    summarize(adjusted_reps=mean(adjusted_reps),
              poststrat_reps=mean(poststrat_reps))

## add adjusted results to output list
output_list <- list(original_incidence=mean(survey_dat$sero_preds),
                    adjusted_incidence=median(seroincidence_dat$adjusted_reps),
                    adjusted_reps=seroincidence_dat$adjusted_reps,
                    adjusted_sensitivity=stan_pars$sens,
                    poststrat_incidence=median(seroincidence_dat$poststrat_reps),
                    poststrat_reps=seroincidence_dat$poststrat_reps,
                    sigma_comm=stan_pars$sigma_comm,
                    sigma_hh=stan_pars$sigma_hh,
                    sens_tv=extract(sens_sample, pars="sens_tv")[[1]],
                    specificity=stan_pars$spec,
                    loc_adjusted_incidence=loc_incidence,
                    seroincidence_dat=seroincidence_dat)
## save output
if(diagnostic=="rf")
    saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,"-",cutoff_type,".rds"))
if(diagnostic=="vib")
    saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,".rds"))

## the seasonality sensitivity analysis is only reported for the 365-day random forest
if(diagnostic=="rf" & days_back==365){
    ## remove objects to free memory
    rm(loc_incidence)
    rm(seroincidence_dat)
    rm(stan_pars)

    #####
    ##### make community estimates with seasonality
    #####
    ## get seasonality for each community
    comm_dat <- survey_dat %>%
        group_by(comm_id_min) %>%
        slice(1) %>%
        select(comm_id_min, daily_probs_100, daily_probs_200, daily_probs_365) %>%
        ungroup

    ## days_back x N_comm matrix of daily infection probabilities
    prob_grid_mtx <- sapply(comm_dat[[paste0("daily_probs_", days_back)]], as.numeric)
    ## extract time-varying sensitivity from Stan
    sens_tv <- extract(sens_sample, pars="sens_tv")[[1]]
    ## combine to get sensitivity by community
    comm_sens <- sens_tv %*% prob_grid_mtx
    comm_sens_mean <- apply(qlogis(comm_sens), 2, mean)
    comm_sens_sd <- apply(qlogis(comm_sens), 2, sd)
    rm(sens_tv)
    rm(comm_sens)

    ## fit a model to find the overall seroincidence across all observations
    seasonality_model <- stan_model("source/random-effects-seasonality.stan")
    # beepr::beep(2)

    time_in <- Sys.time()
    seasonality_est <- sampling(seasonality_model,
                                data=list(N_survey=nrow(survey_dat),
                                          N_comm=max(survey_dat$comm_id_min),
                                          N_hh=max(survey_dat$hh_id_min),
                                          comm=survey_dat$comm_id_min,
                                          hh=survey_dat$hh_id_min,
                                          comm_pop=sampled_comms$pop_both,
                                          z_survey=survey_dat$sero_preds,
                                          K=ncol(X),
                                          X=X,
                                          sens_logit_mean=comm_sens_mean,
                                          sens_logit_sd=comm_sens_sd,
                                          spec_logit_mean=spec_logit_mean,
                                          spec_logit_sd=spec_logit_sd
                                ),
                                pars=c("beta", "gamma_0", "alpha_0", "sigma_comm", "sigma_hh", "alpha_comm", "sens", "spec", "y_survey", "log_lik"),
                                chains=n_chains,
                                iter=n_iter,
                                warmup=n_warmup,
                                control=list(adapt_delta=p_adapt_delta,
                                             max_treedepth=n_max_treedepth),
                                seed=stan_seed,
                                save_warmup=F)
    Sys.time() - time_in
    # beepr::beep(2)
    # loo::extract_log_lik(seasonality_est) %>% loo::loo()
    # loo::extract_log_lik(seasonality_est) %>% loo::loo() %>% plot()
    # print(seasonality_est, pars=c("sens", "spec"))
    # print(seasonality_est, pars=c("beta","gamma_0", "alpha_0", "sigma_comm", "sigma_hh"))

    season_pars <- extract(seasonality_est)
    rm(seasonality_est)
    rm(seasonality_model)

    loc_incidence_season <- season_pars$y_survey %>%
        data.frame() %>%
        gather("sample_id_min", "estimate") %>%
        mutate(sample_id_min=gsub("X", "", sample_id_min) %>%
                   as.character %>% as.numeric) %>%
        group_by(sample_id_min) %>%
        mutate(sim=1:n()) %>%
        left_join(survey_dat %>%
                      select(sample_id_min, community_id, comm_id_min, community_cell, sero_preds)) %>%
        mutate(p_est=estimate) %>%
        group_by(sim, community_id, comm_id_min, comm_id_min, grid_id=community_cell) %>%
        summarize(samp_inc=sum(sero_preds),
                  obs=n(),
                  samp_est=sum(estimate)) %>%
        group_by(rep_ind=sim, community_id, comm_id_min, grid_id, obs) %>%
        transmute(original_incidence=samp_inc/obs,
                  adjusted_reps=samp_est/obs) %>%
        group_by(community_id, comm_id_min, grid_id, obs) %>%
        mutate(adjusted_incidence=median(adjusted_reps)) %>%
        ungroup()

    # ggplot() +
    #     geom_line(aes(x=loc_incidence$adjusted_incidence[1:70], y=loc_incidence$adjusted_incidence[1:70]), linetype="dashed") +
    #     geom_point(aes(x=loc_incidence$adjusted_incidence[1:70], y=loc_incidence_season$adjusted_incidence[1:70])) +
    #     scale_x_continuous("Original adjusted incidence") +
    #     scale_y_continuous("Adjusted incidence with seasonality") +
    #     theme_bw()

    ## post-stratify community results
    registerDoMC(n_chains)
    time_in <- Sys.time()
    poststrat_reps <- foreach(k=1:nrow(loc_incidence_season), .combine=c, .inorder=T)  %dopar% {
        i <- loc_incidence_season$comm_id_min[k]
        j <- loc_incidence_season$rep_ind[k]
        tmp_dat <- analysis_age_props[analysis_age_props$community_id == loc_incidence_season$community_id[k],]
        tmp_mat <- p_comm_mat[analysis_age_props$community_id == loc_incidence_season$community_id[k],]
        tmp_dat$p_est <- NA
        for(m in 1:nrow(tmp_dat)){
            tmp_dat$p_est[m] <- integrate(function(x)
                plogis(qnorm(x, season_pars$alpha_comm[j,i] + season_pars$beta[j,,drop=F] %*% t(tmp_mat[m,,drop=F]),
                             season_pars$sigma_hh[j])), 0,1,
                stop.on.error=F)[[1]]
        }
        sum(tmp_dat$prop*tmp_dat$p_est)
    }
    Sys.time() - time_in

    # beepr::beep(2)
    loc_incidence_season$poststrat_reps <- poststrat_reps
    rm(poststrat_reps)
    loc_incidence_season <- loc_incidence_season %>%
        group_by(community_id) %>%
        mutate(poststrat_incidence=median(poststrat_reps)) %>%
        ungroup

    ggplot() +
        geom_line(aes(x=output_list$loc_adjusted_incidence$poststrat_incidence[1:70],
                      y=output_list$loc_adjusted_incidence$poststrat_incidence[1:70]), linetype="dashed") +
        geom_point(aes(x=output_list$loc_adjusted_incidence$poststrat_incidence[1:70],
                       y=loc_incidence_season$poststrat_incidence[1:70])) +
        scale_x_continuous("Original post-stratified incidence") +
        scale_y_continuous("Post-stratified incidence with seasonality") +
        theme_bw()

    output_list[["loc_season_incidence"]] <- loc_incidence_season
    # rm(loc_incidence_season)
    # rm(season_pars)

    ## save output
    if(diagnostic=="rf")
        saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,"-",cutoff_type,".rds"))
    if(diagnostic=="vib")
        saveRDS(output_list, paste0("generated_data/results-",diagnostic,"-",days_back,".rds"))
}
