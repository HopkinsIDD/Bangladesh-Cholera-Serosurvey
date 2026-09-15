## LOOCV analysis
##

# devtools::install_github('timcdlucas/INLAutils')
library(pacman)
source("source/utils.R")
reload_source()
p_load(INLA,flextable,officer,INLAutils)

## load location data for each community
survey_dat <- readRDS(here("generated_data","dat_with_raster_covs.rds")) %>%
    filter(community_cell==grid_id) %>%
    select(grid_id, community_id,CommLong,CommLat,pop,logpop,friction,
           travel_time,logtravel_time,poverty,
           distance_to_water,altitude,sexF, age_0_9, age_10_19) %>%
    distinct() %>%
    arrange(community_id)

## get community centroid locations
com_locs_btm <- survey_dat %>%
    select(CommLong,CommLat,community_id) %>%
    distinct() %>%
    st_as_sf(coords = c("CommLong","CommLat"),
             crs="+proj=longlat +datum=WGS84 +no_defs") %>%
    transform_to_btm() %>% st_coordinates

## load estimates
adj_rf_ests <- readRDS("generated_data/results-rf-365-youden.rds")$loc_adjusted_incidence %>%
    group_by(community_id, obs) %>%
    summarize(y=floor(mean(poststrat_reps*obs)))

## combine covariate and result data
comm_dat <- survey_dat %>%
    left_join(adj_rf_ests)

## prepare data for INLA step estimation
inla_prep <- prep_inla_est(eqn="logpop + sexF + age_0_9 + age_10_19 + distance_to_water + poverty + friction + travel_time + logtravel_time + altitude",
                           dat=comm_dat,
                           loc_btm=com_locs_btm)

## run backwards selection to find the best model
backward_result <- INLAstep(fam1 = "binomial",
                            dataf=as.data.frame(comm_dat),
                            spde=inla_prep$spde,
                            in_stack=inla_prep$stack,
                            invariant="0 + Intercept +f(spatial.field, model = spde)",
                            direction="backwards",
                            include=which(names(comm_dat) %in%
                                              colnames(inla_prep$model_mtx)),
                            y="y",
                            y2="y",
                            Ntrials=inla_prep$stack_data$Ntrials,
                            num.threads=max(1, parallel::detectCores()-2),
                            thresh=0)
backward_result$best_formula
backward_result$best_model$waic$waic

null_result <- inla(y~0 + Intercept + f(spatial.field, model = inla_prep$spde),
                    data=inla_prep$stack_data,
                    family="binomial",
                    verbose=FALSE,
                    control.predictor=list(
                        A=inla.stack.A(inla_prep$stack),
                        compute=TRUE),
                    control.compute=list(waic=TRUE,config=TRUE),
                    Ntrials=inla_prep$stack_data$Ntrials,
                    num.threads=max(1, parallel::detectCores()-2))
null_result$waic$waic

full_result <- inla(y~0 + Intercept + logpop + sexF + age_0_9 + age_10_19 + distance_to_water + poverty + friction + travel_time + logtravel_time + altitude + f(spatial.field, model = inla_prep$spde),
                    data=inla_prep$stack_data,
                    family="binomial",
                    verbose=FALSE,
                    control.predictor=list(
                        A=inla.stack.A(inla_prep$stack),
                        compute=TRUE),
                    control.compute=list(waic=TRUE,config=TRUE),
                    Ntrials=inla_prep$stack_data$Ntrials,
                    num.threads=max(1, parallel::detectCores()-2))
full_result$waic$waic

## leave one grid-id out of the fitting and predict on it
## save just the mean quantiles for each
## will do with covariates
spde <- inla_prep$spde
rc = NULL
i=1
for(gc in comm_dat$community_id){
    print(i)

    ## remove data from left out community
    tmp_dat <- comm_dat %>%
        mutate(y = ifelse(community_id == gc,NA,y),
               obs = ifelse(community_id == gc,NA,obs))

    ## update the stack
    tmp_est_stack <- inla.stack(
        data=list(y=tmp_dat[,'y'] %>% unlist,
                  Ntrials=tmp_dat[,'obs'] %>% unlist), ## response
        A=list(inla_prep$A_est_multi,1), ## projector matrix
        effects=list(inla_prep$s_index,
                     data.frame(Intercept=rep(1,nrow(tmp_dat)),
                                inla_prep$model_mtx)),
        tag="tmp")

    ## update the stack data
    tmp_est_stack_data <- inla.stack.data(tmp_est_stack)

    ## fit best model
    test <- inla(as.formula(backward_result$best_formula),
                 data=tmp_est_stack_data,
                 family="binomial",
                 verbose=FALSE,
                 control.predictor=list(
                     A=inla.stack.A(tmp_est_stack),link=1,
                     compute=TRUE),
                 control.compute=list(dic=FALSE,config=TRUE),
                 Ntrials=tmp_est_stack_data$Ntrials)

    ## fit null model
    test_null <- inla(y ~ -1 + Intercept + f(spatial.field, model = inla_prep$spde),
                    data=tmp_est_stack_data,
                    family="binomial",
                    verbose=FALSE,
                    control.predictor=list(
                        A=inla.stack.A(tmp_est_stack),link=1,
                        compute=TRUE),
                    control.compute=list(dic=FALSE,config=TRUE),
                    Ntrials=tmp_est_stack_data$Ntrials)

    ## get estimates from full model
    summary_hold_out = test$summary.fitted.values[which(is.na(tmp_dat$y)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]

    ## get estimates from null model
    null_hold_out = test_null$summary.fitted.values[which(is.na(tmp_dat$y)), c("mean", "sd","0.025quant","0.5quant","0.975quant")]

    ## put results into data
    rc = rc %>%
        bind_rows(
            data.frame(left_out_cell = gc,
                       full_mean = summary_hold_out$mean,
                       full_sd = summary_hold_out$sd,
                       full_median = summary_hold_out$`0.5quant`,
                       full_lb = summary_hold_out$`0.025quant`,
                       full_ub = summary_hold_out$`0.975quant`,
                       null_mean = null_hold_out$mean,
                       null_sd = null_hold_out$sd,
                       null_median = null_hold_out$`0.5quant`,
                       null_lb = null_hold_out$`0.025quant`,
                       null_ub = null_hold_out$`0.975quant`,
                       naive_mean = comm_dat %>%
                           filter(community_id!=gc) %>%
                           summarize(mean(y/obs)) %>% unlist,
                       naive_median = comm_dat %>%
                           filter(community_id!=gc) %>%
                           summarize(median(y/obs)) %>% unlist,
                       truth=comm_dat %>%
                           filter(community_id==gc) %>%
                           transmute(est=y/obs) %>% unlist
            ))
    i <- i+1
}

## mean absolute error
model_mae <- (rc$full_mean - rc$truth) %>% abs %>% mean
null_mae <- (rc$null_mean - rc$truth) %>% abs %>% mean
naive_mae <- (rc$naive_mean - rc$truth) %>% abs %>% mean

## relative MAE
model_mae/naive_mae
model_mae/null_mae
null_mae/naive_mae

model_ae <- abs(rc$full_mean-rc$truth)
null_ae <- abs(rc$null_mean-rc$truth)
naive_ae <-abs(rc$naive_mean-rc$truth)
mean(model_ae < naive_ae)
mean(model_ae < null_ae)
mean(null_ae < naive_ae)

## bias
mean(rc$full_mean - rc$truth)
mean(rc$null_mean - rc$truth)

ggplot(rc, aes(x=truth)) +
    geom_abline(intercept = 0, slope = 1,col='grey70') +
    #  geom_point(aes(y=X0.5quant),col='orange') +
    geom_point(aes(y=full_mean),col='steelblue') +
    geom_smooth(aes(y=full_mean),method="lm",col='steelblue',lty=2,se=F)+
    geom_point(aes(y=null_mean),col='red3') +
    geom_smooth(aes(y=null_mean),method="lm",col='red3',lty=2,se=F)+
    #    geom_smooth(aes(y=X0.5quant),method="lm",col='orange')+
    scale_x_log10("proportion seropositive",
                  breaks=c(0.05, 0.1, 0.2, 0.4, 0.8),
                  limits=c(0.05, 0.8)) +
    scale_y_log10("predicted proportion seropositive",
                  breaks=c(0.05, 0.1, 0.2, 0.4, 0.8),
                  limits=c(0.05, 0.8)) +
    theme_bw()

saveRDS(rc, "generated_data/loocv-results.rds")
