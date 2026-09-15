## script for generating smoothed maps for each posterior draw of the
## post-stratified community seroincidence estimates
## usage: Rscript source/sampling_uncert_gridmap.R <rep 1-1000> [main|season]
##   main   - primary estimates (maps saved to generated_data/boot-maps/)
##   season - seasonality sensitivity analysis (generated_data/boot-maps-season/)

options(echo=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)
slurm_rep <- as.numeric(args[1])
variant <- ifelse(length(args) > 1, args[2], "main")
stopifnot(variant %in% c("main", "season"))
map_dir <- ifelse(variant == "main", "boot-maps", "boot-maps-season")
est_name <- ifelse(variant == "main", "loc_adjusted_incidence", "loc_season_incidence")
dir.create(here::here("generated_data", map_dir), showWarnings=FALSE)

library(pacman)
p_load(INLA,conflicted,here)
source(here("source","utils.R"))
reload_source()

time_in <- Sys.time()
print(paste0("Rep: ", slurm_rep, "; ", Sys.time()))
if(find_recent_file(name_start=paste0("boot_map_",slurm_rep,".rds"),
                    path=paste0(here("generated_data", map_dir), "/")) %>% is.na()){
    ## bring in data, which comes with bootstraps
    ## adding an index for replicate to make it easier
    adj_rf_ests <- readRDS(here("generated_data","results-rf-365-youden.rds"))[[est_name]]

    ## get sampled community information
    survey_dat <- readRDS(here("generated_data","dat_with_raster_covs.rds")) %>%
        filter(community_cell==grid_id) %>%
        select(grid_id, community_id,CommLong,CommLat,pop,logpop,friction,
               travel_time,logtravel_time,poverty,
               distance_to_water,altitude,sexF, age_0_9, age_10_19) %>%
        distinct() %>%
        arrange(community_id)

    ## get community centroid locations
    com_locs_btm <- survey_dat %>%
        select(CommLong,CommLat) %>%
        distinct() %>%
        st_as_sf(coords = c("CommLong","CommLat"),
                 crs="+proj=longlat +datum=WGS84 +no_defs") %>%
        transform_to_btm() %>% st_coordinates

    ## choose one set of simulated estimates
    set.seed(1)
    sim_idx <- sample(unique(adj_rf_ests$rep_ind), 1000)[slurm_rep]
    dat_simp <- survey_dat %>%
        ## joining raster data to each replicate
        inner_join(adj_rf_ests %>%
                       filter(rep_ind==sim_idx)) %>%
        mutate(y=floor(poststrat_reps * obs))


    ## bring in grid
    my_grid <- readRDS(here("data","grid_with_covs.rds")) %>%
        arrange(grid_id)

    pred_prep <- prep_inla_pred(eqn="sexF",
                                est_dat=dat_simp,
                                loc_btm=com_locs_btm,
                                pred_dat=my_grid)
    ## run the model
    inla_fit <- inla(y ~ -1 + Intercept + f(spatial.field, model=pred_prep$est_prep$spde),
                     data=pred_prep$est_and_pred_data,
                     family="binomial",
                     quantiles=c(0.025, 0.5, .975),
                     verbose=FALSE,
                     control.predictor=list(
                         A=inla.stack.A(pred_prep$est_and_pred_stack),
                         link=1,
                         compute=FALSE),
                     control.compute=list(dic=FALSE,
                                          config=TRUE,
                                          hyperpar=FALSE,
                                          mlik=FALSE),
                     Ntrials=pred_prep$est_and_pred_data$Ntrials)

        inla_preds <- inla.posterior.sample(result=inla_fit)
        ## get ids for each cell
        id.prd <- inla.stack.index(pred_prep$est_and_pred_stack, "pred")$data
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

    ## samples data and makes data.frame
    saveRDS(rc, file=here("generated_data", map_dir, paste0("boot_map_",slurm_rep,".rds")))
}
time_out <- Sys.time()
time_out-time_in
