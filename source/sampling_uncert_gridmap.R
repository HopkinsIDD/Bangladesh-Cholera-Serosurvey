## script for generating smoothed maps for each bootstrapped dataset

options(echo=TRUE)
Sys.time()

args <- commandArgs(TRUE)
print(args)
slurm_rep <- as.numeric(args[1])

library(pacman)
p_load(INLA,conflicted,here)
setwd(here())
source("source/utils.R")
reload_source()

time_in <- Sys.time()
print(paste0("Rep: ", slurm_rep, "; ", Sys.time()))
if(find_recent_file(paste0("boot_map_",slurm_rep,".rds"),
                    here("generated_data", "boot-maps")) %>% is.na()){
    ## bring in data, which comes with bootstraps
    ## adding an index for replicate to make it easier
    adj_rf_ests <- readRDS("generated_data/results-rf-365-youden.rds")$cell_adjusted_incidence

    ## bring in grid
    my_grid <- readRDS("data/grid_with_covs.rds")

    full_dat <- readRDS("generated_data/serosurvey-with-rf-preds.rds") %>%
        filter(grid_id==community_cell) %>%
        left_join(readRDS("data/community-data.rds")) %>%
        select(community_cell,pop,logpop,distance_to_water,
               grid_easting,grid_northing) %>%
        distinct() %>%
        ## joining raster data to each replicate
        inner_join(adj_rf_ests)

    ## choose a posterior sample
    ## only using 1,000 out of 10,000 samples due to computation time
    ## the estimates for each grid cell do not change much after the first few
    ## hundred samples
    set.seed(1)
    sim_idx <- sample(unique(full_dat$rep_ind), 1000)[slurm_rep]
    dat_simp <- full_dat %>%
        filter(rep_ind==sim_idx) %>%
        mutate(n_tested=obs,
               n_pos_rf_cor=floor(adjusted_reps * obs))

    ## Setting up non-dynamic aspects of SPDE model
    grid_locs <- dat_simp %>%
        select(grid_easting,grid_northing,community_cell) %>%
        st_as_sf(coords=c("grid_easting","grid_northing"),
                 crs="+proj=tmerc +lat_0=0 +lon_0=90 +k=0.9996 +x_0=500000 +y_0=0 +a=6377276.345 +b=6356075.41314024 +towgs84=283.7,735.9,261.1,0,0,0,0 +units=m +no_defs") %>%
        st_coordinates

    ## going to make the convex hull fo estimation aound our household coordinates
    bnd <- inla.nonconvex.hull(grid_locs,convex=-0.1) ## inla mesh segment

    mesh <- inla.mesh.2d(boundary=bnd ,
                         loc=grid_locs,
                         offset=c(-0.05, -0.05),
                         cutoff=2000, # if households are less than 2km apart, builds only a single vertex
                         max.edge=c(30000,50000)
    )

    ## make a grid projector
    grid_cents <- my_grid %>%
        group_by(grid_id) %>%
        st_centroid() %>%
        st_coordinates %>%
        data.frame %>%
        mutate(n=100)

    grid_projector <- inla.mesh.projector(mesh,
                                          loc=grid_cents %>% as.matrix)
    A.pred <- grid_projector$proj$A ## get A matrix

    spde <- inla.spde2.matern(mesh=mesh, alpha=2) # alpha is Fractional operator order
    s_index <- inla.spde.make.index(name="spatial.field",
                                    n.spde=spde$n.spde)

    ## samples data and makes data.frame
    set.seed(slurm_rep)
    boot_map <- make_boot_map(dat_simp,mesh,my_grid,s_index,A.pred)
    saveRDS(boot_map, file=here("generated_data", "boot-maps",
                                paste0("boot_map_",slurm_rep,".rds")))
}
time_out <- Sys.time()
time_out-time_in
