
library(MARS)
library(snowfall)
library(pbapply)
library(dplyr)

source("99-functions-sim.R")

# Data from OM
sims <- readRDS("simulation/sims.rds")

# Estimation model is original model
est <- readRDS("fit/model_spatprior_April2024.rds")

# MARSdata - data from OM
# parameters - parameters for estimation model
# map - map for estimation model
# Example one sim
#MARSdata <- sims[[1]]

#parameters <- as.list(est@SD, what = "Estimate") %>%
#  structure("what" = NULL)
#map <- get_MARSdata(est)@Misc$map
#random <- get_MARSdata(est)@Misc$random

#pars <- list(p = parameters, map = map, random = random)
#fit1 <- fit_sim(MARSdata, pars, nr = 4, ns = 2)

sim_res <- sim_2stock(sims, est)
saveRDS(sim_res, file = "simulation/sim_res.rds")

# Estimation model for eastern area
est_east <- readRDS("fit/modelE_April2024.rds")
sim_east <- sim_1stock(sims, est_east, r = 3:4, s = 1, log_R0 = log(6000))
saveRDS(sim_east, file = "simulation/sim_east.rds")

# Estimation model for western area
est_west <- readRDS("fit/modelW_April2024.rds")
sim_west <- sim_1stock(sims, est_west, r = 1:2, s = 2, log_R0 = log(500))
saveRDS(sim_west, file = "simulation/sim_west.rds")


