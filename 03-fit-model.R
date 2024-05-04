
library(MARS)
library(tidyverse)

dat <- readRDS("data/MARSdata_April2024.rds")
dat <- check_data(dat)

# EBFT recruits in MED
# WBFT recruits in GOM, WATL
log_recdist_rs <- matrix(0, dat@Dmodel@nr, dat@Dmodel@ns)
log_recdist_rs[-4, 1] <- -1000
log_recdist_rs[3:4, 2] <- -1000

map_recdist_rs <- matrix(NA, dat@Dmodel@nr, dat@Dmodel@ns)
map_recdist_rs[2, 2] <- 1

pars <- make_parameters(
  dat,
  start = list(log_recdist_rs = log_recdist_rs, R0_s = c(2000, 20000)),
  map = list(log_recdist_rs = factor(map_recdist_rs)),
  est_mov = "gravity_fixed"
)

saveRDS(pars, file = "data/pars_April2024.rds")
pars <- readRDS("data/pars_April2024.rds")

tictoc::tic()
fit <- fit_MARS(
  dat,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()

saveRDS(fit, file = "fit/model_April2024.rds")

fit <- readRDS("fit/model_April2024.rds")
report(fit, amov = 1, dir = "fit", filename = "preliminary_fit")

