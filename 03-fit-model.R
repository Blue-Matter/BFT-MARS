
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
  start = list(log_recdist_rs = log_recdist_rs, R0_s = c(4000, 400)),
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

r <- MARS::retrospective(fit, cores = 5)
saveRDS(r, file = "fit/ret_April2024.rds")

Rprof <- expand.grid(
  seq(1500, 3000, 250),
  seq(50, 90, 10)
)

p <- profile(
  fit,
  p1 = "R0_s[1]",
  p2 = "R0_s[2]",
  v1 = seq(1000, 5000, 500),
  v2 = seq(30, 90, 10),
  cores = 10
)
saveRDS(p, file = "fit/R0prof_April2024.rds")

p1 <- profile(
  fit,
  p1 = "R0_s[1]",
  v1 = seq(1000, 5000, 250),
  cores = 10
)
saveRDS(p1, file = "fit/R01prof_April2024.rds")


prof <- readRDS("r1prof.rds")
prof[[9]] <- prof[[10]]
prof[[9]][[1]][] <- prof[[9]][[2]][] <- NA
prof2 <- lapply(prof, bind_rows) %>% bind_rows()


R0_1 <- seq(1000, 5000, 250)

p2 <- profile(
  fit,
  p1 = "R0_s[2]",
  v1 = seq(30, 90, 5),
  cores = 10
)
saveRDS(p2, file = "fit/R02prof_April2024.rds")



# Joint profile over R0 for both EBFT/WBFT
p <- readRDS("fit/R0prof_April2024.rds")

plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0")

plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_I_ymi')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_CAL_ymfr')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_SC_ymafr')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'logprior_rdev_ys')


# Profile over R0 for EBFT
p <- readRDS("fit/R01prof_April2024.rds")
plot(p, xlab = "EBFT R0")

plot(p, xlab = "EBFT R0", component = 'loglike_I_ymi')
plot(p, xlab = "EBFT R0", component = 'loglike_CAL_ymfr')
plot(p, xlab = "EBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, xlab = "EBFT R0", component = 'loglike_SC_ymafr')
plot(p, xlab = "EBFT R0", component = 'logprior_rdev_ys')

# Profile over R0 for WBFT
p <- readRDS("fit/R02prof_April2024.rds")
plot(p, xlab = "WBFT R0")

plot(p, xlab = "WBFT R0", component = 'loglike_I_ymi')
plot(p, xlab = "WBFT R0", component = 'loglike_CAL_ymfr')
plot(p, xlab = "WBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, xlab = "WBFT R0", component = 'loglike_SC_ymafr')
plot(p, xlab = "WBFT R0", component = 'logprior_rdev_ys')

