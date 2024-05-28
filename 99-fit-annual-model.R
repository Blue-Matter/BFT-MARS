
library(MARS)
dat_full <- readRDS("data/MARSdata_April2024.rds")

# Two-stock assessment, 4 areas, annual time step
dat <- dat_full

dat@Dmodel@nm <- 1

dat@Dstock@len_ymas <- dat_full@Dstock@len_ymas[, 1, , , drop = FALSE]
dat@Dstock@sdlen_ymas <- dat_full@Dstock@sdlen_ymas[, 1, , , drop = FALSE]
dat@Dstock@swt_ymas <- dat_full@Dstock@swt_ymas[, 1, , , drop = FALSE]

dat@Dfishery@Cobs_ymfr <- array(0, c(58, 1, 18, 4))
dat@Dfishery@Cobs_ymfr[] <- apply(dat_full@Dfishery@Cobs_ymfr, c(1, 3, 4), sum)

dat@Dfishery@CALobs_ymlfr <- array(0, c(58, 1, 15, 18, 4))
dat@Dfishery@CALobs_ymlfr[] <- apply(dat_full@Dfishery@CALobs_ymlfr, c(1, 3:5), sum, na.rm = TRUE)

dat@Dfishery@SC_ymafrs <- dat@Dfishery@SCstdev_ymafrs <- array(0, c(58, 1, 3, 1, 4, 2))
dat@Dfishery@SC_ymafrs[] <- apply(dat_full@Dfishery@SC_ymafrs, c(1, 3:6), sum, na.rm = TRUE)
dat@Dfishery@SCstdev_ymafrs[] <- apply(dat_full@Dfishery@SCstdev_ymafrs, c(1, 3:6), mean, na.rm = TRUE)

dat@Dsurvey@Iobs_ymi <- dat@Dsurvey@Isd_ymi <- array(0, c(58, 1, 28))

dat@Dsurvey@Iobs_ymi[] <- apply(dat_full@Dsurvey@Iobs_ymi, c(1, 3), mean, na.rm = TRUE)
dat@Dsurvey@Isd_ymi[] <- apply(dat_full@Dsurvey@Isd_ymi, c(1, 3), mean, na.rm = TRUE)

dat@Dtag@tag_ymarrs <- array(0, c(1, 1, 1, 4, 4, 2))
dat@Dtag@tag_ymarrs[] <- apply(dat_full@Dtag@tag_ymarrs, c(1, 3:6), sum, na.rm = TRUE)

dat@Dtag@tagN_ymars <- array(0, c(1, 1, 1, 4, 2))
dat@Dtag@tagN_ymars[] <- apply(dat_full@Dtag@tagN_ymars, c(1, 3:5), sum, na.rm = TRUE)

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
  start = list(log_recdist_rs = log_recdist_rs, R0_s = c(20000, 2000)),
  map = list(log_recdist_rs = factor(map_recdist_rs)),
  est_mov = "gravity_fixed"
)

dat@Misc$map <- pars$map
dat@Misc$random <- pars$random

#debug(MARS:::.MARS)
mod <- MARS:::.MARS(pars$p, dat)

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

saveRDS(fit, file = "fit/modelannual_April2024.rds")

fit <- readRDS("fit/modelannual_April2024.rds")
report(fit, amov = 1, dir = "fit", filename = "preliminary_fit_annual")
