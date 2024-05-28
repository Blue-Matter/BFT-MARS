
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

# Show fit to indices (bespoke figure not available in reporting)

# Stock specific index
index <- lapply(22:28, function(i) {
  s <- dat@Dsurvey@samp_irs[i, , ]
  plot_index(fit, i = i) %>%
    mutate(EBFT = colSums(s)[1] > 0)
}) %>%
  bind_rows() %>%
  mutate(stock = ifelse(EBFT, "EBFT", "WBFT"))

g <- ggplot(index, aes(year, obs, shape = stock)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(vars(name), scales = "free_y") +
  geom_line(data = index, colour = "red", aes(year, pred), inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Year", y = "Stock-specific index", shape = "Stock") +
  theme(legend.position = "bottom")
ggsave("fit/figures/index_stock_fit_annual.png", g, height = 6, width = 8)

# Fishery CPUE
cpue <- lapply(1:21, function(i) plot_index(fit, i = i)) %>%
  bind_rows()

g <- ggplot(cpue, aes(year, obs)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(vars(name), scales = "free_y") +
  geom_line(data = cpue, colour = "red", aes(year, pred), inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Year", y = "Fishery CPUE") +
  theme(legend.position = "bottom")
ggsave("fit/figures/index_fishery_fit_annual.png", g, height = 6, width = 10)

png("fit/figures/recruitment_fit_annual.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) plot_R(fit, s = s))
dev.off()

png("fit/figures/spawning_fit_annual.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) plot_S(fit, s = s))
dev.off()

png("fit/figures/bdist_fit_annual.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) {
  plot_B(fit, by = "region", s = s, prop = TRUE)
  title(ifelse(s == 1, "EBFT", "WBFT"), font.main = 1)
})
dev.off()

png("fit/figures/rdist_fit_annual.png", height = 5, width = 5, units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(5, 4, 1, 1))
lapply(1:4, function(r) {
  region <- c("GOM", "WATL", "EATL", "MED")
  plot_B(fit, by = "stock", r = r, prop = TRUE)
  title(region[r], font.main = 1)
})
dev.off()

