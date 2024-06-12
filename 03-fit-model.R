
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
  start = list(log_recdist_rs = log_recdist_rs, R0_s = c(20000, 2000)),
  map = list(log_recdist_rs = factor(map_recdist_rs)),
  est_mov = "gravity_fixed"
)

saveRDS(pars, file = "data/pars_April2024.rds")
pars <- readRDS("data/pars_April2024.rds")

tictoc::tic()
fit_noprior <- fit_MARS(
  dat,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()

saveRDS(fit_noprior, file = "fit/model_noprior_April2024.rds")


# Make fishery sel priors
prior_sel <- lapply(1:dat@Dfishery@nf, function(f) {

  # Uninformative prior for length of full selectivity
  p1 <- paste0("dnorm(p$sel_pf[1, ", f, "], 0, 1.5, log = TRUE)")

  # Ascending limb with lognormal SD = 0.3
  start_p2 <- round(pars$p$sel_pf[2, f], 2)
  p2 <- paste0("dnorm(p$sel_pf[2, ", f, "], ", start_p2, ", 0.4, log = TRUE)")

  # Descending limb with lognormal SD = 0.3
  if (grepl("dome", dat@Dfishery@sel_f[f])) {
    start_p3 <- round(pars$p$sel_pf[3, f], 2)
    p3 <- paste0("dnorm(p$sel_pf[3, ", f, "], ", start_p3, ", 0.4, log = TRUE)")
  } else {
    p3 <- NULL
  }

  c(p1, p2, p3)
}) %>%
  unlist()

prior_recdist <- c("dnorm(p$log_recdist_rs[1, 2], 0, 2, log = TRUE)",
                   "dnorm(p$log_recdist_rs[2, 2], 0, 2, log = TRUE)",
                   "dnorm(p$mov_g_ymars[1, 3, 1, 3, 2], 0, 2, log = TRUE)")
dat@Dmodel@prior <- c(prior_sel, prior_recdist)

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

# Add natal mixing priors
prior_dist <- local({
  d <- readRDS("data/SPrior.rda")
  sapply(1:nrow(d), function(i) {
    val <- d[i, ]

    s <- ifelse(val["Strata"] == 1, 2, 1)

    mov <- paste0("mov_ymarrs[1, ", val["Quarter"], ", 1, , , ", s, "]")
    start <- paste0("recdist_rs[, ", s, "]")

    r <- ifelse(val["Strata"] == 1, 1, 4)

    m <- log(val["Index"]) %>% round(3)
    cv <- val["CV"]

    paste0("calc_eqdist(", mov, ", start = ", start, ")[", r, "] %>% log() %>% dnorm(", m, ", ", cv, ", log = TRUE)")
  })

})

dat@Dmodel@prior <- c(prior_sel, prior_recdist, prior_dist)

tictoc::tic()
fit_spatprior <- fit_MARS(
  dat,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()

saveRDS(fit_spatprior, file = "fit/model_spatprior_April2024.rds")
fit_spatprior <- readRDS(file = "fit/model_spatprior_April2024.rds")
report(fit_spatprior, amov = 1, dir = "fit", filename = "preliminary_fit_spatprior")




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
ggsave("fit/figures/index_stock_fit_all.png", g, height = 6, width = 8)

g <- ggplot(index, aes(year, obs, shape = stock)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(vars(name), scales = "free_y") +
  geom_line(data = index %>% filter(!is.na(obs)), colour = "red", aes(year, pred), inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Year", y = "Stock-specific index", shape = "Stock") +
  theme(legend.position = "bottom")
ggsave("fit/figures/index_stock_fit.png", g, height = 6, width = 8)

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
ggsave("fit/figures/index_fishery_fit_all.png", g, height = 6, width = 10)

g <- ggplot(cpue, aes(year, obs)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  facet_wrap(vars(name), scales = "free_y") +
  geom_line(data = cpue %>% filter(!is.na(obs)), colour = "red", aes(year, pred), inherit.aes = FALSE) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1)) +
  labs(x = "Year", y = "Fishery CPUE") +
  theme(legend.position = "bottom")
ggsave("fit/figures/index_fishery_fit.png", g, height = 6, width = 10)

png("fit/figures/recruitment_fit.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) plot_R(fit, s = s))
dev.off()

png("fit/figures/spawning_fit.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) plot_S(fit, s = s))
dev.off()

png("fit/figures/bdist_fit.png", height = 7, width = 5, units = "in", res = 400)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
lapply(1:2, function(s) {
  plot_B(fit, by = "region", s = s, prop = TRUE)
  title(ifelse(s == 1, "EBFT", "WBFT"), font.main = 1)
})
dev.off()

png("fit/figures/rdist_fit.png", height = 5, width = 5, units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(5, 4, 1, 1))
lapply(1:4, function(r) {
  region <- c("GOM", "WATL", "EATL", "MED")
  plot_B(fit, by = "stock", r = r, prop = TRUE)
  title(region[r], font.main = 1)
})
dev.off()
