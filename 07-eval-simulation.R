


sim2 <- readRDS("simulation/sim_res.rds")
simeast <- readRDS("simulation/sim_east.rds")
simwest <- readRDS("simulation/sim_west.rds")
simdata <- readRDS("simulation/sims.rds")
# OM
fit <- readRDS("fit/model_spatprior_April2024.rds")


gr <- sapply(sim2, function(x) x$gr %>% abs() %>% max())
any(gr > 0.1)

gr <- sapply(simeast, function(x) x$gr %>% abs() %>% max())
any(gr > 0.1)

gr <- sapply(simwest, function(x) x$gr %>% abs() %>% max())
any(gr > 0.1)

get_Isim <- function(simdata, d) {
  year <- seq(min(d@Dlabel@year), max(d@Dlabel@year) + 1, 0.25)

  Iobs <- d@Dsurvey@Iobs_ymi %>%
    MARS:::collapse_yearseason() %>%
    reshape2::melt() %>%
    mutate(Index = d@Dlabel@index[Var2], Year = year[Var1]) %>%
    filter(!is.na(value))

  Isim <- lapply(1:length(simdata), function(x) {
    simdata[[x]]@Dsurvey@Iobs_ymi %>%
      MARS:::collapse_yearseason() %>%
      reshape2::melt() %>%
      mutate(Index = d@Dlabel@index[Var2], Year = year[Var1], Sim = x)
  }) %>%
    bind_rows() %>%
    filter(!is.na(value))

  g <- ggplot(Isim, aes(Year, value)) +
    facet_wrap(vars(Index)) +
    geom_path(alpha = 0.2, aes(group = Sim)) +
    geom_line(data = Iobs, colour = "red") +
    geom_point(data = Iobs, colour = "red") +
    expand_limits(y = 0)
  g
}
g <- get_Isim(simdata, MARS::get_MARSdata(fit))
ggsave("simulation/figures/sim_index.png", g, height = 10, width = 12)

get_Cresid <- function(sim, d) {
  Cdat <- d@Dfishery@Cobs_ymfr
  pos_catch <- Cdat > 0
  Cpred <- apply(sim$report$CB_ymfrs, 1:4, sum)
  Cresid <- Cdat/Cpred
  range(Cresid[pos_catch])
}
Crange <- mapply(get_Cresid, sim = sim2, d = simdata)

get_Catch <- function(x, sim) {
  Cpred <- apply(sim[[x]]$report$CB_ymfrs, 1, sum)

  Cpred
}

Catch <- sapply(1:length(sim2), get_Catch, sim = sim2)
Ctrue <- get_MARSdata(fit)@Dfishery@Cobs_ymfr %>% apply(1, sum)
matplot(Catch, typ = 'l')
lines(Ctrue, col = 2, lty = 2, lwd = 5)

get_SB <- function(x, sim) {
  good_sim <- max(sim[[x]]$report$F_ymfr) > 0.01
  sim[[x]]$report$S_yrs %>%
    apply(c(1, 3), sum) %>%
    reshape2::melt() %>%
    rename(Year = Var1, Stock = Var2, SB = value) %>%
    mutate(Sim = x, good_sim = good_sim)
}

get_SB0 <- function(x, sim) {
  good_sim <- max(sim[[x]]$report$F_ymfr) > 0.01
  sim[[x]]$report$SB0_s %>%
    reshape2::melt() %>%
    rename(SB0 = value) %>%
    mutate(Sim = x, Stock = 1:nrow(.), good_sim = good_sim)
}

get_R <- function(x, sim) {
  good_sim <- max(sim[[x]]$report$F_ymfr) > 0.01
  sim[[x]]$report$R_ys %>%
    reshape2::melt() %>%
    rename(Year = Var1, Stock = Var2, R = value) %>%
    mutate(Sim = x, good_sim = good_sim)
}

get_Rdev <- function(x, sim) {
  good_sim <- max(sim[[x]]$report$F_ymfr) > 0.01
  sim[[x]]$report$Rdev_ys %>%
    reshape2::melt() %>%
    rename(Year = Var1, Stock = Var2, Rdev = value) %>%
    mutate(Sim = x, log_rdev = log(Rdev), good_sim = good_sim)
}

val_2stock <- local({
  SB <- lapply(1:length(sim2), get_SB, sim = sim2) %>%
    bind_rows()

  SB0 <- lapply(1:length(sim2), get_SB0, sim = sim2) %>%
    bind_rows()

  R <- lapply(1:length(sim2), get_R, sim = sim2) %>%
    bind_rows()

  Rdev <- lapply(1:length(sim2), get_Rdev, sim = sim2) %>%
    bind_rows()

  left_join(SB, SB0, by = c("Sim", "Stock", "good_sim")) %>%
    left_join(R, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    left_join(Rdev, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    mutate(`SB/SB0` = SB/SB0) %>%
    select(!SB0) %>%
    #filter(SB < 5e14) %>%
    reshape2::melt(id.vars = c("Year", "Stock", "Sim", "good_sim")) %>%
    mutate(Stock = ifelse(Stock == 1, "EBFT", "WBFT"), Model = "Multi-stock")
})

val_2stock %>% filter(variable == "SB") %>% pull(good_sim) %>% mean()

g <- val_2stock %>%
  filter(good_sim) %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path() +
  facet_grid(vars(variable), vars(Stock), scales = "free_y") +
  expand_limits(y = 0)

val_east <- local({
  SB <- lapply(1:length(sim2), get_SB, sim = simeast) %>%
    bind_rows()

  SB0 <- lapply(1:length(sim2), get_SB0, sim = simeast) %>%
    bind_rows()

  R <- lapply(1:length(sim2), get_R, sim = simeast) %>%
    bind_rows()

  Rdev <- lapply(1:length(sim2), get_Rdev, sim = simeast) %>%
    bind_rows()

  left_join(SB, SB0, by = c("Sim", "Stock", "good_sim")) %>%
    left_join(R, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    left_join(Rdev, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    mutate(`SB/SB0` = SB/SB0) %>%
    select(!SB0) %>%
    #filter(SB < 5e14) %>%
    reshape2::melt(id.vars = c("Year", "Stock", "Sim", "good_sim")) %>%
    mutate(Stock = "EBFT", Model = "Single-stock")
})

val_west <- local({
  SB <- lapply(1:length(sim2), get_SB, sim = simwest) %>%
    bind_rows()

  SB0 <- lapply(1:length(sim2), get_SB0, sim = simwest) %>%
    bind_rows()

  R <- lapply(1:length(sim2), get_R, sim = simwest) %>%
    bind_rows()

  Rdev <- lapply(1:length(sim2), get_Rdev, sim = simwest) %>%
    bind_rows()

  left_join(SB, SB0, by = c("Sim", "Stock", "good_sim")) %>%
    left_join(R, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    left_join(Rdev, by = c("Year", "Sim", "Stock", "good_sim")) %>%
    mutate(`SB/SB0` = SB/SB0) %>%
    select(!SB0) %>%
    #filter(SB < 5e14) %>%
    reshape2::melt(id.vars = c("Year", "Stock", "Sim", "good_sim")) %>%
    mutate(Stock = "WBFT", Model = "Single-stock")
})

val_west %>% filter(variable == "SB") %>% pull(good_sim) %>% mean()

# WBFT
WBFT <- rbind(
  val_west,
  val_2stock
) %>%
  filter(good_sim, Stock == "WBFT")

val_om <-  local({
  SB <- get_SB(1, sim = list(list(report = fit@report)))
  SB0 <- get_SB0(1, sim = list(list(report = fit@report)))

  R <- get_R(1, sim = list(list(report = fit@report)))
  Rdev <- get_Rdev(1, sim = list(list(report = fit@report)))

  left_join(SB, SB0, by = c("Sim", "Stock")) %>%
    left_join(R, by = c("Year", "Sim", "Stock")) %>%
    left_join(Rdev, by = c("Year", "Sim", "Stock")) %>%
    mutate(`SB/SB0` = SB/SB0) %>%
    select(!SB0) %>%
    #filter(SB < 5e14) %>%
    reshape2::melt(id.vars = c("Year", "Stock", "Sim")) %>%
    mutate(Stock = ifelse(Stock == 1, "EBFT", "WBFT"))
})

g1 <- WBFT %>%
  filter(variable == "SB") %>%
  #filter(value < 1e12) %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "WBFT", variable == "SB"), linetype = 2, colour = "red") +
  facet_wrap(vars(Model), scales = "free_y", ncol = 2) +
  expand_limits(y = 0) +
  labs(y = "SB") +
  ggtitle("WBFT")

g2 <- WBFT %>%
  filter(variable == "SB/SB0") %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "WBFT", variable == "SB/SB0"), linetype = 2, colour = "red") +
  facet_wrap(vars(Model), ncol = 2) +
  expand_limits(y = 0) +
  labs(y = "SB/SB0")

g <- ggpubr::ggarrange(g1, g2, ncol = 1)
ggsave("simulation/figures/sim_WBFT.png", g, height = 6, width = 6)


# EBFT
EBFT <- rbind(
  val_east,
  val_2stock
) %>%
  filter(Stock == "EBFT")


  filter(good_sim, Stock == "EBFT")

g1 <- EBFT %>%
  filter(variable == "SB") %>%
  #filter(value < 1e14) %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "EBFT", variable == "SB"), linetype = 2, colour = "red") +
  facet_wrap(vars(Model), scales = "free_y", ncol = 2) +
  expand_limits(y = 0) +
  labs(y = "SB") +
  ggtitle("EBFT")

g2 <- EBFT %>%
  filter(variable == "SB/SB0") %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "EBFT", variable == "SB/SB0"), linetype = 2, colour = "red") +
  facet_wrap(vars(Model), ncol = 2) +
  expand_limits(y = 0) +
  labs(y = "SB/SB0")

g <- ggpubr::ggarrange(g1, g2, ncol = 1)
ggsave("simulation/figures/sim_EBFT.png", g, height = 6, width = 6)

# All state variables
g <- WBFT %>%
  #filter(variable == "SB") %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "WBFT"), linetype = 2, colour = "red") +
  facet_grid(vars(variable), vars(Model), scales = "free_y") +
  expand_limits(y = 0) +
  #labs(y = "SB") +
  ggtitle("WBFT")

g <- EBFT %>%
  #filter(variable == "SB") %>%
  ggplot(aes(Year, value, group = factor(Sim))) +
  geom_path(alpha = 0.5) +
  geom_line(data = val_om %>% filter(Stock == "WBFT"), linetype = 2, colour = "red") +
  facet_grid(vars(variable), vars(Model), scales = "free_y") +
  expand_limits(y = 0) +
  #labs(y = "SB") +
  ggtitle("EBFT")
