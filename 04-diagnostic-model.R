
library(MARS)

fit <- readRDS("fit/model_April2024.rds")

r <- MARS::retrospective(fit, yret = 1:7, cores = 7)
saveRDS(r, file = "fit/ret_April2024.rds")
report(r, dir = "fit", filename = "preliminary_fit_retro")

Rprof <- expand.grid(
  seq(1000, 5000, 500),
  seq(30, 90, 10)
)

p <- profile(
  fit,
  p1 = "R0_s[1]",
  p2 = "R0_s[2]",
  v1 = seq(1500, 5400, 500),
  v2 = seq(35, 95, 10),
  cores = 10
)
saveRDS(p, file = "fit/R0prof_April2024.rds")

p1 <- profile(
  fit,
  p1 = "R0_s[1]",
  v1 = seq(1000, 5000, 250)[13:15],
  cores = 3
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

plot(p, nlevels = 40, xlab = "EBFT R0", ylab = "WBFT R0", main = "Change in objective function")

plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_I_ymi')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_CAL_ymfr')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'loglike_SC_ymafr')
plot(p, levels = seq(0, 200, 10), xlab = "EBFT R0", ylab = "WBFT R0", component = 'logprior_rdev_ys')


# Profile over R0 for EBFT
p <- readRDS("fit/R01prof_April2024.rds")
plot(p, xlab = "EBFT R0", ylab = "Change in objective function")

plot(p, xlab = "EBFT R0", component = 'loglike_I_ymi')
plot(p, xlab = "EBFT R0", component = 'loglike_CAL_ymfr')
plot(p, xlab = "EBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, xlab = "EBFT R0", component = 'loglike_SC_ymafr')
plot(p, xlab = "EBFT R0", component = 'logprior_rdev_ys')

# Profile over R0 for WBFT
p <- readRDS("fit/R02prof_April2024.rds")
plot(p, xlab = "WBFT R0", ylab = "Change in objective function")

plot(p, xlab = "WBFT R0", component = 'loglike_I_ymi')
plot(p, xlab = "WBFT R0", component = 'loglike_CAL_ymfr')
plot(p, xlab = "WBFT R0", component = 'loglike_tag_mov_ymars')
plot(p, xlab = "WBFT R0", component = 'loglike_SC_ymafr')
plot(p, xlab = "WBFT R0", component = 'logprior_rdev_ys')



png("fit/R0_profile.png", height = 8, width = 4, res = 400, units = "in")

par(mfrow = c(3, 1), mar = c(5, 4, 1, 1))
# Joint profile over R0 for both EBFT/WBFT
p <- readRDS("fit/R0prof_April2024.rds")
plot(p, nlevels = 40, xlab = "EBFT R0", ylab = "WBFT R0", main = "Change in objective function")


# Profile over R0 for EBFT
p <- readRDS("fit/R01prof_April2024.rds")
plot(p, xlab = "EBFT R0", ylab = "Change in objective function")

# Profile over R0 for WBFT
p <- readRDS("fit/R02prof_April2024.rds")
plot(p, xlab = "WBFT R0", ylab = "Change in objective function")

dev.off()
