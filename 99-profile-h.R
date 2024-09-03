
library(MARS)

fit <- readRDS("fit/model_April2024.rds")

p <- profile(
  fit,
  p1 = "h_s[1]",
  p2 = "h_s[2]",
  v1 = seq(0.3, 0.95, 0.05),
  v2 = seq(0.3, 0.95, 0.05),
  cores = 4
)
saveRDS(p, file = "fit/hprof_April2024.rds")
p <- readRDS(file = "fit/hprof_April2024.rds")

png("fit/figures/h_profile.png", height = 4, width = 5, res = 400, units = "in")
par(mar = c(5, 4, 1, 1))
plot(p, nlevels = 40, xlab = "EBFT steepness", ylab = "WBFT steepness", main = "Change in objective function")
dev.off()
