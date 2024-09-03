
library(MARS)
library(snowfall)
library(pbapply)

source("99-functions-sim.R")

# Sample from two stock model
fit <- readRDS("fit/model_spatprior_April2024.rds")
fit@obj$retape()
#sim1 <- fit@obj$simulate(fit@obj$env$last.par.best)
#saveRDS(sim1, 'sim1.rds')
#sim1 <- readRDS("sim1.rds")


snowfall::sfInit(TRUE, cpus = 10)
snowfall::sfLibrary(MARS)
snowfall::sfExport(list = c("fit", "sim_fn"))
sims <- pbapply::pblapply(1:50, function(...) sim_fn(fit), cl = snowfall::sfGetCluster())

snowfall::sfStop()

saveRDS(sims, file = "simulation/sims.rds")



