
library(tidyverse)
library(MARS)

# Get data list from ABFTMSE package ----
dat <- local({
  load("data/MARS_input_March24.rda")
  MARSinput
})

# Legends ----
year_df <- data.frame(Year = 1:dat$ny, Real_year = seq(dat$years[1], dat$years[2]))
area_df <- data.frame(Area = 1:dat$nr, AreaName = dat$areas) %>%
  mutate(AreaName = factor(AreaName, levels = AreaName))
fleet_df <- data.frame(Fleet = 1:dat$nf) %>%
  mutate(FleetName = paste0(Fleet, "-", dat$Fleets$name),
         FleetName = factor(FleetName, levels = FleetName))

# Model configuration ----
Dmodel <- new(
  "Dmodel",
  ny = dat$ny, # 58
  nm = 4, #seasons, dat$ns
  na = dat$na, #35
  nl = dat$nl, #15
  nr = dat$nr,
  ns = 2, #stocks, dat$np
  lbin = dat$lenbins, # Length-16
  lmid = dat$lenbins[1:15] + 0.5 * 25,
  Fmax = 3,
  y_phi = 1,
  scale_s = c(1, 10),
  nyinit = 20,
  condition = "catch"
)

# Stock configuration ----
Dstock <- new(
  "Dstock",
  m_spawn = 2,
  m_rec = 2,
  SRR_s = c("BH", "BH"),
  delta_s = c(0, 0)
)

Dstock@len_ymas <- local({
  L2 <- c(318.85, 270.6) # EBFT = Stock 1 !
  A2 <- c(999, 34)

  p <- c(1, -0.12)

  L1 <- c(0, 33)
  A1 <- c(-0.97, 0)

  K <- c(0.093, 0.22)

  #dat$len_age is not season specific?

  sapply(1:Dmodel@ns, function(s) {
    sapply(1:Dmodel@na, function(aa) {
      sapply(1:Dmodel@nm, function(m) {
        tt <- aa + (m - 1)/Dmodel@nm
        ans <- L1[s]^p[s] + (L2[s]^p[s] - L1[s]^p[s]) * (1 - exp(-K[s] * (tt - A1[s])))/(1 - exp(-K[s] * (A2[s] - A1[s])))
        rep(ans^(1/p[s]), Dmodel@ny)
      })
    }, simplify = "array")
  }, simplify = "array")
})

Dstock@sdlen_ymas <- 0.06 * Dstock@len_ymas + 5.84

Dstock@matd_yas <- local({
  #dat$mat # Younger maturing ogive - ns, na, ny
  # Use older ogive below:
  mat_1 <- rep(1, Dmodel@na)
  mat_1[1:8] <- c(0, 0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9)

  mat_2 <- rep(1, Dmodel@na)
  mat_2[1:12] <- c(0, 0, 0, 0, 0, 0, 0.01, 0.04, 0.19, 0.56, 0.88, 0.98)

  matd_yas <- array(NA_real_, c(Dmodel@ny, Dmodel@na, Dmodel@ns))
  matd_yas[, , 1] <- matrix(mat_1, Dmodel@ny, Dmodel@na, byrow = TRUE)
  matd_yas[, , 2] <- matrix(mat_2, Dmodel@ny, Dmodel@na, byrow = TRUE)
  matd_yas
})

Dstock@swt_ymas <- local({
  swt_ymas <- Dstock@len_ymas
  swt_ymas[, , , 1] <- 3.50801e-5 * Dstock@len_ymas[, , , 1] ^ 2.878451
  swt_ymas[, , , 2] <- 1.77054e-5 * Dstock@len_ymas[, , , 2] ^ 3.001252
  swt_ymas
})

Dstock@Md_yas <- dat$Ma %>% # High M scenario
  array(c(Dmodel@ns, Dmodel@na, Dmodel@ny)) %>%
  aperm(3:1)

# Fishery data ----
Dfishery <- new(
  "Dfishery",
  nf = dat$nf
)


# Catch
Dfishery@Cobs_ymfr <- array(0, c(Dmodel@ny, Dmodel@nm, Dfishery@nf, Dmodel@nr))
Dfishery@Cobs_ymfr[dat$Cobs[, c("Year", "Quarter", "Fleet", "Area")]] <- 1e-3 * dat$Cobs[, "Catch"] # Convert kg to tonnes

# Length comp
Dfishery@CALobs_ymlfr <- local({
  CAL <- array(0, c(Dmodel@ny, Dmodel@nm, Dmodel@nl, Dfishery@nf, Dmodel@nr))
  x <- dat$CLobs[dat$CLobs[, "Length_category"] <= dat$nl, ]
  CAL[x[, c("Year", "Subyear", "Length_category", "Fleet", "Area")]] <- x[, "N"]
  CAL
})
Dfishery@fcomp_like <- "lognormal"

Dfishery@sel_f <- ifelse(dat$seltype == 2, "logistic_length", "dome_length")

# Stock composition
SOO_genetic <- dat$SOOobs %>%
  as.data.frame() %>%
  filter(Type == 2, !is.na(r)) %>%
  mutate(EBFT = plogis(probE) * N, WBFT = N - EBFT) %>%
  reshape2::melt(id.vars = c("a", "y", "s", "r", "SE"), measure.vars = c("EBFT", "WBFT")) %>%
  mutate(stock = ifelse(variable == "EBFT", 1, 2), f = 1) %>%
  select(!variable) %>%
  as.matrix()

Dfishery@SC_ymafrs <- array(0, c(Dmodel@ny, Dmodel@nm, 3, 1, Dmodel@nr, Dmodel@ns))
Dfishery@SC_ymafrs[SOO_genetic[, c("y", "s", "a", "f", "r", "stock")]] <- SOO_genetic[, "value"]

Dfishery@SC_aa <- matrix(0, 3, Dmodel@na)
Dfishery@SC_aa[1, 1:4] <- Dfishery@SC_aa[2, 5:8] <- Dfishery@SC_aa[3, 9:Dmodel@na] <- 1

Dfishery@SC_ff <- matrix(1, 1, Dfishery@nf)

Dfishery@SC_like <- "lognormal"
Dfishery@SCstdev_f <- mean(dat$SOOobs[, "SE"])

# Indices/survey ----
nseries <- dat$CPUEobs %>%
  as.data.frame() %>%
  summarise(n = n(), .by = c(Fleet, Area, qNo)) %>%
  arrange(Fleet, Area, qNo)

CPUE_process <- lapply(1:nrow(nseries), function(i) {
  val <- dat$CPUEobs %>%
    as.data.frame() %>%
    filter(Fleet == nseries$Fleet[i], Area == nseries$Area[i], qNo == nseries$qNo[i]) %>%
    as.matrix()

  index <- index_sd <- matrix(NA, Dmodel@ny, Dmodel@nm)
  index[val[, c("Year", "Quarter")]] <- val[, "Index"]
  index_sd[val[, c("Year", "Quarter")]] <- val[, "CV"]

  samp <- matrix(0, Dmodel@nr, Dmodel@ns)
  samp[val[, "Area"], ] <- 1

  ff <- nseries$Fleet[i]
  yy <- year_df$Real_year[val[, "Year"]] %>% range(na.rm = TRUE) %>% paste(collapse = "-")
  name <- paste0(fleet_df$FleetName[nseries$Fleet[i]], ", ", area_df$AreaName[nseries$Area[i]], ", ", yy)

  list(index = index, index_sd = index_sd, samp = samp, name = name, ff = ff)
})

cpue <- sapply(CPUE_process, getElement, "index", simplify = "array")
cpue_sd <- sapply(CPUE_process, getElement, "index_sd", simplify = "array")
cpue_names <- sapply(CPUE_process, getElement, "name")
cpue_samp <- sapply(CPUE_process, getElement, "samp", simplify = "array") %>%
  aperm(c(3, 1, 2))
cpue_f <- sapply(CPUE_process, getElement, "ff")


nsurv <- dat$Iobs %>%
  as.data.frame() %>%
  summarise(n = n(), .by = c(Ino, type, stock, area))

I_process <- lapply(nsurv$Ino, function(i) {
  val <- dat$Iobs %>%
    as.data.frame() %>%
    filter(Ino == i) %>%
    as.matrix()

  index <- index_sd <- matrix(NA, Dmodel@ny, Dmodel@nm)
  index[val[, c("Year", "subyear")]] <- val[, "index"]
  index_sd[val[, c("Year", "subyear")]] <- val[, "CV"]

  samp <- matrix(0, Dmodel@nr, Dmodel@ns)

  if (unique(val[, "type"]) == 2) { # Larval surveys
    samp[val[, c("area", "stock")]] <- 1
  } else {
    samp[val[, "area"], ] <- 1
  }

  ff <- switch(
    unique(val[, "type"]) %>% as.character(),
    "2" = "SB",
    "4" = "15", # VB med, sel RRUSAFS
    "5" = "14"  # combined VB in sGSL, 150 cm+, sel of CANRR
  )

  yy <- year_df$Real_year[val[, "Year"]] %>% range(na.rm = TRUE) %>% paste(collapse = "-")
  name <- paste0(dat$Inames[i], ", ", area_df$AreaName[unique(val[, "area"])], ", ", yy)

  list(index = index, index_sd = index_sd, samp = samp, name = name, ff = ff)
})

index <- sapply(I_process, getElement, "index", simplify = "array")
index_sd <- sapply(I_process, getElement, "index_sd", simplify = "array")
index_names <- sapply(I_process, getElement, "name")
index_samp <- sapply(I_process, getElement, "samp", simplify = "array") %>%
  aperm(c(3, 1, 2))
index_f <- sapply(I_process, getElement, "ff")

Dsurvey <- new(
  "Dsurvey",
  ni = length(c(cpue_f, index_f)),
  Iobs_ymi = abind::abind(cpue, index, along = 3),
  Isd_ymi = abind::abind(cpue_sd, index_sd, along = 3),
  unit_i = rep("B", 28),
  samp_irs = abind::abind(cpue_samp, index_samp, along = 1),
  sel_i = c(cpue_f, index_f),
  delta_i = 0
)

# Tagging data ----
Dtag <- new(
  "Dtag",
  tag_like = "multinomial"
)

PSAT <- dat$PSAT %>%
  as.data.frame() %>%
  mutate(y = 1) %>%
  as.matrix()
Dtag@tag_ymarrs <- array(0, c(1, Dmodel@nm, 3, Dmodel@nr, Dmodel@nr, Dmodel@ns))
Dtag@tag_ymarrs[PSAT[, c("y", "s", "a", "fr", "tr", "p")]] <- PSAT[, "N"]
Dtag@tag_yy <- matrix(1:Dmodel@ny, 1)
Dtag@tag_aa <- Dfishery@SC_aa
Dtag@tagN_ymars <- apply(Dtag@tag_ymarrs, c(1, 2, 3, 4, 6), sum)


## ID possible stock movements ----
mov_seas <- Dtag@tag_ymarrs[1, , , , , ] %>%
  structure(dimnames = list(Season = paste("Season", 1:4), AgeClass = paste("Age class", 1:3), From = area_df$AreaName, To = area_df$AreaName, Stock = c("EBFT", "WBFT"))) %>%
  reshape2::melt()
mov_ann <- Dtag@tag_ymarrs[1, , , , , ] %>%
  structure(dimnames = list(Season = paste("Season", 1:4), AgeClass = paste("Age class", 1:3), From = area_df$AreaName, To = area_df$AreaName, Stock = c("EBFT", "WBFT"))) %>%
  reshape2::melt() %>%
  summarise(value = sum(value), .by = c(AgeClass, From, To, Stock))

## EBFT presence ----
# Age class 1: MED and EATL for each season
# Age class 2: MED and EATL for each season
# Age class 3: MED, WATL, EATL for each season
g <- mov_seas %>%
  filter(Stock == "EBFT") %>%
  ggplot(aes(From, To, label = value)) +
  geom_tile(data = mov_seas %>% filter(Stock == "EBFT", value > 0), alpha = 0.25, aes(fill = value)) +
  geom_text() +
  scale_fill_viridis_c() +
  facet_grid(vars(Season), vars(AgeClass)) +
  labs(fill = "N tags") +
  ggtitle("EBFT")

## WBFT presence ----
# Age class 1-3: EATL, WATL, GOM only for each season
g <- mov_seas %>%
  filter(Stock == "WBFT") %>%
  ggplot(aes(From, To, fill = value, label = value)) +
  geom_tile(data = mov_seas %>% filter(Stock == "WBFT", value > 0), alpha = 0.25, aes(fill = value)) +
  geom_text() +
  scale_fill_viridis_c() +
  facet_grid(vars(Season), vars(AgeClass)) +
  labs(fill = "N tags") +
  ggtitle("WBFT")

## Both stocks by age class
g <- mov_ann %>%
  ggplot(aes(From, To, label = value)) +
  geom_tile(data = mov_ann %>% filter(value > 0), alpha = 0.25, aes(fill = value)) +
  geom_text() +
  scale_fill_viridis_c() +
  facet_grid(vars(Stock), vars(AgeClass)) +
  labs(fill = "N tags")

Dstock@presence_rs <- matrix(FALSE, Dmodel@nr, Dmodel@ns)
Dstock@presence_rs[-1, 1] <- TRUE
Dstock@presence_rs[-4, 2] <- TRUE


# Labels
Dlabel <- new(
  "Dlabel",
  year = year_df$Real_year,
  #season = 1:4,
  #age = 1:Dmodel@na,
  region = area_df$AreaName,
  stock = c("EBFT", "WBFT"),
  fleet = fleet_df$FleetName,
  index = c(cpue_names, index_names)
)

MARSdata <- new(
  "MARSdata",
  Dmodel = Dmodel,
  Dstock = Dstock,
  Dfishery = Dfishery,
  Dsurvey = Dsurvey,
  Dtag = Dtag,
  Dlabel = Dlabel
)

saveRDS(MARSdata, file = 'data/MARSdata_April2024.rds')
