
library(MARS)
dat_full <- readRDS("data/MARSdata_April2024.rds")

# Western assessment, seasonal, one area, one stock
datw <- dat_full

s_west <- 2
datw@Dmodel@nr <- 1
datw@Dmodel@ns <- 1
datw@Dmodel@scale_s <- 1
datw@Dstock@len_ymas <- datw@Dstock@len_ymas[, , , s_west, drop = FALSE]
datw@Dstock@sdlen_ymas <- datw@Dstock@sdlen_ymas[, , , s_west, drop = FALSE]
datw@Dstock@swt_ymas <- datw@Dstock@swt_ymas[, , , s_west, drop = FALSE]
datw@Dstock@Md_yas <- datw@Dstock@Md_yas[, , s_west, drop = FALSE]
datw@Dstock@SRR_s <- "BH"
datw@Dstock@delta_s <- 0
datw@Dstock@presence_rs <- new("matrix")

r_west <- 1:2
datw@Dfishery@Cobs_ymfr <- array(0, c(58, 4, 18, 1))
datw@Dfishery@Cobs_ymfr[] <- apply(dat_full@Dfishery@Cobs_ymfr[, , , r_west], 1:3, sum)
has_catch <- sapply(1:18, function(f) sum(datw@Dfishery@Cobs_ymfr[, , f, ]) > 0)
datw@Dfishery@Cobs_ymfr <- datw@Dfishery@Cobs_ymfr[, , has_catch > 0, , drop = FALSE]
datw@Dfishery@nf <- sum(has_catch)

datw@Dfishery@CALobs_ymlfr <- array(0, c(58, 4, 15, sum(has_catch), 1))
datw@Dfishery@CALobs_ymlfr[] <- apply(dat_full@Dfishery@CALobs_ymlfr[, , , has_catch, r_west], 1:4, sum)
datw@Dfishery@fcomp_like <- "multinomial"
datw@Dfishery@CALN_ymfr <- apply(
  datw@Dfishery@CALobs_ymlfr,
  c(1, 2, 4, 5),
  function(x) ifelse(sum(x) > 0, 40, 0)
)


datw@Dfishery@sel_f <- datw@Dfishery@sel_f[has_catch]

datw@Dfishery@SC_ymafrs <- new("array")

has_index <- sapply(1:28, function(i) sum(dat_full@Dsurvey@samp_irs[i, r_west, s_west]) > 0)

datw@Dsurvey@ni <- sum(has_index)
datw@Dsurvey@Iobs_ymi <- dat_full@Dsurvey@Iobs_ymi[, , has_index]
datw@Dsurvey@Isd_ymi <- dat_full@Dsurvey@Isd_ymi[, , has_index]
datw@Dsurvey@unit_i <- dat_full@Dsurvey@unit_i[has_index]
datw@Dsurvey@samp_irs <- new("array")
datw@Dsurvey@sel_i <- dat_full@Dsurvey@sel_i[has_index]

datw@Dlabel@index[has_index]
datw@Dsurvey@sel_i
datw@Dlabel@fleet[has_catch]

datw@Dsurvey@sel_i[datw@Dsurvey@sel_i == "14"] <- "5"
datw@Dsurvey@sel_i[datw@Dsurvey@sel_i == "15"] <- "6"
datw@Dsurvey@sel_i[datw@Dsurvey@sel_i == "16"] <- "7"
datw@Dsurvey@sel_i[datw@Dsurvey@sel_i == "18"] <- "9"


datw@Dtag <- new("Dtag")

datw@Dlabel@stock <- 1
datw@Dlabel@region <- character(0)
datw@Dlabel@fleet <- dat_full@Dlabel@fleet[has_catch]
datw@Dlabel@index <- dat_full@Dlabel@index[has_index]

datw <- check_data(datw)

pars <- make_parameters(
  datw,
  #map = list(log_rdev_ys = matrix(NA, datw@Dmodel@ny, 1) %>% factor()),
  start = list(R0_s = 500)
)

# Make fishery sel priors
datw@Dmodel@prior <- sapply(1:datw@Dfishery@nf, function(f) {

  # Uninformative prior for length of full selectivity
  p1 <- paste0("dnorm(p$sel_pf[1, ", f, "], 0, 1.75, log = TRUE)")

  # Ascending limb with lognormal SD = 0.3
  start_p2 <- round(pars$p$sel_pf[2, f], 2)
  p2 <- paste0("dnorm(p$sel_pf[2, ", f, "], ", start_p2, ", 0.3, log = TRUE)")

  # Descending limb with lognormal SD = 0.3
  if (grepl("dome", datw@Dfishery@sel_f[f])) {
    start_p3 <- round(pars$p$sel_pf[3, f], 2)
    p3 <- paste0("dnorm(p$sel_pf[3, ", f, "], ", start_p3, ", 0.3, log = TRUE)")
  } else {
    p3 <- NULL
  }

  c(p1, p2, p3)
}) %>%
  unlist()

tictoc::tic()
fit <- fit_MARS(
  datw,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()

saveRDS(fit, file = "fit/modelWmult_April2024.rds")

fit <- readRDS(file = "fit/modelW_April2024.rds")
report(fit, dir = "fit", filename = "preliminary_fit_WATL")

fit@report$loglike
fit@report$loglike_CAL_ymfr %>% apply(3, sum)
fit@report$loglike_I_ymi %>% apply(3, sum)
fit@report$logprior_rdev_ys %>% sum()

obj <- fit@obj
obj$par <- par
opt <- optim(par,  method = "BFGS")
opt <- optimize_RTMB(obj)
SD <- RTMB::sdreport(fit@obj)
