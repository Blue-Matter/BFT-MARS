
library(MARS)
dat_full <- readRDS("data/MARSdata_April2024.rds")

# Western assessment, seasonal, one area, one stock
date <- dat_full #dat_e

s_east <- 1
date@Dmodel@nr <- 1
date@Dmodel@ns <- 1
date@Dmodel@scale_s <- 1
date@Dstock@len_ymas <- date@Dstock@len_ymas[, , , s_east, drop = FALSE]
date@Dstock@sdlen_ymas <- date@Dstock@sdlen_ymas[, , , s_east, drop = FALSE]
date@Dstock@swt_ymas <- date@Dstock@swt_ymas[, , , s_east, drop = FALSE]
date@Dstock@Md_yas <- date@Dstock@Md_yas[, , s_east, drop = FALSE]
date@Dstock@SRR_s <- "BH"
date@Dstock@delta_s <- 0
date@Dstock@presence_rs <- new("matrix")

r_east <- 3:4
date@Dfishery@Cobs_ymfr <- array(0, c(58, 4, 18, 1))
date@Dfishery@Cobs_ymfr[] <- apply(dat_full@Dfishery@Cobs_ymfr[, , , r_east], 1:3, sum)
has_catch <- sapply(1:18, function(f) sum(date@Dfishery@Cobs_ymfr[, , f, ]) > 0)
date@Dfishery@Cobs_ymfr <- date@Dfishery@Cobs_ymfr[, , has_catch > 0, , drop = FALSE]
date@Dfishery@nf <- sum(has_catch)

date@Dfishery@CALobs_ymlfr <- array(0, c(58, 4, 15, sum(has_catch), 1))
date@Dfishery@CALobs_ymlfr[] <- apply(dat_full@Dfishery@CALobs_ymlfr[, , , has_catch, r_east], 1:4, sum)

#datw@Dfishery@fcomp_like <- "multinomial"
#datw@Dfishery@CALN_ymfr <- apply(
#  datw@Dfishery@CALobs_ymlfr,
#  c(1, 2, 4, 5),
#  function(x) ifelse(sum(x) > 0, 40, 0)
#)

date@Dfishery@sel_f <- date@Dfishery@sel_f[has_catch]

date@Dfishery@SC_ymafrs <- new("array")

has_index <- sapply(1:28, function(i) sum(dat_full@Dsurvey@samp_irs[i, r_east, s_east]) > 0)

date@Dsurvey@ni <- sum(has_index)
date@Dsurvey@Iobs_ymi <- dat_full@Dsurvey@Iobs_ymi[, , has_index]
date@Dsurvey@Isd_ymi <- dat_full@Dsurvey@Isd_ymi[, , has_index]
date@Dsurvey@unit_i <- dat_full@Dsurvey@unit_i[has_index]
date@Dsurvey@samp_irs <- new("array")
date@Dsurvey@sel_i <- dat_full@Dsurvey@sel_i[has_index]

date@Dlabel@index[has_index]
date@Dsurvey@sel_i
date@Dlabel@fleet[has_catch]

date@Dsurvey@sel_i[date@Dsurvey@sel_i == "12"] <- "10"
date@Dsurvey@sel_i[date@Dsurvey@sel_i == "13"] <- "11"
date@Dsurvey@sel_i[date@Dsurvey@sel_i == "15"] <- "3" ## Mirror to RRUS?
date@Dsurvey@sel_i[date@Dsurvey@sel_i == "18"] <- "9"


date@Dtag <- new("Dtag")

date@Dlabel@stock <- 1
date@Dlabel@region <- character(0)
date@Dlabel@fleet <- dat_full@Dlabel@fleet[has_catch]
date@Dlabel@index <- dat_full@Dlabel@index[has_index]

date <- check_data(date)

pars <- make_parameters(date, start = list(R0_s = 6000))

# Make fishery sel priors
date@Dmodel@prior <- lapply(1:date@Dfishery@nf, function(f) {

  # Uninformative prior for length of full selectivity
  p1 <- paste0("dnorm(p$sel_pf[1, ", f, "], 0, 1.75, log = TRUE)")

  # Ascending limb with lognormal SD = 0.3
  start_p2 <- round(pars$p$sel_pf[2, f], 2)
  p2 <- paste0("dnorm(p$sel_pf[2, ", f, "], ", start_p2, ", 0.3, log = TRUE)")

  # Descending limb with lognormal SD = 0.3
  if (grepl("dome", date@Dfishery@sel_f[f])) {
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
  date,
  pars$p,
  pars$map,
  pars$random,
  run_model = TRUE,
  do_sd = TRUE
)
tictoc::toc()

saveRDS(fit, file = "fit/modelE_April2024.rds")

fit <- readRDS(file = "fit/modelE_April2024.rds")
report(fit, dir = "fit", filename = "preliminary_fit_EATL")

