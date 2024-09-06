
library(MARS)
library(snowfall)
library(pbapply)

# Simulation function
sim_fn <- function(fit, do_SC = TRUE, do_tag = TRUE) {

  sim1 <- fit@obj$simulate(fit@obj$env$last.par.best)
  MARSdata <- get_MARSdata(fit)
  newdata <- MARSdata

  # Catch
  pos_catch <- MARSdata@Dfishery@Cobs_ymfr > 0
  newdata@Dfishery@Cobs_ymfr[pos_catch] <- MARSdata@Dfishery@Cobs_ymfr[pos_catch] * rlnorm(sum(pos_catch), 0, 0.05)

  # CAL
  CALpred <- apply(fit@report$CN_ymlfrs, 1:5, sum)
  ny <- MARSdata@Dmodel@ny
  nm <- MARSdata@Dmodel@nm
  nf <- MARSdata@Dfishery@nf
  nr <- MARSdata@Dmodel@nr

  solve_N <- function(N = 1, pred, val) {
    est <- N * pred * (1 - pred)
    median(est) - val
  }

  for (y in 1:ny) {
    for (m in 1:nm) {
      for (f in 1:nf) {
        for (r in 1:nr) {
          obs <- MARSdata@Dfishery@CALobs_ymlfr[y, m, , f, r]
          pred <- CALpred[y, m, , f, r]
          if (sum(obs, na.rm = TRUE) && sum(pred)) {

            ppred <- pred/sum(pred)
            pobs <- obs/sum(obs)
            var <- 1/pred
            mvar <- median(var)
            N <- try(uniroot(solve_N, interval = c(0, 500), pred = ppred, val = mvar), silent = TRUE)

            if (is.character(N)) {
              NCAL <- 500
            } else {
              NCAL <- ceiling(N$root)
            }

            # Use multinomial for CAL with variance specified by lognormal var
            newdata@Dfishery@CALobs_ymlfr[y, m, , f, r] <- rmultinom(
              1, NCAL, prob = pred
            )
          }
        }
      }
    }
  }

  # Index
  pos_ind <- !is.na(MARSdata@Dsurvey@Iobs_ymi)
  newdata@Dsurvey@Iobs_ymi[pos_ind] <- sim1$Iobs_ymi[pos_ind]

  # Stock composition
  ns <- MARSdata@Dmodel@ns
  if (do_SC) {

    aa <- dim(MARSdata@Dfishery@SC_ymafrs)[3]
    ff <- dim(MARSdata@Dfishery@SC_ymafrs)[4]

    for (a in 1:aa) {
      for (f in 1:ff) {
        aind <- MARSdata@Dfishery@SC_aa[a, ]
        find <- MARSdata@Dfishery@SC_ff[f, ]

        pred_ymrs <- apply(fit@report$CN_ymafrs[, , aind, find, , , drop = FALSE], c(1, 2, 5, 6), sum)

        for (y in 1:ny) {
          for (m in 1:nm) {
            for (r in 1:nr) {
              obs <- MARSdata@Dfishery@SC_ymafrs[y, m, a, f, r, ]
              pred <- pred_ymrs[y, m, r, ]
              if (sum(obs, na.rm = TRUE) && sum(pred)) {

                ppred <- pred/sum(pred)
                pobs <- obs/sum(obs)
                var <- MARSdata@Dfishery@SCstdev_ymafrs[y, m, a, f, r, ]^2
                mvar <- median(var)
                N <- try(uniroot(solve_N, interval = c(0, 500), pred = ppred, val = mvar), silent = TRUE)

                if (is.character(N)) {
                  NSC <- 500
                } else {
                  NSC <- ceiling(N$root)
                }

                # Use multinomial for SC with variance specified by lognormal var
                newdata@Dfishery@SC_ymafrs[y, m, a, f, r, ] <- rmultinom(
                  1, NSC, prob = pred
                )
              }
            }
          }
        }
      }
    }
  }

  if (do_tag) {
    pos_tag <- !is.na(MARSdata@Dtag@tag_ymarrs)
    newdata@Dtag@tag_ymarrs[] <- sim1$tag_ymarrs[]
  }

  #newdata <- check_data(newdata)

  return(newdata)
}


fit_sim <- function(MARSdata, pars, nr = 1, ns = 1, r, s, mult = FALSE) {

  newdata <- MARSdata

  if (nr == 1 && ns == 1) {

    require(dplyr)

    newdata@Dmodel@nr <- nr
    newdata@Dmodel@ns <- ns
    newdata@Dmodel@scale_s <- rep(1, ns)
    newdata@Dstock@len_ymas <- MARSdata@Dstock@len_ymas[, , , s, drop = FALSE]
    newdata@Dstock@sdlen_ymas <- MARSdata@Dstock@sdlen_ymas[, , , s, drop = FALSE]
    newdata@Dstock@LAK_ymals <- MARSdata@Dstock@LAK_ymals[, , , , s, drop = FALSE]
    newdata@Dstock@fec_yas <- MARSdata@Dstock@fec_yas[, , s, drop = FALSE]
    newdata@Dstock@swt_ymas <- MARSdata@Dstock@swt_ymas[, , , s, drop = FALSE]
    newdata@Dstock@Md_yas <- MARSdata@Dstock@Md_yas[, , s, drop = FALSE]
    newdata@Dstock@SRR_s <- rep("BH", ns)
    newdata@Dstock@delta_s <- rep(0, ns)
    newdata@Dstock@presence_rs <- new("matrix")
    newdata@Dstock@natal_rs <- matrix(1, nr, ns)

    has_catch <- sapply(1:MARSdata@Dfishery@nf, function(f) sum(MARSdata@Dfishery@Cobs_ymfr[, , f, r]) > 0)
    newdata@Dfishery@nf <- sum(has_catch)
    newdata@Dfishery@fwt_ymafs <- new("array")

    newdata@Dfishery@Cobs_ymfr <- array(0, c(MARSdata@Dmodel@ny, MARSdata@Dmodel@nm, newdata@Dfishery@nf, nr))
    newdata@Dfishery@Cobs_ymfr[] <- apply(MARSdata@Dfishery@Cobs_ymfr[, , has_catch, r], 1:3, sum)
    #newdata@Dfishery@Cobs_ymfr <- MARSdata@Dfishery@Cobs_ymfr[, , has_catch > 0, , drop = FALSE]

    newdata@Dfishery@CALobs_ymlfr <- array(0, c(MARSdata@Dmodel@ny, MARSdata@Dmodel@nm, MARSdata@Dmodel@nl, newdata@Dfishery@nf, nr))
    newdata@Dfishery@CALobs_ymlfr[] <- apply(MARSdata@Dfishery@CALobs_ymlfr[, , , has_catch, r], 1:4, sum)

    if (mult) {
      newdata@Dfishery@fcomp_like <- "multinomial"
      newdata@Dfishery@CALN_ymfr <- apply(
        newdata@Dfishery@CALobs_ymlfr,
        c(1, 2, 4, 5),
        function(x) ifelse(sum(x) > 0, 40, 0)
      )
    } else {
      #newdata@Dfishery@fcomp_like <- "multinomial"
      newdata@Dfishery@CALN_ymfr <- new("array")
    }

    newdata@Dfishery@sel_f <- MARSdata@Dfishery@sel_f[has_catch]
    newdata@Dfishery@SC_ymafrs <- new("array")
    newdata@Dfishery@CALtheta_f <- MARSdata@Dfishery@CALtheta_f[has_catch]
    newdata@Dfishery@sel_block_yf <- new("array")
    newdata@Dfishery@Cinit_mfr <- new("array")

    has_index <- sapply(1:MARSdata@Dsurvey@ni, function(i) sum(MARSdata@Dsurvey@samp_irs[i, r, s]) > 0)

    newdata@Dsurvey@ni <- sum(has_index)
    newdata@Dsurvey@Iobs_ymi <- MARSdata@Dsurvey@Iobs_ymi[, , has_index]
    newdata@Dsurvey@Isd_ymi <- MARSdata@Dsurvey@Isd_ymi[, , has_index]
    newdata@Dsurvey@unit_i <- MARSdata@Dsurvey@unit_i[has_index]
    newdata@Dsurvey@samp_irs <- new("array")
    newdata@Dsurvey@sel_i <- MARSdata@Dsurvey@sel_i[has_index]
    newdata@Dsurvey@delta_i <- MARSdata@Dsurvey@delta_i[has_index]

    #date@Dlabel@index[has_index]
    #date@Dsurvey@sel_i
    #date@Dlabel@fleet[has_catch]

    if (all(r == 3:4)) { # Eastern stock
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "12"] <- "10"
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "13"] <- "11"
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "15"] <- "3" ## Mirror to RRUS?
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "18"] <- "9"
    } else { # Western stock
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "14"] <- "5"
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "15"] <- "6"
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "16"] <- "7"
      newdata@Dsurvey@sel_i[newdata@Dsurvey@sel_i == "18"] <- "9"
    }

    newdata@Dtag <- new("Dtag")

    newdata@Dlabel@stock <- ns
    newdata@Dlabel@region <- character(0)
    newdata@Dlabel@fleet <- newdata@Dlabel@fleet[has_catch]
    newdata@Dlabel@index <- newdata@Dlabel@index[has_index]

    # Make fishery sel priors, no spatial targeting priors of course
    newdata@Dmodel@prior <- lapply(1:newdata@Dfishery@nf, function(f) {

      # Uninformative prior for length of full selectivity
      p1 <- paste0("dnorm(p$sel_pf[1, ", f, "], 0, 1.75, log = TRUE)")

      # Ascending limb with lognormal SD = 0.3
      start_p2 <- round(pars$p$sel_pf[2, f], 2)
      p2 <- paste0("dnorm(p$sel_pf[2, ", f, "], ", start_p2, ", 0.3, log = TRUE)")

      # Descending limb with lognormal SD = 0.3
      if (grepl("dome", newdata@Dfishery@sel_f[f])) {
        start_p3 <- round(pars$p$sel_pf[3, f], 2)
        p3 <- paste0("dnorm(p$sel_pf[3, ", f, "], ", start_p3, ", 0.3, log = TRUE)")
      } else {
        p3 <- NULL
      }

      c(p1, p2, p3)
    }) %>%
      unlist()

  }

  newdata <- check_data(newdata)

  fit <- fit_MARS(
    newdata,
    parameters = pars$p,
    map = pars$map,
    random = pars$random,
    run_model = TRUE,
    do_sd = FALSE
  )

  sim_output <- list(
    report = fit@report,
    opt = fit@opt,
    gr = fit@obj$gr(fit@opt$par)
  )

  return(sim_output)
}

sim_2stock <- function(sims, est, cpus = 10) {

  snowfall::sfInit(parallel = TRUE, cpus = cpus)
  on.exit(snowfall::sfStop())

  sfLibrary(MARS)

  parameters <- as.list(est@SD, what = "Estimate") %>%
    structure("what" = NULL)
  map <- get_MARSdata(est)@Misc$map
  random <- get_MARSdata(est)@Misc$random

  pars <- list(p = parameters, map = map, random = random)
  snowfall::sfExport(list = c("pars"))
  output <- pbapply::pblapply(sims, fit_sim, pars = pars, nr = 4, ns = 2, cl = snowfall::sfGetCluster())

  return(output)
}

sim_1stock <- function(sims, est, log_R0, r, s, cpus = 10, mult = FALSE) {

  snowfall::sfInit(parallel = TRUE, cpus = cpus)
  on.exit(snowfall::sfStop())

  sfLibrary(MARS)

  parameters <- as.list(est@SD, what = "Estimate") %>%
    structure("what" = NULL)
  if (!missing(log_R0)) parameters$t_R0_s <- log_R0

  map <- get_MARSdata(est)@Misc$map
  random <- get_MARSdata(est)@Misc$random

  pars <- list(p = parameters, map = map, random = random)
  snowfall::sfExport(list = c("pars", "r", "s"))
  output <- pbapply::pblapply(
    sims, fit_sim, pars = pars, nr = 1, ns = 1, mult = mult,
    r = r, s = s, cl = snowfall::sfGetCluster()
  )

  return(output)
}

