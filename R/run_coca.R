#' COCA Simulation
#'
#' @param case Trial type for stage 2. \code{case = 1} for 4-arm trial comparing AB vs. A vs. B vs. SOC; \code{case = 2} for 3-arm trial comparing AB vs. A (or B) vs. SOC; \code{case = 3} for 2-arm trial comparing AB vs. SOC.
#' @param n.stage1 Sample size for stage 1
#' @param n.stage2 Sample size for stage 2
#' @param Ce,c0 Design cutoffs, obtained using the \code{COCA.calibration} function.
#' @param tox.SOC,tox.A,tox.B True toxicity probabilities for SOC, arm A, and arm B.
#' @param tox.AB A vector of true toxicity probabilities for all combination doses being tested in the trial.
#' @param eff.SOC,eff.A,eff.B True efficacy probabilities for SOC, arm A, and arm B.
#' @param eff.AB A vector of true efficacy probabilities for all combination doses being tested in the trial.
#' @param tox.upper Highest acceptable toxicity rate (\eqn{\phi_{T}})
#' @param eff.lower Lowest acceptable efficacy rate (\eqn{\phi_{E}})
#' @param Cs Probability cutoff in stage 1
#' @param C.f1,C.f2 Probability cutoffs in stage 2
#' @param utility.score Vector of utility score: \code{c(b1, b2, b3, b4)} represents the utility for (toxicity, no efficacy), (no toxicity, no efficacy), (toxicity, efficacy), and (no toxicity, efficacy), respectively.
#' @param rho Correlation between toxicity and efficacy
#' @param period.effect Period effect
#' @param n.simu Number of simulation replicates. The default value \code{n.simu = 10} is used for illustration purposes and is small to reduce computation time. For more accurate results, consider using a larger value, such as 5000.
#' @param seed Random seed
#'
#' @return Returns the operating characteristics of stage 1 (selection and expected sample size) and stage 2 (power, GP, SR, OSR, and expected sample size).
#' @import cli
#' @import tibble
#' @import coda
#' @import isotone
#' @import rjags
#' @import runjags
#' @export
#'
#' @examples
#' # \code{n.simu = 20} is used for illustration purposes. For more accurate results,
#' # consider using a larger value for \code{n.simu}, such as 5000.
#' # Scenario 1 (period effect = 0)
#' \donttest{
#' COCA.getOC(
#'   n.stage2 = 26, Ce = 0.8983, c0 = 0.7,
#'   tox.SOC = 0.10, eff.SOC = 0.25,
#'   tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.25, eff.B = 0.25,
#'   tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.25, 0.25, 0.25),
#'   period.effect = 0, n.simu = 20
#' )}
#' # Scenario 1 (period effect = 0.2)
#' \donttest{
#' COCA.getOC(
#'   n.stage2 = 26, Ce = 0.8983, c0 = 0.7,
#'   tox.SOC = 0.10, eff.SOC = 0.25,
#'   tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.25, eff.B = 0.25,
#'   tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.25, 0.25, 0.25),
#'   period.effect = 0.2, n.simu = 20
#' )}
#' # Scenario 2 (period effect = 0)
#' \donttest{
#' COCA.getOC(
#'   n.stage2 = 26, Ce = 0.8983, c0 = 0.7,
#'   tox.SOC = 0.10, eff.SOC = 0.25,
#'   tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.35, eff.B = 0.35,
#'   tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.55, 0.55, 0.55),
#'   period.effect = 0, n.simu = 20
#' )}
#' # Scenario 2 (period effect = 0.2)
#' \donttest{
#' COCA.getOC(
#'   n.stage2 = 26, Ce = 0.8983, c0 = 0.7,
#'   tox.SOC = 0.10, eff.SOC = 0.25,
#'   tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.35, eff.B = 0.35,
#'   tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.55, 0.55, 0.55),
#'   period.effect = 0.2, n.simu = 20
#' )}
COCA.getOC <- function(case = 1, n.stage1 = 24, n.stage2, Ce, c0,
                       tox.SOC = 0.10, eff.SOC = 0.25,
                       tox.A = 0.25, tox.B = 0.15,
                       eff.A = 0.25, eff.B = 0.25,
                       tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.25, 0.25, 0.25),
                       tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, C.f1 = 0.9, C.f2 = 0.9,
                       utility.score = c(0, 60, 40, 100), rho = 0.2,
                       period.effect = 0, n.simu = 10, seed = 123) {
  runjags.options(
    inits.warning = FALSE, rng.warning = FALSE, silent.jags = TRUE, silent.runjags = TRUE
  )
  Ce <- c(Ce, c0 * Ce, c0 * Ce)

  a <- rep(0.05, 4) # hyperparameters in dirichlet
  n.dose.AB <- 3
  Nmax_p2 <- c(n.stage1, n.stage2) # stage 1 & stage 2
  T_21 <- 2
  A1 <- matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2) # order matrix for pi_tox_hat
  A2 <- matrix(c(1, 2, 1, 3), byrow = TRUE, nrow = 2) # order matrix for prob
  T_22 <- 2
  C_s1 <- C_s2 <- Cs
  C_t <- 0.10
  C.f1.trans <- 1 - C.f1
  C.f2.trans <- 1 - C.f2

  narm_22 <- switch(case,
    4,
    2,
    3
  )
  jags_params_22 <- c("pi", "wt")
  logistic_model <- .logistic_model
  period <- switch(case,
    c(0, 0, 0, 0, 1, 1, 1),
    c(0, 0, 0, 0, 1, 1, 1)[c(1, 4, 5, 6, 7)],
    c(0, 0, 0, 0, 1, 1, 1)[-3]
  )

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = TRUE)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, .dose_standardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = n.dose.AB, byrow = FALSE)
  row.names(dose_std) <- c("A", "B")

  ### simulation settings
  Tox_prob <- tox.AB
  Eff_prob <- eff.AB + period.effect
  Tox_prob_A <- rep(tox.A, n.dose.AB)
  Tox_prob_B <- rep(tox.B, n.dose.AB)
  Eff_prob_A <- rep(eff.A, n.dose.AB)
  Eff_prob_B <- rep(eff.B, n.dose.AB)

  ### Stage I
  multi_prob <- sapply(1:n.dose.AB, function(r) .solve_level(rho, Eff_prob[r], Tox_prob[r]))
  true_utility <- utility.score %*% multi_prob

  narm_21 <- n.dose.AB # #of arms in stage I
  N_21 <- sapply(1:n.simu, function(r) rep(Nmax_p2[1], 3)) # 3xsn matrix
  n_21 <- N_21 / T_21 # sample in each stages
  tox_21 <- switch(narm_21,
    c(Tox_prob, 0, 0),
    c(Tox_prob, 0),
    Tox_prob
  )
  eff_21 <- switch(narm_21,
    c(Eff_prob, 0, 0),
    c(Eff_prob, 0),
    Eff_prob
  )

  Y_21 <- pi_hat_21 <- pi_hat_21_iso <- array(0, dim = c(narm_21, nrow(multi_prob), n.simu))
  piT_hat_21 <- piE_hat_21 <- currn_21 <- matrix(0, nrow = 3, ncol = n.simu) # 3xsn matrix
  proc_21 <- matrix(rep(((tox_21 + eff_21) != 0) * 1, n.simu), nrow = 3) # trial process of stage I
  BCI <- vector(mode = "numeric", n.simu)
  cli_alert("Stage 1: in process")
  cli_progress_bar("Stage 1", total = T_21, clear = FALSE)
  for (t in 1:T_21) {
    cli_progress_update()
    set.seed(seed + 10 * t)
    currn_21 <- currn_21 + n_21 * proc_21 # current sample size
    temp_Y <- sapply(1:narm_21, function(r) rmultinom(n.simu, n_21[1, r], multi_prob[, r]))
    temp_Y <- array(t(temp_Y), dim = dim(Y_21))
    for (i in 1:nrow(multi_prob)) {
      temp_Y[, i, ] <- temp_Y[, i, ] * proc_21
    }
    Y_21 <- Y_21 + temp_Y
    for (i in 1:nrow(multi_prob)) {
      pi_hat_21[, i, ] <- (a[i] + Y_21[, i, ]) / (sum(a) + currn_21)
    }
    piT_hat_21 <- pi_hat_21[, 1, ] + pi_hat_21[, 3, ] # tox rate hat
    w1 <- 1 / apply(piT_hat_21, 1, var)
    piT_hat_21_iso <- sapply(1:n.simu, function(r) activeSet(A1, "LS", weights = w1, y = piT_hat_21[, r])$x)
    piE_hat_21 <- pi_hat_21[, 3, ] + pi_hat_21[, 4, ] # response rate hat
    rho_hat <- (pi_hat_21[, 2, ] * pi_hat_21[, 3, ] - pi_hat_21[, 1, ] * pi_hat_21[, 4, ]) / sqrt(piE_hat_21 * (1 - piE_hat_21) * piT_hat_21_iso * (1 - piT_hat_21_iso))
    for (i in 1:n.dose.AB) {
      pi_hat_21_iso[i, , ] <- sapply(1:n.simu, function(r) .solve_level(rho_hat[i, r], piE_hat_21[i, r], piT_hat_21_iso[i, r]))
    }
    prt_21 <- pbeta(tox.upper, (a[1] + a[3] + Y_21[, 1, ] + Y_21[, 3, ]), (a[2] + a[4] + currn_21 - Y_21[, 1, ] - Y_21[, 3, ]))
    w2 <- 1 / apply(prt_21, 1, var)
    prt_21_iso <- sapply(1:n.simu, function(r) activeSet(A2, "LS", weights = w2, y = prt_21[, r])$x)
    pre_21 <- pbeta(eff.lower, (a[3] + a[4] + Y_21[, 3, ] + Y_21[, 4, ]), (a[1] + a[2] + currn_21 - Y_21[, 3, ] - Y_21[, 4, ]))
    if (t < T_21) { # stage I - interim analyses
      # overly toxic stopping
      proc_21[which(proc_21 != 0 & prt_21_iso < 1 - C_s1)] <- 0
      # ineffective stopping
      proc_21[which(proc_21 != 0 & pre_21 > C_s1)] <- 0
      # if all doses have been terminated
      # terminate the whole trial and no dose is selected
      BCI[which(colSums(proc_21) == 0)] <- -1
    } else if (t == T_21) { # stage I - final analysis
      # admissible set
      # exclude early terminated arms in S by timing 'proc_21'
      Aset <- which(prt_21_iso * proc_21 > 1 - C_s2 & pre_21 * proc_21 < C_s2)
      # deciding the optimal dose
      for (i in 1:nrow(multi_prob)) {
        pi_hat_21_iso[, i, ][-Aset] <- -100
      } # doses not in set A
      utility <- t(sapply(1:n.dose.AB, function(r) utility.score %*% pi_hat_21_iso[r, , ]))

      j_ast1 <- apply(utility, 2, which.max) # 1-dimension indicator
      j_ast1[which(colSums(proc_21) == 0)] <- -1 # early terminated
      j_ast1[which(apply(utility, 2, max) < 0)] <- -1 # no OBD
    }
  }
  cli_alert_info("Stage 1: done")

  ### Stage II (Proof of Concept)
  Eff_prob <- Eff_prob - period.effect

  BCI <- matrix(rep(BCI, (narm_22 - 1)), nrow = (narm_22 - 1), byrow = TRUE)
  Yt_21 <- (Y_21[, 1, ] + Y_21[, 3, ])[, which(j_ast1 > 0)]
  Ye_21 <- (Y_21[, 3, ] + Y_21[, 4, ])[, which(j_ast1 > 0)]
  sn_22 <- sum(j_ast1 > 0)
  N_22 <- sapply(1:sn_22, function(r) rep(Nmax_p2[2], narm_22)) # narm_22xsn matrix
  n_22 <- N_22 / T_22
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]
  tox_22_all <- rbind(rep(tox.SOC, sn_22), Tox_prob_A[j_ast1_tmp], Tox_prob_B[j_ast1_tmp], Tox_prob[j_ast1_tmp])
  tox_22 <- switch(case,
    tox_22_all,
    tox_22_all[c(1, 4), ],
    tox_22_all[-3, ]
  ) # scenario 1-3 C/A/B/AB
  eff_22_all <- rbind(rep(eff.SOC, sn_22), Eff_prob_A[j_ast1_tmp], Eff_prob_B[j_ast1_tmp], Eff_prob[j_ast1_tmp])
  eff_22 <- switch(case,
    eff_22_all,
    eff_22_all[c(1, 4), ],
    eff_22_all[-3, ]
  )

  X1_all <- sapply(1:sn_22, function(r) {
    c(
      dose_level_std["A", 3], dose_std["A", 1], dose_level_std["A", 3],
      dose_std["A", j_ast1_tmp[r]], dose_std["A", ]
    )
  })
  X2_all <- sapply(1:sn_22, function(r) {
    c(
      dose_level_std["B", 3], dose_level_std["B", 3], dose_std["B", 1],
      dose_std["B", j_ast1_tmp[r]], dose_std["B", ]
    )
  })
  X1 <- switch(case,
    X1_all,
    X1_all[c(1, 4, 5, 6, 7), ],
    X1_all[-3, ]
  )
  X2 <- switch(case,
    X2_all,
    X2_all[c(1, 4, 5, 6, 7), ],
    X2_all[-3, ]
  )

  Yt_pre <- sapply(1:sn_22, function(r) c(rep(0, (narm_22 - 1)), (Yt_21[j_ast1_tmp[r], r]))) # phase II(stage I)
  Nt_pre <- sapply(1:sn_22, function(r) c(rep(0, (narm_22 - 1)), (currn_21[, which(j_ast1 > 0)][j_ast1_tmp[r], r])))

  wt <- rep(0, sn_22)
  Yt_22 <- Ye_22 <- matrix(0, nrow = narm_22, ncol = sn_22)
  proc_22 <- matrix(1, nrow = narm_22, ncol = sn_22) # trial process of stage II
  currn_22 <- matrix(0, nrow = narm_22, ncol = sn_22)
  fprob_1 <- fprob_2 <- matrix(NA, nrow = narm_22, ncol = sn_22)
  cli_alert("Stage 2: in process")
  for (t in 1:T_22) {
    set.seed(seed + 10 * t)
    currn_22 <- currn_22 + n_22 * proc_22 # current sample size
    Yt_22 <- Yt_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = tox_22[, r])) * proc_22
    Ye_22 <- Ye_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = eff_22[, r])) * proc_22
    if (t < T_22) {
      cli_progress_bar("Interim", total = sn_22, clear = FALSE)
      for (i in 1:sn_22) {
        cli_progress_update()
        if (proc_22[narm_22, i] == 0) {
          fprob_1[, i] <- fprob_2[, i] <- -1
        } else { # skip if AB has been terminated
          dataYe_22 <- list(
            x1 = X1[, i], x2 = X2[, i], period = period,
            Y = c(Ye_22[, i], Ye_21[, i]),
            N_arms = (narm_22 + n.dose.AB),
            n = c(currn_22[, i], currn_21[, which(j_ast1 > 0)][, i])
          ) # phase II data
          jagsmodel.Ye_22 <- run.jags(
            model = logistic_model, monitor = jags_params_22, data = dataYe_22,
            n.chains = 4, adapt = 2000, burnin = 3000,
            sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
            plots = FALSE, silent.jags = TRUE
          )
          codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
          piE_mcmc_22 <- matrix(NA, nrow = (jagsmodel.Ye_22$sample * length(jagsmodel.Ye_22$mcmc)), ncol = narm_22)
          for (j in 1:narm_22) {
            piE_mcmc_22[, j] <- as.matrix(codasamples.Ye_22[, j])
          }
          # futility stopping prob
          fprob_1[, i] <- c(100, sapply(2:narm_22, function(r) mean(piE_mcmc_22[, r] > piE_mcmc_22[, 1]))) # omit the 1st arm
          if (case %in% c(1, 3)) { # terminate single arm A or B
            fprob_2[, i] <- c(100, sapply(2:(narm_22 - 1), function(r) mean(piE_mcmc_22[, r] > piE_mcmc_22[, narm_22])), 100)
          } else {
            fprob_2[, i] <- c(100, 100)
          }
        }
      }

      # overly toxic stopping
      proc_22[which(proc_22 != 0 & pbeta(tox.upper, (0.1 + Yt_pre + Yt_22), (0.1 + Nt_pre + currn_22 - Yt_pre - Yt_22)) < C_t)] <- 0
      # futility stopping
      proc_22[which(proc_22 != 0 & fprob_1 < C.f1.trans)] <- 0
      proc_22[which(proc_22 != 0 & fprob_2 < C.f2.trans)] <- 0
      # terminate the trial and conclude ineffective
      # to represent the results for early termination of arm AB
      BCI[, which(proc_22[narm_22, ] == 0)] <- -2
    } else if (t == T_22) {
      cli_progress_bar("Final", total = sn_22, clear = FALSE)
      for (i in 1:sn_22) {
        cli_progress_update()
        if (proc_22[narm_22, i] == 0) {
          BCI[, i] <- -2
        } else {
          dataYe_22 <- list(
            x1 = X1[, i], x2 = X2[, i], period = period,
            Y = c(Ye_22[, i], Ye_21[, i]),
            N_arms = (narm_22 + n.dose.AB),
            n = c(currn_22[, i], currn_21[, which(j_ast1 > 0)][, i])
          ) # phase II data

          jagsmodel.Ye_22 <- run.jags(
            model = logistic_model, monitor = jags_params_22, data = dataYe_22,
            n.chains = 4, adapt = 2000, burnin = 3000,
            sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
            plots = FALSE, silent.jags = TRUE
          )
          codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
          sumYe_22 <- summary(codasamples.Ye_22)
          wt[i] <- switch(case,
            sumYe_22$quantiles[8, 3],
            sumYe_22$quantiles[6, 3],
            sumYe_22$quantiles[7, 3]
          )
          piE_mcmc_22 <- matrix(NA, nrow = (jagsmodel.Ye_22$sample * length(jagsmodel.Ye_22$mcmc)), ncol = narm_22)
          for (j in 1:narm_22) {
            piE_mcmc_22[, j] <- as.matrix(codasamples.Ye_22[, j])
          }
          # calculate BCI in the final analysis
          BCI[, i] <- sapply(1:(narm_22 - 1), function(r) mean(piE_mcmc_22[, narm_22] > piE_mcmc_22[, r]))
        }
      }
    }
  }
  cli_alert_info("Stage 2: done")

  j_opt <- which(
    true_utility == max(
      true_utility[which(tox.AB <= tox.upper &
        eff.AB >= eff.lower)]
    )
  )
  ### stage 1 results
  sel <- c()
  sel_all <- mean(sapply(1:n.simu, function(r) j_ast1[r] %in% j_opt)) * 100
  if (min(j_ast1) < 0) {
    sel <- as.vector(table(j_ast1) / n.simu) * 100
  } else {
    sel <- c(0, as.vector(table(j_ast1) / n.simu)) * 100
  }
  EN <- round(mean(colSums(currn_21)), 1)
  stage1_output <- c(sel, sel_all, EN)
  names(stage1_output) <- c("termination", "dose1", "dose2", "dose3", "selection (%)", "EN")

  ### stage 2 results
  BCI <- BCI[, 1:sum(j_ast1 > 0)]
  if (is.null(dim(BCI)[2])) {
    BCI <- matrix(BCI)
  }
  cont.simu <- dim(BCI)[2]
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]

  sel.opt.g <- sapply(1:cont.simu, function(r) j_ast1_tmp[r] %in% j_opt) * 1
  sel.opt.g <- which(sel.opt.g != 0) # correct groups
  for (jj in which(colMeans(BCI) == -2)) {
    currn_22[, jj] <- currn_22[nrow(currn_22), jj]
  }
  currn_22 <- cbind(currn_22, matrix(0, nrow = nrow(currn_22), ncol = (n.simu - cont.simu)))
  avg.n <- mean(colSums(currn_22))
  power <- length(which(BCI[1, ] > Ce[1])) / cont.simu # power
  GP <- length(which(BCI[1, sel.opt.g] > Ce[1])) / cont.simu

  SR <- switch(case,
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2] & BCI[3, ] > Ce[3])) / cont.simu),
    (length(which(BCI[1, ] > Ce[1])) / cont.simu),
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2])) / cont.simu)
  )

  OSR <- switch(case,
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2] & BCI[3, sel.opt.g] > Ce[3])) / cont.simu),
    (length(which(BCI[1, sel.opt.g] > Ce[1])) / cont.simu),
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2])) / cont.simu)
  )

  stage2_output <- data.frame(
    formatC(avg.n, digits = 1, format = "f"),
    formatC(power * 100, digits = 2, format = "f"),
    formatC(GP * 100, digits = 2, format = "f"),
    formatC(SR * 100, digits = 2, format = "f"),
    formatC(OSR * 100, digits = 2, format = "f")
  )
  colnames(stage2_output) <- c("EN", "Power (%)", "GP (%)", "SR (%)", "OSR (%)")

  return(list(
    stage1_output = stage1_output, stage2_output = stage2_output
  ))
}

#' @keywords internal
.solve_level <- function(rho, eff, tox) {
  c <- min(
    (rho * sqrt(eff * (1 - eff) * tox * (1 - tox)) + eff * tox),
    tox, eff
  )
  d <- min((eff - c), (1 - tox))
  a <- min((tox - c), (1 - eff))
  b <- 1 + round((c - eff - tox), 3)
  res <- round(c(a, b, c, d), 3)
  return(res)
}
