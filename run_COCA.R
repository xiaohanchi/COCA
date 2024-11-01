rm(list = ls())

library(isotone)
library(rjags)
library(runjags)
library(doParallel)
library(foreach)

source("scenarios.R")

######################################## FUNCTIONS ########################################
DoseStandardize <- function(d) {
  return(d / max(d))
}

solve.level <- function(rho, eff, tox) {
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

docia_dosage_model <- "
model{
  for(k in 1:N_arms){
    Y[k]~dbin(pi[k],n[k])
    logit(pi[k])<-beta0+beta1*x1[k]+beta2*x2[k]+beta3*x1[k]*x2[k]+beta4*x1[k]*x2[k]*period[k]
  }
  beta0 ~ dnorm(0,0.1)
  beta1 ~ dnorm(0,0.1)
  beta2 ~ dnorm(0,0.1)
  beta3 ~ dnorm(0,0.1)
  beta4 <- wt * dist1 + (1 - wt) * dist0
  dist1 ~ dnorm(0,0.1)
  dist0 ~ dnorm(0,1e4)
  wt ~ dbern(0.50)
}
"

pool_model <- "
model{
  for(k in 1:N_arms){
    Y[k]~dbin(pi[k],n[k])
    logit(pi[k])<-beta0+beta1*x1[k]+beta2*x2[k]+beta3*x1[k]*x2[k]
  }
  beta0 ~ dnorm(0,0.1)
  beta1 ~ dnorm(0,0.1)
  beta2 ~ dnorm(0,0.1)
  beta3 ~ dnorm(0,0.1)
}
"

optpoc_model <- "
model{
  for(k in 1:N_arms){
    Y[k]~dbin(pi[k],n[k])
    logit(pi[k])<-beta0+beta1*x1[k]+beta2*x2[k]+beta3*x1[k]*x2[k]
  }
  beta0 ~ dnorm(0,0.1)
  beta1 ~ dnorm(0,0.1)
  beta2 ~ dnorm(0,0.1)
  beta3 ~ dnorm(0,0.1)
}
"

comb <- function(x, ...) {
  lapply(
    seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
  )
}

get.COCAresult <- function(output.Rdata, scenario.idx) {
  sn_all <- 5000
  Nmax <- c(24, switch(fda_sc,
    26,
    30,
    28
  ))

  Ce <- rep(0, 3)
  Ce[1] <- switch(fda_sc,
    0.8983,
    0.9193,
    0.9068
  )
  Ce[2] <- Ce[3] <- switch(fda_sc,
    0.7 * Ce[1],
    1 * Ce[1],
    0.85 * Ce[1]
  )

  load(output.Rdata)
  q_hat_21 <- piE_hat_21[, which(j_ast1 > 0)]
  sn <- dim(BCI)[2]
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]

  tox.prob <- Tox_prob[, scenario.idx]
  eff.prob <- Eff_prob[, scenario.idx]
  multi.prob <- sapply(1:length(tox.prob), function(r) {
    solve.level(rho, eff.prob[r], tox.prob[r])
  })
  true.utility <- u_score %*% multi.prob
  j_opt <- which(
    true.utility == max(
      true.utility[which(tox.prob <= upper_t &
        eff.prob >= lower_e)]
    )
  )
  sel.opt.g <- sapply(1:sn, function(r) j_ast1_tmp[r] %in% j_opt) * 1
  sel.opt.g <- which(sel.opt.g != 0) # correct groups

  Selection <- table(j_ast1) / sn_all
  CSP <- length(sel.opt.g) / sn_all
  for (jj in which(colMeans(BCI) == -2)) {
    currn_22[, jj] <- currn_22[nrow(currn_22), jj]
  }
  currn_22 <- cbind(currn_22, matrix(0, nrow = nrow(currn_22), ncol = (sn_all - sn)))
  avg.n <- mean(colSums(currn_22))
  # effective + considerable power
  all.power <- length(which(BCI[1, ] > Ce[1])) / sn # power
  opt.power <- length(which(BCI[1, sel.opt.g] > Ce[1])) / length(sel.opt.g) # specific power
  generalized.power <- length(which(BCI[1, sel.opt.g] > Ce[1])) / sn

  # effective power only
  all.effpower <- switch(fda_sc,
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2] & BCI[3, ] > Ce[3])) / sn),
    (length(which(BCI[1, ] > Ce[1])) / sn),
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2])) / sn)
  )
  opt.effpower <- switch(fda_sc,
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2] & BCI[3, sel.opt.g] > Ce[3])) / length(sel.opt.g)),
    (length(which(BCI[1, sel.opt.g] > Ce[1])) / length(sel.opt.g)),
    length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2])) / length(sel.opt.g)
  )

  generalized.effpower <- switch(fda_sc,
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2] & BCI[3, sel.opt.g] > Ce[3])) / sn),
    (length(which(BCI[1, sel.opt.g] > Ce[1])) / sn),
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2])) / sn)
  )


  final.output <- list(
    Selection = Selection * 100,
    CSP = round(CSP * 100, 2),
    EN = round(avg.n, 1),
    Power = round(all.power * 100, 2),
    GP = round(generalized.power * 100, 2),
    SR = round(all.effpower * 100, 2),
    OSR = round(generalized.effpower * 100, 2)
  )
  return(final.output)
}
##################################### MAIN functions: COCA  #######################################
run.COCA.period <- function(sn = 5000, fda_sc = 1, n.s1 = 24, n.s2 = 28, simu_sc, period_eff = 0.2) {
  ### sn: number of simulation replicates
  ### fda_sc: factorial cases 1, 2, or 3 described in the paper
  ### n.s1/n.s2: stage 1/2 sample size per arm
  ### simu_sc: simulation scenario index
  ### period_eff: period effect setting
  Ndoses_ph1 <- 3

  Nmax_p2 <- c(n.s1, n.s2) # stage 1 & stage 2
  T_21 <- 2
  A1 <- matrix(c(2, 1, 3, 1), byrow = T, nrow = 2) # order matrix for pi_tox_hat
  A2 <- matrix(c(1, 2, 1, 3), byrow = T, nrow = 2) # order matrix for prob
  a <- rep(0.05, 4) # hyperparameters in dirichlet

  T_22 <- 2
  Tcontrol <- 0.10 # tox prob for control arm
  Econtrol <- 0.25

  narm_22 <- switch(fda_sc,
    4,
    2,
    3
  )
  jags_params_22 <- c("pi", "wt")
  logistic_model <- docia_dosage_model
  period <- switch(fda_sc,
    c(0, 0, 0, 0, 1, 1, 1),
    c(0, 0, 0, 0, 1, 1, 1)[c(1, 4, 5, 6, 7)],
    c(0, 0, 0, 0, 1, 1, 1)[-3]
  )

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = T)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, DoseStandardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = Ndoses_ph1, byrow = F)
  row.names(dose_std) <- c("A", "B")

  Tox_prob <- Tox_prob[, simu_sc]
  Eff_prob <- Eff_prob[, simu_sc] + period_eff
  Tox_prob_A <- Tox_prob_A[, simu_sc]
  Tox_prob_B <- Tox_prob_B[, simu_sc]
  Eff_prob_A <- Eff_prob_A[, simu_sc]
  Eff_prob_B <- Eff_prob_B[, simu_sc]

  ### Stage I
  multi_prob <- sapply(1:Ndoses_ph1, function(r) solve.level(rho, Eff_prob[r], Tox_prob[r]))
  true_utility <- u_score %*% multi_prob

  narm_21 <- Ndoses_ph1 # #of arms in stage I
  N_21 <- sapply(1:sn, function(r) rep(Nmax_p2[1], 3)) # 3xsn matrix
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

  Y_21 <- pi_hat_21 <- pi_hat_21_iso <- array(0, dim = c(narm_21, nrow(multi_prob), sn))
  piT_hat_21 <- piE_hat_21 <- currn_21 <- matrix(0, nrow = 3, ncol = sn) # 3xsn matrix
  proc_21 <- matrix(rep(((tox_21 + eff_21) != 0) * 1, sn), nrow = 3) # trial process of stage I
  BCI <- vector(mode = "numeric", sn)
  for (t in 1:T_21) {
    set.seed(123 + 10 * t)
    currn_21 <- currn_21 + n_21 * proc_21 # current sample size
    temp_Y <- sapply(1:narm_21, function(r) rmultinom(sn, n_21[1, r], multi_prob[, r]))
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
    piT_hat_21_iso <- sapply(1:sn, function(r) activeSet(A1, "LS", weights = w1, y = piT_hat_21[, r])$x)
    piE_hat_21 <- pi_hat_21[, 3, ] + pi_hat_21[, 4, ] # response rate hat
    rho_hat <- (pi_hat_21[, 2, ] * pi_hat_21[, 3, ] - pi_hat_21[, 1, ] * pi_hat_21[, 4, ]) / sqrt(piE_hat_21 * (1 - piE_hat_21) * piT_hat_21_iso * (1 - piT_hat_21_iso))
    for (i in 1:Ndoses_ph1) {
      pi_hat_21_iso[i, , ] <- sapply(1:sn, function(r) solve.level(rho_hat[i, r], piE_hat_21[i, r], piT_hat_21_iso[i, r]))
    }
    prt_21 <- pbeta(upper_t, (a[1] + a[3] + Y_21[, 1, ] + Y_21[, 3, ]), (a[2] + a[4] + currn_21 - Y_21[, 1, ] - Y_21[, 3, ]))
    w2 <- 1 / apply(prt_21, 1, var)
    prt_21_iso <- sapply(1:sn, function(r) activeSet(A2, "LS", weights = w2, y = prt_21[, r])$x)
    pre_21 <- pbeta(lower_e, (a[3] + a[4] + Y_21[, 3, ] + Y_21[, 4, ]), (a[1] + a[2] + currn_21 - Y_21[, 3, ] - Y_21[, 4, ]))
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
      utility <- t(sapply(1:Ndoses_ph1, function(r) u_score %*% pi_hat_21_iso[r, , ]))

      j_ast1 <- apply(utility, 2, which.max) # 1-dimension indicator
      j_ast1[which(colSums(proc_21) == 0)] <- -1 # early terminated
      j_ast1[which(apply(utility, 2, max) < 0)] <- -1 # no OBD
    }
    # piE_hat_21[proc_21==0]<- -1
  }

  ### Stage II (Proof of Concept)
  Eff_prob <- Eff_prob - period_eff

  BCI <- matrix(rep(BCI, (narm_22 - 1)), nrow = (narm_22 - 1), byrow = T)
  Yt_21 <- (Y_21[, 1, ] + Y_21[, 3, ])[, which(j_ast1 > 0)]
  Ye_21 <- (Y_21[, 3, ] + Y_21[, 4, ])[, which(j_ast1 > 0)]
  sn_22 <- sum(j_ast1 > 0)
  N_22 <- sapply(1:sn_22, function(r) rep(Nmax_p2[2], narm_22)) # narm_22xsn matrix
  n_22 <- N_22 / T_22
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]
  tox_22_all <- rbind(rep(Tcontrol, sn_22), Tox_prob_A[j_ast1_tmp], Tox_prob_B[j_ast1_tmp], Tox_prob[j_ast1_tmp])
  tox_22 <- switch(fda_sc,
    tox_22_all,
    tox_22_all[c(1, 4), ],
    tox_22_all[-3, ]
  ) # scenario 1-3 C/A/B/AB
  eff_22_all <- rbind(rep(Econtrol, sn_22), Eff_prob_A[j_ast1_tmp], Eff_prob_B[j_ast1_tmp], Eff_prob[j_ast1_tmp])
  eff_22 <- switch(fda_sc,
    eff_22_all,
    eff_22_all[c(1, 4), ],
    eff_22_all[-3, ]
  )

  X1_all <- sapply(1:sn_22, function(r) {
    c(
      dose_level_std["A", 3], dose_std["A", 1], dose_level_std["A", 3],
      dose_std["A", j_ast1_tmp[r]], dose_std["A", ]
    )
  }) # dose_level_std["A",3]=0
  X2_all <- sapply(1:sn_22, function(r) {
    c(
      dose_level_std["B", 3], dose_level_std["B", 3], dose_std["B", 1],
      dose_std["B", j_ast1_tmp[r]], dose_std["B", ]
    )
  })
  X1 <- switch(fda_sc,
    X1_all,
    X1_all[c(1, 4, 5, 6, 7), ],
    X1_all[-3, ]
  )
  X2 <- switch(fda_sc,
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

  cl <- makeCluster(num_core - 1)
  registerDoParallel(cl)
  for (t in 1:T_22) {
    set.seed(1234 + 10 * t)
    currn_22 <- currn_22 + n_22 * proc_22 # current sample size
    Yt_22 <- Yt_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = tox_22[, r])) * proc_22
    Ye_22 <- Ye_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = eff_22[, r])) * proc_22

    if (t < T_22) {
      output <- foreach(
        i = 1:sn_22, .packages = c("runjags", "coda"), .combine = "comb", .multicombine = TRUE,
        .init = list(list(), list())
      ) %dopar% {
        if (proc_22[narm_22, i] == 0) {
          fprob_1[, i] <- fprob_2[, i] <- -1
        } else { # skip if AB has been terminated
          dataYe_22 <- list(
            x1 = X1[, i], x2 = X2[, i], period = period,
            Y = c(Ye_22[, i], Ye_21[, i]),
            N_arms = (narm_22 + Ndoses_ph1),
            n = c(currn_22[, i], currn_21[, which(j_ast1 > 0)][, i])
          ) # phase II data
          jagsmodel.Ye_22 <- run.jags(
            model = logistic_model, monitor = jags_params_22, data = dataYe_22,
            n.chains = 4, adapt = 2000, burnin = 5000,
            sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
            plots = FALSE, silent.jags = T
          )
          codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
          piE_mcmc_22 <- matrix(NA, nrow = (jagsmodel.Ye_22$sample * length(jagsmodel.Ye_22$mcmc)), ncol = narm_22)
          for (j in 1:narm_22) {
            piE_mcmc_22[, j] <- as.matrix(codasamples.Ye_22[, j])
          }
          # futility stopping prob
          fprob_1[, i] <- c(100, sapply(2:narm_22, function(r) mean(piE_mcmc_22[, r] > piE_mcmc_22[, 1]))) # omit the 1st arm
          if (fda_sc %in% c(1, 3)) { # terminate single arm A or B
            fprob_2[, i] <- c(100, sapply(2:(narm_22 - 1), function(r) mean(piE_mcmc_22[, r] > piE_mcmc_22[, narm_22])), 100)
          } else {
            fprob_2[, i] <- c(100, 100)
          }
        }
        result <- list(fprob_1[, i], fprob_2[, i])
        return(result)
      }
      fprob_1 <- matrix(unlist(output[[1]]), nrow = narm_22, ncol = sn_22)
      fprob_2 <- matrix(unlist(output[[2]]), nrow = narm_22, ncol = sn_22)
      # overly toxic stopping
      proc_22[which(proc_22 != 0 & pbeta(upper_t, (0.1 + Yt_pre + Yt_22), (0.1 + Nt_pre + currn_22 - Yt_pre - Yt_22)) < C_t)] <- 0
      # futility stopping
      proc_22[which(proc_22 != 0 & fprob_1 < C_f1)] <- 0
      proc_22[which(proc_22 != 0 & fprob_2 < C_f2)] <- 0
      # terminate the trial and conclude ineffective
      # to represent the results for early termination of arm AB
      BCI[, which(proc_22[narm_22, ] == 0)] <- -2
    } else if (t == T_22) {
      output <- foreach(
        i = 1:sn_22, .packages = c("runjags", "coda"), .combine = "comb", .multicombine = TRUE,
        .init = list(list(), list())
      ) %dopar% {
        if (proc_22[narm_22, i] == 0) {
          BCI[, i] <- -2
        } else {
          dataYe_22 <- list(
            x1 = X1[, i], x2 = X2[, i], period = period,
            Y = c(Ye_22[, i], Ye_21[, i]),
            N_arms = (narm_22 + Ndoses_ph1),
            n = c(currn_22[, i], currn_21[, which(j_ast1 > 0)][, i])
          ) # phase II data
          jagsmodel.Ye_22 <- run.jags(
            model = logistic_model, monitor = jags_params_22, data = dataYe_22,
            n.chains = 4, adapt = 2000, burnin = 5000,
            sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
            plots = FALSE, silent.jags = T
          )
          codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
          sumYe_22 <- summary(codasamples.Ye_22)
          wt[i] <- switch(fda_sc,
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
        result <- list(wt[i], BCI[, i])
        return(result)
      }
      wt <- unlist(output[[1]])
      BCI <- matrix(unlist(output[[2]]), nrow = (narm_22 - 1), ncol = sn_22)
    }
  }
  stopCluster(cl)
  save(true_utility, utility, j_ast1, Y_21, currn_21, piE_hat_21,
    currn_22, wt, BCI,
    file = "COCA_output.Rdata"
  )
}


##################################### RUN models #######################################
num_core <- detectCores(logical = F)
upper_t <- 0.35 # overly toxic stopping
lower_e <- 0.25 # ineffective stopping
C_s1 <- 0.85
C_s2 <- 0.85
u_score <- c(0, 60, 40, 100)
rho <- 0.2

C_t <- 0.10
C_f1 <- 0.10
C_f2 <- 0.10

run.COCA.period(fda_sc = 1, n.s2 = 26, simu_sc = 1)
get.COCAresult(output.Rdata = "COCA_output.Rdata", scenario.idx = 1)
