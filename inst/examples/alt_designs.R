rm(list = ls())
library(isotone)
library(Rcpp)
library(RcppArmadillo)
library(RcppEigen)


######################################## FUNCTIONS ########################################
sourceCpp(code='
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP MtxProd(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP expit(Eigen::MatrixXd x) {
  Eigen::MatrixXd result = (1 / (1 + (-x.array()).exp()));
  return Rcpp::wrap(result);
}

')

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

# generate prior samples
get.Beta.prior <- function(n.sample = 1e6, type = 2) {
  set.seed(0)
  # type = 1 for vague prior on beta3; type = 2 for spike and slab prior
  beta0 <- rnorm(n.sample, 0, sqrt(10))
  beta1 <- rnorm(n.sample, 0, sqrt(10))
  beta2 <- rnorm(n.sample, 0, sqrt(10))
  if (type == 1){
    beta3 <- rnorm(n.sample, 0, sqrt(10))
  } else if (type == 2){
    beta3.dist1 <- rnorm(n.sample, 0, sqrt(10))
    beta3.dist0 <- rnorm(n.sample, 0, 0.01)
    beta3.wt <- rbinom(n.sample, size = 1, prob = 0.5)
    beta3 <- beta3.wt * beta3.dist1 + (1 - beta3.wt) * beta3.dist0
  }

  dist1 <- rnorm(n.sample, 0, sqrt(10))
  dist0 <- rnorm(n.sample, 0, 0.01)
  wt <- rbinom(n.sample, size = 1, prob = 0.5)
  beta4 <- wt * dist1 + (1 - wt) * dist0
  Beta_prior <- rbind(beta0, beta1, beta2, beta3, beta4)
  return(Beta_prior)
}


# get posterior sampling
get.post <- function(Beta_prior, X.mtx, pE_prior,
                     logpE_prior0, logpE_prior1, data_Y, data_N) {

  log_likelihood <- (data_Y %*% logpE_prior0) + ((data_N - data_Y) %*% logpE_prior1)
  likelihood <- exp(log_likelihood)
  weights <- likelihood / sum(likelihood)
  idx <- sample(1:ncol(Beta_prior), size = ncol(Beta_prior), replace = TRUE, prob = weights)
  pE_post <- expit(MtxProd(X.mtx, Beta_prior[, idx]))
  return(pE_post)
}

# MTD-Ind main function
mtdind.getOC <- function(case = 1, n.s1 = 24, n.s2 = 28, scenario,
                         u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
                         C_s1 = 0.85, C_s2 = 0.85, C_t = 0.9, C_f1 = 0.9, C_f2 = 0.9,
                         prior.type = 2, seed, sn = 5000) {
  ndose <- 3
  rho <- 0.2
  Nmax_p2 <- c(n.s1, n.s2) # stage 1 & stage 2
  T_22 <- 2
  Tcontrol <- 0.10 # tox prob for control arm
  Econtrol <- 0.25
  narm_21 <- ndose
  fda_sc <- switch(case, 1, 3, 2)
  narm_22 <- switch(fda_sc,
                    4,
                    2,
                    3
  )

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = T)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, DoseStandardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = ndose, byrow = F)
  row.names(dose_std) <- c("A", "B")

  ### simulation settings
  simu_sc <- switch(scenario, 21, 1, 4, 6, 10, 14, 16, 17, 18, 20)
  Tox_prob <- Tox_prob[, simu_sc]
  Eff_prob <- Eff_prob[, simu_sc]
  Tox_prob_A <- Tox_prob_A[, simu_sc]
  Tox_prob_B <- Tox_prob_B[, simu_sc]
  Eff_prob_A <- Eff_prob_A[, simu_sc]
  Eff_prob_B <- Eff_prob_B[, simu_sc]

  N_22 <- sapply(1:sn, function(r) rep(Nmax_p2[2], narm_22)) # narm_22xsn matrix
  n_22 <- N_22 / T_22
  tox_22_all <- rbind(rep(Tcontrol, sn), rep(Tox_prob_A[mtd.dose], sn), rep(Tox_prob_B[mtd.dose], sn), rep(Tox_prob[mtd.dose], sn))
  tox_22 <- switch(fda_sc,
                   tox_22_all,
                   tox_22_all[c(1, 4), ],
                   tox_22_all[-3, ]
  )
  eff_22_all <- rbind(rep(Econtrol, sn), rep(Eff_prob_A[mtd.dose], sn), rep(Eff_prob_B[mtd.dose], sn), rep(Eff_prob[mtd.dose], sn))
  eff_22 <- switch(fda_sc,
                   eff_22_all,
                   eff_22_all[c(1, 4), ],
                   eff_22_all[-3, ]
  )

  X1_all <- c(
    dose_level_std["A", 3], dose_std["A", 1],
    dose_level_std["A", 3], dose_std["A", 1]
  )
  X2_all <- c(
    dose_level_std["B", 3], dose_level_std["B", 3],
    dose_std["B", 1], dose_std["B", 1]
  )

  X1 <- switch(fda_sc,
               X1_all,
               X1_all[c(1, 4)],
               X1_all[-3]
  )
  X2 <- switch(fda_sc,
               X2_all,
               X2_all[c(1, 4)],
               X2_all[-3]
  )

  BCI <- matrix(NA, nrow = (narm_22 - 1), ncol = sn)
  Yt_22 <- Ye_22 <- matrix(0, nrow = narm_22, ncol = sn)
  proc_22 <- matrix(1, nrow = narm_22, ncol = sn) # trial process of stage II
  currn_22 <- matrix(0, nrow = narm_22, ncol = sn)
  fprob_1 <- fprob_2 <- matrix(NA, nrow = narm_22, ncol = sn)

  Beta_prior <- get.Beta.prior(n.sample = 1e6, type = prior.type)
  Beta_prior <- Beta_prior[!rownames(Beta_prior) %in% c("beta4"), ]
  X.mtx.all <- cbind(1, X1, X2, (X1 * X2))
  pE_prior_list <- expit(MtxProd(X.mtx.all, Beta_prior))
  logpE_prior0_list <- log(pE_prior_list)
  logpE_prior1_list <- log(1 - pE_prior_list)

  for (t in 1:T_22) {
    set.seed(seed + 10 * t)
    currn_22 <- currn_22 + n_22 * proc_22 # current sample size
    Yt_22 <- Yt_22 + sapply(1:sn, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = tox_22)) * proc_22
    Ye_22 <- Ye_22 + sapply(1:sn, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = eff_22)) * proc_22
    data_Y <- t(Ye_22)
    data_N <- t(currn_22)
    for(i in 1:sn) {
      if (proc_22[narm_22, i] == 0) {
        fprob_1[, i] <- fprob_2[, i] <- -1
        if (t == T_22) BCI[, i] <- -2
      } else {
        pE_post <- get.post(
          Beta_prior = Beta_prior, X.mtx = X.mtx.all, pE_prior = pE_prior_list,
          logpE_prior0 = logpE_prior0_list, logpE_prior1 = logpE_prior1_list,
          data_Y = data_Y[i, ], data_N = data_N[i, ]
        )
        # futility stopping prob
        fprob_1[, i] <- c(100, sapply(2:narm_22, function(r) mean(pE_post[r, ] > pE_post[1, ])))
        if (fda_sc %in% c(1, 3)) { # terminate single arm A or B
          fprob_2[, i] <- c(100, sapply(2:(narm_22 - 1), function(r) mean(pE_post[r, ] > pE_post[narm_22, ])), 100)
        } else {
          fprob_2[, i] <- c(100, 100)
        }
        if (t == T_22) {
          BCI[, i] <- sapply(1:(narm_22 - 1), function(r) mean(pE_post[narm_22, ] > pE_post[r, ]))
        }
      }
    }
    # overly toxic stopping
    proc_22[which(proc_22 != 0 & pbeta(upper_t, (0.1 + Yt_22), (0.1 + currn_22 - Yt_22)) < (1 - C_t))] <- 0
    # futility stopping
    proc_22[which(proc_22 != 0 & fprob_1 < (1 - C_f1))] <- 0
    proc_22[which(proc_22 != 0 & fprob_2 < (1 - C_f2))] <- 0
    if (t == T_22) {
      BCI[, which(proc_22[narm_22, ] == 0)] <- -2
    }
  }

  output <- list(
    currn_22 = currn_22, BCI = BCI
  )
  return(output)
}

# OBD-Ind main function
obdind.getOC <- function(case = 1, n.s1 = 24, n.s2 = 28, scenario,
                         u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
                         C_s1 = 0.85, C_s2 = 0.85, C_t = 0.9, C_f1 = 0.9, C_f2 = 0.9,
                         period_eff = 0, prior.type = 2, seed, sn = 5000) {
  ndose <- 3
  rho <- 0.2
  Nmax_p2 <- c(n.s1, n.s2) # stage 1 & stage 2
  T_21 <- 2
  A1 <- matrix(c(2, 1, 3, 1), byrow = T, nrow = 2) # order matrix for pi_tox_hat
  A2 <- matrix(c(1, 2, 1, 3), byrow = T, nrow = 2) # order matrix for prob
  a <- rep(0.05, 4) # hyperparameters in dirichlet

  T_22 <- 2
  Tcontrol <- 0.10 # tox prob for control arm
  Econtrol <- 0.25
  fda_sc <- switch(case, 1, 3, 2)
  narm_22 <- switch(fda_sc,
                    4,
                    2,
                    3
  )

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = T)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, DoseStandardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = ndose, byrow = F)
  row.names(dose_std) <- c("A", "B")

  ### simulation settings
  simu_sc <- switch(scenario, 21, 1, 4, 6, 10, 14, 16, 17, 18, 20)
  Tox_prob <- Tox_prob[, simu_sc]
  Eff_prob <- Eff_prob[, simu_sc] + period_eff
  Tox_prob_A <- Tox_prob_A[, simu_sc]
  Tox_prob_B <- Tox_prob_B[, simu_sc]
  Eff_prob_A <- Eff_prob_A[, simu_sc]
  Eff_prob_B <- Eff_prob_B[, simu_sc]

  ### Stage I
  multi_prob <- sapply(1:ndose, function(r) solve.level(rho, Eff_prob[r], Tox_prob[r]))
  true_utility <- u_score %*% multi_prob

  narm_21 <- ndose # #of arms in stage I
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
    for (i in 1:ndose) {
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
    } else if (t == T_21) {
      Aset <- which(prt_21_iso * proc_21 > 1 - C_s2 & pre_21 * proc_21 < C_s2)
      for (i in 1:nrow(multi_prob)) {
        pi_hat_21_iso[, i, ][-Aset] <- -100
      }
      utility <- t(sapply(1:ndose, function(r) u_score %*% pi_hat_21_iso[r, , ]))

      j_ast1 <- apply(utility, 2, which.max) # 1-dimension indicator
      j_ast1[which(colSums(proc_21) == 0)] <- -1 # early terminated
      j_ast1[which(apply(utility, 2, max) < 0)] <- -1 # no OBD
    }
  }

  ### Stage II (Proof of Concept)
  Eff_prob <- Eff_prob - period_eff

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
  )
  eff_22_all <- rbind(rep(Econtrol, sn_22), Eff_prob_A[j_ast1_tmp], Eff_prob_B[j_ast1_tmp], Eff_prob[j_ast1_tmp])
  eff_22 <- switch(fda_sc,
                   eff_22_all,
                   eff_22_all[c(1, 4), ],
                   eff_22_all[-3, ]
  )

  X1_all <- sapply(1:narm_21, function(r) {
    c(
      dose_level_std["A", 3], dose_std["A", 1], dose_level_std["A", 3],
      dose_std["A", r]
    )
  })
  X2_all <- sapply(1:narm_21, function(r) {
    c(
      dose_level_std["B", 3], dose_level_std["B", 3], dose_std["B", 1],
      dose_std["B", r]
    )
  })
  colnames(X1_all) <- colnames(X2_all) <- c("j=1", "j=2", "j=3")
  X1 <- switch(fda_sc,
               X1_all,
               X1_all[c(1, 4), ],
               X1_all[-3, ]
  )
  X2 <- switch(fda_sc,
               X2_all,
               X2_all[c(1, 4), ],
               X2_all[-3, ]
  )

  BCI <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_22)
  Yt_22 <- Ye_22 <- matrix(0, nrow = narm_22, ncol = sn_22)
  proc_22 <- matrix(1, nrow = narm_22, ncol = sn_22) # trial process of stage II
  currn_22 <- matrix(0, nrow = narm_22, ncol = sn_22)
  fprob_1 <- fprob_2 <- matrix(NA, nrow = narm_22, ncol = sn_22)

  Beta_prior <- get.Beta.prior(n.sample = 1e6, type = prior.type)
  Beta_prior <- Beta_prior[!rownames(Beta_prior) %in% c("beta4"), ]
  X.mtx.all <- lapply(1:narm_21, function(r) {
    cbind(1, X1[, r], X2[, r], (X1[, r] * X2[, r]))
  })
  pE_prior_list <- lapply(1:narm_21, function(r) {
    expit(MtxProd(X.mtx.all[[r]], Beta_prior))
  })
  logpE_prior0_list <- lapply(1:narm_21, function(r)
    log(pE_prior_list[[r]])
  )
  logpE_prior1_list <- lapply(1:narm_21, function(r)
    log(1 - pE_prior_list[[r]])
  )


  for (t in 1:T_22) {
    set.seed(seed + 10 * t)
    currn_22 <- currn_22 + n_22 * proc_22 # current sample size
    Yt_22 <- Yt_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = tox_22[, r])) * proc_22
    Ye_22 <- Ye_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = eff_22[, r])) * proc_22
    data_Y <- t(Ye_22)
    data_N <- t(currn_22)

    for(i in 1:sn_22) {
      if (proc_22[narm_22, i] == 0) {
        fprob_1[, i] <- fprob_2[, i] <- -1
        if (t == T_22) BCI[, i] <- -2
      } else {
        pE_prior <- pE_prior_list[[j_ast1_tmp[i]]]
        logpE_prior0 <- logpE_prior0_list[[j_ast1_tmp[i]]]
        logpE_prior1 <- logpE_prior1_list[[j_ast1_tmp[i]]]
        X.mtx <- X.mtx.all[[j_ast1_tmp[i]]]
        pE_post <- get.post(
          Beta_prior = Beta_prior, X.mtx = X.mtx, pE_prior = pE_prior,
          logpE_prior0 = logpE_prior0, logpE_prior1 = logpE_prior1,
          data_Y = data_Y[i, ], data_N = data_N[i, ]
        )
        # futility stopping prob
        fprob_1[, i] <- c(100, sapply(2:narm_22, function(r) mean(pE_post[r, ] > pE_post[1, ])))
        if (fda_sc %in% c(1, 3)) { # terminate single arm A or B
          fprob_2[, i] <- c(100, sapply(2:(narm_22 - 1), function(r) mean(pE_post[r, ] > pE_post[narm_22, ])), 100)
        } else {
          fprob_2[, i] <- c(100, 100)
        }
        if (t == T_22) {
          BCI[, i] <- sapply(1:(narm_22 - 1), function(r) mean(pE_post[narm_22, ] > pE_post[r, ]))
        }
      }
    }
    # overly toxic stopping
    proc_22[which(proc_22 != 0 & pbeta(upper_t, (0.1 + Yt_22), (0.1 + currn_22 - Yt_22)) < (1 - C_t))] <- 0
    # futility stopping
    proc_22[which(proc_22 != 0 & fprob_1 < (1 - C_f1))] <- 0
    proc_22[which(proc_22 != 0 & fprob_2 < (1 - C_f2))] <- 0
    if (t == T_22) {
      BCI[, which(proc_22[narm_22, ] == 0)] <- -2
    }
  }
  output <- list(
    true_utility = true_utility, utility = utility, j_ast1 = j_ast1,
    Y_21 = Y_21, currn_21 = currn_21, piE_hat_21 = piE_hat_21,
    currn_22 = currn_22, BCI = BCI
  )
  return(output)
}

# summarize results
get.results <- function(output, case, scenario, u.score, period.eff = 0,
                        upper_t, lower_e, Ce, c0, method,
                        sn_all = 5000){
  rho <- 0.2
  fda_sc <- switch(case, 1, 3, 2)
  Ce <- c(Ce, (Ce * c0), (Ce * c0))

  BCI <- output$BCI
  if (is.null(dim(BCI)[2])) {
    BCI <- matrix(BCI, nrow = 1)
  }
  sn <- dim(BCI)[2]
  if(method == "MTD-Ind"){
    j_ast1 <- j_ast1_tmp <- rep(1, sn)
  } else if (method %in% c("OBD-Ind", "OBD-Pool")) {
    j_ast1 <- output$j_ast1
    j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]
  }
  sel.OBD <- round((table(factor(j_ast1, levels = c(-1, 1:3)))/sn_all)*100, 1)

  idx <- switch(scenario, 21, 1, 4, 6, 10, 14, 16, 17, 18, 20)
  tox.prob <- Tox_prob[, idx]
  eff.prob <- Eff_prob[, idx] + period.eff
  multi.prob <- sapply(1:length(tox.prob), function(r) {
    solve.level(rho, eff.prob[r], tox.prob[r])
  })
  true.utility <- u.score %*% multi.prob
  j_opt <- which(
    true.utility == max(
      true.utility[which(tox.prob <= upper_t &
                           eff.prob >= lower_e)]
    )
  )
  sel.opt.g <- sapply(1:sn, function(r) j_ast1_tmp[r] %in% j_opt) * 1
  sel.opt.g <- which(sel.opt.g != 0) # correct groups
  CSP <- (length(sel.opt.g) / sn_all) * 100
  currn_22 <- output$currn_22
  for (jj in which(currn_22[nrow(currn_22), ] < max(currn_22))) {
    # comb early stopping
    currn_22[, jj] <- currn_22[nrow(currn_22), jj]
  }
  currn_22 <- cbind(currn_22, matrix(0, nrow = nrow(currn_22), ncol = (sn_all - sn)))
  avg.n <- mean(colSums(currn_22))
  power <- length(which(BCI[1, ] > Ce[1])) / sn # power
  GP <- length(which(BCI[1, sel.opt.g] > Ce[1])) / sn

  # effective power only
  SR <- switch(fda_sc,
               (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2] & BCI[3, ] > Ce[3])) / sn),
               (length(which(BCI[1, ] > Ce[1])) / sn),
               (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2])) / sn)
  )

  OSR <- switch(fda_sc,
                (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2] & BCI[3, sel.opt.g] > Ce[3])) / sn),
                (length(which(BCI[1, sel.opt.g] > Ce[1])) / sn),
                (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2])) / sn)
  )

  results.stage1 <- as.data.frame(t(c(scenario, period.eff, c(sel.OBD), CSP)))
  names(results.stage1) <- c("Scenario", "Period Effect", "sel.none (%)", "sel.dose1 (%)", "sel.dose2 (%)","sel.dose3 (%)", "CSP (%)")

  results.stage2 <- data.frame(
    case, scenario, period.eff, round(avg.n, 1), round(power * 100, 2),
    round(GP * 100, 2), round(SR * 100, 2), round(OSR * 100, 2)
  )
  colnames(results.stage2) <- c("Case", "Scenario", "Period Effect", "EN", "Power (%)", "Generalized Power (%)",  "SR (%)", "OSR (%)")

  return(list(results.stage1 = results.stage1, results.stage2 = results.stage2))
}

############################################ RUN models ##############################################
scenario.path <- system.file("examples", "scenarios.R", package = "COCA")
source(scenario.path)


### MTD-Ind ===================================================================================

# To run the MTD-Ind design for Scenario 1 under Case 1, use:
mtdind.output <- mtdind.getOC(
  case = 1, n.s1 = 24, n.s2 = 40, scenario = 1,
  u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
  seed = 1234, sn = 5000
)

# Note: This code may take approximately 1–2 hours to run.
# For illustration, we provide a pre-generated output file located in `/inst/example_Rdata`.
# Assuming the code above has been successfully run, the function would return:
# (load the saved output as follows)

load(system.file("example_Rdata", "output_mtdind.Rdata", package = "COCA"))
output_mtdind <- output

# This creates an object named `output_mtdind`, which contains the results returned by the function `mtdind.getOC`. To get the results for Scenario 1 under Case 1 of the MTD-Ind design, run:
# (The corresponding Ce and c0 values can be found in Table 1 of the paper)

get.results(
  output = output_mtdind, case = 1, scenario = 1, u.score = c(0, 60, 40, 100),
  period.eff = 0, upper_t = 0.35, lower_e = 0.25, Ce = 0.9049, c0 = 0.81,
  method = "MTD-Ind", sn_all = 5000
  )

# To reproduce other scenarios and cases from the paper, run the above functions with `case` set to 1, 2, or 3, and `scenario` set to 1 through 10.


### OBD-Ind ===================================================================================

# To run the OBD-Ind design for Scenario 1 under Case 1, use:
obdind.output <- obdind.getOC(
  case = 1, n.s1 = 24, n.s2 = 40, scenario = 1,
  u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
  period_eff = 0, seed = 1274, sn = 5000
)

# Note: This code may take approximately 1–2 hours to run.
# For illustration, we provide a pre-generated output file located in `/inst/example_Rdata`.
# Assuming the code above has been successfully run, the function would return:
# (load the saved output as follows)

load(system.file("example_Rdata", "output_obdind.Rdata", package = "COCA"))
output_obdind <- output

# This creates an object named `output_obdind`, which contains the results returned by the function `obdind.getOC`. To get the results for Scenario 1 under Case 1 of the OBD-Ind design, run:
# (The corresponding Ce and c0 values can be found in Table 1 of the paper)

get.results(
  output = output_obdind, case = 1, scenario = 1, u.score = c(0, 60, 40, 100),
  period.eff = 0, upper_t = 0.35, lower_e = 0.25, Ce = 0.9049, c0 = 0.81,
  method = "OBD-Ind", sn_all = 5000
)

# To reproduce other scenarios and cases from the paper, run the above functions with `case` set to 1, 2, or 3, and `scenario` set to 1 through 10.


### OBD-Pool ===================================================================================

# To run the OBD-Ind design for Scenario 1 under Case 1, use:
obdind.output <- obdind.getOC(
  case = 1, n.s1 = 24, n.s2 = 40, scenario = 1,
  u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
  period_eff = 0, seed = 1274, sn = 5000
)

# Note: This code may take approximately 1–2 hours to run.
# For illustration, we provide a pre-generated output file located in `/inst/example_Rdata`.
# Assuming the code above has been successfully run, the function would return:
# (load the saved output as follows)

load(system.file("example_Rdata", "output_obdind.Rdata", package = "COCA"))
output_obdind <- output

# This creates an object named `output_obdind`, which contains the results returned by the function `obdind.getOC`. To get the results for Scenario 1 under Case 1 of the OBD-Ind design, run:
# (The corresponding Ce and c0 values can be found in Table 1 of the paper)

get.results(
  output = output_obdind, case = 1, scenario = 1, u.score = c(0, 60, 40, 100),
  period.eff = 0, upper_t = 0.35, lower_e = 0.25, Ce = 0.9049, c0 = 0.81,
  method = "OBD-Ind", sn_all = 5000
)

# To reproduce other scenarios and cases from the paper, run the above functions with `case` set to 1, 2, or 3, and `scenario` set to 1 through 10.



