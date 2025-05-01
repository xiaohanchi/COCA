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

# main function
mtdind.getOC <- function(case = 1, n.s1 = 24, n.s2 = 28, scenario,
                         u_score = c(0, 60, 40, 100),upper_t = 0.35, lower_e = 0.25,
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

##################################### RUN models #######################################
scenario.path <- system.file("examples", "scenarios.R", package = "COCA")
source(scenario.path)

mtdind.output <- mtdind.getOC(
  case = 1, n.s1 = 24, n.s2 = 40, scenario = 1,
  u_score = c(0, 60, 40, 100), upper_t = 0.35, lower_e = 0.25,
  seed = 1234, sn = 5000
)

# Note: this code may take ~1-2 hours to run
# For illustration, we have provided an example output .Rdata file in `/inst/example_Rdata`
# To get the results, assume we have successfully run the above code, and the function returns:

load(system.file("example_Rdata", "output_mtdind.Rdata", package = "COCA"))



