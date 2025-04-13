
library(readr)
library(filelock)
library(rjags)
library(runjags)

# ====================================FUNCTIONS=================================#
.dose_standardize <- function(d) {
  return(d / max(d))
}

.find_order_stat <- function(order, q) {
  res <- switch(order,
    which.max(q),
    which(q == median(q))[1],
    which.min(q)
  )
  return(res)
}


logistic_model="
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
  wt ~ dbern(0.5)
}
"

.find_klower <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.05) {
  k0 <- seq(k.min, 1, 0.05)
  if (fda.case == 1) {
    typei.eff <- sapply(1:length(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]) & BCI[3, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(typei.eff <= level)])
  } else if (fda.case == 2) {
    k <- 1
  } else if (fda.case == 3) {
    typei.eff <- sapply(1:length(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(typei.eff <= level)])
  }
  if (is.infinite(k)) {
    k <- -1
  }
  return(k)
}

.find_kupper <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.80) {
  k0 <- seq(k.min, 1, 0.05)
  if (fda.case == 1) {
    pw.eff <- sapply(1:length(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]) & BCI[3, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- max(k0[which(pw.eff >= level)])
  } else if (fda.case == 2) {
    k <- 1
  } else if (fda.case == 3) {
    pw.eff <- sapply(1:length(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- max(k0[which(pw.eff >= level)])
  }
  if (is.infinite(k)) {
    k <- -1
  }
  return(k)
}


run.whole <- function(fda.case = 1, n.stage1 = 24, n.stage2,
                      eff.ctrl, eff.A, eff.B, eff.AB, batch.idx, batch.sn = 100) {
  set.seed(1233 + 10 * batch.idx)

  sn_s1 <- batch.sn
  fda_sc <- fda.case
  ndose <- 3
  n_21 <- n.stage1
  n_22 <- n.stage2

  narm_22 <- switch(fda_sc,
    4,
    2,
    3
  )
  jags_params_22 <- c("pi")
  period <- switch(fda_sc,
    c(0, 0, 0, 0, 1, 1, 1),
    c(0, 0, 0, 0, 1, 1, 1)[c(1, 4, 5, 6, 7)],
    c(0, 0, 0, 0, 1, 1, 1)[-3]
  )
  logistic_model <- logistic_model

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = T)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, .dose_standardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = ndose, byrow = F)
  row.names(dose_std) <- c("A", "B")

  Econtrol <- eff.ctrl
  singleA <- eff.A
  singleB <- eff.B
  comb <- eff.AB

  ### stage I
  Ye_21 <- sapply(1:sn_s1, function(r) rbinom(rep(1, ndose), n_21, prob = comb))
  q_hat <- (Ye_21 + 0.1) / (n_21 + 0.2)

  j_ast <- sapply(1:sn_s1, function(r) .find_order_stat(order = 1, q = q_hat[, r]))

  X1_all <- sapply(1:sn_s1, function(r) {
    c(
      dose_level_std["A", 3], dose_std["A", 1], dose_level_std["A", 3],
      dose_std["A", j_ast[r]], dose_std["A", ]
    )
  }) # dose_level_std["A",3]=0
  X2_all <- sapply(1:sn_s1, function(r) {
    c(
      dose_level_std["B", 3], dose_level_std["B", 3], dose_std["B", 1],
      dose_std["B", j_ast[r]], dose_std["B", ]
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

  #### Stage II
  eff_22_all <- rbind(rep(Econtrol, sn_s1), rep(singleA, sn_s1), rep(singleB, sn_s1), rep(comb, sn_s1))
  eff_22 <- switch(fda_sc,
    eff_22_all,
    eff_22_all[c(1, 4), ],
    eff_22_all[-3, ]
  )
  BCI_c_batch <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_s1) # BCI for stage II
  Ye_22 <- sapply(1:sn_s1, function(r) rbinom(rep(1, narm_22), n_22, prob = eff_22[, r]))

  for (i in 1:sn_s1) {
    dataYe_22 <- list(
      x1 = X1[, i], x2 = X2[, i], period = period,
      Y = c(Ye_22[, i], Ye_21[, i]),
      N_arms = (narm_22 + ndose),
      n = c(rep(n_22, narm_22), rep(n_21, ndose))
    ) # phase II data

    jagsmodel.Ye_22 <- run.jags(
      model = logistic_model, monitor = jags_params_22, data = dataYe_22,
      n.chains = 4, adapt = 2000, burnin = 3000,
      sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
      plots = FALSE, silent.jags = T
    )
    codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
    sumYe_22 <- summary(codasamples.Ye_22)
    piE_mcmc_22 <- matrix(NA, nrow = (jagsmodel.Ye_22$sample * length(jagsmodel.Ye_22$mcmc)), ncol = narm_22)
    for (j in 1:narm_22) {
      piE_mcmc_22[, j] <- as.matrix(codasamples.Ye_22[, j])
    }
    # calculate BCI in the final analysis
    BCI_c_batch[, i] <- sapply(1:(narm_22 - 1), function(r) mean(piE_mcmc_22[, narm_22] > piE_mcmc_22[, r]))
  }
  return(t(BCI_c_batch))
}

run.period <- function(fda.case = 1, n.stage1 = 24, n.stage2,
                       eff.ctrl, eff.A, eff.B, eff.AB, period_eff = 0,
                       batch.idx, batch.sn = 100) {
  set.seed(1233 + 10 * batch.idx + 100 * period_eff)

  sn_s1 <- batch.sn
  fda_sc <- fda.case
  ndose <- 3
  n_21 <- n.stage1
  n_22 <- n.stage2

  narm_22 <- switch(fda_sc,
                    4,
                    2,
                    3
  )
  jags_params_22 <- c("pi")
  period <- switch(fda_sc,
                   c(0, 0, 0, 0, 1, 1, 1),
                   c(0, 0, 0, 0, 1, 1, 1)[c(1, 4, 5, 6, 7)],
                   c(0, 0, 0, 0, 1, 1, 1)[-3]
  )
  logistic_model <- logistic_model

  dosage_level <- matrix(c(300, 200, 0, 300, 200, 0), nrow = 2, byrow = T)
  row.names(dosage_level) <- c("A", "B")
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, .dose_standardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = ndose, byrow = F)
  row.names(dose_std) <- c("A", "B")

  Econtrol <- eff.ctrl
  singleA <- eff.A
  singleB <- eff.B
  comb <- eff.AB

  comb_s1 <- comb + period_eff
  #### Stage I
  Ye_21 <- sapply(1:sn_s1, function(r) rbinom(rep(1, ndose), n_21, prob = comb_s1))
  q_hat <- (Ye_21 + 0.1) / (n_21 + 0.2)

  j_ast <- sapply(1:sn_s1, function(r) .find_order_stat(order = 1, q = q_hat[, r]))


  X1_all <- sapply(1:sn_s1, function(r) {
    c(
      dose_level_std["A", 3], dose_std["A", 1], dose_level_std["A", 3],
      dose_std["A", j_ast[r]], dose_std["A", ]
    )
  }) # dose_level_std["A",3]=0
  X2_all <- sapply(1:sn_s1, function(r) {
    c(
      dose_level_std["B", 3], dose_level_std["B", 3], dose_std["B", 1],
      dose_std["B", j_ast[r]], dose_std["B", ]
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

  eff_22_all <- rbind(rep(Econtrol, sn_s1), rep(singleA, sn_s1), rep(singleB, sn_s1), rep(comb, sn_s1))
  eff_22 <- switch(fda_sc,
                   eff_22_all,
                   eff_22_all[c(1, 4), ],
                   eff_22_all[-3, ]
  )
  BCI_c_batch <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_s1) # BCI for stage II

  Ye_22 <- sapply(1:sn_s1, function(r) rbinom(rep(1, narm_22), n_22, prob = eff_22[, r]))

  for (i in 1:sn_s1) {
    dataYe_22 <- list(
      x1 = X1[, i], x2 = X2[, i], period = period,
      Y = c(Ye_22[, i], Ye_21[, i]),
      N_arms = (narm_22 + ndose),
      n = c(rep(n_22, narm_22), rep(n_21, ndose))
    )
    jagsmodel.Ye_22 <- run.jags(
      model = logistic_model, monitor = jags_params_22, data = dataYe_22,
      n.chains = 4, adapt = 2000, burnin = 5000,
      sample = 5000, summarise = FALSE, thin = 1, method = "rjags",
      plots = FALSE, silent.jags = T
    )
    codasamples.Ye_22 <- as.mcmc.list(jagsmodel.Ye_22)
    sumYe_22 <- summary(codasamples.Ye_22)
    piE_mcmc_22 <- matrix(NA, nrow = (jagsmodel.Ye_22$sample * length(jagsmodel.Ye_22$mcmc)), ncol = narm_22)
    for (j in 1:narm_22) {
      piE_mcmc_22[, j] <- as.matrix(codasamples.Ye_22[, j])
    }

    BCI_c_batch[, i] <- sapply(1:(narm_22 - 1), function(r) mean(piE_mcmc_22[, narm_22] > piE_mcmc_22[, r]))
  }
  return(t(BCI_c_batch))

}

### Settings =========
all.config <- expand.grid(
  fda_sc = 1:3,
  n.stage2 = seq(18, 36, 2)
)
lower_e <- 0.25
period_eff <- seq(-0.1, 0.5, 0.05)
# type1.level <- 0.20

### RUN code =======
BCI_null <- run.whole(
  fda.case = all.config$fda_sc[sc00], n.stage2 = all.config$n.stage2[sc00],
  eff.ctrl = lower_e, eff.A = lower_e, eff.B = lower_e, eff.AB = lower_e,
  batch.idx = rep00, batch.sn = 20
)

BCI_alt <- run.whole(
  fda.case = all.config$fda_sc[sc00], n.stage2 = all.config$n.stage2[sc00],
  eff.ctrl = lower_e, eff.A = 0.35, eff.B = 0.35, eff.AB = 0.55,
  batch.idx = rep00, batch.sn = 20
)
BCI_period <- list()
for (pp in 1:length(period_eff)) {
  BCI_period[[pp]] <- run.period(
    fda.case = all.config$fda_sc[sc00], n.stage2 = all.config$n.stage2[sc00],
    eff.ctrl = lower_e, eff.A = lower_e, eff.B = lower_e, eff.AB = lower_e,
    period_eff = period_eff[pp],
    batch.idx = rep00, batch.sn = 20
  )
}

