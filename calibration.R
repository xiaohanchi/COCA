

rm(list = ls())
library(rjags)
library(runjags)
library(foreach)
library(doParallel)

#### ==================================== FUNCTIONS ====================================
DoseStandardize <- function(d) {
  return(d / max(d))
}

find.order.stat <- function(order, q) {
  # q[which(q<0)]<-NA
  res <- switch(order,
    which.max(q),
    which(q == median(q))[1],
    which.min(q)
  )
  # if(length(res)==0){res<- -1}else if(is.na(res)){res<- -1}
  return(res)
}

logistic_model <- "
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

foreach.comb <- function(x, ...) {
  lapply(
    seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
  )
}

find.k.lower <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.05) {
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

find.k.upper <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.80) {
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


run.noperiod <- function(fda.case = 1, n.stage1 = 24, n.stage2, 
                         eff.ctrl, eff.A, eff.B, eff.AB, sn = 1e4) {
  ### fda.case = 1, 3, or 2 for cases 1, 2, or 3 described in the paper
  ### n.stage1: stage 1 sample size / arm
  ### n.stage2: stage 2 sample size / arm
  ### eff.ctrl/eff.A/eff.B/eff.AB: efficacy setting for control/agent A (if exists)/agent B (if exists)/AB
  ### sn: number of simulation replicates
  
  set.seed(1233)

  sn_s1 <- sn
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
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, DoseStandardize))
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

  num_core <- detectCores(logical = F)
  cl <- makeCluster(num_core - 1)
  registerDoParallel(cl)
  ### stage I
  Ye_21 <- sapply(1:sn_s1, function(r) rbinom(rep(1, ndose), n_21, prob = comb))
  q_hat <- (Ye_21 + 0.1) / (n_21 + 0.2)

  j_ast <- sapply(1:sn_s1, function(r) find.order.stat(order = 1, q = q_hat[, r]))

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
  BCI_c <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_s1) # BCI for stage II

  Ye_22 <- sapply(1:sn_s1, function(r) rbinom(rep(1, narm_22), n_22, prob = eff_22[, r]))
  BCI_c <- foreach(i = 1:sn_s1, .packages = c("runjags", "coda"), .combine = "cbind") %dopar% {
    dataYe_22 <- list(
      x1 = X1[, i], x2 = X2[, i], period = period,
      Y = c(Ye_22[, i], Ye_21[, i]),
      N_arms = (narm_22 + ndose),
      n = c(rep(n_22, narm_22), rep(n_21, ndose))
    ) # phase II data

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
    # calculate BCI in the final analysis
    BCI_c[, i] <- sapply(1:(narm_22 - 1), function(r) mean(piE_mcmc_22[, narm_22] > piE_mcmc_22[, r]))
    return(BCI_c[, i])
  }
  stopCluster(cl)
  return(BCI_c)
}


run.period <- function(fda.case = 1, n.stage1 = 24, n.stage2, type1.level = 0.20, lower_e, sn = 1e4){
  ### fda.case = 1, 3, or 2 for cases 1, 2, or 3 described in the paper
  ### n.stage1: stage 1 sample size / arm
  ### n.stage2: stage 2 sample size / arm
  ### type1.level: nomial type I error rate under non-aero period effect
  ### lower_e: null efficacy setting 
  ### sn: number of simulation replicates
  
  ndose <- 3
  fda_sc <- fda.case
  sn_s1 <- sn
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
  dose_level_std <- t(apply(dosage_level, MARGIN = 1, DoseStandardize))
  dose_std <- matrix(c(
    dose_level_std["A", 1], dose_level_std["B", 1],
    dose_level_std["A", 1], dose_level_std["B", 2],
    dose_level_std["A", 2], dose_level_std["B", 1]
  ), ncol = ndose, byrow = F)
  row.names(dose_std) <- c("A", "B")
  
  Econtrol <- singleA <- singleB <- comb <- lower_e
  period_eff <- seq(-0.1, 0.5, 0.05)
  
  type1.period <- Ce1_p <- c()
  num_core <- detectCores(logical = F)
  cl <- makeCluster(num_core - 1)
  registerDoParallel(cl)
  
  for (pp in 1:length(period_eff)) {
    set.seed(1233 + 10 * pp)
    
    comb_s1 <- comb + period_eff[pp]
    #### Stage I
    Ye_21 <- sapply(1:sn_s1, function(r) rbinom(rep(1, ndose), n_21, prob = comb_s1))
    q_hat <- (Ye_21 + 0.1) / (n_21 + 0.2)
    
    j_ast <- sapply(1:sn_s1, function(r) find.order.stat(order = 1, q = q_hat[, r]))
    
    
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
    BCI_c <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_s1) # BCI for stage II
    
    Ye_22 <- sapply(1:sn_s1, function(r) rbinom(rep(1, narm_22), n_22, prob = eff_22[, r]))
    BCI_c <- foreach(i = 1:sn_s1, .packages = c("runjags", "coda"), .combine = "cbind") %dopar% {
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
      
      BCI_c[, i] <- sapply(1:(narm_22 - 1), function(r) mean(piE_mcmc_22[, narm_22] > piE_mcmc_22[, r]))
      return(BCI_c[, i])
    }
    
    Ce1_p[pp] <- round(quantile(BCI_c[1, ], (1 - type1.level)), digits = 4)
    # type1.period[pp] <- length(which(BCI_c[1, ] > Ce1_0)) / dim(BCI_c)[2] # type I error (period) using cutoff_0
  }
  stopCluster(cl)
  return(Ce1_p)
}

one.search <- function(fda_sc, n.stage2, target_pw, lower_e){
  ### fda_sc: factorial cases 1, 2, or 3 described in the paper
  ### n.stage2: stage 2 sample size/arm
  ### target_pw: target power under zero period effect
  ### lower_e: null efficacy setting
  
  ########## Period Effect = 0 
  
  ### calibrate Ce: type I error
  BCI_null <- run.noperiod(
    fda.case = fda_sc, n.stage2 = n.stage2,
    eff.ctrl = lower_e, eff.A = lower_e, eff.B = lower_e, eff.AB = lower_e
  )
  Ce1_0 <- round(quantile(BCI_null[1, ], 0.9), digits = 4)
  
  ### power
  BCI_alt <- run.noperiod(
    fda.case = fda_sc, n.stage2 = n.stage2,
    eff.ctrl = lower_e, eff.A = 0.35, eff.B = 0.35, eff.AB = 0.55
  )
  
  power <- length(which(BCI_alt[1, ] > Ce1_0)) / dim(BCI_alt)[2] # power
  
  if (power < target_pw) {
    stop(paste0("Search stopped: Target power (", (target_pw * 100),"%) not reached."))
  }
  
  ########## Period Effect > 0 
  Ce1_p <- run.period(fda.case = fda_sc, n.stage2 = n.stage2, type1.level = 0.20, lower_e = lower_e)
  
  power_p <- length(which(BCI_alt[1, ] > max(Ce1_p))) / dim(BCI_alt)[2] # power under cutoff from period effect scenarios
  if (power_p < target_pw) {
    stop("Search stopped: Target power (90%) not reached under non-zero period effect cutoff.")
  }
  
  Ce1 <- max(Ce1_0, max(Ce1_p))
  
  Ce.k.lower <- find.k.lower(fda.case = fda_sc, k.min = 0.1, Ce.1 = Ce1, BCI = BCI_null, level = 0.05)
  if (Ce.k.lower == -1) {
    stop("Search stopped: Successful rate under H0 cannot be controlled.")
  }
  
  Ce.k <- find.k.upper(fda.case = fda_sc, k.min = Ce.k.lower, Ce.1 = Ce1, BCI = BCI_alt, level = 0.80)
  if (Ce.k == -1) {
    stop("Search stopped: Target successful rate (80%) not reached.")
  }
  return(
    list(n.stage2 = n.stage2, Ce1 = Ce1, k = Ce.k, 
         type.I.error = (length(which(BCI_null[1, ] > Ce1)) / dim(BCI_null)[2]))
  )
}

all.search <- function(fda.case = 1, n.start, n.stop, target_pw){
  ### fda.case = 1, 3, or 2 for cases 1, 2, or 3 described in the paper
  ### n.start: initial search point for the sample size in Stage 2
  ### n.stop: final search point for the sample size in Stage 2
  ### target_pw: target power under zero period effect
  for(n.stage2 in seq(n.start, n.stop, 2)){
    output.try <- try(one.search(
      fda_sc = fda.case, n.stage2 = n.stage2, target_pw = target_pw, lower_e = 0.25
    ), silent = TRUE)
    if (!("try-error" %in% class(output.try))) {
      break
    }
  }
  return(output.try)
}


#### RUN Calibration ========

all.search(n.start = 16, n.stop = 40, target_pw = 0.90)

