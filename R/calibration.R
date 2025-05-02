#' Calibrate COCA parameters
#'
#' Computes the optimal design cutoffs (\eqn{C_{e1}} and \eqn{c_0}) and the corresponding power across a range of stage 2 sample sizes for a given configuration.
#'
#' @details
#' For each candidate value in `n.stage2`, the function computes the calibrated cutoffs \eqn{C_{e1}} and \eqn{c_0}, along with power and type I error rates. It returns the full results across all sample sizes considered, as well as the optimal configuration that first meets the target power and validity criteria.
#'
#' @param case Trial type for stage 2. \code{case = 1} for 4-arm trial comparing AB vs. A vs. B vs. SOC; \code{case = 2} for 3-arm trial comparing AB vs. A (or B) vs. SOC; \code{case = 3} for 2-arm trial comparing AB vs. SOC.
#' @param n.comb.dose Number of combination arms in stage 1
#' @param n.stage1 Sample size for stage 1
#' @param n.stage2 A numeric vector of candidate stage 2 sample sizes to evaluate during calibration.
#' @param dosage.ctrl Dosage level of the control arm for stage 2. For an SOC control, use \code{c(A = 0, B = 0)}. If one of the single agents is used as the SOC (e.g., drug A 300 mg), use \code{c(A = 300, B = 0)}.
#' @param dosage.singleA Dosage level of drug A in the single arm for stage 2.
#' @param dosage.singleB Dosage level of drug B in the single arm for stage 2.
#' @param dosage.comb A named list specifying the dosage levels of drugs A and B across combination arms in stage 1. For example, \code{list(A = c(300, 300, 200), B = c(300, 200, 300))} defines three dose combinations:
#'   dose1 = (A = 300, B = 300), dose2 = (A = 300, B = 200), and dose3 = (A = 200, B = 300).
#'   All dosage values can be on any scale but must use the same scale across \code{dosage.singleA}, \code{dosage.singleB}, and \code{dosage.comb}.
#' @param eff.null Unpromising efficacy rate (\eqn{\widetilde{q}_1}) in the global null hypothesis
#' @param eff.alt.SOC Efficacy rate of the SOC arm (\eqn{\widetilde{q}_{SOC}}) in the alternative hypothesis
#' @param eff.alt.A Efficacy rate of single arm A (\eqn{\widetilde{q}_{A}}) in the alternative hypothesis. This argument is ignored if arm A is not included.
#' @param eff.alt.B Efficacy rate of single arm B (\eqn{\widetilde{q}_{B}}) in the alternative hypothesis. This argument is ignored if arm B is not included.
#' @param eff.alt.AB Efficacy rate of combination arm AB (\eqn{\widetilde{q}_{AB}}) in the alternative hypothesis.
#' @param period.effect Vector of possible period effects
#' @param alpha.level Type I error level (\eqn{\alpha_0}) under no period effect assumption
#' @param alpha.max Maximum type I error level (\eqn{\alpha_{\max}}) under non-zero period effect
#' @param fsr.level False success rate level (\eqn{\gamma_0}) under no period effect assumption
#' @param power.target Target power (\eqn{\phi_0}) under no period effect assumption
#' @param tsr.target Target true success rate (\eqn{\zeta_0}) under no period effect assumption
#' @param prior.sample Number of prior draws in each simulation
#' @param seed Random seed
#' @param n.simu Number of simulation replicates
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{full.results}{A tibble data frame containing the calibrated cutoffs \eqn{C_{e1}} and \eqn{c_0}, and the corresponding power and type I error rates for all tried sample sizes.}
#'   \item{optimal.config}{A tibble data frame with the same metrics for the first sample size meeting the target criteria.}
#' }
#'
#' @examples
#'
#' \donttest{
#' COCA.calibration(
#'   case = 1, n.stage1 = 24, n.stage2 = seq(20, 26, 2),
#'   dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
#'   dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
#'   eff.null = 0.25, eff.alt.SOC = 0.25, eff.alt.A = 0.35,
#'   eff.alt.B = 0.35, eff.alt.AB = 0.55, period.effect = c(0.1, 0.2, 0.3),
#'   alpha.level = 0.10, alpha.max = 0.20, fsr.level = 0.05,
#'   power.target = 0.90, tsr.target = 0.80,
#'   prior.sample = 1e4, seed = 123, n.simu = 1000
#' )
#' }
#'

#' @import stats
#' @import pbapply
#' @import cli
#' @import tibble
#' @import coda

#' @export
#'
COCA.calibration <- function(
    case, n.stage1 = 24, n.stage2 = c(),
    dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 0, dosage.singleB = 0,
    dosage.comb = list(A = c(), B = c()), eff.null, eff.alt.SOC, eff.alt.A = 0,
    eff.alt.B = 0, eff.alt.AB, period.effect = c(0.1, 0.2, 0.3),
    alpha.level = 0.10, alpha.max = 0.20, fsr.level = 0.05, power.target = 0.90, tsr.target = 0.80,
    prior.sample = 1e6, seed = 123, n.simu = 1e4){

  for (ii in seq_along(n.stage2)) {
    cli_alert(paste0("n.stage2 = ", n.stage2[ii], ": in process"))
    output.tmp <- one.calibration(
      case = case, n.stage1 = n.stage1, n.stage2 = n.stage2[ii],
      dosage.ctrl = dosage.ctrl, dosage.singleA = dosage.singleA, dosage.singleB = dosage.singleB,
      dosage.comb = dosage.comb, eff.null = eff.null, eff.alt.SOC = eff.alt.SOC, eff.alt.A = eff.alt.A,
      eff.alt.B = eff.alt.B, eff.alt.AB = eff.alt.AB, period.effect = period.effect,
      alpha.level = alpha.level, alpha.max = alpha.max, fsr.level = fsr.level, tsr.target = tsr.target,
      prior.sample = prior.sample, seed = seed, n.simu = n.simu
    )

    if(ii == 1){
      output <- output.tmp
    } else {
      output <- rbind(output, output.tmp)
    }

    if(output$Power[ii] >= power.target & output$c0[ii] > 0){
      optimal <- output.tmp
      break
    } else if (ii == length(n.stage2)) {
      optimal <- c()
      warning(paste0("Target power or true success rate not reached. Please consider increasing the stage 2 sample size by exploring the range (", max(n.stage2), ",  ", (max(n.stage2) + 10), ")."))
    }
  }

  return(list(full.results = output, optimal.config = optimal))

}

#' @keywords internal
one.calibration <- function(
    case, n.stage1, n.stage2,
    dosage.ctrl, dosage.singleA, dosage.singleB, dosage.comb,
    eff.null, eff.alt.SOC, eff.alt.A, eff.alt.B, eff.alt.AB, period.effect,
    alpha.level, alpha.max, fsr.level, tsr.target, prior.sample, seed, n.simu) {
  # Check input
  if (!case %in% c(1, 2, 3)) stop("'case' must be one of: 1, 2, or 3.")
  if (case == 1 & (dosage.singleA == 0 | dosage.singleB == 0)) {
    stop("In case 1, both single arms must have dosages greater than zero.")
  }
  if (case == 2 & all(dosage.ctrl == 0) & (dosage.singleA * dosage.singleB != 0)) {
    stop("In case 2, please set the dosage of the unavailable single arm (A or B) to zero.")
  }
  if (case == 2 & all(dosage.ctrl == 0) & (dosage.singleA == 0 & dosage.singleB == 0)) {
    stop("In case 2, at least one of the single arms must have a dosage greater than zero.")
  }
  if (case == 2 & any(dosage.ctrl != 0)) {
    if ((dosage.ctrl["A"] != 0 & dosage.singleA != 0) | (dosage.ctrl["B"] != 0 & dosage.singleB != 0)) {
      stop("In case 2, where SOC is one of the single agents, please set the dosage for that single arm to zero. ")
    }
  }

  if (!is.numeric(n.stage1) || length(n.stage1) != 1 || n.stage1 <= 0 || n.stage1 != as.integer(n.stage1)) {
    stop("'n.stage1' must be a positive integer.")
  }
  if (!is.numeric(n.stage2) || length(n.stage2) != 1 || n.stage2 <= 0 || n.stage2 != as.integer(n.stage2)) {
    stop("'n.stage2' must be a positive integer.")
  }

  if (dosage.ctrl["A"] * dosage.ctrl["B"] != 0) {
    stop("dosage.ctrl: Both dosages for A and B in the control arm cannot be non-zero; at least one of them must be zero.")
  }
  if (any(c(dosage.ctrl, dosage.singleA, dosage.singleB, do.call(c, dosage.comb)) < 0)) {
    stop("All dosages must be greater than 0.")
  }

  for (param_name in c("eff.null", "eff.alt.SOC", "eff.alt.A", "eff.alt.B", "eff.alt.AB")) {
    param_value <- get(param_name)
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(sprintf("'%s' must be between 0 and 1.", param_name))
    }
  }

  if (case == 1 & (eff.alt.SOC > eff.alt.A || eff.alt.SOC > eff.alt.B)) {
    stop(sprintf("'eff.alt.SOC' should not exceed 'eff.alt.A' or 'eff.alt.B'."))
  }
  if (case == 2 & (eff.alt.SOC > max(eff.alt.A, eff.alt.B))) {
    stop(sprintf("'eff.alt.SOC' should not exceed 'eff.alt.A' or 'eff.alt.B'."))
  }
  if (eff.alt.A > eff.alt.AB || eff.alt.B > eff.alt.AB) {
    stop(sprintf("'eff.alt.A' or 'eff.alt.B' should not exceed 'eff.alt.AB'."))
  }

  if (!is.numeric(period.effect) || any(period.effect < 0)) stop("'period.effect' must be a numeric vector with non-negative values.")
  if (any(eff.null + period.effect > 1)) {
    stop(sprintf("Each 'eff.null + period.effect' must be <= 1. Found max = %s", max(eff.null + period.effect)))
  }

  for (param_name in c("alpha.level", "alpha.max", "fsr.level", "tsr.target")) {
    param_value <- get(param_name)
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(sprintf("'%s' must be between 0 and 1.", param_name))
    }
  }
  if (alpha.max < alpha.level) stop(sprintf("'alpha.max' must be greater than or equal to 'alpha.level'."))

  # Main function
  summary_tab <- tibble(
    case = case, n.stage2 = n.stage2,
    Ce1 = NA, c0 = NA, Power = NA, TypeI = NA
  )
  n.comb.dose <- length(dosage.comb[["A"]])
  BCI_null <- run.whole(
    fda.case = case, n.comb.dose = n.comb.dose, n.stage1 = n.stage1, n.stage2 = n.stage2,
    dosage.ctrl = dosage.ctrl, dosage.singleA = dosage.singleA,
    dosage.singleB = dosage.singleB, dosage.comb = dosage.comb,
    eff.ctrl = eff.null, eff.A = eff.null, eff.B = eff.null, eff.AB = eff.null,
    prior.sample = prior.sample, batch.idx = seed, batch.sn = n.simu
  )

  BCI_alt <- run.whole(
    fda.case = case, n.comb.dose = n.comb.dose, n.stage1 = n.stage1, n.stage2 = n.stage2,
    dosage.ctrl = dosage.ctrl, dosage.singleA = dosage.singleA,
    dosage.singleB = dosage.singleB, dosage.comb = dosage.comb,
    eff.ctrl = eff.alt.SOC, eff.A = eff.alt.A, eff.B = eff.alt.B, eff.AB = eff.alt.AB,
    prior.sample = prior.sample, batch.idx = seed, batch.sn = n.simu
  )

  Ce1_0 <- round(quantile(BCI_null[1, ], (1 - alpha.level)), digits = 4)
  BCI_period <- list()
  Ce1_p <- c()
  for (pp in seq_along(period.effect)) {
    BCI_period[[pp]] <- run.whole(
      fda.case = case, n.comb.dose = n.comb.dose, n.stage1 = n.stage1, n.stage2 = n.stage2,
      dosage.ctrl = dosage.ctrl, dosage.singleA = dosage.singleA,
      dosage.singleB = dosage.singleB, dosage.comb = dosage.comb,
      eff.ctrl = eff.null, eff.A = eff.null, eff.B = eff.null, eff.AB = eff.null,
      period_eff = period.effect[pp],
      prior.sample = prior.sample, batch.idx = seed, batch.sn = n.simu
    )
    Ce1_p[pp] <- round(quantile(BCI_period[[pp]][1, ], (1 - alpha.max)), digits = 4)
  }

  Ce1 <- max(Ce1_0, max(Ce1_p))
  Ce.k.lower <- .find_klower(
    fda.case = case,
    k.min = 0.1, Ce.1 = Ce1, BCI = BCI_null, level = fsr.level
  ) %>% suppressWarnings()
  if (Ce.k.lower == -1) {
    Ce.k <- -1
  } else {
    Ce.k <- .find_kupper(
      fda.case = case,
      k.min = Ce.k.lower, Ce.1 = Ce1, BCI = BCI_alt, level = tsr.target
    ) %>% suppressWarnings()
  }

  summary_tab$Ce1 <- Ce1
  summary_tab$c0 <- Ce.k
  summary_tab$Power <- round((length(which(BCI_alt[1, ] > Ce1)) / dim(BCI_alt)[2]), 4)
  summary_tab$TypeI <- round((length(which(BCI_null[1, ] > Ce1)) / dim(BCI_null)[2]), 4)
  return(summary_tab)
}


#' @keywords internal
.dose_standardize <- function(d) {
  return(d / max(d))
}

#' @keywords internal
.find_order_stat <- function(q) {
  which.max(q)
}

#' @keywords internal
.get_Beta_prior <- function(n.sample = 1e6, control = TRUE, type = 2, seed = 0) {
  set.seed(seed)
  # type = 1 for vague prior on beta3; type = 2 for spike and slab prior
  if (control) {
    beta0 <- rnorm(n.sample, 0, sqrt(10))
  } else {
    beta0 <- rnorm(n.sample, 0, 0.000001)
  }
  beta1 <- rnorm(n.sample, 0, sqrt(10))
  beta2 <- rnorm(n.sample, 0, sqrt(10))
  if (type == 1) {
    beta3 <- rnorm(n.sample, 0, sqrt(10))
  } else if (type == 2) {
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

#' @keywords internal
.get_post <- function(Beta_prior, X.mtx, pE_prior,
                      logpE_prior0, logpE_prior1, data_Y, data_N) {
  log_likelihood <- (data_Y %*% logpE_prior0) + ((data_N - data_Y) %*% logpE_prior1)
  likelihood <- exp(log_likelihood)
  weights <- likelihood / sum(likelihood)
  idx <- sample(1:ncol(Beta_prior), size = ncol(Beta_prior), replace = TRUE, prob = weights)
  pE_post <- expit(MtxProd(X.mtx, Beta_prior[, idx]))
  return(pE_post)
}



#' @keywords internal
.find_klower <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.05) {
  k0 <- seq(k.min, 1, 0.01)
  if (fda.case == 1) {
    typei.eff <- sapply(seq_along(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]) & BCI[3, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(typei.eff <= level)])
  } else if (fda.case == 2) {
    typei.eff <- sapply(seq_along(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(typei.eff <= level)])
  } else if (fda.case == 3) {
    k <- 1
  }
  if (is.infinite(k)) {
    k <- -1
  }
  return(k)
}

#' @keywords internal
.find_kupper <- function(fda.case = 1, k.min, BCI, Ce.1, level = 0.80) {
  k0 <- seq(k.min, 1, 0.01)
  if (fda.case == 1) {
    pw.eff <- sapply(seq_along(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]) & BCI[3, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(pw.eff >= level)])
  } else if (fda.case == 2) {
    pw.eff <- sapply(seq_along(k0), function(r) {
      length(which(BCI[1, ] > Ce.1 & BCI[2, ] > (Ce.1 * k0[r]))) / (dim(BCI)[2])
    })
    k <- min(k0[which(pw.eff >= level)])
  } else if (fda.case == 3) {
    k <- 1
  }
  if (is.infinite(k)) {
    k <- -1
  }
  return(k)
}

#' @keywords internal
run.whole <- function(fda.case, n.comb.dose, n.stage1, n.stage2,
                      dosage.ctrl, dosage.singleA = 0, dosage.singleB = 0, dosage.comb,
                      eff.ctrl, eff.A, eff.B, eff.AB, period_eff = 0,
                      prior.sample, batch.idx, batch.sn = 100) {
  set.seed(1233 + 10 * batch.idx + 100 * period_eff)

  sn_s1 <- batch.sn
  ndose <- n.comb.dose
  n_21 <- n.stage1
  n_22 <- n.stage2
  narm_22 <- switch(fda.case,
    4,
    3,
    2
  )

  if (fda.case == 1) {
    trial.arm <- 1:4
  } else if (fda.case == 2) {
    trial.arm <- c(1, ifelse(dosage.singleB == 0, 2, 3), 4)
  } else if (fda.case == 3) {
    trial.arm <- c(1, 4)
  }
  period <- c(rep(0, 4)[trial.arm], rep(1, ndose))
  prior.control <- ifelse(all(dosage.ctrl == 0), TRUE, FALSE)

  dosage.comb <- do.call(rbind, dosage.comb)
  dose_std <- t(apply(
    cbind(dosage.ctrl, c(dosage.singleA, dosage.singleB), dosage.comb),
    MARGIN = 1, .dose_standardize
  ))
  row.names(dose_std) <- c("A", "B")
  dose_std_ctrl <- dose_std[, 1]
  dose_std_single <- dose_std[, 2]
  dose_std_comb <- dose_std[, -(1:2)]

  Econtrol <- eff.ctrl
  singleA <- eff.A
  singleB <- eff.B
  comb <- eff.AB
  comb_s1 <- comb + period_eff

  #### Stage I
  Ye_21 <- sapply(1:sn_s1, function(r) rbinom(rep(1, ndose), n_21, prob = comb_s1))
  q_hat <- (Ye_21 + 0.1) / (n_21 + 0.2)
  j_ast <- sapply(1:sn_s1, function(r) .find_order_stat(q = q_hat[, r]))

  X1_all <- sapply(1:ndose, function(r) {
    c(dose_std_ctrl["A"], dose_std_single["A"], 0, dose_std_comb["A", r], dose_std_comb["A", ])
  })
  X2_all <- sapply(1:ndose, function(r) {
    c(dose_std_ctrl["B"], 0, dose_std_single["B"], dose_std_comb["B", r], dose_std_comb["B", ])
  })
  colnames(X1_all) <- colnames(X2_all) <- paste0("j=", 1:ndose)
  X1 <- X1_all[c(trial.arm, 5:(ndose + 4)), ]
  X2 <- X2_all[c(trial.arm, 5:(ndose + 4)), ]

  eff_22_all <- rbind(rep(Econtrol, sn_s1), rep(singleA, sn_s1), rep(singleB, sn_s1), rep(comb, sn_s1))
  eff_22 <- eff_22_all[trial.arm, ]
  BCI_c_batch <- matrix(NA, nrow = (narm_22 - 1), ncol = sn_s1) # BCI for stage II
  Ye_22 <- sapply(1:sn_s1, function(r) rbinom(rep(1, narm_22), n_22, prob = eff_22[, r]))

  X.mtx.all <- lapply(1:ndose, function(r) {
    cbind(1, X1[, r], X2[, r], (X1[, r] * X2[, r]), (X1[, r] * X2[, r] * period))
  })
  Beta_prior <- .get_Beta_prior(n.sample = prior.sample, control = prior.control, type = 2)
  pE_prior_list <- lapply(1:ndose, function(r) {
    expit(MtxProd(X.mtx.all[[r]], Beta_prior))
  })
  logpE_prior0_list <- lapply(1:ndose, function(r) {
    log(pE_prior_list[[r]])
  })
  logpE_prior1_list <- lapply(1:ndose, function(r) {
    log(1 - pE_prior_list[[r]])
  })
  data_Y <- t(rbind(Ye_22, Ye_21))
  data_N <- c(rep(n_22, narm_22), rep(n_21, ndose))

  BCI_c_batch <- pbsapply(1:sn_s1, function(i) {
    pE_prior <- pE_prior_list[[j_ast[i]]]
    logpE_prior0 <- logpE_prior0_list[[j_ast[i]]]
    logpE_prior1 <- logpE_prior1_list[[j_ast[i]]]
    X.mtx <- X.mtx.all[[j_ast[i]]]
    pE_post <- .get_post(
      Beta_prior = Beta_prior, X.mtx = X.mtx, pE_prior = pE_prior,
      logpE_prior0 = logpE_prior0, logpE_prior1 = logpE_prior1,
      data_Y = data_Y[i, ], data_N = data_N
    )
    output <- sapply(1:(narm_22 - 1), function(r) mean(pE_post[narm_22, ] > pE_post[r, ]))
    return(output)
  })
  return(BCI_c_batch)
}
