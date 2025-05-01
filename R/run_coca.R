#' Run the COCA Design
#'
#' Run the COCA design and get its operating characteristics.
#'
#' @param case Trial type for stage 2. \code{case = 1} for 4-arm trial comparing AB vs. A vs. B vs. SOC; \code{case = 2} for 3-arm trial comparing AB vs. A (or B) vs. SOC; \code{case = 3} for 2-arm trial comparing AB vs. SOC.
#' @param n.stage1 Sample size for stage 1
#' @param n.stage2 Sample size for stage 2
#' @param Ce,c0 Design cutoffs, obtained using the \code{COCA.calibration} function.
#' @param nlook.stage1 The total number of planned analyses in stage 1, including all interim analyses and the final analysis (e.g., \code{nlook.stage1 = 2} for one interim and one final analysis in stage 2).
#' @param nlook.stage2 The total number of planned analyses in stage 2, including all interim analyses and the final analysis (e.g., \code{nlook.stage2 = 2} for one interim and one final analysis in stage 2).
#' @param dosage.ctrl Dosage level of the control arm for stage 2. For an SOC control, use \code{c(A = 0, B = 0)}.If one of the single agents is used as the SOC (e.g., drug A 300 mg), use \code{c(A = 300, B = 0)}.
#' @param dosage.singleA Dosage level of drug A in the single arm for stage 2.
#' @param dosage.singleB Dosage level of drug B in the single arm for stage 2.
#' @param dosage.comb A named list specifying the dosage levels of drugs A and B across combination arms in stage 1.
#'   For example, \code{list(A = c(300, 300, 200), B = c(300, 200, 300))} defines three dose combinations:
#'   dose1 = (A = 300, B = 300), dose2 = (A = 300, B = 200), and dose3 = (A = 200, B = 300).
#'   All dosage values can be on any scale but must use the same scale across \code{dosage.singleA}, \code{dosage.singleB}, and \code{dosage.comb}.
#' @param tox.SOC,tox.A,tox.B True toxicity probabilities for SOC, arm A, and arm B.
#' @param tox.AB A vector of true toxicity probabilities for all combination doses being tested in the trial.
#' @param eff.SOC,eff.A,eff.B True efficacy probabilities for SOC, arm A, and arm B.
#' @param eff.AB.s1 A vector of true efficacy probabilities for all combination doses in stage 1.
#' @param eff.AB.s2 A vector of true efficacy probabilities for all combination doses in stage 2.
#' @param tox.isomat Matrix with 2 columns that contains isotonicity conditions for the toxicity order among the J combination doses. The format should follow the structure expected by the `activeSet()` function in the `isotone` package. For details, refer to the [isotone package documentation](https://cran.r-project.org/web/packages/isotone/index.html).
#' @param tox.upper Highest acceptable toxicity rate (\eqn{\phi_{T}})
#' @param eff.lower Lowest acceptable efficacy rate (\eqn{\phi_{E}})
#' @param Cs,Ct Probability cutoffs in stage 1
#' @param C.f1,C.f2 Probability cutoffs in stage 2
#' @param utility.score Vector of utility score: \code{c(b1, b2, b3, b4)} represents the utility for (toxicity, no efficacy), (no toxicity, no efficacy), (toxicity, efficacy), and (no toxicity, efficacy), respectively.
#' @param rho Correlation between toxicity and efficacy
#' @param period.effect Period effect
#' @param n.simu Number of simulation replicates. The default value \code{n.simu = 10} is used for illustration purposes and is small to reduce computation time. For more accurate results, consider using a larger value, such as 5000.
#' @param prior.sample Number of prior draws in each simulation
#' @param seed Random seed
#' @description This function uses `activeSet()` from the `isotone` package to perform isotonic regression.
#'
#' @return Returns the operating characteristics of stage 1 (selection and expected sample size) and stage 2 (power, GP, SR, OSR, and expected sample size).
#'
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
#'
#' # Scenario 1 (period effect = 0)
#' \donttest{
#' COCA.getOC(
#'   case = 1, n.stage1 = 24, n.stage2 = 26, Ce = 0.9152, c0 = 0.66,
#'   dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
#'   dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
#'   tox.SOC = 0.10, eff.SOC = 0.25, tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.25, eff.B = 0.25, tox.AB = c(0.30, 0.30, 0.15),
#'   eff.AB.s1 = c(0.25, 0.25, 0.25), eff.AB.s2 = c(0.25, 0.25, 0.25),
#'   tox.isomat = matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2),
#'   tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
#'   utility.score = c(0, 60, 40, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
#' )
#' }
#'
#' # Scenario 1 (period effect = 0.2)
#' \donttest{
#' COCA.getOC(
#'   case = 1, n.stage1 = 24, n.stage2 = 26, Ce = 0.9152, c0 = 0.66,
#'   dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
#'   dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
#'   tox.SOC = 0.10, eff.SOC = 0.25, tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.25, eff.B = 0.25, tox.AB = c(0.30, 0.30, 0.15),
#'   eff.AB.s1 = c(0.45, 0.45, 0.45), eff.AB.s2 = c(0.25, 0.25, 0.25),
#'   tox.isomat = matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2),
#'   tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
#'   utility.score = c(0, 60, 40, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
#' )
#' }
#'
#' # Scenario 2 (period effect = 0)
#' \donttest{
#' COCA.getOC(
#'   case = 1, n.stage1 = 24, n.stage2 = 26, Ce = 0.9152, c0 = 0.66,
#'   dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
#'   dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
#'   tox.SOC = 0.10, eff.SOC = 0.25, tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.35, eff.B = 0.35, tox.AB = c(0.30, 0.30, 0.15),
#'   eff.AB.s1 = c(0.55, 0.55, 0.55), eff.AB.s2 = c(0.55, 0.55, 0.55),
#'   tox.isomat = matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2),
#'   tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
#'   utility.score = c(0, 60, 40, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
#' )
#' }
#'
#' # Scenario 2 (period effect = 0.2)
#' \donttest{
#' COCA.getOC(
#'   case = 1, n.stage1 = 24, n.stage2 = 26, Ce = 0.9152, c0 = 0.66,
#'   dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
#'   dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
#'   tox.SOC = 0.10, eff.SOC = 0.25, tox.A = 0.25, tox.B = 0.15,
#'   eff.A = 0.35, eff.B = 0.35, tox.AB = c(0.30, 0.30, 0.15),
#'   eff.AB.s1 = c(0.75, 0.75, 0.75), eff.AB.s2 = c(0.55, 0.55, 0.55),
#'   tox.isomat = matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2),
#'   tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
#'   utility.score = c(0, 60, 40, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
#' )
#' }
COCA.getOC <- function(case = 1, n.stage1 = 24, n.stage2, Ce, c0, nlook.stage1 = 2, nlook.stage2 = 2,
                       dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 0, dosage.singleB = 0,
                       dosage.comb = list(A = c(), B = c()),
                       tox.SOC, eff.SOC, tox.A = 0, tox.B = 0, eff.A = 0, eff.B = 0,
                       tox.AB = c(), eff.AB.s1 = c(), eff.AB.s2 = c(), tox.isomat,
                       tox.upper, eff.lower, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
                       utility.score = c(0, 60, 40, 100), rho = 0.2,
                       prior.sample = 1e5, seed = 123, n.simu = 5000) {
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

  if (!is.numeric(nlook.stage1) || length(nlook.stage1) != 1 || nlook.stage1 <= 0 || nlook.stage1 != as.integer(nlook.stage1)) {
    stop("'nlook.stage1' must be a positive integer.")
  }
  if (!is.numeric(nlook.stage2) || length(nlook.stage2) != 1 || nlook.stage2 <= 0 || nlook.stage2 != as.integer(nlook.stage2)) {
    stop("'nlook.stage2' must be a positive integer.")
  }

  if (n.stage1 %% nlook.stage1 != 0) {
    stop("The stage 1 total sample size ('n.stage1') must be an integer multiple of the number of planned analyses ('nlook.stage1') to allow evenly spaced interim analyses, each conducted after every n.stage1 / nlook.stage1 patients.")
  }

  if (n.stage2 %% nlook.stage2 != 0) {
    stop("The stage 2 total sample size ('n.stage2') must be an integer multiple of the number of planned analyses ('nlook.stage2') to allow evenly spaced interim analyses, each conducted after every n.stage2 / nlook.stage2 patients.")
  }

  if (dosage.ctrl["A"] * dosage.ctrl["B"] != 0) {
    stop("dosage.ctrl: Both dosages for A and B in the control arm cannot be non-zero; at least one of them must be zero.")
  }
  if (any(c(dosage.ctrl, dosage.singleA, dosage.singleB, do.call(c, dosage.comb)) < 0)) {
    stop("All dosages must be greater than 0.")
  }

  for (param_name in c("tox.SOC", "eff.SOC", "tox.A", "tox.B", "eff.A", "eff.B", "tox.upper", "eff.lower", "Ce", "Cs", "Ct", "C.f1", "C.f2")) {
    param_value <- get(param_name)
    if (!is.numeric(param_value) || param_value < 0 || param_value > 1) {
      stop(sprintf("'%s' must be between 0 and 1.", param_name))
    }
  }

  if (!(length(tox.AB) == length(eff.AB.s1) && length(eff.AB.s1) == length(eff.AB.s2))) {
    stop(sprintf("'tox.AB', 'eff.AB.s1', and 'eff.AB.s2' must have the same length."))
  }

  for (param_name in c("tox.AB", "eff.AB.s1", "eff.AB.s2")) {
    param_value <- get(param_name)
    if (!is.numeric(param_value) || any(param_value < 0) || any(param_value > 1)) {
      stop(sprintf("'%s' must be between 0 and 1.", param_name))
    }
  }

  if (length(utility.score) != 4) stop(sprintf("'utility.score' must have length 4."))
  if (!is.numeric(rho) || length(rho) != 1 || rho < -1 || rho > 1) {
    stop(sprintf("'rho' must be a numeric value in [-1, 1]."))
  }


  # Main function
  Ce <- c(Ce, c0 * Ce, c0 * Ce)
  n.dose.AB <- length(tox.AB)
  Nmax_p2 <- c(n.stage1, n.stage2) # stage 1 & stage 2
  T_21 <- nlook.stage1
  T_22 <- nlook.stage2
  if (nrow(tox.isomat) == 1) {
    A1 <- rbind(tox.isomat, tox.isomat)
  } else {
    A1 <- tox.isomat
  }
  A2 <- A1[, c(2, 1)]
  C_s1 <- C_s2 <- Cs

  a <- rep(0.05, 4) # hyperparameters in dirichlet

  narm_22 <- switch(case,
    4,
    3,
    2
  )
  if (case == 1) {
    trial.arm <- 1:4
  } else if (case == 2) {
    trial.arm <- c(1, ifelse(dosage.singleB == 0, 2, 3), 4)
  } else if (case == 3) {
    trial.arm <- c(1, 4)
  }
  period <- c(rep(0, 4)[trial.arm], rep(1, n.dose.AB))
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

  ### simulation settings
  Tox_prob <- tox.AB
  Eff_prob <- eff.AB.s1

  ### Stage I
  multi_prob <- sapply(1:n.dose.AB, function(r) .solve_level(rho, Eff_prob[r], Tox_prob[r]))
  true_utility <- utility.score %*% multi_prob
  N_21 <- sapply(1:n.simu, function(r) rep(Nmax_p2[1], n.dose.AB))
  n_21 <- N_21 / T_21
  Y_21 <- pi_hat_21 <- pi_hat_21_iso <- array(0, dim = c(n.dose.AB, nrow(multi_prob), n.simu))
  piT_hat_21 <- piE_hat_21 <- currn_21 <- matrix(0, nrow = n.dose.AB, ncol = n.simu)
  proc_21 <- matrix(rep(((Tox_prob + Eff_prob) != 0) * 1, n.simu), nrow = n.dose.AB)
  BCI <- vector(mode = "numeric", n.simu)
  cli_alert("Stage 1: in process")
  for (t in 1:T_21) {
    set.seed(seed + 10 * t)
    currn_21 <- currn_21 + n_21 * proc_21 # current sample size
    temp_Y <- sapply(1:n.dose.AB, function(r) rmultinom(n.simu, n_21[1, r], multi_prob[, r]))
    temp_Y <- array(t(temp_Y), dim = dim(Y_21))
    for (i in 1:nrow(multi_prob)) {
      temp_Y[, i, ] <- temp_Y[, i, ] * proc_21
    }
    Y_21 <- Y_21 + temp_Y
    for (i in 1:nrow(multi_prob)) {
      pi_hat_21[, i, ] <- (a[i] + Y_21[, i, ]) / (sum(a) + currn_21)
    }
    piT_hat_21 <- pi_hat_21[, 1, ] + pi_hat_21[, 3, ]
    w1 <- 1 / apply(piT_hat_21, 1, var)
    piT_hat_21_iso <- sapply(1:n.simu, function(r) activeSet(A1, "LS", weights = w1, y = piT_hat_21[, r])$x)
    piE_hat_21 <- pi_hat_21[, 3, ] + pi_hat_21[, 4, ]
    rho_hat <- (pi_hat_21[, 2, ] * pi_hat_21[, 3, ] - pi_hat_21[, 1, ] * pi_hat_21[, 4, ]) / sqrt(piE_hat_21 * (1 - piE_hat_21) * piT_hat_21_iso * (1 - piT_hat_21_iso))
    for (i in 1:n.dose.AB) {
      pi_hat_21_iso[i, , ] <- sapply(1:n.simu, function(r) .solve_level(rho_hat[i, r], piE_hat_21[i, r], piT_hat_21_iso[i, r]))
    }
    prt_21 <- pbeta(tox.upper, (a[1] + a[3] + Y_21[, 1, ] + Y_21[, 3, ]), (a[2] + a[4] + currn_21 - Y_21[, 1, ] - Y_21[, 3, ]))
    w2 <- 1 / apply(prt_21, 1, var)
    prt_21_iso <- sapply(1:n.simu, function(r) activeSet(A2, "LS", weights = w2, y = prt_21[, r])$x)
    pre_21 <- pbeta(eff.lower, (a[3] + a[4] + Y_21[, 3, ] + Y_21[, 4, ]), (a[1] + a[2] + currn_21 - Y_21[, 3, ] - Y_21[, 4, ]))
    if (t < T_21) {
      proc_21[which(proc_21 != 0 & prt_21_iso < 1 - C_s1)] <- 0
      proc_21[which(proc_21 != 0 & pre_21 > C_s1)] <- 0
    } else if (t == T_21) {
      Aset <- which(prt_21_iso * proc_21 > 1 - C_s2 & pre_21 * proc_21 < C_s2)
      for (i in 1:nrow(multi_prob)) {
        pi_hat_21_iso[, i, ][-Aset] <- -100
      }
      utility <- t(sapply(1:n.dose.AB, function(r) utility.score %*% pi_hat_21_iso[r, , ]))
      j_ast1 <- apply(utility, 2, which.max)
      j_ast1[which(colSums(proc_21) == 0)] <- -1
      j_ast1[which(apply(utility, 2, max) < 0)] <- -1
    }
  }
  cli_alert_info("Stage 1: done")

  ### Stage II (Proof of Concept)
  Eff_prob <- eff.AB.s2

  Yt_21 <- (Y_21[, 1, ] + Y_21[, 3, ])[, which(j_ast1 > 0)]
  Ye_21 <- (Y_21[, 3, ] + Y_21[, 4, ])[, which(j_ast1 > 0)]
  sn_22 <- sum(j_ast1 > 0)
  N_22 <- sapply(1:sn_22, function(r) rep(Nmax_p2[2], narm_22)) # narm_22xsn matrix
  n_22 <- N_22 / T_22
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]
  tox_22_all <- rbind(
    rep(tox.SOC, sn_22), rep(tox.A, sn_22), rep(tox.B, sn_22), Tox_prob[j_ast1_tmp]
  )
  tox_22 <- tox_22_all[trial.arm, ]
  eff_22_all <- rbind(
    rep(eff.SOC, sn_22), rep(eff.A, sn_22), rep(eff.B, sn_22), Eff_prob[j_ast1_tmp]
  )
  eff_22 <- eff_22_all[trial.arm, ]

  X1_all <- sapply(1:n.dose.AB, function(r) {
    c(dose_std_ctrl["A"], dose_std_single["A"], 0, dose_std_comb["A", r], dose_std_comb["A", ])
  })
  X2_all <- sapply(1:n.dose.AB, function(r) {
    c(dose_std_ctrl["B"], 0, dose_std_single["B"], dose_std_comb["B", r], dose_std_comb["B", ])
  })
  colnames(X1_all) <- colnames(X2_all) <- paste0("j=", 1:n.dose.AB)
  X1 <- X1_all[c(trial.arm, 5:(n.dose.AB + 4)), ]
  X2 <- X2_all[c(trial.arm, 5:(n.dose.AB + 4)), ]

  Yt_pre <- sapply(1:sn_22, function(r) c(rep(0, (narm_22 - 1)), (Yt_21[j_ast1_tmp[r], r])))
  Nt_pre <- sapply(1:sn_22, function(r) c(rep(0, (narm_22 - 1)), (currn_21[, which(j_ast1 > 0)][j_ast1_tmp[r], r])))

  Yt_22 <- Ye_22 <- matrix(0, nrow = narm_22, ncol = sn_22)
  proc_22 <- matrix(1, nrow = narm_22, ncol = sn_22)
  currn_22 <- matrix(0, nrow = narm_22, ncol = sn_22)

  X.mtx.all <- lapply(1:n.dose.AB, function(r) {
    cbind(1, X1[, r], X2[, r], (X1[, r] * X2[, r]), (X1[, r] * X2[, r] * period))
  })
  cli_alert("Stage 2: in process")
  for (t in 1:T_22) {
    set.seed(seed + 10 * t)
    currn_22 <- currn_22 + n_22 * proc_22
    Yt_22 <- Yt_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = tox_22[, r])) * proc_22
    Ye_22 <- Ye_22 + sapply(1:sn_22, function(r) rbinom(rep(1, narm_22), n_22[, r], prob = eff_22[, r])) * proc_22
    data_Y <- t(rbind(Ye_22, Ye_21))
    data_N <- t(rbind(currn_22, currn_21[, which(j_ast1 > 0)]))
    fprob_all <- pbsapply(1:sn_22, function(i) {
      if (proc_22[narm_22, i] == 0) {
        fprob_1 <- fprob_2 <- rep(-1, narm_22)
        if (t == T_22) BCI <- rep(-2, (narm_22 - 1))
      } else {
        Beta_prior <- .get_Beta_prior(
          n.sample = prior.sample, control = prior.control, type = 2, seed = i
          )
        pE_prior <- expit(MtxProd(X.mtx.all[[j_ast1_tmp[i]]], Beta_prior))
        logpE_prior0 <- log(pE_prior)
        logpE_prior1 <- log(1 - pE_prior)
        X.mtx <- X.mtx.all[[j_ast1_tmp[i]]]
        pE_post <- .get_post(
          Beta_prior = Beta_prior, X.mtx = X.mtx, pE_prior = pE_prior,
          logpE_prior0 = logpE_prior0, logpE_prior1 = logpE_prior1,
          data_Y = data_Y[i, ], data_N = data_N[i, ]
        )
        # futility stopping prob
        fprob_1 <- c(100, sapply(2:narm_22, function(r) mean(pE_post[r, ] > pE_post[1, ])))
        if (case %in% c(1, 2)) { # terminate single arm A or B
          fprob_2 <- c(100, sapply(2:(narm_22 - 1), function(r) mean(pE_post[r, ] > pE_post[narm_22, ])), 100)
        } else {
          fprob_2 <- c(100, 100)
        }
        if (t == T_22) {
          BCI <- sapply(1:(narm_22 - 1), function(r) mean(pE_post[narm_22, ] > pE_post[r, ]))
        }
      }
      output <- ifelse(t == T_22, c(fprob_1, fprob_2, BCI), c(fprob_1, fprob_2))
      if(t == T_22) {
        return(c(fprob_1, fprob_2, BCI))
      } else {
        return(c(fprob_1, fprob_2))
      }
    })
    fprob_1 <- fprob_all[1:narm_22, ]
    fprob_2 <- fprob_all[(narm_22 + 1):(2*narm_22), ]
    if (t == T_22) BCI <- fprob_all[-(1:(2*narm_22)), ]

    proc_22[which(proc_22 != 0 & pbeta(tox.upper, (0.1 + Yt_pre + Yt_22), (0.1 + Nt_pre + currn_22 - Yt_pre - Yt_22)) < (1 - Ct))] <- 0
    proc_22[which(proc_22 != 0 & fprob_1 < (1 - C.f1))] <- 0
    proc_22[which(proc_22 != 0 & fprob_2 < (1 - C.f2))] <- 0
    if (t == T_22) {
      BCI[, which(proc_22[narm_22, ] == 0)] <- -2
    }
  }
  cli_alert_info("Stage 2: done")

  j_opt <- which(
    true_utility == max(
      true_utility[which(tox.AB <= tox.upper &
        eff.AB.s1 >= eff.lower)]
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
  names(stage1_output) <- c("termination (%)", paste0("dose", 1:n.dose.AB, " (%)"), "selection (%)", "EN")

  ### stage 2 results
  BCI <- BCI[, 1:sum(j_ast1 > 0)]
  if (is.null(dim(BCI)[2])) BCI <- matrix(BCI, nrow = 1)
  cont.simu <- dim(BCI)[2]
  j_ast1_tmp <- j_ast1[which(j_ast1 > 0)]

  sel.opt.g <- sapply(1:cont.simu, function(r) j_ast1_tmp[r] %in% j_opt) * 1
  sel.opt.g <- which(sel.opt.g != 0)
  for (jj in which(currn_22[nrow(currn_22), ] < n.stage2)) {
    currn_22[, jj] <- currn_22[nrow(currn_22), jj]
  }
  currn_22 <- cbind(currn_22, matrix(0, nrow = nrow(currn_22), ncol = (n.simu - cont.simu)))
  avg.n <- mean(colSums(currn_22))
  power <- length(which(BCI[1, ] > Ce[1])) / cont.simu
  GP <- length(which(BCI[1, sel.opt.g] > Ce[1])) / cont.simu

  SR <- switch(case,
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2] & BCI[3, ] > Ce[3])) / cont.simu),
    (length(which(BCI[1, ] > Ce[1] & BCI[2, ] > Ce[2])) / cont.simu),
    (length(which(BCI[1, ] > Ce[1])) / cont.simu)
  )

  OSR <- switch(case,
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2] & BCI[3, sel.opt.g] > Ce[3])) / cont.simu),
    (length(which(BCI[1, sel.opt.g] > Ce[1] & BCI[2, sel.opt.g] > Ce[2])) / cont.simu),
    (length(which(BCI[1, sel.opt.g] > Ce[1])) / cont.simu)
  )

  stage2_output <- data.frame(
    formatC(power * 100, digits = 2, format = "f"),
    formatC(GP * 100, digits = 2, format = "f"),
    formatC(SR * 100, digits = 2, format = "f"),
    formatC(OSR * 100, digits = 2, format = "f"),
    formatC(avg.n, digits = 1, format = "f")
  )
  colnames(stage2_output) <- c("Power (%)", "GP (%)", "SR (%)", "OSR (%)", "EN")

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
