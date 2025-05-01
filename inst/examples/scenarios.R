Nscenario <- 21
ndose <- 3
mtd.dose <- 1

Tox_prob <- Eff_prob <- matrix(NA, nrow = ndose, ncol = Nscenario)
Tox_prob_A <- Eff_prob_A <- matrix(NA, nrow = ndose, ncol = Nscenario)
Tox_prob_B <- Eff_prob_B <- matrix(NA, nrow = ndose, ncol = Nscenario)

Tox_prob_A <- matrix(0.25, nrow = ndose, ncol = Nscenario)
Tox_prob_B <- matrix(0.15, nrow = ndose, ncol = Nscenario)

Tcontrol <- 0.10 # tox prob for control arm
Econtrol <- 0.25
# scenario 1
Tox_prob[, 1] <- c(0.30, 0.30, 0.15)
Eff_prob[, 1] <- c(0.55, 0.55, 0.55)

Eff_prob_A[, 1] <- rep(0.35, ndose)
Eff_prob_B[, 1] <- rep(0.35, ndose)

# scenario 2
Tox_prob[, 2] <- c(0.30, 0.30, 0.30)
Eff_prob[, 2] <- c(0.55, 0.55, 0.35)

Eff_prob_A[, 2] <- rep(0.40, ndose)
Eff_prob_B[, 2] <- rep(0.40, ndose)

# scenario 3
Tox_prob[, 3] <- c(0.30, 0.30, 0.30)
Eff_prob[, 3] <- c(0.55, 0.45, 0.35)

Eff_prob_A[, 3] <- rep(0.35, ndose)
Eff_prob_B[, 3] <- rep(0.35, ndose)

# scenario 4
Tox_prob[, 4] <- c(0.30, 0.30, 0.30)
Eff_prob[, 4] <- c(0.55, 0.35, 0.35)

Eff_prob_A[, 4] <- rep(0.35, ndose)
Eff_prob_B[, 4] <- rep(0.35, ndose)

# scenario 5
Tox_prob[, 5] <- c(0.30, 0.25, 0.15)
Eff_prob[, 5] <- c(0.55, 0.55, 0.55)

Eff_prob_A[, 5] <- rep(0.40, ndose)
Eff_prob_B[, 5] <- rep(0.40, ndose)

# scenario 6
Tox_prob[, 6] <- c(0.30, 0.25, 0.15)
Eff_prob[, 6] <- c(0.55, 0.55, 0.50)

Eff_prob_A[, 6] <- rep(0.35, ndose)
Eff_prob_B[, 6] <- rep(0.30, ndose)

# scenario 7
Tox_prob[, 7] <- c(0.30, 0.25, 0.15)
Eff_prob[, 7] <- c(0.55, 0.50, 0.45)

Eff_prob_A[, 7] <- rep(0.35, ndose)
Eff_prob_B[, 7] <- rep(0.30, ndose)

# scenario 8
Tox_prob[, 8] <- c(0.30, 0.25, 0.15)
Eff_prob[, 8] <- c(0.55, 0.45, 0.45)

Eff_prob_A[, 8] <- rep(0.25, ndose)
Eff_prob_B[, 8] <- rep(0.25, ndose)

# scenario 9
Tox_prob[, 9] <- c(0.30, 0.25, 0.20)
Eff_prob[, 9] <- c(0.55, 0.55, 0.55)

Eff_prob_A[, 9] <- rep(0.35, ndose)
Eff_prob_B[, 9] <- rep(0.25, ndose)

# scenario 10
Tox_prob[, 10] <- c(0.30, 0.25, 0.20)
Eff_prob[, 10] <- c(0.55, 0.55, 0.35)

Eff_prob_A[, 10] <- rep(0.40, ndose)
Eff_prob_B[, 10] <- rep(0.25, ndose)

# scenario 11
Tox_prob[, 11] <- c(0.30, 0.25, 0.20)
Eff_prob[, 11] <- c(0.55, 0.45, 0.35)

Eff_prob_A[, 11] <- rep(0.30, ndose)
Eff_prob_B[, 11] <- rep(0.25, ndose)

# scenario 12
Tox_prob[, 12] <- c(0.30, 0.25, 0.20)
Eff_prob[, 12] <- c(0.55, 0.30, 0.30)

Eff_prob_A[, 12] <- rep(0.35, ndose)
Eff_prob_B[, 12] <- rep(0.35, ndose)

# scenario 13
Tox_prob[, 13] <- c(0.30, 0.20, 0.15)
Eff_prob[, 13] <- c(0.55, 0.55, 0.55)

Eff_prob_A[, 13] <- rep(0.40, ndose)
Eff_prob_B[, 13] <- rep(0.25, ndose)


# scenario 14
Tox_prob[, 14] <- c(0.30, 0.20, 0.15)
Eff_prob[, 14] <- c(0.55, 0.55, 0.35)

Eff_prob_A[, 14] <- rep(0.45, ndose)
Eff_prob_B[, 14] <- rep(0.25, ndose)

# scenario 15
Tox_prob[, 15] <- c(0.30, 0.20, 0.15)
Eff_prob[, 15] <- c(0.55, 0.45, 0.45)

Eff_prob_A[, 15] <- rep(0.35, ndose)
Eff_prob_B[, 15] <- rep(0.35, ndose)

# scenario 16
Tox_prob[, 16] <- c(0.45, 0.30, 0.30)
Eff_prob[, 16] <- c(0.70, 0.45, 0.35)

Eff_prob_A[, 16] <- rep(0.35, ndose)
Eff_prob_B[, 16] <- rep(0.25, ndose)

# scenario 17
Tox_prob[, 17] <- c(0.45, 0.30, 0.30)
Eff_prob[, 17] <- c(0.60, 0.55, 0.35)

Eff_prob_A[, 17] <- rep(0.30, ndose)
Eff_prob_B[, 17] <- rep(0.20, ndose)

# scenario 18
Tox_prob[, 18] <- c(0.30, 0.25, 0.10)
Eff_prob[, 18] <- c(0.30, 0.30, 0.10)

Eff_prob_A[, 18] <- rep(0.20, ndose)
Eff_prob_B[, 18] <- rep(0.15, ndose)

# scenario 19
Tox_prob[, 19] <- c(0.30, 0.25, 0.10)
Eff_prob[, 19] <- c(0.55, 0.40, 0.10)

Eff_prob_A[, 19] <- rep(0.40, ndose)
Eff_prob_B[, 19] <- rep(0.35, ndose)

# scenario 20
Tox_prob[, 20] <- c(0.45, 0.40, 0.40)
Eff_prob[, 20] <- c(0.15, 0.10, 0.10)

Eff_prob_A[, 20] <- rep(0.10, ndose)
Eff_prob_B[, 20] <- rep(0.10, ndose)

# scenario 21
Tox_prob[, 21] <- c(0.30, 0.30, 0.15)
Eff_prob[, 21] <- c(0.25, 0.25, 0.25)

Eff_prob_A[, 21] <- rep(0.25, ndose)
Eff_prob_B[, 21] <- rep(0.25, ndose)
