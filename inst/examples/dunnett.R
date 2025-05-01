rm(list = ls())

dunnett.getOC <- function(p, n, n.simu){
  set.seed(1)
  d1 <- 1.62
  d2 <- 1.04
  count1 <- count2 <- count3 <- 0
  n.simu <- 200000
  for (i in 1:n.simu) {
    y <- rbinom(4, n, p)
    phat1 <- sum(y[-2]) / (n * 3)
    phat2 <- sum(y[-1]) / (n * 3)
    s1 <- sqrt(phat1 * (1 - phat1))
    s2 <- sqrt(phat2 * (1 - phat2))
    t11 <- (y[3] / n - y[1] / n) - d1 * s1 * sqrt(1 / n + 1 / n)
    t12 <- (y[4] / n - y[1] / n) - d1 * s1 * sqrt(1 / n + 1 / n)
    if (y[4] > y[3]) {
      count3 <- count3 + 1
    }
    t21 <- (y[3] / n - y[2] / n) - d2 * s2 * sqrt(1 / n + 1 / n)
    t22 <- (y[4] / n - y[2] / n) - d2 * s2 * sqrt(1 / n + 1 / n)
    if (t11 > 0 | t12 > 0) {
      count1 <- count1 + 1
      if (t21 > 0 | t22 > 0) {
        count2 <- count2 + 1
      }
    }
  }

  power <- round(count1 / n.simu * 100, 2)
  SR <- round(count2 / n.simu * 100, 2)

  output <- data.frame(power, SR)
  names(output) <- c("Power (%)", "SR (%)")
  return(output)
}


# To get the operating characteristics of  Dunnett's test under the null hypothesis, run:
dunnett.getOC(p = c(0.07,0.07,0.07,0.07), n = 46)
#   Power (%) SR (%)
#       9.94   4.99

# Operating characteristics under the alternative hypothesis:
dunnett.getOC(p = c(0.07, 0.12, 0.25, 0.25), n = 46)
#   Power (%) SR (%)
#       91.05  80.41


# Operating characteristics under observed rates:
dunnett.getOC(p = c(0.07, 0.11, 0.10, 0.24), n = 46)
#   Power (%) SR (%)
#       79.46  66.33
