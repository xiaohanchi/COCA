
<!-- README.md is generated from README.Rmd. Please edit that file -->

# COCA: Combination Dose Optimization and Component Contribution Assessment

<!-- badges: start -->
<!-- badges: end -->

This work introduces **COCA**, a two-stage randomized phase II design
that seamlessly integrates combination dose optimization with component
contribution assessment. In stage 1, the optimal combination dose is
determined by maximizing the risk–benefit tradeoff across multiple
candidate combination doses. In stage 2, a multi-arm randomized phase is
initiated to evaluate the contribution of each component within the
combination therapy.

<img src="./man/figures/flowchart_v0.png" width="80%"/>

## Installation

JAGS is required for implementing Bayesian MCMC sampling. Please
download JAGS from <https://sourceforge.net/projects/mcmc-jags/>.

To install the COCA package from GitHub, run:

``` r
install.packages("devtools")
devtools::install_github("xiaohanchi/COCA")
```

## Usage

Load COCA package:

``` r
library(COCA)
```

For a specific stage 2 sample size (e.g., 20), get the calibrated design
cutoffs and power:

``` r
COCA.calibration(
  case = 1, n.stage1 = 24, n.stage2 = 20, eff.null = 0.25,
  eff.alt.SOC = 0.25, eff.alt.A = 0.35, eff.alt.B = 0.35, eff.alt.AB = 0.55,
  period.effect = c(0.1, 0.2, 0.3),
  alpha.level = 0.10, alpha.max = 0.20, fsr.level = 0.05, tsr.level = 0.80, 
  n.simu = 100
)
```

- `n.simu` is set to 100 for illustration. For more accurate
  calibration, consider using a larger value, such as 10000, though this
  may require additional computation time.
- If the power does not reach the target, increase `n.stage2` and repeat
  the process.

Once the optimal `n.stage2` is found, run simulations to get the
operating characteristics of the COCA design with the calibrated
configurations:

``` r
# E.g., scenario 1 (period effect = 0)
COCA.getOC(
  n.stage2 = 26, Ce = 0.8983, c0 = 0.7,
  tox.SOC = 0.10, eff.SOC = 0.25,
  tox.A = 0.25, tox.B = 0.15,
  eff.A = 0.25, eff.B = 0.25,
  tox.AB = c(0.30, 0.30, 0.15), eff.AB = c(0.25, 0.25, 0.25),
  period.effect = 0, n.simu = 100
)
```

- Again, `n.simu` is set to 100 for illustration. For more accurate
  simulation, consider using a larger value, such as 10000.

## Authors and Reference

Chi, X., Lin, R.<sup>\*</sup>, Yuan, Y.<sup>\*</sup> (2025+). COCA: A
Randomized Bayesian Design Integrating Dose Optimization and Component
Contribution Assessment for Combination Therapies. Under Minor Revision
in *Biometrics*.
