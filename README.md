
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

<img src="./man/figures/flowchart_v0.png" width="80%">

## Installation

JAGS is required for implementing Bayesian MCMC sampling. Please
download JAGS from <https://sourceforge.net/projects/mcmc-jags/>.

To install the COCA package from GitHub, run:

``` r
install.packages("devtools")
devtools::install_github("xiaohanchi/COCA")
```

## Usage

Get the calibrated design cutoffs and stage 1 sample size:

``` r
# library(COCA)
## basic example code
```

Run simulation to get the operating characteristics of the COCA design
with the calibrated configurations:

``` r
# library(COCA)
## basic example code
```

## Authors and Reference

Chi, X., Lin, R.<sup>\*</sup>, Yuan, Y.<sup>\*</sup> (2025+). COCA: A
Randomized Bayesian Design Integrating Dose Optimization and Component
Contribution Assessment for Combination Therapies. Under Minor Revision
in *Biometrics*.
