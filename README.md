
<!-- README.md is generated from README.Rmd. Please edit that file -->

# COCA: Combination Dose Optimization and Component Contribution Assessment

<!-- badges: start -->
<!-- badges: end -->

### Overview

**COCA** is a two-stage randomized phase II design that seamlessly
integrates combination dose optimization with component contribution
assessment. In stage 1, the optimal combination dose is determined by
maximizing the risk–benefit tradeoff across multiple candidate
combination doses. In stage 2, a multi-arm randomized phase is initiated
to evaluate the contribution of each component within the combination
therapy.

<p align="center">
<img src="./man/figures/flowchart_v0.png" width="80%"/>
</p>

### Installation

To install the COCA package from GitHub, run:

``` r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("xiaohanchi/COCA")
```

### Usage

Load COCA package:

``` r
library(COCA)
```

To search for the optimal stage 2 sample size between 20 and 30 and
obtain the corresponding design parameters under `case = 1` and the
settings described in the ‘Numerical Studies’ section of paper \[1\],
run:

``` r
# Estimated runtime: ~5 minutes on a MacBook Air (M1, 16 GB RAM)
COCA.calibration(
  case = 1, n.stage1 = 24, n.stage2 = seq(20, 30, 2), 
  dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300, 
  dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
  eff.null = 0.25, eff.alt.SOC = 0.25, eff.alt.A = 0.35, 
  eff.alt.B = 0.35, eff.alt.AB = 0.55, period.effect = c(0.1, 0.2, 0.3), 
  alpha.level = 0.10, alpha.max = 0.20, fsr.level = 0.05, power.target = 0.90, tsr.target = 0.80,
  prior.sample = 1e5, seed = 123, n.simu = 1000
)
```

If the power or success rate does not reach the target, increase
`n.stage2` and rerun the calibration. **To fully reproduce the results
in our paper, please use the following inputs:
`period.effect = seq(-0.1, 0.5, 0.05)`, `prior.sample = 1e6`, and
`n.simu = 1e4`, while keeping all other inputs unchanged. For other
trial cases, please use `case = 2` or `case = 3` as appropriate.**

> **NOTE:** The above example is a simplified version with fewer prior
> samples and simulations to facilitate user testing. Please note that
> Monte Carlo error may be non-negligible if the number of simulations
> is not sufficiently large.  
> The full calibration procedure is time-consuming, as it requires
> identifying a suitable sample size `n.stage2` while controlling the
> maximum type I error inflation across a range of potential period
> effects. To reproduce our results, finer prior sampling and a larger
> number of simulation replications are required (e.g., increase
> `prior.sample` to $10^6$ and setting `n.simu` to $10^4$). On a MacBook
> Air (Apple M1 chip, 16 GB RAM), running 10,000 replicated studies for
> a given sample size and true probability scenario (e.g., $H_0$ or
> $H_1$) takes approximately 20 to 30 minutes. Since our calibration
> involves exploring a range of candidate sample sizes across multiple
> true scenarios, the entire process typically takes 12 to 16 hours. In
> our paper, we used a high-performance computing cluster (Seadragon,
> the MD Anderson research cluster) for the calibration. With 200
> parallel jobs, the complete process takes approximately 20 to 30
> minutes on average.

Once the optimal `n.stage2` is found, run simulations to get the
operating characteristics of the COCA design with the calibrated
configurations:

``` r
# E.g., scenario 1 (period effect = 0)
# Estimated runtime: ~13 minutes on a MacBook Air (M1, 16 GB RAM)
COCA.getOC(
  case = 1, n.stage1 = 24, n.stage2 = 26, Ce = 0.9152, c0 = 0.66, nlook.stage1 = 2, 
  nlook.stage2 = 2, dosage.ctrl = c(A = 0, B = 0), dosage.singleA = 300, dosage.singleB = 300,
  dosage.comb = list(A = c(300, 300, 200), B = c(300, 200, 300)),
  tox.SOC = 0.10, eff.SOC = 0.25, tox.A = 0.25, tox.B = 0.15,
  eff.A = 0.25, eff.B = 0.25, tox.AB = c(0.30, 0.30, 0.15),
  eff.AB.s1 = c(0.25, 0.25, 0.25), eff.AB.s2 = c(0.25, 0.25, 0.25),
  tox.isomat = matrix(c(2, 1, 3, 1), byrow = TRUE, nrow = 2),
  tox.upper = 0.35, eff.lower = 0.25, Cs = 0.85, Ct = 0.9, C.f1 = 0.9, C.f2 = 0.9,
  utility.score = c(0, 60, 40, 100), rho = 0.2, prior.sample = 1e5, seed = 1254, n.simu = 5000
)
```

`prior.sample` is set to $10^5$ for illustration. For more accurate
results, we recommend using `prior.sample = 1e6`. **To fully reproduce
the results in our paper, please use `prior.sample = 1e6`, while keeping
all other inputs unchanged. For other simulation settings, replace the
inputs as appropriate, following the guidelines in the ‘Numerical
Studies’ section of paper \[1\].**

##### \* Competing Approaches

We provide R code to implement the competing approaches (MTD-Ind,
OBD-Ind, and OBD-Pool) mentioned in the paper<sup>\[1\]</sup>. To access
it, use the following R code:

``` r
open_example("alt_designs.R")
```

Detailed instructions for reproducing the results in paper \[1\] are
provided at the end of the script, along with some example output
`.Rdata` files in `inst/example_Rdata`.

### Examples: Redesigning the NCT02519348 Trial

This example provides a step-by-step tutorial on redesigning the
NCT02519348 trial<sup>\[2\]</sup> using COCA. This trial involves four
arms: T 300 mg $\times$ 1 dose plus D 1500 mg (T300+D), T 75 mg $\times$
4 doses plus D 1500 mg (T75+D), D 1500 mg monotherapy, and T 750mg
monotherapy.

##### 1. Preparation

To get started, we need to specify the appropriate input arguments.
Since this trial lacks a control arm, we assume the T monotherapy arm as
the standard of care (SOC) control for illustration. Therefore, this
trial falls into `case = 2` (S0C(A) vs. B vs AB). We use the log scale
dosage as the dosage input for each arm:

``` r
dosage.ctrl = c(A = log(750), B = 0)
dosage.singleA = 0 # Set to zero since single agent A is considered as the control arm
dosage.singleB = log(1500) 
dosage.comb = list(A = c(log(300), log(75)), B = c(log(1500), log(1500)))
```

The null hypothesis is $H_0: q_{21}=q_{22}=q_{23}=0.07$, so
`eff.null = 0.07`. The alternative hypothesis is
$H_1: q_{21}=0.07, q_{22}=0.12, q_{23}=0.25$, so:

``` r
eff.alt.SOC = 0.07
eff.alt.B = 0.12
eff.alt.AB = 0.25
```

We also hope to control the maximum type I error under potential period
effects within the range of 0 to 0.1, so we specify
`period.effect = seq(0, 0.1, 0.02)`. For other configurations, please
refer to the ‘Redesign NCT02519348 Trial’ section of paper \[1\]. With
all these arguments in place, we can proceed to calibrate our design
parameters.

##### 2. Calibration

The `COCA.calibration` function will return the optimal stage 2 sample
size along with the corresponding power and calibrated design cutoffs
($C_{e1}$ and $c_0$). If the target power or TSR is not achieved under
the current settings, the function will issue a warning suggesting an
increase in the stage 2 sample size. We would like to assume
`n.stage1 = 24` and search for the optimal stage 2 sample size in the
range of 32 to 40:

``` r
# Estimated runtime: ~13 hours on a MacBook Air (M1, 16 GB RAM)
COCA.calibration(
    case = 2, n.stage1 = 24, n.stage2 = seq(32, 40, 2), 
    dosage.ctrl = dosage.ctrl, dosage.singleA = 0, 
    dosage.singleB = dosage.singleB,  dosage.comb = dosage.comb,
    eff.null = 0.07, eff.alt.SOC = 0.07, eff.alt.B = 0.12, eff.alt.AB = 0.25, 
    period.effect = seq(0, 0.1, 0.02), alpha.level = 0.10, alpha.max = 0.20, 
    fsr.level = 0.05, power.target = 0.90, tsr.target = 0.80, 
    prior.sample = 1e6, seed = 123, n.simu = 1e4
  )
```

    #> # A tibble: 1 × 6
    #>    case n.stage2   Ce1    c0 Power TypeI
    #>   <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl>
    #> 1     2       38 0.771  0.81 0.914   0.1

(This code may take a while to run…) For illustration, consider using a
smaller number of prior draws (e.g., `prior.sample = 1e5`) and
simulation replicates (e.g., `n.simu = 1000`) to reduce runtime. Once
completed, we obtained the optimal stage 2 sample size of 38, along with
design cutoffs $C_{e1}=0.7710$ and $c_0=0.81$, as reported in our paper.

##### 3. Run COCA Design

The `COCA.getOC` function is used to obtain the operating
characteristics of stage 1 (selection and expected sample size) and
stage 2 (power, GP, SR, OSR, and expected sample size) of COCA. So far,
we have obtained the optimal design parameters: `n.stage2 = 38`,
`Ce = 0.7710`, and `c0 = 0.81`. To assess the design performance, we
still need to specify the true efficacy and toxicity rates. Using the
observed trial outcomes: the toxicity rates in the four arms (T vs. D
vs. T300+D vs. T75+D) were 24.6%, 10.9%, 17.6%, and 14.6%, and efficacy
rates were 7.2%, 10.6%, 24.0%, and 9.5%, respectively. Therefore:

``` r
tox.SOC = 0.246
tox.B = 0.109
tox.AB = c(0.176, 0.146)

eff.SOC = 0.072
eff.B = 0.106
eff.AB.s1 = c(0.240, 0.095) # assume no period effect
eff.AB.s2 = c(0.240, 0.095)
```

We assume the toxicity ordering between the two combination doses is
T300+D $\geq$ T75+D, so

``` r
tox.isomat = matrix(c(2, 1), byrow = T, nrow = 1)
```

We use a utility score of `utility.score = c(0, 40, 60, 100)` in stage 1
dose optimization. For additional configurations, please refer to the
‘Redesign NCT02519348 Trial’ section of paper \[1\].

First, let’s assume no period effect between stages 1 and 2 and obtain
the design operating characteristics:

``` r
# Estimated runtime: ~10 minutes on a MacBook Air (M1, 16 GB RAM)
COCA.getOC(
  case = 2, n.stage1 = 24, n.stage2 = 38, Ce = 0.7710, c0 = 0.81, nlook.stage1 = 2, 
  nlook.stage2 = 2, dosage.ctrl = dosage.ctrl, dosage.singleA = 0, 
  dosage.singleB = dosage.singleB,  dosage.comb = dosage.comb,
  tox.SOC = tox.SOC, eff.SOC = eff.SOC, tox.B = tox.B, eff.B = eff.B, 
  tox.AB = tox.AB, eff.AB.s1 = eff.AB.s1, eff.AB.s2 = eff.AB.s2, 
  tox.isomat = matrix(c(2, 1), byrow = T, nrow = 1), 
  tox.upper = 0.30, eff.lower = 0.07, Cs = 0.80, Ct = 0.9, C.f1 = 0.90, C.f2 = 0.90, 
  utility.score = c(0, 40, 60, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
  )
```

    #> $stage1_output
    #> termination (%)       dose1 (%)       dose2 (%)   selection (%)              EN 
    #>            1.56           83.34           15.10           83.34           43.80 
    #> 
    #> $stage2_output
    #> Power (%)    GP (%)    SR (%)   OSR (%)        EN 
    #>     76.94     73.47     68.24     66.72    104.20

In stage 1, we have an 83.34% chance of selecting the correct
combination dose as the OBD, with an average sample size of 43.80
patients. In stage 2, the power of our design is 76.94%, the GP is
73.47%, the SR is 68.24%, and the OSR is 66.72%, with an average sample
size of 104.2 patients. In total, the trial requires an average of 148.0
patients (43.8 in stage 1 and 104.2 in stage 2). For illustration, here
we set `prior.sample = 1e5` (i.e., $10^5$ prior draws per simulation) to
balance computational speed and result accuracy. This code may take
approximately 10–20 minutes to run, depending on your system
specifications. **To reproduce the results in our paper, please use
`prior.sample = 1e6`, though this will require additional computation
time.**

If the ORRs of the combinations in stage 1 are 5% higher than the stage
2 rates (i.e., period effect = 0.05), run:

``` r
# Estimated runtime: ~9 minutes on a MacBook Air (M1, 16 GB RAM)
eff.AB.s1 <- c(0.240, 0.095) + 0.05
COCA.getOC(
  case = 2, n.stage1 = 24, n.stage2 = 38, Ce = 0.7710, c0 = 0.81, nlook.stage1 = 2, 
  nlook.stage2 = 2, dosage.ctrl = dosage.ctrl, dosage.singleA = 0, 
  dosage.singleB = dosage.singleB,  dosage.comb = dosage.comb,
  tox.SOC = tox.SOC, eff.SOC = eff.SOC, tox.B = tox.B, eff.B = eff.B, 
  tox.AB = tox.AB, eff.AB.s1 = eff.AB.s1, eff.AB.s2 = eff.AB.s2, 
  tox.isomat = matrix(c(2, 1), byrow = T, nrow = 1), 
  tox.upper = 0.30, eff.lower = 0.07, Cs = 0.80, Ct = 0.9, C.f1 = 0.90, C.f2 = 0.90, 
  utility.score = c(0, 40, 60, 100), rho = 0.2, prior.sample = 1e5, n.simu = 5000
  )
```

    #> $stage1_output
    #> termination (%)       dose1 (%)       dose2 (%)   selection (%)              EN 
    #>            0.36           81.38           18.26           81.38           45.87 
    #> 
    #> $stage2_output
    #> Power (%)    GP (%)    SR (%)   OSR (%)        EN 
    #>     77.72     73.12     70.75     68.04    105.10

##### \* Comparison with Dunnett’s test

This package also provides the code to reproduce the Dunnett’s test
approach described in Web Appendix C of \[1\]. To access it, run:

``` r
open_example("dunnett.R")
```

Detailed instructions for reproducing the results in paper \[1\] are
provided at the end of the script.

### References

\[1\]. Chi, X., Lin, R.<sup>\*</sup>, Yuan, Y.<sup>\*</sup> (2025+).
COCA: A Randomized Bayesian Design Integrating Dose Optimization and
Component Contribution Assessment for Combination Therapies.  
\[2\]. Kelley, R. K., Sangro, B., Harris, W., et al. Safety, Efficacy,
and Pharmacodynamics of Tremelimumab Plus Durvalumab for Patients With
Unresectable Hepatocellular Carcinoma: Randomized Expansion of a Phase
I/II Study. *Journal of Clinical Oncology*. 2021;39(27):2991-3001.
<doi:10.1200/JCO.20.03555>
