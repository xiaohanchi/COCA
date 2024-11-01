# COCA
R code for implementing the COCA design. 

### Session Information
R version 4.3.1 was used for the simulations. Bayesian MCMC sampling were done by rjags (Download JAGS from
https://sourceforge.net/projects/mcmc-jags/). 


### Repository Contents
* scenarios.R: simulation scenarios
* calibration.R: get calibrated stage 2 sample size and design cutoffs. More detail is commented out at the beginning of each function. 
* run_COCA.R: run COCA design and get the summarized results. More detail is commented out at the beginning of each function. 
An example to run the code is provided at the end of each file. 

### Authors and Reference
Chi, X., Lin, R.<sup>\*</sup>, Yuan, Y.<sup>\*</sup> (2024+). COCA: A Randomized Bayesian Design Integrating Dose Optimization and Component Contribution Assessment for Combination Therapies. Under Revision in _Biometrics_.
