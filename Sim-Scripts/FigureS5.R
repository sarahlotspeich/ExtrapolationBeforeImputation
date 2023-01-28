# /////////////////////////////////////////////////////////////////////////
# Run simulation results for Figure S5 ////////////////////////////////////
# Compare CMI based on estimated survival function and trapezoidal rule ///
# w/ different extrapolation methods for Weibull X ////////////////////////
# /////////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd)
library(imputeCensRd) # To impute censored covariates 

# Data generation function based on Atem et al. (2017)'s "independent censoring" censoringing
## We note that they generate data from a Weibull with shape = 1 and scale varied, 
## But this is equivalent to an exponential (used here).
generate_AtemEtAl2017 = function(n, censoring = "light") {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  x = rweibull(n = n, shape = 0.75, scale = 0.25)  # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(test = censoring == "light", 
             yes = 0.5, ## ~ 12%
             no = ifelse(test = censoring == "heavy", 
                         yes = 2.9, ## ~ 41%
                         no = 20) ## ~ 78%
  ) # Rate parameter for censoring
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}

# Set the number of replicates per setting
reps = 1 ## We used a total of 1000, but see NOTES below

# Choose seed 
sim_seed = 114

tic("start")
# Loop over different censoring rates: light, heavy, extra heavy
for (censoring in c("light", "heavy", "extra_heavy")) {
  # And different sample sizes n = 100, 500, 1000, 2000
  for (n in c(100, 500, 1000, 2000)){
    # For reproducibility
    set.seed(sim_seed) 
    
    # Create dataframe to save results for setting
    sett_res = data.frame(sim = paste(sim_seed, 1:reps, sep = "-"), censoring, n, perc_censored = NA, 
                          alpha_d = NA, beta_d = NA, gamma_d = NA, 
                          alpha_e = NA, beta_e = NA, gamma_e = NA, 
                          alpha_w = NA, beta_w = NA, gamma_w = NA
    )
    
    # Loop over replicates 
    for (r in 1:reps) {
      # Generate data
      dat = generate_AtemEtAl2017(n = n, 
                                  censoring = censoring)
      
      # Save % censored
      sett_res$perc_censored[r] = 1 - mean(dat$d)
      
      # Loop over the three extrapolation methods
      for (extrap in c("d", "e", "w")) {
        ## Create imputed dataset
        imp_dat = cmi_sp(W = "w", Delta = "d", Z = "z", data = dat, 
                         trapezoidal_rule = TRUE, ## approximate integral using adaptive quadrature
                         surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                         surv_beyond = extrap) ## and extrapolated using the Weibull extension
        ## Check that imputation was successful 
        if (imp_dat$code) {
          ## Fit model to imputed dataset
          fit = lm(y ~ imp + z, data = imp_dat$imputed_data)
          sett_res[r, paste0(c("alpha_", "beta_", "gamma_"), extrap)] <- fit$coefficients
        }
      }
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0(censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}
toc()
# //////////////////////////////////////////////////////////////////////
# NOTES: When using the estimated survival function, we need to  ///////
# to fit an imputation model or use an extrapolation method. This //////
# makes the simulations for this table slower than Table 1. It took ////
# ~1 second to run 1 replication per setting MacBook Pro (M1) with 16GB/
# RAM. Based on this, it would take ~19 minutes to run 1000 replications
# per setting. We parallelized instead, using sim_seed = 114-133 and ///
# running reps = 50 replications per seed. (See FigureS5_Parallel.R) ///
# //////////////////////////////////////////////////////////////////////