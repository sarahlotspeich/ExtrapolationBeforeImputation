# /////////////////////////////////////////////////////////////////////////
# Run simulation results for Table S1 /////////////////////////////////////
# Compare full cohort analysis to CMI based on estimated survival /////////
# function and adaptive quadrature vs. trapezoidal rule for Log-Normal X //
# /////////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd)
library(imputeCensRd) # To impute censored covariates 

# Data generation function 
## This setting is new; it was created for this manuscript to consider a non-Weibull X
generate_data = function(n, censoring = "light") {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  x = rweibull(n = n, shape = 0.75, scale = 0.25)  # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(censoring == "light", 2, 
             ifelse(censoring == "heavy", 0.35, 0.05)) # Rate parameter for censoring
  c = rweibull(n = n, shape = 1, scale = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}

# Set the number of replicates per setting
reps = 50 ## Run 50 each across 20 jobs in parallel for 1000 total per setting

# Choose seed 
## Get job ID
args = commandArgs(TRUE) 
## Add job ID to seed to get unique seeds for each set
sim_seed = 114 + as.integer(args) 

# Loop over different censoring rates: light, heavy, extra heavy
for (censoring in c("light", "heavy", "extra_heavy")) {
  # And different sample sizes n = 100, 500, 1000, 2000
  for (n in c(100, 500, 1000, 2000)){
    # For reproducibility
    set.seed(sim_seed) 
    
    # Create dataframe to save results for setting
    sett_res = data.frame(sim = paste(sim_seed, 1:reps, sep = "-"), censoring, n, perc_censored = NA, 
                          alpha_fc = NA, beta_fc = NA, gamma_fc = NA, 
                          alpha_aq = NA, beta_aq = NA, gamma_aq = NA, 
                          alpha_tr = NA, beta_tr = NA, gamma_tr = NA
    )
    
    # Loop over replicates 
    for (r in 1:reps) {
      # Generate data
      dat = generate_AtemEtAl2017(n = n, 
                                  censoring = censoring)
      
      # Save % censored
      sett_res$perc_censored[r] = 1 - mean(dat$d)
      
      # Method 1: Full cohort analysis
      fit_fc = lm(y ~ x + z, data = dat)
      sett_res[r, c("alpha_fc", "beta_fc", "gamma_fc")] <- fit_fc$coefficients
      
      # Method 2: CMI using adaptive quadrature 
      ## Create imputed dataset
      imp_dat = cmi_sp(W = "w", Delta = "d", Z = "z", data = dat, 
                       trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      ## Check that imputation was successful 
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_aq = lm(y ~ imp + z, data = imp_dat$imputed_data)
        sett_res[r, c("alpha_aq", "beta_aq", "gamma_aq")] <- fit_aq$coefficients
      }
      
      # Method 3: CMI using trapezoidal rule
      imp_dat = cmi_sp(W = "w", Delta = "d", Z = "z", data = dat, 
                       trapezoidal_rule = TRUE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      ## Check that imputation was successful 
      ## (it always is when using trueSURV)
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_tr = lm(y ~ imp + z, data = imp_dat$imputed_data)
        sett_res[r, c("alpha_tr", "beta_tr", "gamma_tr")] <- fit_tr$coefficients
      }
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0(censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}