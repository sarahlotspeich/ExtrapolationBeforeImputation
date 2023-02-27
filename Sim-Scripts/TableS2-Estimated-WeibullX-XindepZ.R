# //////////////////////////////////////////////////////////////////////
# Run simulation results for Table S2 ///////////////////////////////////
# Compare full cohort analysis to CMI based on estimated survival / ////
# function and adaptive quadrature vs. trapezoidal rule for Weibull X //
# where X was generated independently of Z /////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) # To impute censored covariates 

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/generate_data.R")

# Set the number of replicates per setting
reps = 1 ## We used a total of 1000, but see NOTES below

# Choose seed 
sim_seed = 114

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
      dat = generate_data(n = n, ## Sample size
                          censoring = censoring, ## Censoring setting
                          distX = "weibull", ## Distribution for X
                          XdepZ = FALSE) ## Since FALSE, assume that X is independent of Z
      
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
                file = paste0("TableS2_", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: When using the estimated survival function, we need to  ///////
# to fit an imputation model or use an extrapolation method. This //////
# makes the simulations for this table slower than Table 1. It took ~4 /
# minutes to run 1 replication per setting MacBook Pro (M1) with 16GB //
# RAM. Based on this, it would take ~63 hours to run 1000 replications /
# per setting. We parallelized instead, using sim_seed = 114-133 and ///
# running reps = 50 replications per seed. /////////////////////////////
# //////////////////////////////////////////////////////////////////////