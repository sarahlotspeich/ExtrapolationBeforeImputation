# /////////////////////////////////////////////////////////////////////////
# Run simulation results for Figure S4 ////////////////////////////////////
# Compare CMI based on estimated survival function and adaptive ///////////
# quadrature w/ different extrapolation methods for Log-Normal X //////////
# /////////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) # To impute censored covariates 

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/hybridCMI/main/generate_data.R")

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
                          alpha_d = NA, beta_d = NA, gamma_d = NA, 
                          alpha_e = NA, beta_e = NA, gamma_e = NA, 
                          alpha_w = NA, beta_w = NA, gamma_w = NA
    )
    
    # Loop over replicates 
    for (r in 1:reps) {
      # Generate data
      dat = generate_data(n = n, ## Sample size
                          censoring = censoring, ## Censoring setting
                          distX = "lognormal", ## Distribution for X
                          XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
      
      # Save % censored
      sett_res$perc_censored[r] = 1 - mean(dat$d)
      
      # Loop over the three extrapolation methods
      for (extrap in c("d", "e", "w")) {
        ## Create imputed dataset
        imp_dat = cmi_sp(imputation_model = Surv(time = w, event = d) ~ z, 
                         data = dat, 
                         trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
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
                file = paste0("FigureS4_", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: When using the estimated survival function, we need to  ///////
# to fit an imputation model or use an extrapolation method. This //////
# makes the simulations for this table slower than Table 1. It took ~30/
# minutes to run 1 replication per setting MacBook Pro (M1) with 16GB //
# RAM. Based on this, it would take ~519 hours to run 1000 replications/
# per setting. We parallelized instead, using sim_seed = 114-133 and ///
# running reps = 50 replications per seed. /////////////////////////////
# //////////////////////////////////////////////////////////////////////
