# //////////////////////////////////////////////////////////////////////
# Run simulation results for Table 1 ///////////////////////////////////
# Compare full cohort analysis to "gold standard" CMI based on true ////
# function and adaptive quadrature vs. trapezoidal rule for Weibull X //
# //////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) # To impute censored covariates 

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/hybridCMI/main/generate_data.R")

# Write a function for the true survival function used to generate Weibull X 
trueSURV = function(q, z) {
  pweibull(q = q, 
           shape = 0.75, 
           scale = 0.25 + 0.25 * z,
           lower.tail = FALSE)
}

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
                          XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
      
      # Save % censored
      sett_res$perc_censored[r] = 1 - mean(dat$d)
      
      # Method 1: Full cohort analysis
      fit_fc = lm(y ~ x + z, data = dat)
      sett_res[r, c("alpha_fc", "beta_fc", "gamma_fc")] <- fit_fc$coefficients
      
      # Method 2: CMI using adaptive quadrature 
      ## Create imputed dataset
      imp_dat = cmi_custom(W = "w", Delta = "d", Z = "z", data = dat, 
                           useSURV = trueSURV, trapezoidal_rule = FALSE)
      ## Check that imputation was successful 
      ## (it always is when using trueSURV)
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_aq = lm(y ~ imp + z, data = imp_dat$imputed_data)
        sett_res[r, c("alpha_aq", "beta_aq", "gamma_aq")] <- fit_aq$coefficients
      }
      
      # Method 3: CMI using trapezoidal rule
      imp_dat = cmi_custom(W = "w", Delta = "d", Z = "z", data = dat, 
                           useSURV = trueSURV, trapezoidal_rule = TRUE)
      ## Check that imputation was successful 
      ## (it always is when using trueSURV)
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_tr = lm(y ~ imp + z, data = imp_dat$imputed_data)
        sett_res[r, c("alpha_tr", "beta_tr", "gamma_tr")] <- fit_tr$coefficients
      }
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0("Table1_", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: When using the true survival function, there is no need  //////
# to fit an imputation model or use an extrapolation method. This //////
# makes the simulations for this table quicker than other tables. //////
# It took <1 second to run 1 replication per setting MacBook Pro (M1) //
# with 16GB RAM. Based on this, it would take ~11 minutes to run 1000 //
# replications per setting. ////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
