# //////////////////////////////////////////////////////////////////////
# Run simulation results for Table 1 ///////////////////////////////////
# Simulation results for Weibull X from the full cohort analysis and ///
# conditional mean imputation (CMI) approaches /////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) ## To impute censored covariates 
library(sandwich) ## To get sandwich covariance

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/hybridCMI/main/generate_data.R")

# Set the number of replicates per setting
reps = 50 ## We used a total of 1000, but see NOTES below

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
                          alpha_extrap = NA, beta_extrap = NA, gamma_extrap = NA, 
                          se_alpha_extrap = NA, se_beta_extrap = NA, se_gamma_extrap = NA, 
                          alpha_non_extrap = NA, beta_non_extrap = NA, gamma_non_extrap = NA,
                          se_alpha_non_extrap = NA, se_beta_non_extrap = NA, se_gamma_non_extrap = NA
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
      
      # Method 2: Extrapolated CMI
      ## Create imputed dataset
      imp_dat = cmi_sp(imputation_model = Surv(time = w, event = d) ~ z, ## imputation model specification
                       data = dat, ## dataset
                       trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      
      ## Check that imputation was successful 
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_extrap = lm(y ~ imp + z, 
                        data = imp_dat$imputed_data)
        sett_res[r, c("alpha_extrap", "beta_extrap", "gamma_extrap")] = fit_extrap$coefficients
        
        ## Estimate sandwich covariance
        sandcov = sandwich(x = fit_extrap)
        sett_res[r, c("se_alpha_extrap", "se_beta_extrap", "se_gamma_extrap")] = sqrt(diag(sandcov))
      }
      
      # Method 3: CMI using trapezoidal rule
      imp_dat = cmi_sp(imputation_model = Surv(time = w, event = d) ~ z, ## imputation model specification
                       data = dat, ## dataset
                       trapezoidal_rule = TRUE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      ## Check that imputation was successful 
      ## (it always is when using trueSURV)
      if (imp_dat$code) {
        ## Fit model to imputed dataset
        fit_non_extrap = lm(y ~ imp + z, 
                            data = imp_dat$imputed_data)
        sett_res[r, c("alpha_non_extrap", "beta_non_extrap", "gamma_non_extrap")] = fit_non_extrap$coefficients
        
        ## Estimate sandwich covariance
        sandcov = sandwich(x = fit_non_extrap)
        sett_res[r, c("se_alpha_non_extrap", "se_beta_non_extrap", "se_gamma_non_extrap")] = sqrt(diag(sandcov))
      }
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0("Rev-Table1/", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: It took ~3 minutes to run 1 replication per setting MacBook ///
# Pro (M1) with 16GB RAM. Based on this, it would take ~50 hours to ////
# run 1000 replications per setting. We parallelized instead, using ////
# sim_seed = 114-133 and running reps = 50 replications per seed. //////
# //////////////////////////////////////////////////////////////////////
