# //////////////////////////////////////////////////////////////////////
# Run simulation results for revised Table 1 ///////////////////////////
# Simulation results for Weibull X from the full cohort analysis and ///
# conditional mean imputation (CMI) approaches /////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) ## To impute censored covariates 

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/hybridCMI/main/generate_data.R")

# Set the number of replicates per setting
reps = 50 ## We used a total of 1000, but see NOTES below

# Choose seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
sim_seed = 11422 + as.integer(args)

# Set censoring rate
censoring = "heavy" ## other options: "light" or "extra_heavy"

# Loop over different sample sizes n = 100, 500, 2000
for (n in c(100, 500, 2000)){
  # For reproducibility
  set.seed(sim_seed) 
  
  # Create dataframe to save results for setting
  sett_res = data.frame(sim = paste(sim_seed, 1:reps, sep = "-"), censoring, n, perc_censored = NA, 
                        alpha_fc = NA, beta_fc = NA, gamma_fc = NA, 
                        alpha_extrap = NA, beta_extrap = NA, gamma_extrap = NA, 
                        se_alpha_extrap = NA, se_beta_extrap = NA, se_gamma_extrap = NA, 
                        alpha_extrap_no = NA, beta_extrap_no = NA, gamma_extrap_no = NA,
                        se_alpha_extrap_no = NA, se_beta_extrap_no = NA, se_gamma_extrap_no = NA,
                        time_extrap = NA, time_extrap_no = NA
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
    sett_res[r, c("alpha_fc", "beta_fc", "gamma_fc")] = fit_fc$coefficients
    
    # Method 2: Extrapolated CMI (B = 20 imputations via resampling coefficients)
    time_imp = system.time(
      fit_extrap <- csi_w_boot(imputation_model = Surv(time = w, event = d) ~ z, ## imputation model specification
                               analysis_model = y ~ imp + z, ## analysis model specification 
                               data = dat, ## dataset
                               trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
                               surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                               surv_beyond = "w", ## and extrapolated using the Weibull extension
                               B = 500) ## number of bootstraps
    )
    sett_res[r, "time_extrap"] = time_imp["elapsed"]
    sett_res[r, c("alpha_extrap", "beta_extrap", "gamma_extrap")] = fit_extrap$Est
    sett_res[r, c("se_alpha_extrap", "se_beta_extrap", "se_gamma_extrap")] = fit_extrap$SE
    
    # Method 3: CMI using trapezoidal rule
    time_imp = system.time(
      fit_extrap_no <- csi_w_boot(imputation_model = Surv(time = w, event = d) ~ z, ## imputation model specification
                                  analysis_model = y ~ imp + z, ## analysis model specification 
                                  data = dat, ## dataset
                                  trapezoidal_rule = TRUE, ## approximate integral using trapezoidal rule
                                  surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                                  surv_beyond = "w", ## and extrapolated using the Weibull extension
                                  B = 500) ## number of bootstraps
    )
    sett_res[r, "time_extrap_no"] = time_imp["elapsed"]
    sett_res[r, c("alpha_extrap_no", "beta_extrap_no", "gamma_extrap_no")] = fit_extrap_no$Est
    sett_res[r, c("se_alpha_extrap_no", "se_beta_extrap_no", "se_gamma_extrap_no")] = fit_extrap_no$SE
    
    # Save results
    write.csv(x = sett_res, 
              file = paste0("Table1-Rev-CSB/", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
              row.names = F)
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: It took ~2 minutes to run 1 replication per setting MacBook ///
# Pro (M1) with 16GB RAM. We parallelized, using sim_seed = 114-133 ////
# and running reps = 50 replications per seed. /////////////////////////
# //////////////////////////////////////////////////////////////////////
