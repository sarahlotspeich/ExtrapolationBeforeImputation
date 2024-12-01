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
source_url("https://raw.githubusercontent.com/sarahlotspeich/extrapolationbeforeimputation/main/generate_data.R")

# Set the number of replicates per setting
reps = 100 ## We used a total of 1000, but see NOTES below

# Choose seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
sim_seed = 11422 + as.integer(args)

# Set censoring rate
censoring = "heavy" ## other options: "light" or "extra_heavy"

# Number of bootstraps for CSI-B standard errors
B = 1000

# Set distribution for X | Z
distX = "gamma"

# Loop over different sample sizes n = 100, 500, 2000
for (n in c(100, 500, 2000)){
  # For reproducibility
  set.seed(sim_seed) 
  
  # Create dataframe to save results for setting
  sett_res = data.frame(sim = paste(sim_seed, 1:reps, sep = "-"), censoring, n, distX, perc_censored = NA, 
                        alpha_fc = NA, beta_fc = NA, gamma_fc = NA, 
                        alpha_weibull = NA, beta_weibull = NA, gamma_weibull = NA, 
                        se_alpha_weibull = NA, se_beta_weibull = NA, se_gamma_weibull = NA)
  
  # Loop over replicates 
  for (r in 1:reps) {
    # Generate data
    dat = generate_data(n = n, ## Sample size
                        censoring = censoring, ## Censoring setting
                        distX = distX, ## Distribution for X
                        XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
    
    # Save % censored
    sett_res$perc_censored[r] = 1 - mean(dat$d)
    
    # Method 1: Full cohort analysis
    fit_fc = lm(y ~ x + z, data = dat)
    sett_res[r, c("alpha_fc", "beta_fc", "gamma_fc")] = fit_fc$coefficients
    
    # Method 2: Extrapolated CMI using Weibull extension 
    ## Create imputed dataset
    imp_dat = cmi_sp(imputation_model = Surv(time = w, event = d) ~ z,
                     data = dat, 
                     trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
                     surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                     surv_beyond = "w") ## and extrapolated using Weibull extension
    
    ## Check that imputation was successful 
    if (imp_dat$code) {
      if (B == 0) {
        ## Fit model to imputed dataset
        fit_aq = lm(y ~ imp + z, 
                    data = imp_dat$imputed_data)
        sett_res[r, est_cols] = fit_aq$coefficients
        sett_res[r, se_cols] = sqrt(diag(vcov(fit_aq)))
      } else {
        ## Fit model to B bootstrapped resamples from the imputed dataset
        all_fits = do.call(what = rbind, 
                           args = replicate(n = B,
                                            expr = coefficients(summary(lm(y ~ imp + z,
                                                                           data = imp_dat$imputed_data[sample(x = 1:n, 
                                                                                                              size = n, 
                                                                                                              replace = TRUE), ])
                                            )
                                            ),
                                            simplify = FALSE))
        
        ## Pool estimates and standard errors using Rubin's Rules
        pooled_est = rowsum(x = all_fits[, "Estimate"], group = rep(1:3, times = B)) / B
        sett_res[r, c("alpha_weibull", "beta_weibull", "gamma_weibull")] = pooled_est
        VW = rowsum(x = all_fits[, "Std. Error"] ^ 2, group = rep(1:3, times = B)) / B
        VB = rowsum(x = (all_fits[, "Estimate"] - pooled_est[rep(1:3, times = B)]) ^ 2, group = rep(1:3, times = B)) / (B - 1)
        pooled_var = VW + VB + (VB / B)
        sett_res[r, c("se_alpha_weibull", "se_beta_weibull", "se_gamma_weibull")] = sqrt(pooled_var)
      }
    }

    # Save results
    write.csv(x = sett_res,
              file = paste0(distX, "_", censoring, "_n", n, "_seed", sim_seed, ".csv"),
              row.names = F)
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: It took ~2 minutes to run 1 replication per setting MacBook ///
# Pro (M1) with 16GB RAM. We parallelized, using sim_seed = 114-133 ////
# and running reps = 50 replications per seed. /////////////////////////
# //////////////////////////////////////////////////////////////////////
