# //////////////////////////////////////////////////////////////////////
# Run simulation results for Table 1 ///////////////////////////////////
# Compare full cohort analysis to CMI based on estimated survival / ////
# function and adaptive quadrature vs. trapezoidal rule for Weibull X //
# //////////////////////////////////////////////////////////////////////

# Load packages
## Run once: install.packages("devtools")
## Run once: devtools::install_github("sarahlotspeich/imputeCensRd")
library(imputeCensRd) # To impute censored covariates 

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/generate_data.R")

# Set the number of replicates per setting
reps = 100 ## We used a total of 1000, but see NOTES below

# Choose seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
start_sim_seed = 11422 + 100 * 0 #as.integer(args)

# Number of bootstraps for CSI-B standard errors
B = 1000

# Loop over different censoring rates: light, heavy, extra heavy
for (censoring in c("light", "heavy", "extra_heavy")) {
  # And different sample sizes n = 100, 500, 1000, 2000
  for (n in c(100, 500, 2000)){
    # Create dataframe to save results for setting
    sett_res = data.frame(sim = start_sim_seed + 1:reps, 
                          censoring, n, perc_censored = NA, 
                          alpha_fc = NA, beta_fc = NA, gamma_fc = NA, 
                          alpha_aq = NA, beta_aq = NA, gamma_aq = NA, 
                          se_alpha_aq = NA, se_beta_aq = NA, se_gamma_aq = NA, 
                          alpha_tr = NA, beta_tr = NA, gamma_tr = NA, 
                          se_alpha_tr = NA, se_beta_tr = NA, se_gamma_tr = NA, 
                          time_aq = NA, conv_msg_aq = NA, time_tr = NA)
    
    # Loop over replicates 
    for (r in 1:reps) {
      # For reproducibility
      set.seed(sett_res$sim[r]) 
      
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
      
      # Method 2: CMI using adaptive quadrature 
      ## Create imputed dataset
      time_aq = system.time(
        imp_dat <- cmi_sp(imputation_model = Surv(time = w, event = d) ~ z,
                       data = dat, 
                       trapezoidal_rule = FALSE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      )
      sett_res[r, "time_aq"] = as.numeric(time_aq["elapsed"]) ## save computing time
      
      ## Check that imputation was successful 
      if (imp_dat$code) {
        if (B == 0) {
          ## Fit model to imputed dataset
          fit_aq = lm(y ~ imp + z, 
                      data = imp_dat$imputed_data)
          sett_res[r, c("alpha_aq", "beta_aq", "gamma_aq")] = fit_aq$coefficients
          sett_res[r, c("se_alpha_aq", "se_beta_aq", "se_gamma_aq")] = sqrt(diag(vcov(fit_aq)))
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
          sett_res[r, c("alpha_aq", "beta_aq", "gamma_aq")] = pooled_est
          VW = rowsum(x = all_fits[, "Std. Error"] ^ 2, group = rep(1:3, times = B)) / B
          VB = rowsum(x = (all_fits[, "Estimate"] - pooled_est[rep(1:3, times = B)]) ^ 2, group = rep(1:3, times = B)) / (B - 1)
          pooled_var = VW + VB + (VB / B)
          sett_res[r, c("se_alpha_aq", "se_beta_aq", "se_gamma_aq")] = sqrt(pooled_var)
        }
      } else {
        sett_res$conv_msg_aq[r] = paste(unique(imp_dat$imputed_data$imp_msg), collapse = ", ")
      }
      
      # Method 3: CMI using trapezoidal rule
      time_tr = system.time(
        imp_dat <- cmi_sp(imputation_model = Surv(time = w, event = d) ~ z,
                       data = dat, 
                       trapezoidal_rule = TRUE, ## approximate integral using adaptive quadrature
                       surv_between = "cf", ## Breslow's estimator is carry-forward interpolated
                       surv_beyond = "w") ## and extrapolated using the Weibull extension
      )
      sett_res[r, "time_tr"] = as.numeric(time_tr["elapsed"]) ## save computing time
      
      ## Check that imputation was successful 
      if (imp_dat$code) {
        if (B == 0) {
          ## Fit model to imputed dataset
          fit_tr = lm(y ~ imp + z, 
                      data = imp_dat$imputed_data)
          sett_res[r, c("alpha_tr", "beta_tr", "gamma_tr")] = fit_tr$coefficients
          sett_res[r, c("se_alpha_tr", "se_beta_tr", "se_gamma_tr")] = sqrt(diag(vcov(fit_tr)))
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
          sett_res[r, c("alpha_tr", "beta_tr", "gamma_tr")] = pooled_est
          VW = rowsum(x = all_fits[, "Std. Error"] ^ 2, group = rep(1:3, times = B)) / B
          VB = rowsum(x = (all_fits[, "Estimate"] - pooled_est[rep(1:3, times = B)]) ^ 2, group = rep(1:3, times = B)) / (B - 1)
          pooled_var = VW + VB + (VB / B)
          sett_res[r, c("se_alpha_tr", "se_beta_tr", "se_gamma_tr")] = sqrt(pooled_var)
        }
      }
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0("rev_table1_", censoring, "_n", n, "_seed", start_sim_seed, ".csv"), 
                row.names = F)
    }
  }
}

# //////////////////////////////////////////////////////////////////////
# NOTES: When using the estimated survival function, we need to  ///////
# to fit an imputation model or use an extrapolation method. This //////
# makes the simulations for this table slower than Table 1. It took ~3 /
# minutes to run 1 replication per setting MacBook Pro (M1) with 16GB //
# RAM. Based on this, it would take ~50 hours to run 1000 replications /
# per setting. We parallelized instead, using sim_seed = 114-133 and ///
# running reps = 50 replications per seed. /////////////////////////////
# //////////////////////////////////////////////////////////////////////
