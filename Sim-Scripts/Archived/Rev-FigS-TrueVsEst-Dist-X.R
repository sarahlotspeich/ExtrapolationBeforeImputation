# //////////////////////////////////////////////////////////////////////
# Run simulation results for Table 1 ///////////////////////////////////
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
reps = 100

# Choose seed 
sim_seed = 114

# Loop over different censoring rates: light, heavy, extra heavy
for (censoring in c("light", "heavy", "extra_heavy")) {
  # And different sample sizes n = 100, 500, 1000, 2000
  for (n in c(100, 500, 2000)){
    # For reproducibility
    set.seed(sim_seed) 
    
    # Create dataframe to save results for setting
    sett_res = data.frame(sim = paste(sim_seed, 1:reps, sep = "-"), censoring, n, perc_censored = NA, 
                          distX = "Weibull", XdepZ = TRUE, alpha = NA, lambda = NA)
    
    # Loop over replicates 
    for (r in 1:reps) {
      # Generate data
      dat = generate_data(n = n, ## Sample size
                          censoring = censoring, ## Censoring setting
                          distX = "weibull", ## Distribution for X
                          XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
      
      # Save % censored
      sett_res$perc_censored[r] = 1 - mean(dat$d)
      
      # Estimate Weibull parameters for the extension
      alpha_lambda = get_weib_params(imputation_model = Surv(time = w, event = d) ~ z, ## imputation model specification
                                     data = dat) ## dataset
      
      # Estimate the hazard ratio
      imp_fit = coxph(formula = Surv(time = w, event = d) ~ z, 
                      data = dat)
      
      # Use closed-form to compute the true conditional means
      alpha = 0.75
      lambda = 1 / ((0.25 + 0.25 * dat[, "z"]) ^ 0.75)
      ## Save quantities for use in formula
      inside_exp = lambda * dat[, "w"] ^ alpha ## inside exp() for Weibull survival function
      gamma_surv = pgamma(q = inside_exp,
                          shape = 1 / alpha,
                          scale = 1,
                          lower.tail = FALSE) ## survival function of a gamma
      dat$cm_true = dat[, "w"] * exp(- inside_exp) +
        gamma(1 / alpha) / (alpha * lambda ^ (1 / alpha)) * gamma_surv ## start with numerator
      dat$cm_true = dat$cm_true / exp(- inside_exp) ## divide by denominator
      
      # Calculate true conditional means 
      to_integrate = function(t, z) {
        pweibull(q = t, 
                 shape = 0.75, 
                 scale = 0.25 + 0.25 * z, 
                 lower.tail = FALSE)
      }
      int_surv = sapply(X = 1:n, 
                        FUN = function(i) {
                          tryCatch(expr = integrate(f = to_integrate, 
                                                    lower = dat[i, "w"], 
                                                    upper = Inf, 
                                                    subdivisions = 2000, 
                                                    z = dat[i, "z"])$value, 
                                   error = function(e) return(NA))}
      )
      dat$cm_aq = dat[, "w"] + int_surv / to_integrate(t = dat[, "w"], z = dat[, "z"])
      
      # Calculate estimated conditional means
      alpha_hat = alpha_lambda[1]
      lambda_tilde = alpha_lambda[2] * exp(imp_fit$coeff * dat$z)
      ## Save quantities for use in formula
      int_surv = exp(lambda_tilde * dat[, "w"] ^ alpha_hat) * 
        expint::gammainc(a = 1 / alpha_hat,
                         x = lambda_tilde * dat[, "w"] ^ alpha_hat)
      dat$cm_est = dat[, "w"] + int_surv / (alpha_hat * lambda_tilde ^ (1 / alpha_hat))
      
      # Save estimated Weibull parameters
      sett_res[r, c("alpha", "lambda")] = alpha_lambda
      
      # Save results
      write.csv(x = sett_res, 
                file = paste0("~/Documents/ExtrapolationBeforeImputation/Figure-Data/CMs_", censoring, "_n", n, "_seed", sim_seed, ".csv"), 
                row.names = F)
    }
  }
}