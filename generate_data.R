# //////////////////////////////////////////////////////////////////////
# Data generation function based on  ///////////////////////////////////
# Atem et al. (2017)'s "independent censoring" setting  ////////////////
# //////////////////////////////////////////////////////////////////////

generate_AtemEtAl2017 = function(n, censoring = "light", independent = TRUE) {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  x = rweibull(n = n, shape = 0.75, scale = 0.25)  # To-be-censored covariate
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  q = ifelse(censoring == "light", 2, 
             ifelse(censoring == "moderate", 0.35, 0.05)) # Rate parameter for censoring
  c = rweibull(n = n, shape = 1, scale = q) # Random censoring mechanism
  w <- pmin(x, c) # Observed covariate value
  d <- as.numeric(x <= c) # "Event" indicator
  dat = data.frame(z, w, y, d) # Construct data set
  return(dat)
}