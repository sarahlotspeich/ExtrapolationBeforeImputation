# /////////////////////////////////////////////////////////////////////////
# Data generation function for all simulations ////////////////////////////
# /////////////////////////////////////////////////////////////////////////

generate_data = function(n, censoring = "light", distX = "weibull", XdepZ = TRUE) {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  if (XdepZ) {
    if (distX == "weibull") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25 + 0.25 * z)  # To-be-censored covariate    
    } else if (distX == "lognormal") {
      x = rlnorm(n = n, meanlog = 0 + 0.05 * z, sdlog = 0.5) # To-be-censored covariate
    }
  } else {
    if (distX == "weibull") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25)  # To-be-censored covariate  
    } else if (distX == "lognormal") {
      x = rlnorm(n = n, meanlog = 0, sdlog = 0.5) # To-be-censored covariate
    }
  }
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  if (distX == "weibull") {
    q = ifelse(test = censoring == "light", 
               yes = 0.5, ## ~ 12%
               no = ifelse(test = censoring == "heavy", 
                           yes = 2.9, ## ~ 41%
                           no = 20) ## ~ 78%
    ) # Rate parameter for censoring
  } else if (distX == "lognormal") {
    q = ifelse(test = censoring == "light", 
               yes = 0.2, # ~ 20% 
               no = ifelse(test = censoring == "heavy", 
                           yes = 0.4, # ~35%
                           no = 1.67) # ~79%
    ) # Rate parameter for censoring
  }
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}