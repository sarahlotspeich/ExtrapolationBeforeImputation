# /////////////////////////////////////////////////////////////////////////
# Data generation function for all simulations ////////////////////////////
# /////////////////////////////////////////////////////////////////////////

generate_data = function(n, censoring = "light", distX = "weibull", XdepZ = TRUE, contZ = FALSE) {
  # Uncensored covariate
  if (contZ) {
    z = rlnorm(n = n, meanlog = 0, sdlog = 1 / 2) ## continuous Z
  } else {
    z = rbinom(n = n, size = 1, prob = 0.5) ## binary Z 
  }
  # To-be-censored covariate
  if (XdepZ) {
    if (distX == "weibull") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25 + 0.25 * z)
    } else if (distX == "lognormal") {
      x = rlnorm(n = n, meanlog = 0 + 0.05 * z, sdlog = 0.5)
    } else if (distX == "normal") {
      x = truncnorm::rtruncnorm(n = n, a = 0, b = Inf, mean = 1.5 + 0.5 * z, sd = 0.5)
    } else if (distX == "gamma") {
      x = rgamma(n = n, shape = 1.5, scale = 2 + 0.5 * z)
    } else if (distX == "t") {
      x = rt(n = n, df = 3 + z)
    }
  } else {
    if (distX == "weibull") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25)  
    } else if (distX == "lognormal") {
      x = rlnorm(n = n, meanlog = 0, sdlog = 0.5) 
    } else if (distX == "normal") {
      x = rnorm(n = n, mean = 1.5, sd = 0.5)
    } else if (distX == "gamma") {
      x = rgamma(n = n, shape = 1.5, scale = 2) 
    } else if (distX == "t") {
      x = rt(n = n, df = 3)
    }
  }
  # Random errors
  e = rnorm(n = n, mean = 0, sd = 1) 
  # Continuous outcome
  y = 1 + 0.5 * x + 0.25 * z + e 
  if (contZ) {
    # Rate parameter for censoring
    if (distX == "weibull") {
      q = ifelse(test = censoring == "light", 
                 yes = 0.23, ## ~ 12%
                 no = ifelse(test = censoring == "heavy", 
                             yes = 2, ## ~ 49%
                             no = 10) ## ~ 78%
      ) 
    } else if (distX == "lognormal") {
      q = ifelse(test = censoring == "light",
                 yes = 0.11, # ~ 20%
                 no = ifelse(test = censoring == "heavy",
                             yes = 0.62, # ~49%
                             no = 1.55) # ~78%
      )
    } 
    # else if (distX == "normal") {
    #   q = ifelse(test = censoring == "light",
    #              yes = 0.13, # ~ 20%
    #              no = ifelse(test = censoring == "heavy",
    #                          yes = 0.4, # ~49%
    #                          no = 1) # ~79%
    #   )
    # } else if (distX == "gamma") {
    #   q = ifelse(test = censoring == "light",
    #              yes = 0.07, # ~ 20%
    #              no = ifelse(test = censoring == "heavy",
    #                          yes = 0.27, # ~50%
    #                          no = 0.8) # ~79%
    #   )
    # } else if (distX == "t") {
    #   q = ifelse(test = censoring == "light",
    #              yes = 0.6, # ~ 20%
    #              no = ifelse(test = censoring == "heavy",
    #                          yes = 100, # ~50%
    #                          no = 10000) # ~47%
    #   )
    # }
  } else {
    # Rate parameter for censoring
    if (distX == "weibull") {
      q = ifelse(test = censoring == "light", 
                 yes = 0.5, ## ~ 12%
                 no = ifelse(test = censoring == "heavy", 
                             yes = 2.9, ## ~ 49%
                             no = 20) ## ~ 78%
      ) 
    } else if (distX == "lognormal") {
      q = ifelse(test = censoring == "light",
                 yes = 0.2, # ~ 20%
                 no = ifelse(test = censoring == "heavy",
                             yes = 0.65, # ~50%
                             no = 1.67) # ~79%
      )
    } else if (distX == "normal") {
      q = ifelse(test = censoring == "light",
                 yes = 0.13, # ~ 20%
                 no = ifelse(test = censoring == "heavy",
                             yes = 0.4, # ~49%
                             no = 1) # ~79%
      )
    } else if (distX == "gamma") {
      q = ifelse(test = censoring == "light",
                 yes = 0.07, # ~ 20%
                 no = ifelse(test = censoring == "heavy",
                             yes = 0.27, # ~50%
                             no = 0.8) # ~79%
      )
    } else if (distX == "t") {
      q = ifelse(test = censoring == "light",
                 yes = 0.6, # ~ 20%
                 no = ifelse(test = censoring == "heavy",
                             yes = 100, # ~50%
                             no = 10000) # ~47%
      )
    }
  } 
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}