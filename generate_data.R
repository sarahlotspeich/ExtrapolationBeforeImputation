generate_data = function(n, censoring = "light", distX = "weibull", XZ = "indep") {
  z = rbinom(n = n, size = 1, prob = 0.5) # Uncensored covariate
  if (distX == "weibull") {
    if (XZ == "indep") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25)  # To-be-censored covariate
    } else if (XZ == "PH") {
      x = rweibull(n = n, shape = 0.75, scale = 0.25 + 0.15 * z)  # To-be-censored covariate
    } else if (XZ == "NPH") {
      x = rweibull(n = n, shape = 0.75 - 0.25 * z, scale = 0.25)  # To-be-censored covariate
    }
    
    # Rate parameter for censoring
    q = ifelse(test = censoring == "light", 
               yes = 0.5, ## ~ 12%
               no = ifelse(test = censoring == "heavy", 
                           yes = 2.9, ## ~ 41%
                           no = 20) ## ~ 78%
    ) 
  } else if (distX == "log-normal") {
    if (XZ == "indep") {
      x = rlnorm(n = n, meanlog = 0, sdlog = 0.5) # To-be-censored covariate
    } else if (XZ == "PH") {
      x = rlnorm(n = n, meanlog = 0 + 0.05 * z, sdlog = 0.5) # To-be-censored covariate
    } else if (XZ == "NPH") {
      x = rlnorm(n = n, meanlog = 0 + 0.25 * z, sdlog = 0.5) # To-be-censored covariate
    }
    
    # Rate parameter for censoring
    q = ifelse(test = censoring == "light", 
               yes = 0.2, # ~ 20% 
               no = ifelse(test = censoring == "heavy", 
                           yes = 0.4, # ~35%
                           no = 1.67) # ~79%
    ) 
  }
  e = rnorm(n = n, mean = 0, sd = 1) # Random errors
  y = 1 + 0.5 * x + 0.25 * z + e # Continuous outcome
  c = rexp(n = n, rate = q) # Random censoring mechanism
  w = pmin(x, c) # Observed covariate value
  d = as.numeric(x <= c) # "Event" indicator
  dat = data.frame(x, z, w, y, d) # Construct data set
  return(dat)
}

summary(replicate(n = 1000, 
               expr = mean(generate_data(n = 10000, 
                                         censoring = "light", 
                                         distX = "weibull", 
                                         XZ = "indep")$d == 0)))
# 

summary(replicate(n = 1000, 
                  expr = mean(generate_data(n = 100, 
                                            censoring = "extra heavy", 
                                            distX = "weibull", 
                                            XZ = "PH")$d == 0)))


summary(replicate(n = 1000, 
                  expr = mean(generate_data(n = 10000, 
                                            censoring = "light", 
                                            distX = "weibull", 
                                            XZ = "NPH")$d == 0)))
