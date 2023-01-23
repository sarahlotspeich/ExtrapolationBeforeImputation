# //////////////////////////////////////////////////////////////////////
# Replicate Figure S1 in Supplementary Materials  //////////////////////
# Illustration of the four extrapolation methods for a step ////////////
# survival function $\widehat{S}(t)$ in simulated data...." ////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # For data wrangling
library(tidyr) # To collapse wide data
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(survival) # To fit step survival function 

# Be reproducible
set.seed(114)

# Generate data 
rateC = 20  # Rate parameter for censoring ("heavy")
x = rweibull(n = 1000, shape = 0.25, scale = 0.75) # To-be-censored covariate
c = rexp(n = 1000, rate = rateC) # Censoring variable
w = pmin(x, c) # Observed covariate 
d = as.numeric(x <= c) # "Event" indicator
data = data.frame(w, d) # Construct dataset
data = data |> 
  arrange(w)

# Fit Kaplan-Meier survival
fit = survfit(formula = Surv(time = w, event = d) ~ 1, 
              data = data, 
              timefix = FALSE) # need timefix = FALSE to avoid rounding issues

# Save \widetilde{X}: the largest uncensored covariate value 
# (beyond which survival needs to be extrapolated)
Xtilde = max(data$w[data$d == 1])

# Create index for rows needing extrapolation 
needs_extrap = which(data$d == 0 & data$w > Xtilde)

# //////////////////////////////////////////////////////////////////////
# Extrapolation method: Carry forward (default) ////////////////////////
# //////////////////////////////////////////////////////////////////////
## Create dataframe with times and survival from it
surv_df = with(fit, 
               data.frame(w = time, cf_surv = surv)) # Note: survfit() carries forward ("cf") by default beyond last "event"
## Merge  survival estimates into data
data = merge(x = data, 
             y = surv_df, 
             all.x = TRUE, 
             sort = TRUE)

# //////////////////////////////////////////////////////////////////////
# Extrapolation method: Immediate drop-off /////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Create variable for immediate drop-off survival 
data$drop_surv = data$cf_surv # start with carry forward 
data$drop_surv[needs_extrap] = 0 # replace extrapolated rows 

# //////////////////////////////////////////////////////////////////////
# Extrapolation method: Exponential extension //////////////////////////
# //////////////////////////////////////////////////////////////////////

# This function (included in imputeCensRd, but not exported) 
# Is one of the helper functions implementing the  extrapolation methods
extrap_surv_beyond <- function(x, t, surv, surv_beyond, weibull_params = NULL) {
  if (surv_beyond == "cf") { # Carry forward 
    before <- which(t <= x)
    surv[max(before)]
  } else if (surv_beyond == "d") { # Immediate drop-off
    0 
  } else if (surv_beyond == "e") { # Exponential extension
    exp(x * log(surv[length(surv)]) / t[length(t)])
  } else if (surv_beyond == "w") { # Weibull extension
    alpha_hat <- weibull_params[1]
    lambda_hat <- weibull_params[2]
    exp(- lambda_hat * x ^ alpha_hat)
  }
}
data$exp_surv = data$cf_surv # start with carry forward 
data$exp_surv[needs_extrap] = sapply(X = data$w[needs_extrap], # replace extrapolated rows 
                                     FUN = extrap_surv_beyond, # using extrap_surv_beyond()
                                     t = data$w[data$d == 1], # subset to event times
                                     surv = data$cf_surv[data$d == 1], # subset to survival at event times
                                     surv_beyond = "e") 

# //////////////////////////////////////////////////////////////////////
# Extrapolation method: Weibull extension //////////////////////////////
# //////////////////////////////////////////////////////////////////////

# For Weibull extension, need to find the constrained MLEs first 
S_Xtilde = min(data$cf_surv[data$d == 1]) # Save survival at largest uncensored covariate value (Xtilde)

# These functions (included in imputeCensRd, but not exported) 
# Help find the constrained MLEs for the Weibull extension
## Function 1: the "usual" log-likelihood for a Weibull
weibull_loglik = function(alpha, lambda, t, I_event) {
  n1 <- sum(I_event)
  ll <- - lambda * sum(t ^ alpha)
  ll <- ll + (alpha - 1) * sum(I_event * log(t)) 
  ll <- ll + n1 * log(lambda)
  ll <- ll + n1 * log(alpha)
  return(- ll)
}
## Fnuction 2: the constrained log-likelihood for a Weibull 
constr_weibull_loglik = function(alpha, t, I_event, Xtilde, rho) {
  lambda <- - log(rho) / (Xtilde ^ exp(alpha))
  n1 <- sum(I_event)
  weibull_loglik(alpha = alpha, lambda = lambda, t = t, I_event = I_event)
}
## Function 3: obtains the constrained Weibull MLEs
constr_weibull_mle = function(t, I_event, Xtilde, rho, alpha0, tol = 1E-4, max_iter = 1E3) {
  t[t == 0] = 1E-4 # Replace t = 0 with arbitrarily small but nonzero value to avoid errors
  suppressWarnings(
    nlm_res <- nlm(f = constr_weibull_loglik, 
                  p = alpha0, 
                  t = t, 
                  I_event = I_event, 
                  Xtilde = Xtilde, 
                  rho = rho)
  )
  conv = nlm_res$code <= 2 & nlm_res$iterations > 0 
  if (conv) {
    alpha1 = nlm_res$estimate
    lambda1 = - log(rho) / (Xtilde ^ alpha1)
    return(c(alpha1, lambda1))  
  } else {
    return(c(NA, NA))
  }
}

weibull_params = constr_weibull_mle(t = data$w, 
                                    I_event = data$d, 
                                    Xtilde = Xtilde, 
                                    rho = S_Xtilde, 
                                    alpha0 = 1E-4) # initial value for alpha (must be > 0)

# Then we can use weibull_params in extrap_surv_beyond() to complete Weibull extension
data$weib_surv = data$cf_surv # start with carry forward 
data$weib_surv[needs_extrap] = sapply(X = data$w[needs_extrap], # replace extrapolated rows 
                                      FUN = extrap_surv_beyond, # using extrap_surv_beyond()
                                      t = data$w[data$d == 1], # subset to event times
                                      surv = data$cf_surv[data$d == 1], # subset to survival at event times
                                      surv_beyond = "w", 
                                      weibull_params = weibull_params) 

# //////////////////////////////////////////////////////////////////////
# Create plot //////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
max_W = max(data$w)
data |>
  tidyr::gather(key = "extrapolation", value = "surv", -c(1:2)) |> 
  dplyr::mutate(extrapolation = factor(extrapolation, 
                                       levels = c("cf_surv", "drop_surv", "exp_surv", "weib_surv"),
                                       labels = c("Carry forward", "Immediate drop-off", "Exponential extension", "Weibull extension")),
                surv = as.numeric(surv)) |>
  ggplot(aes(x = w, y = surv, linetype = extrapolation)) + 
  geom_rect(aes(xmin = Xtilde, xmax = max_W, ymin = 0, ymax = 1), 
            fill = "aliceblue", col = "aliceblue", alpha = 0.3) + 
  geom_line() + 
  theme_minimal(base_size = 12) + theme(legend.position = "top") + 
  xlab(label = TeX("$t$")) + ylab(label = TeX("$\\hat{S}(t)$")) + 
  scale_linetype(name = "Extrapolation Method:")
  
# Save as 1000 wide x 500 tall