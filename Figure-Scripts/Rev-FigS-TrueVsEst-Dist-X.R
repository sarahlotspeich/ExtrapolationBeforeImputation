# //////////////////////////////////////////////////////////////////////
# Run simulation results for Figure S# /////////////////////////////////

# Load packages
library(ggplot2) ## for pretty plots
library(dplyr) ## for data wrangling


res = do.call(what = rbind, 
              args = lapply(X = paste0("~/Documents/hybridCMI/Figure-Data/", 
                                       list.files("~/Documents/hybridCMI/Figure-Data/", pattern = "distX")), 
                            FUN = read.csv))


res_summ = res |> 
  group_by(n, censoring, distX, XdepZ) |> 
  summarize(avg_alpha = mean(alpha, na.rm = TRUE), 
            lambda = mean(lambda, na.rm = TRUE))
  

data.frame(x = rweibull(n = 10000, shape = 0.75, scale = 0.25)) |> 
  ggplot(aes(x = x)) + 
  stat_function(fun = pweibull, 
                args = list(shape = 0.75, scale = 0.25)) + 

alpha_hat <- weibull_params[1]
lambda_hat <- weibull_params[2]
exp(- lambda_hat * x ^ alpha_hat)