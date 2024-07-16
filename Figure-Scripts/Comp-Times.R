# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption begins "Simulation results for $\pmb{\alpha}$ (the coefficient/
# on censored Weibull $X$) from the full cohort analysis and ///////////
# conditional mean imputation (CMI) approaches." ///////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(ggplot2) # To make pretty graphs
library(tidyr) # To "gather" data

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Figure-Data/data_comp_times.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# Get average comp time per setting 
res_summ = res |> 
  group_by(n, censoring) |> 
  summarize(mean_time_extrap = mean(time_extrap), 
            median_time_extrap = median(time_extrap), 
            mean_time_extrap_no = mean(time_extrap_no),
            median_time_extrap_no = median(time_extrap_no))

# Make plot 
res_summ |> 
  gather(key = "avg_method", value = "time", -c(1:2)) |> 
  mutate(summ = ifelse(test = grepl("mean", avg_method), 
                       "Mean", 
                       "Median"),
         method = ifelse(test = grepl("extrap_no", avg_method), 
                         yes = "Non-Extrapolated CMI", 
                         no = "Extrapolated CMI")) |> 
  ggplot(aes(x = n, y = time, color = method, linetype = summ)) + 
  geom_line() + 
  facet_wrap(~censoring) + 
  labs(x = "Sample Size")
