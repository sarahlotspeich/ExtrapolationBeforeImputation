# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption begins "Simulation results for Weibull $X$ from the full /////
# cohort analysis and imputation approaches using the true survival ////
# function and adaptive quadrature versus the trapezoidal rule..." /////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # For data wrangling
library(tidyr) # To gather wide tables
library(kableExtra) # To format pretty tables

# //////////////////////////////////////////////////////////////////////
# Read in simulation results from GitHub ///////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Table-Data/data_Table1.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with 100 reps per seed, per setting)
## This information is captured in the "sim" variable which is of the form seed-replicate. 

# Calculate average % censoring per censoring setting
res |> 
  group_by(censoring) |> 
  summarize(avg_perc_censored = mean(perc_censored))

# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Unpivot results to get rows for each parameter (rather than columns)
res_summ_long = res |> 
  dplyr::select(-perc_censored) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "est", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", "Adaptive Quadrature", "Trapezoidal Rule")),
         param = sub("_.*", "", param_calc),
         censoring = factor(x = censoring,
                            levels = c("light", "heavy", "extra_heavy"), 
                            labels = c("Light", "Heavy", "Extra Heavy"))) |> 
  group_by(censoring, n, calc, param) |> 
  summarize(bias = mean(est - ifelse(test = param == "alpha", 
                                     yes = 1,
                                     no = ifelse(test = param == "beta", 
                                                 yes = 0.5, 
                                                 no = 0.25
                                                 ))),
            se = sd(est))

# Then pivot them back out by method 
res_summ_wide = res_summ_long |> 
  pivot_wider(names_from = calc, 
              values_from = c("bias", "se")) |> 
  arrange(param, censoring) |> 
  mutate(`re_Adaptive Quadrature` = `se_Full Cohort` ^ 2 / `se_Adaptive Quadrature` ^ 2,
         `re_Trapezoidal Rule` = `se_Full Cohort` ^ 2 / `se_Trapezoidal Rule` ^ 2,
         mid1 = "", mid2 = "", mid3 = "") |> 
  dplyr::select(censoring, n, param, mid1,
                ends_with("Cohort"), mid2, 
                ends_with("Quadrature"), mid3, 
                ends_with("Rule")) 

# //////////////////////////////////////////////////////////////////////
# Format table for export to LaTex /////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Write function to add "padded" zeros and wrap with $$ for consistency 
format_num = function(num) {
  paste0("$", format(round(num, 3), nsmall = 3), "$")
}

# Format res_summ_wide for LaTex
res_summ_wide |> 
  dplyr::select(-param) |> # delete param column - ordering is alpha, beta, gamma
  mutate_if(.predicate = is.numeric, .funs = format_num) |> 
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = "llrcccccccrccccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters