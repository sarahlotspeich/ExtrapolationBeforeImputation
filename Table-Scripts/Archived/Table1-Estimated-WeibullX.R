# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption begins "Simulation results for Weibull $X$ from the full /////
# cohort analysis and imputation approaches using the estimated survival
# function and adaptive quadrature versus the trapezoidal rule..." /////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To gather wide tables
library(kableExtra) # To format pretty tables
library(ggplot2) # To make pretty plots

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table1.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# //////////////////////////////////////////////////////////////////////
# Get convergence numbers for footnote /////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  summarize(reps_na_aq = sum(is.na(beta_aq)),
            reps_na_tr = sum(is.na(beta_tr))
  ) ## Just 24 replicates out of 12,000 did not converge

res |> 
  group_by(censoring, n) |> 
  summarize(reps_na_aq = sum(is.na(beta_aq)),
            reps_na_tr = sum(is.na(beta_tr))
            ) |> 
  arrange(desc(reps_na_aq)) ## Converged in $\geq 99.4\%$ of replicates per setting

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
                            labels = c("Light", "Heavy", "Extra Heavy")),
         truth = ifelse(test = param == "alpha", 
                        yes = 1,
                        no = ifelse(test = param == "beta", 
                                    yes = 0.5, 
                                    no = 0.25)
         ) 
  ) |> 
  group_by(censoring, n, calc, param, truth) |> 
  summarize(bias = mean(est - truth, na.rm = TRUE), # Exclude small number of replicates where Weibull didn't converge
            se = sd(est, na.rm = TRUE)) |> # Exclude small number of replicates where Weibull didn't converge
  mutate(perc_bias = paste0("($", format(round(bias / truth * 100, 2), nsmall = 2), "$)"),
         bias = paste0("$", format(round(bias, 3), nsmall = 3), "$")
  ) |> 
  ungroup() |> 
  select(param, censoring, n, calc, bias, perc_bias, se)

# Then pivot them back out by method 
res_summ_wide = res_summ_long |> 
  pivot_wider(names_from = calc, 
              values_from = c("bias", "perc_bias", "se")) |> 
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
  mutate_at(.vars = c("se_Full Cohort", 
                      "se_Adaptive Quadrature", "re_Adaptive Quadrature", 
                      "se_Trapezoidal Rule", "re_Trapezoidal Rule"), .funs = format_num) |>
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = "llrcccccccrccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
