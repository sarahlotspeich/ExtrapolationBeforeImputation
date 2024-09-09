# //////////////////////////////////////////////////////////////////////
# Replicate Table 2 ////////////////////////////////////////////////////
# Caption begins "Simulation results" /////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To gather wide tables
library(kableExtra) # To format pretty tables
library(ggplot2) # To make pretty plots

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table2_rev.csv") |> 
  dplyr::filter(!(dist == "GAMMA" & n == 2000), ## exclude original runs of gamma with n = 2000 (lots of roundoff errors)
                !(dist == "LOGNORMAL" & n == 2000)) ## exclude original runs of log-normal with n = 2000 (lots of roundoff errors)
res = res |> 
  dplyr::bind_rows(
    read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table2_rev_absol1E2.csv")
  )
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# Standardize the conv_msg_aq 
res = res |> 
  mutate(conv_msg_aq = sub(pattern = ", ", replacement = "", x = conv_msg_aq), 
         div_err = grepl(pattern = "Divergent", x = conv_msg_aq), 
         round_err = grepl(pattern = "Roundoff", x = conv_msg_aq))

# //////////////////////////////////////////////////////////////////////
# Get convergence numbers for footnote /////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  summarize(reps_na_aq = sum(is.na(beta_aq)),
            reps_na_tr = sum(is.na(beta_tr))) 
## 3 replicates of Non-Extrapolated CMI out of 9,000 (0.03%) did not converge 
### (all presumably due to Cox model)
## 3 replicates of Extrapolated CMI out of 9,000 (0.03%) did not converge 
### (all presumably due to Cox model)

res |> 
  group_by(dist, censoring, n) |> 
  dplyr::summarize(reps_na_aq = sum(is.na(beta_aq)),
                   reps_na_tr = sum(is.na(beta_tr))
  ) |> 
  arrange(desc(reps_na_aq)) 
## Extrapolated CMI was successful in $\geq 99.9\%$ of replicates per setting
## Non-Extrapolated CMI was successful in $\geq 99.9\%$ of replicates per setting

# //////////////////////////////////////////////////////////////////////////////
# Add coverage indicators for confidence intervals /////////////////////////////
# //////////////////////////////////////////////////////////////////////////////
res = res |> 
  mutate(cov_beta_aq = (beta_aq - 1.96 * se_beta_aq) <= 0.5 & 0.5 <= (beta_aq + 1.96 * se_beta_aq),
         cov_beta_tr = (beta_tr - 1.96 * se_beta_tr) <= 0.5 & 0.5 <= (beta_tr + 1.96 * se_beta_tr))

# //////////////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////
## Bias and empirical standard errors ------------------------------------------
bias_ese_long = res |> 
  dplyr::select(sim, dist, n, beta_fc, beta_aq, beta_tr) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "est", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", 
                                  "Adaptive Quadrature", 
                                  "Trapezoidal Rule")),
         truth = 0.5) |> 
  group_by(dist, n, calc, truth) |> 
  dplyr::summarize(bias = mean(est - truth, na.rm = TRUE), # Exclude small number of replicates where Weibull didn't converge
                   ese = sd(est, na.rm = TRUE)) |> # Exclude small number of replicates where Weibull didn't converge
  mutate(perc_bias = paste0("($", format(round(bias / truth * 100, 2), nsmall = 2), "$)"),
         bias = paste0("$", format(round(bias, 3), nsmall = 3), "$")
  ) |> 
  ungroup() |> 
  select(dist, n, calc, bias, perc_bias, ese)

## Empirical coverage probabilities --------------------------------------------
cp_long = res |> 
  dplyr::select(sim, dist, n, cov_beta_aq, cov_beta_tr) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "cov", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", "Adaptive Quadrature", "Trapezoidal Rule"))) |> 
  group_by(dist, n, calc) |> 
  dplyr::summarize(cp = mean(cov, na.rm = TRUE)) # Exclude small number of replicates where Weibull didn't converge

# Unpivot results to get rows for each parameter (rather than columns)
ase_long = res |> 
  dplyr::select(sim, dist, n, se_beta_aq, se_beta_tr) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "se", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", "Adaptive Quadrature", "Trapezoidal Rule"))) |> 
  group_by(dist, n, calc) |> 
  dplyr::summarize(ase = mean(se, na.rm = TRUE)) |> # Exclude small number of replicates where Weibull didn't converge
  ungroup() |> 
  select(dist, n, calc, ase)

# Merge bias, ESE, and ASE 
res_summ_long = bias_ese_long |> 
  full_join(ase_long) |> 
  full_join(cp_long)

# Then pivot them back out by method 
res_summ_wide = res_summ_long |> 
  pivot_wider(names_from = calc, 
              values_from = c("bias", "perc_bias", "ese", "ase", "cp")) |> 
  arrange(dist) |> 
  mutate(`re_Adaptive Quadrature` = `ese_Full Cohort` ^ 2 / `ese_Adaptive Quadrature` ^ 2,
         `re_Trapezoidal Rule` = `ese_Full Cohort` ^ 2 / `ese_Trapezoidal Rule` ^ 2,
         mid1 = "", mid2 = "", mid3 = "") |> 
  dplyr::select(n, dist, mid1,
                ends_with("Cohort"), mid2, 
                ends_with("Quadrature"), mid3, 
                ends_with("Rule")) |> 
  dplyr::select(-`ase_Full Cohort`, -`cp_Full Cohort`)

# //////////////////////////////////////////////////////////////////////
# Format table for export to LaTex /////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Write function to add "padded" zeros and wrap with $$ for consistency 
format_num = function(num) {
  paste0("$", format(round(num, 3), nsmall = 3), "$")
}

# Format res_summ_wide for LaTex
res_summ_wide |>
  dplyr::select(-dist) |> ## delete param column - ordering is alpha, beta, gamma
  mutate_at(.vars = c("ese_Full Cohort", 
                      "ese_Adaptive Quadrature", "ase_Adaptive Quadrature", "cp_Adaptive Quadrature",  "re_Adaptive Quadrature", 
                      "ese_Trapezoidal Rule", "ase_Trapezoidal Rule", "cp_Trapezoidal Rule", "re_Trapezoidal Rule"),
            .funs = format_num) |>
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = "llrcccccccrccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
