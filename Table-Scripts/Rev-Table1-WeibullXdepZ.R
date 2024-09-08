# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption "Simulation results for Weibull $X$ dependent on $Z$ from ////
# the full cohort analysis (i.e., where all $n$ observations had ///////
# uncensored $X$) and conditional mean imputation (CMI) approaches." ///
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To gather wide tables
library(kableExtra) # To format pretty tables
library(ggplot2) # To make pretty plots

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table1_rev.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# Standardize the conv_msg_aq 
res = res |> 
  mutate(conv_msg_aq = sub(pattern = ", ", replacement = "", x = conv_msg_aq), 
         div_err = grepl(pattern = "Divergent", x = conv_msg_aq), 
         round_err = grepl(pattern = "Roundoff", x = conv_msg_aq))

# //////////////////////////////////////////////////////////////////////
# Summarize computing times ////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
## Overall average per-replicate computing times:
res |> 
  filter(!(round_err | div_err)) |> 
  select(starts_with("time")) |> 
  summarize_all(mean)
### 24.6 seconds (extrapolated) vs. 0.12 (non-extrapolated)

## By setting average per-replicate computing times:
res |> 
  group_by(n, censoring) |> 
  select(starts_with("time")) |> 
  summarize_all(mean)
### For extrapolated CMI, computing time increased with sample size, 
### but for all sample sizes extra heavy < light < heavy

# //////////////////////////////////////////////////////////////////////
# Get convergence numbers for footnote /////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  summarize(reps_na_aq = sum(is.na(beta_aq)),
            reps_na_tr = sum(is.na(beta_tr))) 
## 16 replicates of Non-Extrapolated CMI out of 9,000 (0.2%) did not converge 
### (all presumably due to Cox model)
## 44 replicates of Extrapolated CMI out of 9,000 (0.5%) did not converge 
### (presumably 16 due to the Cox model; 
### the other 28 due to roundoff or potential divergence error)

res |> 
  group_by(censoring, n) |> 
  dplyr::summarize(reps_na_aq = sum(is.na(beta_aq)),
                   reps_na_tr = sum(is.na(beta_tr))
  ) |> 
  arrange(desc(reps_na_aq)) 
## Extrapolated CMI was successful in $\geq 96.8\%$ of replicates per setting
## Non-Extrapolated CMI was successful in $\geq 99.4\%$ of replicates per setting

## The roundoff and potential divergence errors were all in n = 100, extra heavy censoring setting
res |> 
  filter(div_err | round_err) |> 
  group_by(censoring, n, div_err, round_err) |> 
  summarize(num = n())

# //////////////////////////////////////////////////////////////////////////////
# Add coverage indicators for confidence intervals /////////////////////////////
# //////////////////////////////////////////////////////////////////////////////
res = res |> 
  mutate(cov_alpha_aq = (alpha_aq - 1.96 * se_alpha_aq) <= 1 & 1 <= (alpha_aq + 1.96 * se_alpha_aq),
         cov_beta_aq = (beta_aq - 1.96 * se_beta_aq) <= 0.5 & 0.5 <= (beta_aq + 1.96 * se_beta_aq),
         cov_gamma_aq = (gamma_aq - 1.96 * se_gamma_aq) <= 0.25 & 0.25 <= (gamma_aq + 1.96 * se_gamma_aq),
         cov_alpha_tr = (alpha_tr - 1.96 * se_alpha_tr) <= 1 & 1 <= (alpha_tr + 1.96 * se_alpha_tr),
         cov_beta_tr = (beta_tr - 1.96 * se_beta_tr) <= 0.5 & 0.5 <= (beta_tr + 1.96 * se_beta_tr),
         cov_gamma_tr = (gamma_tr - 1.96 * se_gamma_tr) <= 0.25 & 0.25 <= (gamma_tr + 1.96 * se_gamma_tr)
         )

# //////////////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////////////
# //////////////////////////////////////////////////////////////////////////////
## Bias and empirical standard errors ------------------------------------------
bias_ese_long = res |> 
  dplyr::select(-perc_censored, -starts_with(c("se_", "cov", "time")), -ends_with("err"), -conv_msg_aq) |> # use package prefix to avoid conflict with MASS::select
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
                                    yes = 0.50, 
                                    no = 0.25)
         ) 
  ) |> 
  group_by(censoring, n, calc, param, truth) |> 
  dplyr::summarize(bias = mean(est - truth, na.rm = TRUE), # Exclude small number of replicates where Weibull didn't converge
            ese = sd(est, na.rm = TRUE)) |> # Exclude small number of replicates where Weibull didn't converge
  mutate(perc_bias = paste0("($", format(round(bias / truth * 100, 2), nsmall = 2), "$)"),
         bias = paste0("$", format(round(bias, 3), nsmall = 3), "$")
  ) |> 
  ungroup() |> 
  select(param, censoring, n, calc, bias, perc_bias, ese)

## Empirical coverage probabilities --------------------------------------------
cp_long = res |> 
  dplyr::select(sim, censoring, n, starts_with("cov")) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "cov", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", "Adaptive Quadrature", "Trapezoidal Rule")),
         param = sub("_.*", "", sub(".*cov_", "", param_calc)),
         censoring = factor(x = censoring,
                            levels = c("light", "heavy", "extra_heavy"), 
                            labels = c("Light", "Heavy", "Extra Heavy"))
  ) |> 
  group_by(censoring, n, calc, param) |> 
  dplyr::summarize(cp = mean(cov, na.rm = TRUE)) # Exclude small number of replicates where Weibull didn't converge

# Unpivot results to get rows for each parameter (rather than columns)
ase_long = res |> 
  dplyr::select(sim, censoring, n, starts_with("se")) |> # use package prefix to avoid conflict with MASS::select
  gather(key = "param_calc", value = "se", -c(1:3)) |> 
  mutate(calc = gsub(pattern = ".*_", replacement = "", x = param_calc),
         calc = factor(x = calc, 
                       levels = c("fc", "aq", "tr"),
                       labels = c("Full Cohort", "Adaptive Quadrature", "Trapezoidal Rule")),
         param = sub("_.*", "", sub(".*se_", "", param_calc)),
         censoring = factor(x = censoring,
                            levels = c("light", "heavy", "extra_heavy"), 
                            labels = c("Light", "Heavy", "Extra Heavy"))
  ) |> 
  group_by(censoring, n, calc, param) |> 
  dplyr::summarize(ase = mean(se, na.rm = TRUE)) |> # Exclude small number of replicates where Weibull didn't converge
  ungroup() |> 
  select(param, censoring, n, calc, ase)

# Merge bias, ESE, and ASE 
res_summ_long = bias_ese_long |> 
  full_join(ase_long) |> 
  full_join(cp_long)

# Then pivot them back out by method 
res_summ_wide = res_summ_long |> 
  pivot_wider(names_from = calc, 
              values_from = c("bias", "perc_bias", "ese", "ase", "cp")) |> 
  arrange(param, censoring) |> 
  mutate(`re_Adaptive Quadrature` = `ese_Full Cohort` ^ 2 / `ese_Adaptive Quadrature` ^ 2,
         `re_Trapezoidal Rule` = `ese_Full Cohort` ^ 2 / `ese_Trapezoidal Rule` ^ 2,
         mid1 = "", mid2 = "", mid3 = "") |> 
  dplyr::select(censoring, n, param, mid1,
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
  dplyr::select(-param) |> ## delete param column - ordering is alpha, beta, gamma
  mutate_at(.vars = c("ese_Full Cohort", 
                      "ese_Adaptive Quadrature", "ase_Adaptive Quadrature", "cp_Adaptive Quadrature",  "re_Adaptive Quadrature", 
                      "ese_Trapezoidal Rule", "ase_Trapezoidal Rule", "cp_Trapezoidal Rule", "re_Trapezoidal Rule"),
            .funs = format_num) |>
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = "llrcccccccrccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
