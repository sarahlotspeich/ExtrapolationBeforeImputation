# //////////////////////////////////////////////////////////////////////
# Replicate Table 1 ////////////////////////////////////////////////////
# Caption begins "Simulation results for $\pmb{\beta}$ (the coefficient/
# on censored Weibull $X$) from the full cohort analysis and ///////////
# conditional mean imputation (CMI) approaches." ///////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(kableExtra) # To format pretty tables

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table1_R1_CSB.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# //////////////////////////////////////////////////////////////////////
# Get convergence numbers for footnote /////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  summarize(reps_na_extrap = sum(is.na(beta_extrap)),
            reps_na_extrap_no = sum(is.na(beta_extrap_no))
  ) ## 13 replicates out of 9,000 did not converge (~ 0.1%)

res |> 
  group_by(censoring, n) |> 
  summarize(reps_na_extrap = sum(is.na(beta_extrap)),
            reps_na_extrap_no = sum(is.na(beta_extrap_no))
            ) |> 
  arrange(desc(reps_na_extrap)) ## Converged in 99.7%+ of replicates per setting

# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting //////////////////////////////
# //////////////////////////////////////////////////////////////////////
## Exclude small number of replicates where Weibull didn't converge (above)
res_summ = res |> 
  group_by(censoring, n) |> 
  summarize(bias_fc = mean(beta_fc - 0.5, na.rm = TRUE), 
            ese_fc = sd(beta_fc, na.rm = TRUE),
            bias_extrap = mean(beta_extrap - 0.5, na.rm = TRUE), 
            ese_extrap = sd(beta_extrap, na.rm = TRUE), 
            ase_extrap = mean(se_beta_extrap, na.rm = TRUE), 
            cp_extrap = mean((beta_extrap - 1.96 * se_beta_extrap) <= 0.5 & 
                               0.5 <= (beta_extrap + 1.96 * se_beta_extrap), 
                             na.rm = TRUE),
            bias_extrap_no = mean(beta_extrap_no - 0.5, na.rm = TRUE), 
            ese_extrap_no = sd(beta_extrap_no, na.rm = TRUE), 
            ase_extrap_no = mean(se_beta_extrap_no, na.rm = TRUE), 
            cp_extrap_no = mean((beta_extrap_no - 1.96 * se_beta_extrap_no) <= 0.5 & 
                               0.5 <= (beta_extrap_no + 1.96 * se_beta_extrap_no), 
                             na.rm = TRUE)
            ) |> 
  mutate(perc_bias_fc = paste0("($", format(round(bias_fc / 0.5 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap = paste0("($", format(round(bias_extrap / 0.5 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap_no = paste0("($", format(round(bias_extrap_no / 0.5 * 100, 2), nsmall = 2), "$)"),
         censoring = factor(x = censoring, 
                            levels = c("light", "heavy", "extra_heavy"), 
                            labels = c("Light", "Heavy", "Extra Heavy")),
         re_extrap = ese_fc ^ 2 / ese_extrap ^ 2,
         re_extrap_no = ese_fc ^ 2 / ese_extrap_no ^ 2,
         mid1 = "", mid2 = "", mid3 = "") |> 
  select(censoring, n, starts_with(c("bias", "perc_bias", "ese", "ase", "cp", "re")), everything()) |> 
  arrange(censoring, n) |> 
  select(censoring, n, mid1,
         ends_with("_fc"), mid2, 
         ends_with("_extrap"), mid3, 
         ends_with("_extrap_no")) 

# //////////////////////////////////////////////////////////////////////
# Format table for export to LaTex /////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Write function to add "padded" zeros and wrap with $$ for consistency 
format_num = function(num) {
  paste0("$", format(round(num, 3), nsmall = 3), "$")
}

# Format res_summ_wide for LaTex
res_summ |> 
  mutate_at(.vars = c("bias_fc", "ese_fc",
                      "bias_extrap", "ese_extrap", "ase_extrap", "cp_extrap", "re_extrap",
                      "bias_extrap_no", "ese_extrap_no", "ase_extrap_no", "cp_extrap_no", "re_extrap_no"), 
            .funs = format_num) |>
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        align = "llrcccccccrccccc") |> 
  kable_styling() 
## Note: For visual reasons, the \addlinespace were manually deleted in LaTex
## And a \multicolumn used to separate the three parameters
