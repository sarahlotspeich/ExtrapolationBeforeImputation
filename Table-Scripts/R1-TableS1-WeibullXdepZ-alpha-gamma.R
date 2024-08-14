# //////////////////////////////////////////////////////////////////////
# Replicate Table S1 ///////////////////////////////////////////////////
# Caption begins "Simulation results for $\pmb{\alpha}$ (the coefficient/
# on censored Weibull $X$) from the full cohort analysis and ///////////
# conditional mean imputation (CMI) approaches." ///////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(kableExtra) # To format pretty tables

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ExtrapolationBeforeImputation/main/Table-Data/data_Table1_rev.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with ~100 reps per seed, per setting)

# //////////////////////////////////////////////////////////////////////
# Get convergence numbers for footnote /////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  summarize(reps_na_extrap = sum(is.na(alpha_extrap)),
            reps_na_extrap_no = sum(is.na(alpha_extrap_no))
  ) ## No replicates out of 9,000 did not converge

# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting (alpha) //////////////////////
# //////////////////////////////////////////////////////////////////////
## Exclude small number of replicates where Weibull didn't converge (above)
res_summ = res |> 
  group_by(censoring, n) |> 
  summarize(bias_fc = mean(alpha_fc - 1, na.rm = TRUE), 
            ese_fc = sd(alpha_fc, na.rm = TRUE),
            bias_extrap = mean(alpha_extrap - 1, na.rm = TRUE), 
            ese_extrap = sd(alpha_extrap, na.rm = TRUE), 
            ase_extrap = mean(se_alpha_extrap, na.rm = TRUE), 
            cp_extrap = mean((alpha_extrap - 1.96 * se_alpha_extrap) <= 1 & 
                               1 <= (alpha_extrap + 1.96 * se_alpha_extrap), 
                             na.rm = TRUE),
            bias_extrap_no = mean(alpha_extrap_no - 1, na.rm = TRUE), 
            ese_extrap_no = sd(alpha_extrap_no, na.rm = TRUE), 
            ase_extrap_no = mean(se_alpha_extrap_no, na.rm = TRUE), 
            cp_extrap_no = mean((alpha_extrap_no - 1.96 * se_alpha_extrap_no) <= 1 & 
                               1 <= (alpha_extrap_no + 1.96 * se_alpha_extrap_no), 
                             na.rm = TRUE)
            ) |> 
  mutate(perc_bias_fc = paste0("($", format(round(bias_fc / 1 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap = paste0("($", format(round(bias_extrap / 1 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap_no = paste0("($", format(round(bias_extrap_no / 1 * 100, 2), nsmall = 2), "$)"),
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


# //////////////////////////////////////////////////////////////////////
# Summarize simulation results by setting (gamma) //////////////////////
# //////////////////////////////////////////////////////////////////////
## Exclude small number of replicates where Weibull didn't converge (above)
res_summ = res |> 
  group_by(censoring, n) |> 
  summarize(bias_fc = mean(gamma_fc - 0.25, na.rm = TRUE), 
            ese_fc = sd(gamma_fc, na.rm = TRUE),
            bias_extrap = mean(gamma_extrap - 0.25, na.rm = TRUE), 
            ese_extrap = sd(gamma_extrap, na.rm = TRUE), 
            ase_extrap = mean(se_gamma_extrap, na.rm = TRUE), 
            cp_extrap = mean((gamma_extrap - 1.96 * se_gamma_extrap) <= 0.25 & 
                               0.25 <= (gamma_extrap + 1.96 * se_gamma_extrap), 
                             na.rm = TRUE),
            bias_extrap_no = mean(gamma_extrap_no - 0.25, na.rm = TRUE), 
            ese_extrap_no = sd(gamma_extrap_no, na.rm = TRUE), 
            ase_extrap_no = mean(se_gamma_extrap_no, na.rm = TRUE), 
            cp_extrap_no = mean((gamma_extrap_no - 1.96 * se_gamma_extrap_no) <= 0.25 & 
                                  0.25 <= (gamma_extrap_no + 1.96 * se_gamma_extrap_no), 
                                na.rm = TRUE)
  ) |> 
  mutate(perc_bias_fc = paste0("($", format(round(bias_fc / 0.25 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap = paste0("($", format(round(bias_extrap / 0.25 * 100, 2), nsmall = 2), "$)"),
         perc_bias_extrap_no = paste0("($", format(round(bias_extrap_no / 0.25 * 100, 2), nsmall = 2), "$)"),
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
