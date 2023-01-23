# //////////////////////////////////////////////////////////////////////
# Replicate Figure S4 in Supplementary Materials  //////////////////////
# Caption begins "Interpolating Breslow's estimator $\widehat{S}_0(t)$ /
# between uncensored values with either of the two interpolation ///////
# methods offered similar bias and efficiency ..." /////////////////////
# //////////////////////////////////////////////////////////////////////
  
# Load packages
library(dplyr) # For data wrangling
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(wesanderson) # For colors

# //////////////////////////////////////////////////////////////////////
# Read in simulation results from GitHub ///////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Fig-Data/sims_FigureS4.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with 100 reps per seed, per setting)
## This information is captured in the "sim" variable which is of the form seed-replicate. 

# //////////////////////////////////////////////////////////////////////
# Create plot //////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |>
  dplyr::mutate(censoring = factor(x = censoring, 
                                   levels = c("Light", "Moderate", "Heavy"),
                                   labels = c("Light Censoring", 
                                              "Moderate Censoring",
                                              "Heavy Censoring")),
                n = factor(x = n, 
                           levels = c(100, 500, 1000, 2000),
                           labels = paste("n =", c(100, 500, 1000, 2000))),
                interp = factor(x = interpolate, 
                                levels = c("carry-forward", "mean"),
                                labels = c("Carry forward", "Mean"))) |>
  ggplot(aes(x = factor(n), 
             y = beta_aq, 
             fill = interp)) +
  geom_hline(yintercept = 0.5, linetype = 2) + 
  geom_boxplot(alpha = 0.7) + 
  facet_grid(rows = vars(censoring), scales = "free") +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") +
  scale_fill_manual(values = wes_palette("Zissou1", n = 5, type = "discrete")[c(1,3,5)],
                    name = TeX("Method of Interpolation for $\\hat{S}_{0}(t)$")) +
  xlab("Sample Size") +
  ylab(TeX("Parameter Estimate $\\hat{\\beta}$")) 
## Note: 44 rows will be removed because beta_aq is NA 
## These are the instances where the Weibull extension did not converge