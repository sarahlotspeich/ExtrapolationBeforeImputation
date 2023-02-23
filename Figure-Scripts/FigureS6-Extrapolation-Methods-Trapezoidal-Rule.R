# //////////////////////////////////////////////////////////////////////
# Replicate Figure S6 in Supplementary Materials  //////////////////////
# Caption begins "Extrapolating Breslow's estimator $\widehat{S}_0(t)$ /
# beyond ... with any of the three extrapolation methods offered ///////
# similar bias and efficiency ... with the trapezoidal rule" ///////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # For data wrangling
library(tidyr) # To transform data
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(wesanderson) # For colors

# //////////////////////////////////////////////////////////////////////
# Read in simulation results from GitHub ///////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Read in simulation results 
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Figure-Data/data_FigureS6.csv")
## Note: Simulations were run in parallel on random seeds 114-123 (with 100 reps per seed, per setting)
## This information is captured in the "sim" variable which is of the form "seed-replicate." 

# Focus on estimates of beta for figure 
## and transform from wide (one column per extrapolation method per replicate)
## to long (one row per extrapolation method per replicate)
beta_res = res |> 
  select(sim, censoring, n, starts_with("beta")) |> 
  gather(key = "extrap", value = "beta", -c(1:3))

# //////////////////////////////////////////////////////////////////////
# Create plot //////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
beta_res |>
  dplyr::mutate(censoring = factor(censoring, levels = c("light", "heavy", "extra_heavy"),
                                   labels = c("Light Censoring", 
                                              "Heavy Censoring",
                                              "Extra Heavy Censoring")),
                n = factor(x = n, 
                           levels = c(100, 500, 1000, 2000),
                           labels = paste("n =", c(100, 500, 1000, 2000))),
                extrap = factor(x = extrap, 
                                levels = c("beta_d", "beta_e", "beta_w"),
                                labels = c("Immediate drop-off", "Exponential extension", "Weibull extension"))
  ) |>
  ggplot(aes(x = n, y = beta, fill = extrap)) +
  geom_hline(yintercept = 0.5, linetype = 2) + 
  geom_boxplot(alpha = 0.7) + 
  facet_grid(rows = vars(censoring), scales = "free") +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") +
  scale_fill_manual(values = wes_palette("Zissou1", n = 5, type = "discrete")[c(1,3,5)],
                    name = "Extrapolation Method:") +
  xlab("Sample Size") +
  ylab(TeX("Parameter Estimate $\\hat{\\beta}$")) 

# Save as 10" wide x 10" tall
ggsave("FigureS6.png", width = 10, height = 10, units = "in")

## Note: 24 rows will be removed because beta is NA 
beta_res |> 
  group_by(extrap, censoring, n) |> 
  summarize(reps_na = sum(is.na(beta))) |> 
  arrange(desc(reps_na))
## These are the instances where the Weibull extension did not converge
## There were <= 6 instances per setting (<=0.6%)