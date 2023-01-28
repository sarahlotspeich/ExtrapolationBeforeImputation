# //////////////////////////////////////////////////////////////////////
# Replicate Figure S6 in Supplementary Materials  //////////////////////
# Caption begins "Due to the Weibull distribution's skewness, //////////
# higher censoring rates led to smaller values of $W_{(n)}$ ////////////
# (the maximum of the observed covariate), which led to worse //////////
# performance (i.e., higher bias) when calculating the conditional /////
# mean with the trapezoidal rule... //////////////////////////////////// 
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # For data wrangling
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(wesanderson) # For colors
library(ggpubr) # To panel plots A and B together

# Be reproducible
set.seed(721)

# //////////////////////////////////////////////////////////////////////
# Create plot A: Distribution of $W_{(n)}$ (the maximum observed ///////
# covariate) when X was generated from a Weibull distribution  /////////
# //////////////////////////////////////////////////////////////////////

# For n = 1000, simulate 1000 datasets from each censoring setting
# And save the largest value of W = min(X, C) for each 
sims = 1000
Wn = vector(length = 3 * sims)
for (cens in 1:3) {
  rateC = c(0.5, 2.9, 20)[cens] # Rate parameter for censoring
  for (s in 1:sims) {
    x = rweibull(n = 1000, shape = 0.25, scale = 0.75) # True covariate
    c = rexp(n = 1000, rate = rateC) # Censoring variable
    w = pmin(x, c) # Observed covariate 
    Wn[(cens - 1) * sims + s] = max(w) # Save maximum observed covariate
  }
}

# Add column with factor for censoring rate
Wn_df = data.frame(Censoring = rep(x = c("Light", "Moderate", "Heavy"), 
                                   each = 1000), 
                   Wn) |> 
  mutate(Censoring = factor(Censoring, 
                            levels = c("Light", "Moderate", "Heavy")))

# Create plot 
x_lower = 0
x_upper = 22
plot_a = ggplot() +
  xlim(x_lower, x_upper) +
  geom_density(data = Wn_df, mapping = aes(x = Wn, fill = Censoring), alpha = 0.5) + 
  xlab(TeX("$W_{(n)}$")) +
  ylab("Density") +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") + 
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) 

# //////////////////////////////////////////////////////////////////////
# Create plot B: Distribution of $W_{(n)}$ (the maximum observed ///////
# covariate) when X was generated from a log-normal distribution  //////
# //////////////////////////////////////////////////////////////////////

# For n = 1000, simulate 1000 datasets from each censoring setting
# And save the largest value of W = min(X, C) for each 
sims = 1000
Wn = vector(length = 3 * sims)
for (cens in 1:3) {
  rateC = c(0.2, 0.4, 1.67)[cens] # Rate parameter for censoring
  for (s in 1:sims) {
    x = rlnorm(n = 1000, meanlog = 0, sdlog = 0.5) # True covariate
    c = rexp(n = 1000, rate = rateC) # Censoring variable
    w = pmin(x, c) # Observed covariate 
    Wn[(cens - 1) * sims + s] = max(w) # Save maximum observed covariate
  }
}

# Add column with factor for censoring rate
Wn_df = data.frame(Censoring = rep(x = c("Light", "Moderate", "Heavy"), 
                                   each = 1000), 
                   Wn) |> 
  mutate(Censoring = factor(Censoring, 
                            levels = c("Light", "Moderate", "Heavy")))

# Create plot 
x_lower = 0
x_upper = 11
plot_b = ggplot() +
  xlim(x_lower, x_upper) +
  geom_density(data = Wn_df, mapping = aes(x = Wn, fill = Censoring), alpha = 0.5) + 
  xlab(TeX("$W_{(n)}$")) +
  ylab("Density") +
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") + 
  scale_fill_manual(values = wes_palette("Zissou1")[c(1,3,5)]) + 
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) 

# //////////////////////////////////////////////////////////////////////
# Put plots A and B side-by-side  //////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
ggpubr::ggarrange(plot_a, plot_b, common.legend = TRUE, ncol = 1, labels = "AUTO")