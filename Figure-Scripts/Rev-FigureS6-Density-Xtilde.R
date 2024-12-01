# //////////////////////////////////////////////////////////////////////
# Replicate Figure S8 in Supplementary Materials  //////////////////////
# Caption begins "Due to the Weibull distribution's skewness, //////////
# higher censoring rates led to smaller values of $W_{(n)}$ ..." ///////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # For data wrangling
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(wesanderson) # For colors
library(ggpubr) # To panel plots A and B together

# Be reproducible
set.seed(721)

# Load data generating function generate_data() from GitHub 
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/hybridCMI/main/generate_data.R")

# //////////////////////////////////////////////////////////////////////
# Create plot A: Distribution of $W_{(n)}$ (the maximum observed ///////
# covariate) when X was generated from a Weibull distribution  /////////
# //////////////////////////////////////////////////////////////////////

# For n = 1000, simulate 1000 datasets from each censoring setting
# And save the largest value of W = min(X, C) for each 
sims = 1000
Wn = vector(length = 3 * sims)
count = 1 # counter to index the Wn vector
for (censoring in c("light", "heavy", "extra heavy")) {
  for (s in 1:sims) {
    # Generate data
    dat = generate_data(n = 1000, ## Sample size
                        censoring = censoring, ## Censoring setting
                        distX = "weibull", ## Distribution for X
                        XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
    Wn[(count - 1) * sims + s] = max(dat$w) # Save maximum observed covariate
  }
  count = count + 1 # increment for next censoring 
}

# Add column with factor for censoring rate
Wn_df = data.frame(Censoring = rep(x = c("Light", "Heavy", "Extra Heavy"), 
                                   each = 1000), 
                   Wn) |> 
  mutate(Censoring = factor(Censoring, 
                            levels = c("Light", "Heavy", "Extra Heavy")))

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
count = 1 # counter to index the Wn vector
for (censoring in c("light", "heavy", "extra heavy")) {
  for (s in 1:sims) {
    # Generate data
    dat = generate_data(n = 1000, ## Sample size
                        censoring = censoring, ## Censoring setting
                        distX = "lognormal", ## Distribution for X
                        XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
    Wn[(count - 1) * sims + s] = max(dat$w) # Save maximum observed covariate
  }
  count = count + 1 # increment for next censoring 
}

# Add column with factor for censoring rate
Wn_df = data.frame(Censoring = rep(x = c("Light", "Heavy", "Extra Heavy"), 
                                   each = 1000), 
                   Wn) |> 
  mutate(Censoring = factor(Censoring, 
                            levels = c("Light", "Heavy", "Extra Heavy")))

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
# Create plot B: Distribution of $W_{(n)}$ (the maximum observed ///////
# covariate) when X was generated from a gamma distribution  ///////////
# //////////////////////////////////////////////////////////////////////

# For n = 1000, simulate 1000 datasets from each censoring setting
# And save the largest value of W = min(X, C) for each 
sims = 1000
Wn = vector(length = 3 * sims)
count = 1 # counter to index the Wn vector
for (censoring in c("light", "heavy", "extra heavy")) {
  for (s in 1:sims) {
    # Generate data
    dat = generate_data(n = 1000, ## Sample size
                        censoring = censoring, ## Censoring setting
                        distX = "gamma", ## Distribution for X
                        XdepZ = TRUE) ## Since TRUE, assume that X depends on Z
    Wn[(count - 1) * sims + s] = max(dat$w) # Save maximum observed covariate
  }
  count = count + 1 # increment for next censoring 
}

# Add column with factor for censoring rate
Wn_df = data.frame(Censoring = rep(x = c("Light", "Heavy", "Extra Heavy"), 
                                   each = 1000), 
                   Wn) |> 
  mutate(Censoring = factor(Censoring, 
                            levels = c("Light", "Heavy", "Extra Heavy")))

# Create plot 
x_lower = 0
x_upper = 31
plot_c = ggplot() +
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
ggpubr::ggarrange(plot_a, plot_b, plot_c, common.legend = TRUE, ncol = 1, labels = "AUTO")

# Save as 10" wide x 10" tall
ggsave("Rev_FigureS6.png", width = 10, height = 10, units = "in")
