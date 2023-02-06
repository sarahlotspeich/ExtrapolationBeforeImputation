# //////////////////////////////////////////////////////////////////////
# Replicate Figure S2 in Supplementary Materials  //////////////////////
# Caption begins "We explored light (~17%), heavy (~49%), and extra ////
# heavy (~82%) censoring in Weibull $X$, induced by generating ... /////
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
res = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Table-Data/data_Table1_PH.csv")

# Calculate average % censoring per censoring setting
res |> 
  group_by(censoring) |> 
  summarize(min_perc_censored = min(perc_censored),
            avg_perc_censored = mean(perc_censored),
            max_perc_censored = max(perc_censored))
## Light 17% (5-29%)
## Heavy 49% (32-67%)
## Extra heavy 82% (65-92%)

# //////////////////////////////////////////////////////////////////////
# Create plot //////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res |> 
  mutate(censoring = factor(censoring, levels = c("light", "heavy", "extra_heavy"),
                            labels = c("Light Censoring", 
                                       "Heavy Censoring",
                                       "Extra Heavy Censoring")),
         n = factor(x = n, 
                    levels = c(100, 500, 1000, 2000),
                    labels = paste("n =", c(100, 500, 1000, 2000)))
  ) |> 
  ggplot(aes(x = n, y = perc_censored, fill = censoring)) + 
  geom_boxplot(alpha = 0.7) + 
  theme_minimal(base_size = 14) +
  geom_hline(yintercept = 0.17, linetype = 2) + 
  geom_hline(yintercept = 0.49, linetype = 2) + 
  geom_hline(yintercept = 0.815, linetype = 2) + 
  xlab("Sample Size") + 
  ylab("Percent Censored Observations") + 
  scale_fill_manual(values = wes_palette("Zissou1", n = 5, type = "discrete")[c(1,3,5)], name = "") +
  theme(legend.position = "top") + 
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.17, 0.49, 0.815, 1))

# Save as 10" wide x 6" tall
ggsave("FigureS2.png", width = 10, height = 6, units = "in")