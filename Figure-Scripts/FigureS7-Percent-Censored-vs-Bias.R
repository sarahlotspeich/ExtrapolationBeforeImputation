# //////////////////////////////////////////////////////////////////////
# Replicate Figure S7 in Supplementary Materials  //////////////////////
# Caption begins "Due to the Weibull distribution's skewness, //////////
# higher censoring rates led to smaller values of $W_{(n)}$ ..." ///////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To transform data
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots
library(scales) # To get % labels for plots
library(wesanderson) # To get fun colors

# //////////////////////////////////////////////////////////////////////
# Read in simulation results  //////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

## (A) Using *true* survival function
res_true = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Table-Data/data_Table1_PH.csv")

# Calculate average % censoring per censoring setting
res_true |> 
  group_by(censoring) |> 
  summarize(min_perc_censored = min(perc_censored),
            avg_perc_censored = mean(perc_censored),
            max_perc_censored = max(perc_censored))
## Light 17% (5-29%)
## Heavy 49% (32-67%)
## Extra heavy 82% (65-92%)

## (B) Using *estimated* survival function
res_est = read.csv(file = "https://raw.githubusercontent.com/sarahlotspeich/ItsIntegral/main/Table-Data/data_Table2_PH.csv")

# Calculate average % censoring per censoring setting
res_est |> 
  group_by(censoring) |> 
  summarize(min_perc_censored = min(perc_censored),
            avg_perc_censored = mean(perc_censored),
            max_perc_censored = max(perc_censored))
## Light 17% (5-32%)
## Heavy 49% (33-63%)
## Extra heavy 82% (68-95%)

# //////////////////////////////////////////////////////////////////////
# Reshape simulation results ///////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Unpivot results to get rows for each parameter (rather than columns)
res_summ_long = res_est |> 
  gather(key = "param_calc", value = "est", -c(1:4)) |> 
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
         ),
         bias = est - truth,
         surv = "Estimated"
  ) |> 
  bind_rows(
    res_true |> 
      gather(key = "param_calc", value = "est", -c(1:4)) |> 
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
             ),
             bias = est - truth,
             surv = "Truth"
      )
  )

# //////////////////////////////////////////////////////////////////////
# Create plot //////////////////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////
res_summ_long |> 
  filter(calc != "Full Cohort") |> 
  mutate(surv = factor(surv, levels = c("Truth", "Estimated"), 
                       labels = c(TeX("True $S(x|z)$"), 
                                  TeX("Estimated $\\widehat{S}(x|z)$"))),
         param = factor(param, levels = c("alpha", "beta", "gamma"), 
                        labels = c(TeX("$\\hat{\\alpha}$: Intercept"),
                                   TeX("$\\hat{\\beta}$: Coefficient on $X$"),
                                   TeX("$\\hat{\\gamma}$: Coefficient on $Z$")
                                   )),
         ) |> 
  ggplot(aes(x = perc_censored, y = bias / truth, col = calc, fill = calc)) + # / truth
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(cols = vars(surv), rows = vars(param), scales = "free", labeller = label_parsed) + 
  xlab("Percent Censored") + 
  ylab("Percent Bias") + 
  scale_color_manual(values = wes_palette("Zissou1", n = 5, type = "discrete")[c(1,5)],
                     name = "Method") + 
  scale_fill_manual(values = wes_palette("Zissou1", n = 5, type = "discrete")[c(1,5)],
                     name = "Method") + 
  theme_minimal(base_size = 14) + 
  theme(legend.position = "top") + 
  scale_x_continuous(labels = percent, breaks = c(0, 0.17, 0.49, 0.815, 1)) + 
  scale_y_continuous(labels = percent)

# Save as 10" wide x 10" tall
ggsave("FigureS7.png", width = 10, height = 10, units = "in")