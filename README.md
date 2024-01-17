# A hybrid approach to impute heavily censored covariates with their conditional means

This repository contains R code and simulation data to reproduce results from the [manuscript](https://arxiv.org/abs/2209.04716) by Lotspeich and Garcia (2022+). 

For the `imputeCensRd` package, which implements the conditional mean imputation approaches from the paper, can be found in its own repo [here](https://github.com/sarahlotspeich/imputeCensRd). 

Each of the "Script (Run Simulations)" files is coded to run 1 replication of each setting for demonstration. Per the NOTES at the bottom of the scripts, some more time-intensive simulations were run in parallel.

## Tables 

**Table 1.** Simulation results for Weibull $X$ from the full cohort analysis and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations)](Sim-Scripts/Table1-Estimated-WeibullX.R)
  - [Script (Make Table)](Table-Scripts/Table1-Estimated-WeibullX.R)
  - [Data (Simulation Results)](Table-Data/data_Table1.csv)  

**Table S1.** Simulation results for Weibull $X$ independent of $Z$ from the full cohort analysis (i.e., where all $n$ observations had uncensored $X$) and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations)](Sim-Scripts/TableS1-Estimated-WeibullX-XindepZ.R)
  - [Script (Make Table)](Table-Scripts/TableS1-Estimated-WeibullX-XindepZ.R)
  - [Data (Simulation Results)](Table-Data/data_TableS1.csv)  

**Table S2.** Simulation results for log-normal $X$ from the full cohort analysis (i.e., where all $n$ observations had uncensored $X$) and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations)](Sim-Scripts/TableS2-Estimated-LogNormalX.R)
  - [Script (Make Table)](Table-Scripts/TableS2-Estimated-LogNormalX.R)
  - [Data (Simulation Results)](Table-Data/data_TableS2.csv)  

## Figures 

**Figure S1.** Illustration of the four extrapolation methods for a step survival function $\widehat{S}(t)$ in simulated data.

  - [Script (Make Figure)](Figure-Scripts/FigureS1-Illustrate-Extrapolation-Methods.R)

**Figure S2.** We explored light ($\sim 17\%$), heavy ($\sim 49\%$), and extra heavy ($\sim 82\%$) censoring in Weibull $X$, induced by generating $C$ from an exponential distribution with rates $= 0.5$, $2.9$, and $20$, respectively.

  - [Script (Make Figure)](Figure-Scripts/FigureS2-Percent-Censored.R)
  
**Figure S3.** With Weibull $X$, extrapolating Breslow's estimator $\widehat{S}_0(t)$ beyond the largest uncensored value $\widetilde{X}$ with the Weibull extension offered the lowest bias and best efficiency for $\hat{\beta}$ in conditional mean imputation with adaptive quadrature.

  - [Script (Run Simulations)](Sim-Scripts/FigureS3-Extrapolation-Methods-Weibull.R)
  - [Script (Make Figure)](Figure-Scripts/FigureS3-Extrapolation-Methods-Weibull.R)
  - [Data (Simulation Results)](Figure-Data/data_FigureS3.csv)  

**Figure S4.** With log-normal $X$, extrapolating Breslow's estimator $\widehat{S}_0(t)$ beyond the largest uncensored value $\widetilde{X}$ with any of the three extrapolation methods offered similar bias and efficiency for $\hat{\beta}$ in conditional mean imputation with adaptive quadrature.

  - [Script (Run Simulations)](Sim-Scripts/FigureS4-Extrapolation-Methods-Log-Normal.R)
  - [Script (Make Figure)](Figure-Scripts/FigureS4-Extrapolation-Methods-Log-Normal.R)
  - [Data (Simulation Results)](Figure-Data/data_FigureS4.csv)  

**Figure S5.** Interpolating Breslow's estimator $\widehat{S}_0(t)$ between uncensored values with either of the two interpolation methods offered similar bias and efficiency for $\hat{\beta}$ in conditional mean imputation with adaptive quadrature. 

  - [Script (Run Simulations)](Sim-Scripts/FigureS5-Interpolation-Methods.R)
  - [Script (Make Figure)](Figure-Scripts/FigureS5-Interpolation-Methods.R)
  - [Data (Simulation Results)](Figure-Data/data_FigureS5.csv)  

**Figure S6.** Extrapolating Breslow's estimator $\widehat{S}_0(t)$ beyond the largest uncensored value $\widetilde{X}$ with any of the three extrapolation methods offered similar bias and efficiency for $\hat{\beta}$ in conditional mean imputation with the trapezoidal rule.

  - [Script (Run Simulations)](Sim-Scripts/FigureS6-Extrapolation-Methods-Trapezoidal-Rule.R)
  - [Script (Make Figure)](Figure-Scripts/FigureS6-Extrapolation-Methods-Trapezoidal-Rule.R)
  - [Data (Simulation Results)](Figure-Data/data_FigureS6.csv)  

**Figure S7.** Due to the Weibull distribution's skewness, higher censoring rates led to smaller values of $W_{(n)}$ (the maximum of the observed covariate), which led to worse performance (i.e., higher bias) when calculating the conditional mean with the trapezoidal rule.

  - [Script (Make Figure)](Figure-Scripts/FigureS8-Weibull-vs-Log-Normal.R) 
