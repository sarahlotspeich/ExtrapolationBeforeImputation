# Extrapolation before imputation reduces bias when imputing censored covariates

This repository contains R code and simulation data to reproduce results from the [manuscript](https://arxiv.org/abs/2209.04716) by Lotspeich and Garcia (2022+). 

For the `imputeCensRd` package, which implements the conditional mean imputation approaches from the paper, can be found in its own repo [here](https://github.com/sarahlotspeich/imputeCensRd). 

Each of the "Script (Run Simulations)" files is coded to run 1 replication of each setting for demonstration. Per the NOTES at the bottom of the scripts, some more time-intensive simulations were run in parallel.

## Tables 

**Table 1.** Simulation results for Weibull $X$ dependent on $Z$ from the full cohort analysis (i.e., where all $n$ observations had uncensored $X$) and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations)](Sim-Scripts/Rev-Table1-WeibullXdepZ.R)
  - [Script (Make Table)](Table-Scripts/Rev-Table1-WeibullXdepZ.R)
  - [Data (Simulation Results)](Table-Data/data_Table1_rev.csv)  

**Table 2.** Simulation results for $\hat{\beta}$ under various distributions of $X$ dependent on $Z$ from the full cohort analysis (i.e., where all $n$ observations had uncensored $X$) and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations with Gamma Distribution)](Sim-Scripts/Rev-Table2-Misspec-Gamma.R)
  - [Script (Run Simulations with Log-Normal Distribution)](Sim-Scripts/Rev-Table2-Misspec-LogNormal.R)
  - [Script (Run Simulations with Normal Distribution)](Sim-Scripts/Rev-Table2-Misspec-Normal.R)
  - [Script (Make Table)](Table-Scripts/Rev-Table2-Misspec.R)
  - [Data (Simulation Results - Smaller Absolute Tolerance)](Table-Data/data_Table2_rev.csv)
  - [Data (Simulation Results - Larger Absolute Tolerance)](Table-Data/data_Table2_rev_absol1E2.csv)  

**Table S1.** Simulation results for Weibull $X$ independent of $Z$ from the full cohort analysis (i.e., where all $n$ observations had uncensored $X$) and conditional mean imputation (CMI) approaches.

  - [Script (Run Simulations)](Sim-Scripts/Rev-TableS1-WeibullXindepZ.R)
  - [Script (Make Table)](Table-Scripts/Rev-TableS1-WeibullXindepZ.R)
  - [Data (Simulation Results)](Table-Data/data_TableS1_rev.csv)  

## Figures 

**Figure S1.** Interpolating Breslow's estimator $\widehat{S}_0(t)$ between uncensored values with either of the two interpolation methods offered similar bias and efficiency for $\hat{\beta}$ in extrapolated conditional mean imputation. Between uncensored values, $\widehat{S}_0(\cdot)$ was either be carried forward from the last uncensored value or taken as the mean of the uncensored values immediately before and after. The dashed line denotes the true parameter value, $\beta = 0.5$. Extrapolated CMI using either imputation method was successful in all but 44 replications out of 9000; these few replications encountered errors with numerical integration or non-convergence with the Cox model. Data were simulated following Section 3.1.

  - [Script (Run Simulations)](Sim-Scripts/Rev-FigureS1-Interpolation-Methods.R)
  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS1-Interpolation-Methods.R)
  - [Data (Simulation Results)](Figure-Data/rev_data_figureS1.csv)  

**Figure S2.** Illustration of the four extrapolation methods for a step survival function $\widehat{S}(t)$ in simulated data. The shaded area represents values of $t > \widetilde{X}$ (the largest uncensored value), where extrapolation is needed.

  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS2-Illustrate-Extrapolation-Methods.R)

**Figure S3.** Extrapolating Breslow's estimator beyond the largest uncensored value $\widetilde{X}$ to the overall largest value $X_{(n)}$ with any of the three extrapolation methods offered similar bias and efficiency for $\hat{\beta}$ in non-extrapolated conditional mean imputation using the trapezoidal rule. The dashed line denotes the true parameter value, $\beta = 0.5$. Data were simulated following Section 3.1.

  - [Script (Run Simulations)](Sim-Scripts/Rev-FigureS3-Extrapolation-Methods-Trapezoidal-Rule.R)
  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS3-Extrapolation-Methods-Trapezoidal-Rule.R)
  - [Data (Simulation Results)](Figure-Data/rev_data_figureS3.csv) 

**Figure S4.** In our simulations with Weibull $X$, we explored light ($\sim 17\%$), heavy ($\sim 49\%$), and extra heavy ($\sim 78\%$) censoring in $X$ by generating $C$ from an exponential distribution with rates $= 0.23$, $2$, and $10$, respectively.

  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS4-Percent-Censored.R)

**Figure S5.** Comparison of the probability density functions (**A**) and hazard functions (**B**) of the different distributions considered for the censored covariate $X$ in Sections 3.2 and 3.3.

  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS5-DistX.R)

**Figure S6.** Due to the Weibull distribution's skewness, higher censoring rates led to smaller values of $W_{(n)}$ (the maximum of the observed covariate used as the integral's upper bound by the trapezoidal rule). With smaller values of $W_{(n)}$, the trapezoidal rule was cutting off more of the survival function, leading to worse performance (i.e., higher bias) with non-extrapolated conditional mean imputation. **A**, **B**, and **C** are the empirical densities of $W_{(n)}$ when $X$ was generated from a Weibull, log-normal, and gamma distribution, respectively, under light, heavy, or extra heavy censoring.

  - [Script (Make Figure)](Figure-Scripts/Rev-FigureS6-Density-Xtilde.R) 
