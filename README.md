# It’s integral: Replacing the trapezoidal rule to remove bias and correctly impute censored covariates with their conditional means

This repository contains R code and simulation data to reproduce results from the manuscript by Lotspeich and Garcia (2022+).

## `Sim-Scripts/`: Scripts to reproduce simulations 

`R` scripts to reproduce all simulated results can be found in the `Sim-Scripts/` subdirectory. 

  -  `Table1.R`: 
  -  `Table2.R`:
  -  `TableS1.R`:
  -  `FigureS1.R`: 
  -  `FigureS2.R`:
  -  `FigureS3.R`:
  -  `FigureS4.R`:
  -  `FigureS5.R`:
  -  `FigureS6.R`:

## `Table-Scripts/`: Scripts to create tables

`R` scripts to reproduce the tables of simulated results can be found in the `Table-Scripts/` subdirectory. 

  -  `Table1.R`: 
  -  `Table2.R`:
  -  `TableS1.R`:

## `Table-Data/`: Results from simulations to be displayed in tables 

Many simulations were run in parallel for more efficient computation. Individual `.csv` files from each "batch" of simulations (i.e., 50 replicates with a particular random seed) can be found in the named subdirectories for each table (e.g., `Table2/`). To build the tables above, these individual files have been combined (i.e., stacked using `rbind()`) into a single file per table, named like `data_Table#.csv`. 

  -  `data_Table1.csv`: 
  -  `data_Table2.csv`:
  -  `data_TableS1.csv`:

## `Figure-Scripts/`: Scripts to create figures

`R` scripts to reproduce the figures of simulated results can be found in the `Figure-Scripts/` subdirectory. 

  -  `FigureS1.R`: 
  -  `FigureS2.R`:
  -  `FigureS3.R`:
  -  `FigureS4.R`:
  -  `FigureS5.R`:
  -  `FigureS6.R`:

## `Figure-Data/`: Results from simulations to be displayed in figures

Many simulations were run in parallel for more efficient computation. Individual `.csv` files from each "batch" of simulations (i.e., 50 replicates with a particular random seed) can be found in the named subdirectories for each table (e.g., `FigureS2/`). To build the tables above, these individual files have been combined (i.e., stacked using `rbind()`) into a single file per figure, named like `data_Figure#.csv`. 

  -  `data_FigureS1.csv`: 
  -  `data_FigureS2.csv`:
  -  `data_FigureS3.csv`:
  -  `data_FigureS4.csv`:
  -  `data_FigureS5.csv`:
  -  `data_FigureS6.csv`:
