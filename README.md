# Reproducibility package for "Do Professionals' Inflation Forecasts Incorporate the Beliefs of Others? A Functional Data Approach"
**Date Assembled:** April 27, 2026

**Journal:** International Journal of Forecasting (IJF)

## Authors and Contact Information
* **Michael P. Clements**, ICMA Centre, Henley Business School, University of Reading, UK.
* **Shixuan Wang** (Corresponding Author), Department of Economics, University of Reading, UK. Email: shixuan.wang@reading.ac.uk

## Computing Environment
* **Operating System:** Windows 11
* **Software Language:** R version 4.3.0, MATLAB R2022a (only for Figure 1)
* **Required R Packages:** `readxl`, `zoo`, `pracma`, `e1071`, `fdadensity`, `fdapace`, `plm`, `ggplot2`, `viridis`, `hrbrthemes`, `ftsa`, `lmtest`, `tidyverse`
* **Hardware Setup:** Intel Core i7, 16GB RAM

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```R
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Repository Structure
The main contents of the repository are the following:

* `Data/`: folder of data files and calculated log quantile densities (LQDs)
* `Plots/`: folder of generated plots as EPS or PDF files
* `Tables/`: folder of generated tables as csv files
* `renv/` and `renv.lock`: Dependency management files for the R environment.
* `Main.R`: R script to create all figures (except for Figure 1) and all tables
* `myFunctions.R`: R script that contains all supporting user-defined functions
* `Plot_Den_QLD.m`: MATLAB script for create Figure 1.
* `README.md`: Readme file that gives all relevant information to run the reproducibility package.

## Instructions
All file paths are relative to the root of the reproducibility package. Please set your working directory accordingly, or open the `.Rproj` file using RStudio. Please follow the instructions below to run different scripts to generate relevant figures and tables.

* `Plot_Den_QLD.m`: to generate `Figure 1` (Note that the densities and LQDs are obtained in the R script `Main.R`, and Matlab is used to generate the 3-D plot of Figure 1)
* `Main.R`: to generate `Figures 2-6` and `Table 1-19`

## Data availability and provenance

### Survey of Professional Forecasters
The forecast data is from the US Survey of Professional Forecasters (SPF). The individual-level data are made available on the website of the Federal Reserve Bank at Philadelphia (https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters). It is important that the respondent to a given survey will know the responses made to the survey in the previous quarter, as well as the fact that individuals can be tracked through time even though they remain anonymous.

We use the individual expectations from 1981:Q3 - 2022:Q2 (inclusive):

* **CPI Inflation Rate (CPI) Headline**: annualized percentage points, seasonally adjusted, based on quarterly average index level,
* **Civilian Unemployment Rate (UNEMP)**: percentage points, seasonally adjusted, quarterly average,

which are located at `Data/Individual_CPI_2022Q2.xlsx` and `Individual_UNEMP_2022Q2.xlsx`, respectively. 

Additionally, long-term inflation forecasts from 1981:Q3 to 1991:Q3 from Blue Chip Economic Indicators are located at `Data/Additional-CPIE10.xlsx`, which is also sourced from the SPF website.

### External Variables

The dataset of prepared external varaible  used for this study is located in the `Data/External_Variables`, including:

* **VIX**: The raw quarterly VIX data (1990:Q1 – 2022:Q2), originating from the Federal Reserve Economic Data (FRED) database (https://fred.stlouisfed.org).
* **Recession**: The business cycle expansions and contractions dates used to define "Normal" and "Abnormal" times are sourced from the National Bureau of Economic Research (NBER) (https://www.nber.org/research/data/us-business-cycle-expansions-and-contractions).
* **Inflation**: The CPI inflation rate headline (CPIAUCSL), originating from the Federal Reserve Economic Data (FRED) database (https://fred.stlouisfed.org).
* **Oil**: The data for crude oil prices (WTISPLC), originating from the Federal Reserve Economic Data (FRED) database (https://fred.stlouisfed.org).
* **LabourCost**: The data for unit labor costs (ULCNFB), originating from the Federal Reserve Economic Data (FRED) database (https://fred.stlouisfed.org). 
