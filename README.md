# DOCM_validation
# Diabetes, Obesity, and Cardiovascular Disease Microsimulation (DOC-M) Model
## Project Overview
The DOCM model forecasts US trends in the prevalence of obesity, diabetes, and cardiovascular disease, as well as all-cause mortality, quality-adjusted life years, and health care costs. Results are stratified by race/ethnicity in order to highlight health disparities. 

## Getting Started
### Required R Packages
The following packages are required to run the code. You can install packages in R by running `install.packages('package_name')`.
* survey
* svMisc
* psych
* gdata
* dplyr
* plyr
* data.table
* foreach
* doParallel
* abind

### Model Settings
The 01_DOCM_model.R script is the main script for the model. Under "0.5 Creating Working Directory," you can set your working directory. You can also specify a seed, the number of simulations/probabilistic samples (`n.sim`), the number of cycles (`n.cycle`), the number of individuals to sample (`n.sample`), and the intervention type. These can be set via command line arguments (in that order), or manually under "0.6 Manually set modeling choices or read command-line arguments." For forecasting trends without any intervention, use the "No_Policy" setting. The model can be adapted for predicting policy impacts.

## Contact
David D. Kim, PhD
Assistant Professor of Medicine
Biological Sciences Division and the College
University of Chicago
ddk@uchicago.edu
