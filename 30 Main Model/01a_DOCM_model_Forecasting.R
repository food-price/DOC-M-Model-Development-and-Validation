#Obesity, Diabetes, and Cardiovascular Disease Microsimulation (ODC-M) Model 
#The FOOD-PRICE project (https://food-price.org/)
#Primary Author: David D. Kim, Tufts Medical Center
#Contact: DKim3@tuftsmedicalcenter.org 
#Main Purpose: 
#1) To develop a microsimulation model to link individual risk factors with cardiometabolic diseases, including diabets and CVD
#2) To generate estimations of dietary intakes and changes in risk factors and disease risk from dietary changes

#Date: May 21, 2020

#Updates: Fixing the mortality beyond age 101 so that the model runs really long enough. 
#Updates (12/20/2019): Adding DM adjustment factors for non-whites (2.4 for NHBW & HW and 1.5 for NHBM & HM) [Source: Brancati et al, JAMA, 2000 & Narayan et al, JAMA, 2003]
#Updates (01/16/2020): Complete the full CEA Model with the dietary and policy modules
#Updates (02/12/2020): Full probablistic models
#Updates (03/12/2020): Improving Efficiency of models
#Updates (04/10/2020): Summarizing Outputs
#Updates (05/21/2020): Modified for use on HPC: seed now read from CL args, removed Windows-specific components
#                       Improved file I/O speed by replacing read.csv/write.csv with data.table::fread/fwrite
#Updates (07/16/2020): Extract race-stratified outputs

#################
# 0.Preparation #
#################

# 0.1 remove all data from memory
rm(list = ls())
#memory.size(max=T) # this is Windows-specific and will not work on the HPC running Linux

# 0.2 Start the clock
ptm <- proc.time()

# 0.3 Install/Read R Packages
#install.packages("abind")
library(survey)
library(svMisc)
library(psych)
library(gdata)
library(dplyr)
library(plyr)
library(data.table)   # needed for fread/fwrite
library(foreach)
library(doParallel)
library(abind)

# 0.4 Create Functions 

# 0.4-1 To estimate gamma/beta parameters based on Mean and SD
estGammaParams <- function(mu, var) {
  beta <- var / mu
  alpha <- mu / beta
  return(params = list(alpha = alpha, beta = beta))
}

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

calc_nsims_rgamma <- function(n.sim, mu, se) {
  gamma_par <- estGammaParams(mu, se^2)
  return(rgamma(n.sim, gamma_par$alpha, 1/gamma_par$beta))
}

calc_nsims_rbeta <- function(n.sim, mu, se) {
  beta_par <- estBetaParams(mu, se^2)
  return(rbeta(n=n.sim, beta_par$alpha, beta_par$beta))
}

# 0.4-2 Converting multi-year risk (e.g., 10-year ASCVD risk) to annual probabilities
Multi_yr_Risk_to_annual_prob <- function(time, risk) {
  annual_rate <- -(1/time)*log(1-risk)
  annual_prob <- 1 - exp(-annual_rate)
}

# 0.5 Creating Working Directory 
setwd("C:/Users/dkim23/Box/Working Manuscripts/2021/Model Development - ODC-M/ODC-M Model 2022 06 21")
#setwd("/cluster/tufts/kimlab/lwang18/ODC-M_Validation")

#setwd("C:\\Users\\lwang18\\Documents\\GitHub\\ODC-M_Validation")

# Source other scriptsqa
source("30 Main Model/1a - Diabetes_risk_prediction_FHS (categorical points).R")
source("30 Main Model/1b - ASCVD_risk_calculator (calibrated_v2).R")
source("30 Main Model/1c - FHS Subsequent_CVD_risk_calculator.R")
source("30 Main Model/2a - HrQOL estimator for US general population.R")
source("30 Main Model/2b - HCE estimator for US general population.R")

# 0.6 Manually set modeling choices or read command-line arguments
args <- commandArgs(trailingOnly = TRUE)  # get all arguments after script.R
# If no arguments are read from command line, set modeling choices manually. Otherwise, read modeling choices from command line.
if (length(args) == 0) {
  seed <- 1234
  n.sim <- 1000 #Number of probabilistic samplings
  n.cycle <- 15 #Number of Cycle Length: How long does the model run (i.e., analytic time horizon)
  n.sample <- "ALL" #Number of individuals to be selected from the full sample; if full sample, enter "ALL"
  intervention <- "No_Policy" #SELECT between "Policy" or "No_Policy"
} else {
  # expecting 4 arguments
  if (length(args) != 4) {
    stop("ERROR: Invalid command line arguments", call. = FALSE)
  }
  seed <- as.numeric(args[1]) # extracting first arg as seed and attempt to cast to number
  n.sim <- as.numeric(args[2])
  n.cycle <- as.numeric(args[3])
  intervention <- args[4]
  # Assume entire sample
  n.sample = "ALL"

  # check that modeling choices were set
  if (is.na(seed)) {
    stop("ERROR: missing seed", call. = FALSE)
  }
  if (is.na(n.sim)) {
    stop("ERROR: missing n.sim", call. = FALSE)
  }
  if (is.na(n.cycle)) {
    stop("ERROR: missing n.cycle", call. = FALSE)
  }
  if (!intervention %in% c("No_Policy", "Policy")) {
    stop("ERROR: invalid intervention", call. = FALSE)
  }
}

# 0.7 Complete Modeling Choices
beta_cost <- beta_QALY <- 0.03 #Annual discounting rate of costs and QALYs
set.seed(seed)

################################################
# 1 Importing Necessary Data Inputs + Cleaning #
################################################
print('Importing data')
# 1.1 Read in master input file
NHANES<- fread("00 Input Data/NHANES0102_10_Imp_New.csv", stringsAsFactors = TRUE, data.table = FALSE)

# 1.2 Select only necessary variables from master input file
variables <- c("seqn", "wtint2yr", "wtmec2yr", "sdmvpsu", "sdmvstra", 
               "Age", "Female", "Race", "CVD_history", "Diabetes",
               "Total_Chol","HDL", "SBP", "DBP", "HPT_Txt", "Smoking",
               "DM_family", "BMI", "Glucose", "Trig",
               "ssb","added_sugar", "sodium", "kcal", "sfat")


NHANES <- NHANES[variables]

# 1.3 Select the target starting population to model (Age 25-79) + remove observations with missing 
NHANES_age25_79 <- na.omit(subset(NHANES, Age >= 25 & Age < 80)) 
NHANES_age25_79$Subject_ID <- c(1:nrow(NHANES_age25_79))

# 1.4 Define the analytic dataset to carry it forward
if (n.sample == "ALL"){
  data_for_analysis <- NHANES_age25_79
} else {
  random_sample <- sample(1:nrow(NHANES_age25_79),n.sample,replace=F)
  data_for_analysis <- NHANES_age25_79[random_sample,]
}
data_for_analysis <- data_for_analysis[order(data_for_analysis$Subject_ID), ]
data_for_analysis$SEQN<-data_for_analysis$seqn
# 1.5 Important other input data   

# 1.5-1 Time-varying policy-effect size
policy_effect <- fread("00 Input Data/Policy_effect.csv", stringsAsFactors = TRUE, data.table = FALSE)

if (intervention == "Policy") {
  policy_effect_ssb <- t(mapply(calc_nsims_rbeta, n.sim, policy_effect$ssb, policy_effect$ssb.se))
  policy_effect_added_sugar <- t(mapply(calc_nsims_rbeta, n.sim, policy_effect$added_sugar, policy_effect$added_sugar.se))
  policy_effect_sodium <- t(mapply(calc_nsims_rbeta, n.sim, policy_effect$sodium, policy_effect$sodium.se))
  policy_effect_kcal <- t(mapply(calc_nsims_rbeta, n.sim, policy_effect$kcal, policy_effect$kcal.se))
  policy_effect_sfat <- t(mapply(calc_nsims_rbeta, n.sim, policy_effect$sfat, policy_effect$sfat.se))
} else {
  policy_effect_ssb <- policy_effect_added_sugar <- policy_effect_sodium <- policy_effect_kcal <- policy_effect_sfat <- matrix(0, nrow=nrow(policy_effect), ncol=n.sim)
}

# 1.5-2 Age-specific relative risk estimates between dietary intake and disease outcomes  
RR_diet_disease <- fread("00 Input Data/final_diet_RRs_MainAnalysis v11.csv", stringsAsFactors = TRUE, data.table = FALSE)
RR_diet_disease <- subset(RR_diet_disease, minage>=35)


# 1.5-3 Health-state/Event specific mortality data  

#1.5-3a: Re-calibration
Adj_Non_CVD_DM_mortality_Male <- Adj_Non_CVD_DM_mortality_Female <- 1.034
Adj_Non_CVD_DM_mortality_NHWM <- Adj_Non_CVD_DM_mortality_NHWF <- 0.993
Adj_Non_CVD_DM_mortality_NHBM <- Adj_Non_CVD_DM_mortality_NHBF <- 1.121*1.087628
Adj_Non_CVD_DM_mortality_HM <- Adj_Non_CVD_DM_mortality_HF <- 1.445*1.049

Adj_CVD_mortality_Male <- Adj_CVD_mortality_Female <- 0.647
Adj_CVD_mortality_NHWM <- Adj_CVD_mortality_NHWF <- 0.8
Adj_CVD_mortality_NHBM <- Adj_CVD_mortality_NHBF <- 0.705
Adj_CVD_mortality_HM <- Adj_CVD_mortality_HF <- 1.090*1.103*1.169

Adj_DM_mortality_Male <- Adj_DM_mortality_Female <- 0.715
Adj_DM_mortality_NHWM <- Adj_DM_mortality_NHWF <- 0.469
Adj_DM_mortality_NHBM <- Adj_DM_mortality_NHBF <- 1.237
Adj_DM_mortality_HM <- Adj_DM_mortality_HF <- 1.438*1.122*1.437

Non_CVD_DM_mortality <- fread("00 Input Data/Non_DM_IHD_Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ corrected).csv", stringsAsFactors = TRUE, data.table = FALSE)
Non_CVD_DM_mortality_Male <- Adj_Non_CVD_DM_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Male, Non_CVD_DM_mortality$Male_SE))
Non_CVD_DM_mortality_Female <- Adj_Non_CVD_DM_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Female, Non_CVD_DM_mortality$Female_SE))
Non_CVD_DM_mortality_NHWM <- Adj_Non_CVD_DM_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWM, Non_CVD_DM_mortality$NHWM_SE))
Non_CVD_DM_mortality_NHWF <- Adj_Non_CVD_DM_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWF, Non_CVD_DM_mortality$NHWF_SE))
Non_CVD_DM_mortality_NHBM <- Adj_Non_CVD_DM_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBM, Non_CVD_DM_mortality$NHBM_SE))
Non_CVD_DM_mortality_NHBF <- Adj_Non_CVD_DM_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBF, Non_CVD_DM_mortality$NHBF_SE))
Non_CVD_DM_mortality_HM <- Adj_Non_CVD_DM_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HM, Non_CVD_DM_mortality$HM_SE))
Non_CVD_DM_mortality_HF <- Adj_Non_CVD_DM_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HF, Non_CVD_DM_mortality$HF_SE))

stroke_mortality <- fread("00 Input Data/Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
stroke_mortality_Male <- Adj_CVD_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Male, stroke_mortality$Male_SE))
stroke_mortality_Female <- Adj_CVD_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Female, stroke_mortality$Female_SE))
stroke_mortality_NHWM <- Adj_CVD_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWM, stroke_mortality$NHWM_SE))
stroke_mortality_NHWF <- Adj_CVD_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWF, stroke_mortality$NHWF_SE))
stroke_mortality_NHBM <- Adj_CVD_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBM, stroke_mortality$NHBM_SE))
stroke_mortality_NHBF <- Adj_CVD_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBF, stroke_mortality$NHBF_SE))
stroke_mortality_HM <- Adj_CVD_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HM, stroke_mortality$HM_SE))
stroke_mortality_HF <- Adj_CVD_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HF, stroke_mortality$HF_SE))
  
CHD_mortality <- fread("00 Input Data/IHD_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
CHD_mortality_Male <- Adj_CVD_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Male, CHD_mortality$Male_SE))
CHD_mortality_Female <- Adj_CVD_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Female, CHD_mortality$Female_SE))
CHD_mortality_NHWM <- Adj_CVD_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWM, CHD_mortality$NHWM_SE))
CHD_mortality_NHWF <- Adj_CVD_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWF, CHD_mortality$NHWF_SE))
CHD_mortality_NHBM <- Adj_CVD_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBM, CHD_mortality$NHBM_SE))
CHD_mortality_NHBF <- Adj_CVD_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBF, CHD_mortality$NHBF_SE))
CHD_mortality_HM <- Adj_CVD_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HM, CHD_mortality$HM_SE))
CHD_mortality_HF <- Adj_CVD_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HF, CHD_mortality$HF_SE))

DM_mortality <- fread("00 Input Data/DM_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
DM_mortality_Male <- Adj_DM_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Male, DM_mortality$Male_SE))
DM_mortality_Female <- Adj_DM_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Female, DM_mortality$Female_SE))
DM_mortality_NHWM <- Adj_DM_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWM, DM_mortality$NHWM_SE))
DM_mortality_NHWF <- Adj_DM_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWF, DM_mortality$NHWF_SE))
DM_mortality_NHBM <- Adj_DM_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBM, DM_mortality$NHBM_SE))
DM_mortality_NHBF <- Adj_DM_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBF, DM_mortality$NHBF_SE))
DM_mortality_HM <- Adj_DM_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HM, DM_mortality$HM_SE))
DM_mortality_HF <- Adj_DM_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HF, DM_mortality$HF_SE))


# 1.5-4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 
risk_factor_trend <- fread("00 Input Data/Risk_Factors_Trends_Summary (Final)_calibrated.csv", stringsAsFactors = TRUE, data.table = FALSE)
risk_factor_trend_sim <- matrix(NA, nrow=nrow(risk_factor_trend), ncol=n.sim)
colnames(risk_factor_trend_sim) <- paste("simulation_", 1:n.sim, sep = " ")

risk_factor_trend_sim <- t(mapply(rnorm, n=n.sim, risk_factor_trend$APC, risk_factor_trend$APC_SE))

for (i in 1:nrow(risk_factor_trend)) {
  risk_factor_trend_sim[i,] <- as.matrix(rnorm(n.sim, risk_factor_trend$APC[i], risk_factor_trend$APC_SE[i]))
}

risk_factor_trend_sim <- cbind(risk_factor_trend, risk_factor_trend_sim)

# 1.5-5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke)
Prop_CHD <- fread("00 Input Data/Prop_CHD.csv", stringsAsFactors = FALSE, data.table = FALSE)

#1.5-6 Health Care Expenditure Model
HCE_parameter_sim <- matrix(NA, nrow=nrow(HCE_parameters), ncol=n.sim)
rownames(HCE_parameter_sim) <- rownames(HCE_parameters)
colnames(HCE_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HCE_parameters)) {
  HCE_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HCE_parameters$Beta[i], HCE_parameters$SE[i]))
}

#1.5-7 HrQOL Prediction
HRQOL_parameter_sim <- matrix(NA, nrow=nrow(HRQOL_parameters), ncol=n.sim)
rownames(HRQOL_parameter_sim) <- rownames(HRQOL_parameters)
colnames(HRQOL_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HRQOL_parameters)) {
  HRQOL_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HRQOL_parameters$Beta[i], HRQOL_parameters$SE[i]))
}

# 1.6 DATA CLEANING

data_for_analysis$initial_H[data_for_analysis$CVD_history == 0 & data_for_analysis$Diabetes == 0] <- "No CVD, No Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 0 & data_for_analysis$Diabetes == 1] <- "No CVD, With Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 1 & data_for_analysis$Diabetes == 0] <- "CVD History, No Diabetes"
data_for_analysis$initial_H[data_for_analysis$CVD_history == 1 & data_for_analysis$Diabetes == 1] <- "CVD History, With Diabetes"

data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 1] <- "NHWM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 1] <- "NHWF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 2] <- "NHBM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 2] <- "NHBF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 3] <- "HM"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 3] <- "HF"
data_for_analysis$DEMO[data_for_analysis$Female == 0 & data_for_analysis$Race == 4] <- "Male"
data_for_analysis$DEMO[data_for_analysis$Female == 1 & data_for_analysis$Race == 4] <- "Female"

data_for_analysis$Age_cycle <- NA

data_for_analysis$DM_parent <- data_for_analysis$DM_family

data_for_analysis$BMI_cat[data_for_analysis$BMI < 18.5] <- "Underweight"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 18.5 & data_for_analysis$BMI < 25] <- "Normal"
data_for_analysis$BMI_cat[data_for_analysis$BMI >=25 & data_for_analysis$BMI < 30] <- "Overweight"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 30 & data_for_analysis$BMI < 35] <- "Obesity I"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 35 & data_for_analysis$BMI < 40] <- "Obesity II"
data_for_analysis$BMI_cat[data_for_analysis$BMI >= 40] <- "Obesity III"

data_for_analysis$Obesity[data_for_analysis$BMI < 30] <- 0
data_for_analysis$Obesity[data_for_analysis$BMI >= 30] <- 1

# 1.7 DM risk adjustment for non-whites
data_for_analysis$risk_adjustment.DM <- ifelse(data_for_analysis$DEMO %in% c("NHWM", "NHWF", "Female", "Male"), 1.0,
                                               ifelse(data_for_analysis$DEMO %in% c("NHBM", "HM"), 1.5, 2.4))

# 1.7 Creating dupilicate (counter-factual) observations to predict the effect under policy vs. no policy
data_for_analysis$source <- intervention

# 1.8 Diet-disease data inputs 

# 1.8.1 Indirect effect on disease
sugar_bmi_low <- subset(RR_diet_disease, outcome == 'BMI_low' & riskfactor == 'sugar', select = c(minage, logRR.perchange.exp, se.perchange.exp))
sugar_bmi_high <- subset(RR_diet_disease, outcome == 'BMI_high' & riskfactor == 'sugar', select = c(minage, logRR.perchange.exp, se.perchange.exp))

ssb_bmi_low <- subset(RR_diet_disease, outcome == 'BMI_low' & riskfactor == 'ssb', select = c(minage, logRR.perchange.exp, se.perchange.exp))
ssb_bmi_high <- subset(RR_diet_disease, outcome == 'BMI_high' & riskfactor == 'ssb', select = c(minage, logRR.perchange.exp, se.perchange.exp))

na_sbp_main <- subset(RR_diet_disease, outcome == 'SBP_main' & riskfactor == 'sodium', select = c(minage, logRR.perchange.exp, se.perchange.exp))
na_sbp_black <- subset(RR_diet_disease, outcome == 'SBP_black' & riskfactor == 'sodium', select = c(minage, logRR.perchange.exp, se.perchange.exp))
na_sbp_hpt <- subset(RR_diet_disease, outcome == 'SBP_hpt' & riskfactor == 'sodium', select = c(minage, logRR.perchange.exp, se.perchange.exp))

# 1.8.2 Direct effect on disease

RR_sugar_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'sugar', select = c(minage, logRR.perchange, se.perchange))
RR_sugar_t2d <- subset(RR_diet_disease, outcome == 'DIAB' & riskfactor == 'sugar', select = c(minage, logRR.perchange, se.perchange))

RR_ssb_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'ssb', select = c(minage, logRR.perchange, se.perchange))
RR_ssb_t2d <- subset(RR_diet_disease, outcome == 'DIAB' & riskfactor == 'ssb', select = c(minage, logRR.perchange, se.perchange))

RR_sfat_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'sfat', select = c(minage, logRR.perchange, se.perchange))

RR_bmi_ihd <- subset(RR_diet_disease, outcome == 'IHD' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))
RR_bmi_tstk <- subset(RR_diet_disease, outcome == 'TSTK' & riskfactor == 'bmi', select = c(minage, logRR.perchange, se.perchange))

# Simulate direct effect of changes in dietary factor and BMI on disease 
# SSB -> CHD & T2d; BMI -> CHD & Stroke; Sat.Fat -> CHD; Sugar -> CHD & T2D
# Set number of diet/riskfactor-disease pairs
# Random draws of dietary effect on risk factors
n.beta <- 7

# Prepare matrix
random_logrr <- matrix(NA, nrow=nrow(RR_ssb_ihd), ncol=n.sim)
rownames(random_logrr) <- c("35-44 y", "45-54 y", "55-64 y", "65-74 y","75+ y")
colnames(random_logrr) <- paste("simulation_", 1:n.sim, sep = " ")

random_logrr_SSB_IHD <- random_logrr_SSB_T2D <- random_logrr_BMI_IHD <- random_logrr_BMI_TSTK <-
  random_logrr_SFAT_IHD <- random_logrr_Sugar_IHD <- random_logrr_Sugar_T2D <- random_logrr

for (g in 1:nrow(random_logrr)) {
  random_logrr_SSB_IHD[g,] <- rnorm(n=n.sim, mean = RR_ssb_ihd$logRR.perchange[g], sd = RR_ssb_ihd$se.perchange[g])
  random_logrr_SSB_T2D[g,] <- rnorm(n=n.sim, mean = RR_ssb_t2d$logRR.perchange[g], sd = RR_ssb_t2d$se.perchange[g])
  random_logrr_BMI_IHD[g,] <- rnorm(n=n.sim, mean = RR_bmi_ihd$logRR.perchange[g], sd = RR_bmi_ihd$se.perchange[g])
  random_logrr_BMI_TSTK[g,] <- rnorm(n=n.sim, mean = RR_bmi_tstk$logRR.perchange[g], sd = RR_bmi_tstk$se.perchange[g])
  random_logrr_SFAT_IHD[g,] <- rnorm(n=n.sim, mean = RR_sfat_ihd$logRR.perchange[g], sd = RR_sfat_ihd$se.perchange[g])
  random_logrr_Sugar_IHD[g,] <- rnorm(n=n.sim, mean = RR_sugar_ihd$logRR.perchange[g], sd = RR_sugar_ihd$se.perchange[g])
  random_logrr_Sugar_T2D[g,] <- rnorm(n=n.sim, mean = RR_sugar_t2d$logRR.perchange[g], sd = RR_sugar_t2d$se.perchange[g])
}



# Prepare matrix
random_beta <- matrix(NA, nrow=nrow(ssb_bmi_low), ncol=n.sim)
rownames(random_beta) <- c("35-44 y", "45-54 y", "55-64 y", "65-74 y","75+ y")
colnames(random_beta) <- paste("simulation_", 1:n.sim, sep = " ")

random_beta_Sugar_BMI_Low <- random_beta_Sugar_BMI_high <- random_beta_SSB_BMI_Low <- random_beta_SSB_BMI_high <-
  random_beta_Na_SBP_main <-  random_beta_Na_SBP_black <- random_beta_Na_SBP_hpt <- random_beta

# Populate matrix
for (g in 1:nrow(ssb_bmi_low)) {
  random_beta_Sugar_BMI_Low[g,] <- rnorm(n=n.sim, mean = sugar_bmi_low$logRR.perchange[g], sd = sugar_bmi_low$logRR.perchange[g])
  random_beta_Sugar_BMI_high[g,] <- rnorm(n=n.sim, mean = sugar_bmi_high$logRR.perchange[g], sd = sugar_bmi_high$se.perchange[g])
  random_beta_SSB_BMI_Low[g,] <- rnorm(n=n.sim, mean = ssb_bmi_low$logRR.perchange[g], sd = ssb_bmi_low$se.perchange[g])
  random_beta_SSB_BMI_high[g,] <- rnorm(n=n.sim, mean = ssb_bmi_high$logRR.perchange[g], sd = ssb_bmi_high$se.perchange[g])
  random_beta_Na_SBP_main[g,] <- rnorm(n=n.sim, mean = na_sbp_main$logRR.perchange[g], sd = na_sbp_main$se.perchange[g])
  random_beta_Na_SBP_black[g,] <- rnorm(n=n.sim, mean = na_sbp_black$logRR.perchange[g], sd = na_sbp_black$se.perchange[g])
  random_beta_Na_SBP_hpt[g,] <- rnorm(n=n.sim, mean = na_sbp_hpt$logRR.perchange[g], sd = na_sbp_hpt$se.perchange[g])
}

#############################################
# 2 Describing the baseline characteristics #
#############################################
  #OMITTED - SEE THE PREVIOUS VERSION

#################################################################################################################################
# 3 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
#################################################################################################################################
print('part 3')
p3_start <- proc.time()
# 3.1 FHS 8-year Diabetes Risk Prediction (With BMI)
variable_for_raw.input <- names(data_for_analysis)
raw.input.data <- data_for_analysis
data_for_analysis$DM_risk_8yr <- calc_DM_risk(raw.input.data)
data_for_analysis$DM_prob <- Multi_yr_Risk_to_annual_prob(time=8, risk=data_for_analysis$DM_risk_8yr)

# 3.2 ACC/AHA ASCVD 10-year Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
data_for_analysis$CVD_prob <- Multi_yr_Risk_to_annual_prob(time=10, risk=data_for_analysis$ASCVD_Risk_10yr)

# 3.3 FHS 2-year CVD recurrent Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
data_for_analysis$CVD_recurrent_prob <- Multi_yr_Risk_to_annual_prob(time=2, risk=data_for_analysis$CVD_Recurrent_risk_2yr)

# 3.4 Individual HrQOL prediction
raw.input.data <- data_for_analysis
data_for_analysis$HRQOL_scores <- calc_HRQOL(raw.input.data, HRQOL_parameters[,1])

# 3.5 Individual HCE prediction
raw.input.data <- data_for_analysis
data_for_analysis$HCE_predict <- calc_HCE(raw.input.data, HCE_parameters[,1])

p3_end <- proc.time()
p3_end - p3_start

###############################################################################################################################################################
# 4 Microsimulation Model - Predicting individual-levle changes in diet, risk factors, disease outcomes, health-related quality of life, and healthcare costs # #
###############################################################################################################################################################
print('part 4')
# 4.1 Modeling Input
initial_H <- data_for_analysis$initial_H #Vector of Initial states for individuals
n.individual <- nrow(data_for_analysis) # number of individuals in the model
name.health.state <- c("No CVD, No Diabetes", "No CVD, With Diabetes", "First Stroke", "First CHD w/o RVSC", "First CHD with RVSC",
                       "CVD History, No Diabetes", "CVD History, With Diabetes", "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC", 
                       "DM_Death", "Stroke_Death", "CHD_Death", "Non_DM_Non_CVD_Death") 
n.health.state <- length(name.health.state) # number of health state that individuals can transition over time
name.death.states <- c("DM_Death", "Stroke_Death", "CHD_Death", "Non_DM_Non_CVD_Death")

#Health care costs associated with clinical ASCVD events  
#Source: https://www.cms.gov/Research-Statistics-Data-and-Systems/Statistics-Trends-and-Reports/Medicare-Provider-Charge-Data/Inpatient2017.html
c_CHD <- 10034.34 # Medicare Inpatient Prospective Payment Hospitals across six DRGs for MI [280-285] (Average total payments)
c_CHD_sim <- calc_nsims_rgamma(n.sim, c_CHD, c_CHD*0.2)

c_stroke <- 15994.49 # Medicare Inpatient Prospective Payment Hospitals across three DRGs for stroke [061-063] (Average total payments)
c_stroke_sim <- calc_nsims_rgamma(n.sim, c_stroke, c_stroke*0.2)

#Disutility associated with major clinical events
#Source: Davies et al. Health and Quality of Life Outcomes (2015) 13:90.DOI 10.1186/s12955-015-0266-9

u_CHD <- -0.055 #Average disutility among MI (-0.06) and unstable angina (-0.05)
u_CHD_sim <- -calc_nsims_rbeta(n.sim, -u_CHD, -u_CHD*0.2)

u_stroke <- -0.3 
u_stroke_sim <- -calc_nsims_rbeta(n.sim, -u_stroke, -u_stroke*0.2)

#RVSC-related Input
p.RVSC <- 0.67 #Probabilities of receiving RVSC after CVD event, Source: Benjamin et al, Circulation, 2018, Table 13-1[# of Stroke Discharge] & 18-1[# of CHD Hospital Discharge] & 24-2[# of inpatient RVSC]
p.RVSC_sim <- calc_nsims_rbeta(n.sim, p.RVSC, p.RVSC*0.2)

prop_PCI <- 0.711 #Proportions of PCI (including PCI with stent) among all RVSC (PI & CABG) Source: Heart Disease and Stroke Statistics (Benjamin et al., Circulation, 2018) Table 24-2
prop_PCI_sim <- calc_nsims_rbeta(n.sim, prop_PCI, prop_PCI*0.2)

c_CABG <- 44538.01 # Medicare Inpatient Prospective Payment Hospitals across six DRGs for CABG[231-236] (Average total payments)
c_CABG_sim <- calc_nsims_rgamma(n.sim, c_CABG, c_CABG*0.2)

c_PCI <- 18476.82 # Medicare Inpatient Prospective Payment Hospitals across nine DRGs for PCI[034-036; 246-251] (Average total payments)
c_PCI_sim <- calc_nsims_rgamma(n.sim, c_PCI, c_PCI*0.2)

c_RVSC <- prop_PCI*c_PCI + (1-prop_PCI)*c_CABG
c_RVSC_sim <- prop_PCI_sim*c_PCI_sim + (1-prop_PCI_sim)*c_CABG_sim


mortality_CABG <- 0.0178 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_CABG_sim <- calc_nsims_rbeta(n.sim, mortality_CABG, mortality_CABG*0.2)

mortality_PCI <- 0.0207 #In-hospital death rate, Source: Benjamin et al, Circulation, 2018, Table 24-1
mortality_PCI_sim <- calc_nsims_rbeta(n.sim, mortality_PCI, mortality_PCI*0.2)

p.death.RVSC <- prop_PCI*mortality_PCI + (1-prop_PCI)*mortality_CABG
p.death.RVSC_sim <- prop_PCI_sim*mortality_PCI_sim + (1-prop_PCI_sim)*mortality_CABG_sim

#########################################################################
# 5 Main Model - PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING #
#########################################################################
print('running model')

update_risk_factor <- function(sim_out_t, risk_factor, predictor, s) {
  # Updates risk factor values from time t to time t+1
  # References risk_factor_trend_sim for APC trends in Total_Chol, HDL, SBP/DBP, BMI, Trig, and Glucose
  risk_factor_vals <- case_when(
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHWM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHWM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHWM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHWM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHWF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHWF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHWF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHWF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHBM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHBM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHBM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHBM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "NHBF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "NHBF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "NHBF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "NHBF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "HM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "HM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "HM" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "HM"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "HF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "HF"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "HF" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "HF"),s+5]/100)), 
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "Male" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "Male"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "Male" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "Male"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) < 45 & sim_out_t[,"DEMO"] == "Female" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 1 & risk_factor_trend_sim$DEMO == "Female"),s+5]/100)),
    as.numeric(sim_out_t[,"Age_cycle"]) >= 45 & sim_out_t[,"DEMO"] == "Female" ~ as.numeric(sim_out_t[,risk_factor])*(1+(risk_factor_trend_sim[which(risk_factor_trend_sim$Predictors==predictor & risk_factor_trend_sim$Age_40_59 == 0 & risk_factor_trend_sim$DEMO == "Female"),s+5]/100)),
    TRUE ~ as.numeric(NA))
  
  return(risk_factor_vals)
}

run_sim <- function(s) {
  # Runs single simulation
  # Returns array with variables for each individual at each cycle
#  key_variables <- c("seqn", "Age", "Age_cycle", "DEMO", "Female", "Race", "Total_Chol","HDL", "SBP", "DBP", "BMI", "BMI_cat", "Obesity", "HPT_Txt", "Trig", "Glucose", 
 #                    "Smoking", "Diabetes", "CVD_history", "risk_adjustment.DM", "DM_parent", "ssb", "added_sugar", "sodium", "kcal" , "sfat", "DM_prob", "CVD_prob", "CVD_recurrent_prob", "HRQOL_scores", "HCE_predict")  
  #all_variables <- c(key_variables, "state", "cost_disc", "effect_disc")
  
  
  fixedvar<-c("Female", "Race", "DM_parent",  "HPT_Txt",  "Smoking")
  Input_variables <- c("SEQN","DEMO", "Total_Chol","HDL", "SBP", "DBP", "BMI", "BMI_cat", "Obesity", "Trig", "Glucose", 
                       "Diabetes", "CVD_history","DM_prob", "ssb", "added_sugar", "sodium", "kcal" , "sfat", "CVD_prob", "CVD_recurrent_prob")
  
  Output_variables<- c("HRQOL_scores", "Obesity",  "Diabetes", "CVD_history", "HCE_predict","state", "HCE_disc", "effect_disc")
  
  datain_t<-array(NA, dim=c(n.individual, length(Input_variables)+1),
                  dimnames = list(data_for_analysis$Subject_ID,
                                  c("Age_cycle", Input_variables)))
  datain_t1<-array(NA, dim=c(n.individual, length(Input_variables)+1),
                   dimnames = list(data_for_analysis$Subject_ID,
                                   c("Age_cycle", Input_variables)))
  
  # Initialize array to hold output from single simulation
  sim_out<-array(NA, dim=c(n.individual, length(Output_variables), n.cycle+1),
                 dimnames = list(data_for_analysis$Subject_ID,
                                 Output_variables, paste("cycle", 0:n.cycle, sep = " ")))
  
  # Set initial state (cycle 0)
  datain_t[,Input_variables] <- as.matrix(data_for_analysis[,Input_variables])
  datain_t[,"Age_cycle"]=data_for_analysis[,"Age"]
  sim_out[,c( "HRQOL_scores", "HCE_predict","Diabetes","CVD_history","Obesity"),1]<-as.matrix(data_for_analysis[,c( "HRQOL_scores", "HCE_predict", "Diabetes","CVD_history","Obesity")])
  sim_out[,"state",1] <- initial_H
  sim_out[,"HCE_disc",1] <- data_for_analysis$HCE_predict
  sim_out[,"effect_disc",1] <- data_for_analysis$HRQOL_scores
  datain_t[,"CVD_recurrent_prob"]<-0
  
  
  for (t in 1:n.cycle) {
    
    
    #Non-time varying data inputs: carry it over from the baseline data  
    datain_t1[,c("SEQN","DEMO")]<-datain_t[,c("SEQN","DEMO")]
    #Time-varying data inputs
    datain_t1[,"Age_cycle"] <- data_for_analysis[,"Age"] + t
    datain_t1[,"ssb"] <- data_for_analysis[,"ssb"]*(1-policy_effect_ssb[t,s])
    datain_t1[,"added_sugar"] <- data_for_analysis[,"added_sugar"]*(1-policy_effect_added_sugar[t,s])
    datain_t1[,"sodium"] <- data_for_analysis[,"sodium"]*(1-policy_effect_sodium[t,s])
    datain_t1[,"kcal"] <- data_for_analysis[,"kcal"]*(1-policy_effect_kcal[t,s])
    datain_t1[,"sfat"] <- data_for_analysis[,"sfat"]*(1-policy_effect_sfat[t,s])
    
    # 5.4 Updating Risk Factors over time: Applying secular trends based on age, gender, R/E
    datain_t1[,"BMI"] <- update_risk_factor(datain_t, "BMI", "BMI", s)
    
    # Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
    # Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
    
    datain_t1[,"BMI"] <- ifelse(as.numeric(datain_t1[,"BMI"]) < 12, 12,
                                ifelse(as.numeric(datain_t1[,"BMI"])>70, 70, as.numeric(datain_t1[,"BMI"])))
    
    
    # Simulating Age-specific indirect effect of dietary change on disease through modifying risk factors (sugar & SSB -> BMI; Sodium -> SSB) 
    datain_t1[,"BMI"]  <- ifelse(as.numeric(datain_t1[,"BMI"] ) <25,
                                 as.numeric(datain_t1[,"BMI"] )+random_beta_Sugar_BMI_Low[1,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))+random_beta_SSB_BMI_Low[1,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"])),
                                 as.numeric(datain_t1[,"BMI"] )+random_beta_Sugar_BMI_high[1,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))+random_beta_SSB_BMI_high[1,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"])))
    
    # sim_out[,"BMI_cat",t+1] <- case_when(
    #   as.numeric(datain_t1[,"BMI"]) < 18.5 ~ "Underweight",
    #   as.numeric(datain_t1[,"BMI"]) < 25 ~ "Normal",
    #   as.numeric(datain_t1[,"BMI"]) < 30 ~ "Overweight",
    #   as.numeric(datain_t1[,"BMI"]) < 35 ~ "Obesity I",
    #   as.numeric(datain_t1[,"BMI"]) < 40 ~ "Obesity II",
    #   as.numeric(datain_t1[,"BMI"]) >= 40 ~ "Obesity III"
    # )
    # 
    
    datain_t1[,"Total_Chol"] <- update_risk_factor(datain_t, "Total_Chol", "Total_Chol", s)
    datain_t1[,"HDL"] <- update_risk_factor(datain_t, "HDL", "HDL", s)
 
    datain_t1[,"DBP"] <- update_risk_factor(datain_t, "DBP", "SBP", s)
    datain_t1[,"Trig"] <- update_risk_factor(datain_t, "Trig", "Trig", s)
    datain_t1[,"Glucose"] <- update_risk_factor(datain_t, "Glucose", "Glucose", s)
    
    datain_t1[,"SBP"] <- update_risk_factor(datain_t, "SBP", "SBP", s)
    # Simulate additional SBP changes due to dietary change dependent on BMI status. Take into consideration intake change of sodium. Allows for age-specific main effects of sodium
    datain_t1[,"SBP"] <- case_when(
      as.numeric(datain_t[,"Age_cycle"]) < 45 ~ as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_main[1,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3),
      as.numeric(datain_t[,"Age_cycle"]) < 55 ~ as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_main[2,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3),
      as.numeric(datain_t[,"Age_cycle"]) < 65 ~ as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_main[3,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3),
      as.numeric(datain_t[,"Age_cycle"]) < 75 ~ as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_main[4,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3),
      TRUE ~  as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_main[5,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3)
    )
    
    # Additional SBP effect of sodium intake in Blacks
    datain_t1[,"SBP"] <- ifelse(data_for_analysis$Race == 2, as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_black[1,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3),  as.numeric(datain_t1[,"SBP"]))
    datain_t1[,"SBP"] <- ifelse((as.numeric(datain_t[,"SBP"]) >= 140 | as.numeric(datain_t[,"DBP"]) >= 90), as.numeric(datain_t1[,"SBP"])+random_beta_Na_SBP_hpt[1,s]*((as.numeric(datain_t1[,"sodium"]) - as.numeric(data_for_analysis[,"sodium"]))/2.3), as.numeric(datain_t1[,"SBP"]))
    
    
    
    # SSB-IHD/CHD
    RR_diff_ssb_ihd <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_SSB_IHD[1,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_SSB_IHD[2,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_SSB_IHD[3,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_SSB_IHD[4,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      TRUE ~ exp(random_logrr_SSB_IHD[5,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"])))
    )
    
    # SSB-T2DM
    RR_diff_ssb_t2d <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_SSB_T2D[1,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_SSB_T2D[2,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_SSB_T2D[3,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_SSB_T2D[4,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"]))),
      TRUE ~ exp(random_logrr_SSB_T2D[5,s]*(as.numeric(datain_t1[,"ssb"]) - as.numeric(data_for_analysis[,"ssb"])))
    )
    
    # BMI-IHD/CHD
    RR_diff_bmi_ihd <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_IHD[1,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5), 
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_IHD[2,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_IHD[3,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_IHD[4,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),  
      TRUE ~ exp(random_logrr_BMI_IHD[5,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5)
    )
    
    # BMI-Total stroke
    RR_diff_bmi_tstk <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_BMI_TSTK[1,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5), 
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_BMI_TSTK[2,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_BMI_TSTK[3,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_BMI_TSTK[4,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5),
      TRUE ~ exp(random_logrr_BMI_TSTK[5,s]*(as.numeric(datain_t1[,"BMI"]) - as.numeric(data_for_analysis[,"BMI"]))/5)
    )
    
    # Sat.Fat-IHD/CHD
    RR_diff_sfat_ihd <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_SFAT_IHD[1,s]*(as.numeric(datain_t1[,"sfat"]) - as.numeric(data_for_analysis[,"sfat"]))/5),
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_SFAT_IHD[2,s]*(as.numeric(datain_t1[,"sfat"]) - as.numeric(data_for_analysis[,"sfat"]))/5),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_SFAT_IHD[3,s]*(as.numeric(datain_t1[,"sfat"]) - as.numeric(data_for_analysis[,"sfat"]))/5),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_SFAT_IHD[4,s]*(as.numeric(datain_t1[,"sfat"]) - as.numeric(data_for_analysis[,"sfat"]))/5),
      TRUE ~ exp(random_logrr_SFAT_IHD[5,s]*(as.numeric(datain_t1[,"sfat"]) - as.numeric(data_for_analysis[,"sfat"]))/5)
    )
    
    # Sugar-IHD/CHD
    RR_diff_sugar_ihd <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_Sugar_IHD[1,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_Sugar_IHD[2,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_Sugar_IHD[3,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_Sugar_IHD[4,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      TRUE ~ exp(random_logrr_Sugar_IHD[5,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"])))
    )
    
    # Sugar-T2DM
    RR_diff_sugar_t2d <- case_when(
      datain_t[,"Age_cycle"] < 45 ~ exp(random_logrr_Sugar_T2D[1,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 55 ~ exp(random_logrr_Sugar_T2D[2,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 65 ~ exp(random_logrr_Sugar_T2D[3,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      datain_t[,"Age_cycle"] < 75 ~ exp(random_logrr_Sugar_T2D[4,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"]))),
      TRUE ~ exp(random_logrr_Sugar_T2D[5,s]*(as.numeric(datain_t1[,"added_sugar"]) - as.numeric(data_for_analysis[,"added_sugar"])))
    )
    
    # Combine direct and indirect effects by multiplying RRs
    RR_diff_total_CHD <- RR_diff_bmi_ihd*RR_diff_ssb_ihd*RR_diff_sfat_ihd*RR_diff_sugar_ihd
    RR_diff_total_DM <- RR_diff_ssb_t2d*RR_diff_sugar_t2d
    RR_diff_total_Stroke <- RR_diff_bmi_tstk
    
    #Updating the disease risk
    
    fixedvar<-c("Female", "Race", "DM_parent",  "HPT_Txt",  "Smoking")
    keyinvars<-c("SEQN","Age_cycle","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history")
    
    #Updating the disease risk
    raw.input.data<-matrix(as.numeric(datain_t[,keyinvars]), ncol=length(keyinvars))
    colnames(raw.input.data) <- c("seqn","Age","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history")
    raw.input.data <- as.data.frame(cbind(raw.input.data, data_for_analysis[,fixedvar]))
    
    if (t%/%2 > 0 & t%%2 == 0)  {
      CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
      datain_t1[,"CVD_recurrent_prob"] <- Multi_yr_Risk_to_annual_prob(time=2, risk=CVD_Recurrent_risk_2yr)
    } else {
      datain_t1[,"CVD_recurrent_prob"] <-datain_t[,"CVD_recurrent_prob"]
    }
    
    if (t%/%8 > 0 & t%%8 == 0)  {
      DM_risk_8yr <- calc_DM_risk(raw.input.data)
      datain_t1[,"DM_prob"] <- Multi_yr_Risk_to_annual_prob(time=8, risk=DM_risk_8yr)
    } else {
      datain_t1[,"DM_prob"] <-datain_t[,"DM_prob"]
    }
    
    if (t%/%10 > 0 & t%%10 == 0) {
      ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
      datain_t1[,"CVD_prob"] <- Multi_yr_Risk_to_annual_prob(time=10, risk=ASCVD_Risk_10yr)
    } else {
      datain_t1[,"CVD_prob"] <-datain_t[,"CVD_prob"]
    }
    
    #Defining transition probabilities
    # Subtract 19 from age because indexing starts at age 20
    p.death.DM <- case_when(
      data_for_analysis[,"DEMO"] == "Male" ~ DM_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "Female" ~ DM_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWM" ~ DM_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWF" ~ DM_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBM" ~ DM_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBF" ~ DM_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HM" ~ DM_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HF" ~ DM_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
    )
    
    p.death.CHD <- case_when(
      data_for_analysis[,"DEMO"] == "Male" ~ CHD_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "Female" ~ CHD_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWM" ~ CHD_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWF" ~ CHD_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBM" ~ CHD_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBF" ~ CHD_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HM" ~ CHD_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HF" ~ CHD_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
    )
    
    p.death.Stroke <- case_when(
      data_for_analysis[,"DEMO"] == "Male" ~ stroke_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "Female" ~ stroke_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWM" ~ stroke_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWF" ~ stroke_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBM" ~ stroke_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBF" ~ stroke_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HM" ~ stroke_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HF" ~ stroke_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
    )
    
    p.death <- case_when(
      data_for_analysis[,"DEMO"] == "Male" ~ Non_CVD_DM_mortality_Male[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "Female" ~ Non_CVD_DM_mortality_Female[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWM" ~ Non_CVD_DM_mortality_NHWM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHWF" ~ Non_CVD_DM_mortality_NHWF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBM" ~ Non_CVD_DM_mortality_NHBM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "NHBF" ~ Non_CVD_DM_mortality_NHBF[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HM" ~ Non_CVD_DM_mortality_HM[as.numeric(datain_t[,"Age_cycle"])-19,s],
      data_for_analysis[,"DEMO"] == "HF" ~ Non_CVD_DM_mortality_HF[as.numeric(datain_t[,"Age_cycle"])-19,s]
    )
    
    
    #Disaggregating ASCVD
    p.CHD_first <- as.numeric(Prop_CHD[1,data_for_analysis[,"DEMO"]]) #Prob. of CHD event (Angina, MI, Fatal MI, Fatal CHD) Source: Benjamin et al, Circulation, 2018 Table 13-1 & 18-1, and 18-2
    p.CHD_recurrent <- as.numeric(Prop_CHD[2,data_for_analysis[,"DEMO"]])
    
    
    #Markov State #1: "No CVD, No Diabetes"
    p.H.2.DM <- as.numeric(datain_t[,"DM_prob"])*as.numeric(data_for_analysis[,"risk_adjustment.DM"])
    p.H.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(datain_t[,"CVD_prob"])
    p.H.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(datain_t[,"CVD_prob"])
    p.H.2.H <- ifelse((p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death) > 1, 0,
                      1 - (p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death))
    
    #Markov State #2: "No CVD, With Diabetes
    p.DM.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(datain_t[,"CVD_prob"])
    p.DM.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(datain_t[,"CVD_prob"])
    sum.p.DM.2.death <- p.death + p.death.DM
    p.DM.2.DM <- ifelse((p.DM.2.initial_CHD + p.DM.2.initial_Stroke + sum.p.DM.2.death) > 1, 0,
                        1 - (p.DM.2.initial_CHD + p.DM.2.initial_Stroke + sum.p.DM.2.death))
    
    #Markov State #3: "First Stroke"
    sum.p.initial_Stroke.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.Stroke+p.death.DM,
                                           p.death+p.death.Stroke)
    p.initial_Stroke.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.initial_Stroke.2.death)
    p.initial_Stroke.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 
                                        1 - sum.p.initial_Stroke.2.death, 0)
    
    #Markov State #4: "First CHD w/o RVSC" 
    sum.p.initial_CHD_No_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.CHD+p.death.DM,
                                                p.death+p.death.CHD)
    p.initial_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.initial_CHD_No_RVSC.2.death)
    p.initial_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 
                                             1 - sum.p.initial_CHD_No_RVSC.2.death, 0)
    
    #Markov State #5: "First CHD with RVSC"
    sum.p.initial_CHD_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s],
                                             p.death+p.death.CHD+p.death.RVSC_sim[s])
    p.initial_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.initial_CHD_RVSC.2.death)
    p.initial_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 
                                          1 - sum.p.initial_CHD_RVSC.2.death, 0)
    
    #Markov State #6: "CVD History, No Diabetes"
    p.CVD.2.DM <- as.numeric(datain_t[,"DM_prob"])*as.numeric(data_for_analysis[,"risk_adjustment.DM"])
    p.CVD.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(datain_t[,"CVD_recurrent_prob"])
    p.CVD.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(datain_t[,"CVD_recurrent_prob"])
    sum.p.CVD.2.death <- p.death+p.death.CHD+p.death.Stroke
    p.CVD.2.CVD <- ifelse((p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + sum.p.CVD.2.death) > 1, 0,
                          1 - (p.CVD.2.DM+p.CVD.2.Sub_CHD+p.CVD.2.Sub_Stroke+sum.p.CVD.2.death))
    
    #Markov State #7: "CVD History, With Diabetes"
    p.CVD_DM.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(datain_t[,"CVD_recurrent_prob"])
    p.CVD_DM.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(datain_t[,"CVD_recurrent_prob"])
    sum.p.CVD_DM.2.death <- p.death+p.death.CHD+p.death.Stroke+p.death.DM
    p.CVD_DM.2.CVD_DM <- ifelse((sum.p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke) > 1, 0,
                                1 - (p.CVD_DM.2.Sub_CHD+p.CVD_DM.2.Sub_Stroke+sum.p.CVD_DM.2.death))
    
    #Markov State #8: "Subsequent Stroke"
    sum.p.sub_Stroke.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.Stroke+p.death.DM,
                                       p.death+p.death.Stroke)
    p.sub_Stroke.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.sub_Stroke.2.death)
    p.sub_Stroke.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - sum.p.sub_Stroke.2.death, 0)
    
    #Markov State #9: "Subsequent CHD w/o RVSC" 
    sum.p.sub_CHD_No_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.CHD+p.death.DM,
                                            p.death+p.death.CHD)
    p.sub_CHD_No_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.sub_CHD_No_RVSC.2.death)
    p.sub_CHD_No_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - sum.p.sub_CHD_No_RVSC.2.death, 0)
    
    #Markov State #10: "Subsequent CHD with RVSC"
    sum.p.sub_CHD_RVSC.2.death <- ifelse(as.numeric(sim_out[,"Diabetes",t]) == 1, p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s],
                                         p.death+p.death.CHD+p.death.RVSC_sim[s])
    p.sub_CHD_RVSC.2.CVD_No_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 0, 1 - sum.p.sub_CHD_RVSC.2.death)
    p.sub_CHD_RVSC.2.CVD_DM <- ifelse(as.numeric(sim_out[,"Diabetes",t]) ==1, 1 - sum.p.sub_CHD_RVSC.2.death, 0)
    
    #Assign transition probablities
    p.transition <- array(NA, dim=c(n.individual, n.health.state),
                          dimnames = list(data_for_analysis$Subject_ID,name.health.state))
    p.transition[,"No CVD, No Diabetes"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.H, rep(0, n.individual))
    p.transition[,"No CVD, With Diabetes"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.DM, 
                                                     ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.DM.2.DM, rep(0, n.individual)))
    p.transition[,"First Stroke"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.H.2.initial_Stroke, 
                                            ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.DM.2.initial_Stroke, rep(0, n.individual)))
    p.transition[,"First CHD w/o RVSC"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", (1-p.RVSC)*p.H.2.initial_CHD, 
                                                  ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", (1-p.RVSC)*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"First CHD with RVSC"] <- ifelse(sim_out[,"state",t] == "No CVD, No Diabetes", p.RVSC*p.H.2.initial_CHD, 
                                                   ifelse(sim_out[,"state",t] == "No CVD, With Diabetes", p.RVSC*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"CVD History, No Diabetes"] <- case_when(
      sim_out[,"state",t] == "First Stroke" ~ p.initial_Stroke.2.CVD_No_DM, 
      sim_out[,"state",t] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out[,"state",t] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_No_DM, 
      sim_out[,"state",t] == "CVD History, No Diabetes" ~ p.CVD.2.CVD, 
      sim_out[,"state",t] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_No_DM, 
      sim_out[,"state",t] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out[,"state",t] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_No_DM, 
      TRUE ~ rep(0, n.individual)
    )
    p.transition[,"CVD History, With Diabetes"] <- case_when(
      sim_out[,"state",t] == "First Stroke" ~ p.initial_Stroke.2.CVD_DM, 
      sim_out[,"state",t] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_DM, 
      sim_out[,"state",t] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_DM, 
      sim_out[,"state",t] == "CVD History, No Diabetes" ~ p.CVD.2.DM, 
      sim_out[,"state",t] == "CVD History, With Diabetes" ~ p.CVD_DM.2.CVD_DM, 
      sim_out[,"state",t] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_DM, 
      sim_out[,"state",t] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_DM, 
      sim_out[,"state",t] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_DM, 
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Subsequent Stroke"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", p.CVD.2.Sub_Stroke,
                                                 ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", p.CVD_DM.2.Sub_Stroke, rep(0, n.individual)))
    p.transition[,"Subsequent CHD w/o RVSC"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", (1-p.RVSC)*p.CVD.2.Sub_CHD, 
                                                       ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", (1-p.RVSC)*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    p.transition[,"Subsequent CHD with RVSC"] <- ifelse(sim_out[,"state",t] == "CVD History, No Diabetes", p.RVSC*p.CVD.2.Sub_CHD, 
                                                        ifelse(sim_out[,"state",t] == "CVD History, With Diabetes", p.RVSC*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    p.transition[,"DM_Death"] <- case_when(
      sim_out[,"state",t] == "DM_Death" ~ rep(1, n.individual),
      sim_out[,"state",t] %in% c("Non_DM_Non_CVD_Death", "Stroke_Death", "CHD_Death") ~ rep(0, n.individual),
      sim_out[,"Diabetes",t] == 1 ~ p.death.DM, 
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Stroke_Death"] <- case_when(
      sim_out[,"state",t] %in% 
        c("First Stroke", "Subsequent Stroke", "CVD History, No Diabetes", "CVD History, With Diabetes") ~ p.death.Stroke, 
      sim_out[,"state",t] == "Stroke_Death" ~ rep(1, n.individual),
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"CHD_Death"] <- case_when(
      sim_out[,"state",t] %in%
        c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes") ~ p.death.CHD,
      sim_out[,"state",t] %in% c("First CHD with RVSC", "Subsequent CHD with RVSC") ~ p.death.CHD+p.death.RVSC_sim[s],
      sim_out[,"state",t] == "CHD_Death" ~ rep(1, n.individual),
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Non_DM_Non_CVD_Death"] <- case_when(
      sim_out[,"state",t] %in% c("DM_Death", "Stroke_Death", "CHD_Death") ~ rep(0, n.individual),
      sim_out[,"state",t] == "Non_DM_Non_CVD_Death" ~ rep(1, n.individual),
      TRUE ~ p.death
    )
    # Check that all transition probabilities add to 1 (rounded to 3 digits)
    # if (sum(round(rowSums(p.transition),3)!=1) != 0) {
    #  print("Transition probabilities do not add to 1")
    #}
    # set.seed(seed)
    # Transition to the next health state 
    sim_out[,"state",t+1] <- apply(p.transition, 1, function(x) sample(name.health.state, 1, prob = x))
    
    #Updating Disease history status
    sim_out[,"Diabetes",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("No CVD, With Diabetes", "CVD History, With Diabetes"), 1, 
                                       ifelse(sim_out[,"state",t+1] %in% name.death.states, NA, 
                                              ifelse(as.numeric(datain_t1[,"Glucose"]) > 126, 1, sim_out[,"Diabetes",t])))
    
    sim_out[,"CVD_history",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "First CHD w/o RVSC", "First CHD with RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes",
                                                                       "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), 1, 
                                          ifelse(sim_out[,"state",t+1] %in% name.death.states, NA, sim_out[,"CVD_history",t]))
    
    datain_t1[,"BMI"] <- ifelse(sim_out[,"state",t+1] %in% name.death.states, NA, datain_t1[,"BMI"])
    sim_out[,"Obesity",t+1] <- ifelse(as.numeric(datain_t1[,"BMI"]) < 30, 0, 
                                      ifelse((as.numeric(datain_t1[,"BMI"]) >= 30), 1, NA))
    
    datain_t1[,"CVD_history"] <-sim_out[,"CVD_history",t+1]
    datain_t1[,"Diabetes"] <-sim_out[,"Diabetes",t+1]
    datain_t1[,"Obesity"] <-sim_out[,"Obesity",t+1]

    keyvars2 <- c("SEQN","Age_cycle","SBP", "DBP",  "BMI")
    raw.input.data<-cbind(matrix(as.numeric(datain_t1[,keyvars2]), ncol=length(keyvars2)),as.numeric(sim_out[,"Diabetes",t+1]), as.numeric(sim_out[,"CVD_history",t+1]) )
    colnames(raw.input.data)<-c("seqn","Age","SBP", "DBP",  "BMI","Diabetes","CVD_history")
    raw.input.data <- as.data.frame(cbind(raw.input.data, data_for_analysis[,c("Female", "Race","HPT_Txt")]))
    
    # Update QALYs and costs
    sim_out[,"HRQOL_scores",t+1] <- calc_HRQOL(raw.input.data, HRQOL_parameter_sim[,s])
    sim_out[,"HRQOL_scores",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_stroke_sim[s], 
                                           ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "First CHD with RVSC", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HRQOL_scores",t+1])+u_CHD_sim[s],
                                                  ifelse(sim_out[,"state",t+1] %in% name.death.states, 0, sim_out[,"HRQOL_scores",t+1])))
    
    #ppp: for sensitivity analysis, add economic effect of HA1c
    
    
    sim_out[,"HCE_predict",t+1] <- calc_HCE(raw.input.data, HCE_parameter_sim[,s])
    sim_out[,"HCE_predict",t+1] <- ifelse(sim_out[,"state",t+1] %in% c("First Stroke", "Subsequent Stroke"), as.numeric(sim_out[,"HCE_predict",t+1])+c_stroke_sim[s], 
                                          ifelse(sim_out[,"state",t+1] %in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s], 
                                                 ifelse(sim_out[,"state",t+1] %in% c("First CHD with RVSC", "Subsequent CHD with RVSC"), as.numeric(sim_out[,"HCE_predict",t+1])+c_CHD_sim[s]+c_RVSC_sim[s],
                                                        ifelse(sim_out[,"state",t+1] %in% name.death.states, NA, sim_out[,"HCE_predict",t+1]))))
    
    # Discount effects and costs
    sim_out[,"effect_disc",t+1] <- as.numeric(sim_out[,"HRQOL_scores",t+1])/((1+beta_QALY)^(t-1))
    sim_out[,"HCE_disc",t+1] <- as.numeric(sim_out[,"HCE_predict",t+1])/((1+beta_cost)^(t-1))
    
datain_t=datain_t1
  }
  
  # Subset output to include variables of interest
 out_variables <- c("Obesity", "Diabetes", "CVD_history", "state")
 sim_out <- sim_out[,out_variables,]
  
options(survey.lonely.psu = "adjust")
design <- svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~wtmec2yr, nest=TRUE, data=data_for_analysis)
design_NHW <- subset(design, (DEMO=="NHWM" | DEMO=="NHWF"))
design_NHB <- subset(design, (DEMO=="NHBM" | DEMO=="NHBF"))
design_Hisp <- subset(design, (DEMO=="HM" | DEMO=="HF"))

# Subset sim_out array
m.state <- array(sim_out[,"state",1:(n.cycle+1)], dim=c(n.individual, n.cycle+1),
                 dimnames = list(paste("ind", 1:n.individual, sep = " "), paste("cycle", 0:n.cycle, sep = " ")))

#m.alive <- array(!(m.state %in% name.death.states), dim = c(n.individual,n.cycle), 
 #                    dimnames = list(data_for_analysis$Subject_ID, paste("cycle", 1:n.cycle, sep = " ")))

# Define relevant variables for output
#vars <- c("Obesity", "Diabetes", "CVD_history", "HCE_predict", "effect_disc", "cost_disc", "HRQOL_scores", 
 #         "All-cause mortality","DM mortality", "CVD mortality", "Non CMD mortality" ,"Incident Stroke", "Incident CHD","race")
vars <- c("Obesity", "Diabetes", "CVD_history", "All-cause mortality", "DM mortality", "CVD mortality", "Non CMD mortality", "race")

Vars<- c("Obesity", "Diabetes", "CVD_history")


races <- c("all", "NHW", "NHB", "Hisp")

stat<-c("mean", "var")
aa<-(n.cycle+1)*4
mean_out <- array(NA, dim = c(length(vars), aa, 2), 
                  dimnames = list(vars, rep(paste("cycle", 0:n.cycle, sep = " "), 4), stat))

for (race in races) {
  print(race)

  # Initialize arrays to hold output

  race_num = switch(race, "NHW" = 1, "NHB" = 2, "Hisp" = 3)
  str=switch(race, "all"=1, "NHW" = 2, "NHB" = 3, "Hisp" = 4)
  race_design = switch(race,
                       "NHW" = design_NHW,
                       "NHB" = design_NHB,
                       "Hisp" = design_Hisp)
    
for (t in 1:(n.cycle+1)) {
    # Calculate weighted means and variances (SE^2) across individuals (new dimensions: cycle x simulation)

t2<-t+(n.cycle+1)*(str-1) 
if (race == "all") {

allcause_mort<-m.state[,t]=="DM_Death"|m.state[,t]=="CHD_Death"|m.state[,t]=="Stroke_Death"|m.state[,t]=="Non_DM_Non_CVD_Death"
dm_mort<-m.state[,t]=="DM_Death"
cvd_mort<-m.state[,t]=="CHD_Death"|m.state[,t]=="Stroke_Death"
noncmd_mort<- m.state[,t]=="Non_DM_Non_CVD_Death"

for (var in Vars) {
      mean_out[var,t2, "mean"] <- svymean(as.numeric(sim_out[,var,t]), design, deff=F, na.rm=T)
      mean_out [var,t2, "var"] <- SE(svymean(as.numeric(sim_out[,var,t]), design, deff=F, na.rm=T))^2
}

mean_out["All-cause mortality",t2, "mean"] <- svymean(allcause_mort, design, deff=F, na.rm=T)
mean_out["All-cause mortality",t2, "var"] <- SE(svymean(allcause_mort, design, deff=F, na.rm=T))^2

mean_out["DM mortality",t2, "mean"] <- svymean(dm_mort, design, deff=F, na.rm=T)
mean_out["DM mortality",t2, "var"] <- SE(svymean(dm_mort, design, deff=F, na.rm=T))^2

mean_out["Non CMD mortality",t2, "mean"] <- svymean(noncmd_mort, design, deff=F, na.rm=T)
mean_out["Non CMD mortality",t2, "var"] <- SE(svymean(noncmd_mort, design, deff=F, na.rm=T))^2


mean_out["CVD mortality",t2, "mean"] <- svymean(cvd_mort, design, deff=F, na.rm=T)
mean_out["CVD mortality",t2, "var"] <- SE(svymean(cvd_mort, design, deff=F, na.rm=T))^2

mean_out["race",t2, "mean"] <- race
mean_out["race",t2, "var"] <- race

#mean_out["Incident CHD",t2, "mean"] <- svymean(chd, design, deff=F, na.rm=T)
#mean_out["Incident CHD",t2, "var"] <- SE(svymean(chd, design, deff=F, na.rm=T))^2
#mean_out["Incident Stroke",t2, "mean"] <- svymean(stroke, design, deff=F, na.rm=T)
#mean_out["Incident Stroke",t2, "var"] <- SE(svymean(stroke, design, deff=F, na.rm=T))^2

} else {

cvd_mort<- m.state[which(data_for_analysis$Race==race_num),t]=="CHD_Death"|m.state[which(data_for_analysis$Race==race_num),t]=="Stroke_Death"

dm_mort<- m.state[which(data_for_analysis$Race==race_num),t]=="DM_Death"

allcause_mort<- m.state[which(data_for_analysis$Race==race_num),t]=="DM_Death"|m.state[which(data_for_analysis$Race==race_num),t]=="CHD_Death"|m.state[which(data_for_analysis$Race==race_num),t]=="Stroke_Death"|m.state[which(data_for_analysis$Race==race_num),t]=="Non_DM_Non_CVD_Death"

noncmd_mort<- m.state[which(data_for_analysis$Race==race_num),t]=="Non_DM_Non_CVD_Death"
      
#stroke<-m.state[which(data_for_analysis$Race==race_num),t]=="First Stroke"

#chd<-m.state[which(data_for_analysis$Race==race_num),t]=="First CHD with RVSC"| m.state[which(data_for_analysis$Race==race_num),t]=="First CHD w/o RVSC"


for (var in Vars) {
mean_out[var,t2, "mean"] <- svymean(as.numeric(sim_out[which(data_for_analysis$Race==race_num),var,t]), race_design, deff=F, na.rm=T)
mean_out[var,t2, "var"] <- SE(svymean(as.numeric(sim_out[which(data_for_analysis$Race==race_num),var,t]), race_design, deff=F, na.rm=T))^2
}


mean_out["All-cause mortality",t2, "mean"] <- svymean(allcause_mort, race_design, deff=F, na.rm=T)
mean_out["All-cause mortality",t2, "var"] <- SE(svymean(allcause_mort, race_design, deff=F, na.rm=T))^2

mean_out["DM mortality",t2, "mean"] <- svymean(dm_mort, race_design, deff=F, na.rm=T)
mean_out["DM mortality",t2, "var"] <- SE(svymean(dm_mort, race_design, race_deff=F, na.rm=T))^2

mean_out["Non CMD mortality",t2, "mean"] <- svymean(noncmd_mort, race_design, deff=F, na.rm=T)
mean_out["Non CMD mortality",t2, "var"] <- SE(svymean(noncmd_mort, race_design, deff=F, na.rm=T))^2

#mean_out["Incident Stroke",t2, "mean"] <- svymean(stroke, race_design, deff=F, na.rm=T)
#mean_out["Incident Stroke",t2, "var"] <- SE(svymean(stroke, race_design, deff=F, na.rm=T))^2

mean_out["CVD mortality",t2, "mean"] <- svymean(cvd_mort, race_design, deff=F, na.rm=T)
mean_out["CVD mortality",t2, "var"] <- SE(svymean(cvd_mort, race_design, deff=F, na.rm=T))^2

#mean_out["Incident CHD",t2, "mean"] <- svymean(chd, race_design, deff=F, na.rm=T)
#mean_out["Incident CHD",t2, "var"] <- SE(svymean(chd, race_design, deff=F, na.rm=T))^2


mean_out["race",t2, "mean"] <- race
mean_out["race",t2, "var"] <- race

    }
  }

}


return (mean_out)
}

# Detect system type
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
    if (os == "windows")
      os <- "windows"
  } else { ## if we still don't know
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

# Check if we are on Windows or Mac using our function.
cluster_type <- if (get_os() == "windows") {"PSOCK"} else {"FORK"}

no_cores <- detectCores() - 1
s = 1 # initialize iteration variable before forking
cl<-makeCluster(no_cores, type=cluster_type) # Make the cluster
# clusterEvalQ(cl)
registerDoParallel(cl)

acomb <- function(...) abind(..., along = 4)
sim_out <- foreach(s=1:n.sim, .combine = 'acomb', .verbose = T) %do% {
  run_sim(s)
}

stopCluster(cl)

proc.time() - ptm

#########################################################################
# 6 Summarize and save output                                           #
#########################################################################
print("Summarizing and saving output")
summ_start = proc.time()


  saveRDS(sim_out, file = paste("40 Results/sim_out_1000", intervention, "SEED", seed, Sys.Date(), ".rda", sep = "_"))


print("Time to summarize and save output:")
proc.time() - summ_start

print("Total time:")
proc.time() - ptm
