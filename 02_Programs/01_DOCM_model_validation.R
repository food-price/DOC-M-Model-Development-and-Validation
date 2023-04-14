#Diabetes-Obesity, CVD Microsimulation (DOC-M) Model 
#The FOOD-PRICE project (https://food-price.org/)
#model validation
#Primary Author: David D. Kim, University of Chicago (ddk@uchicago.edu)
#Secondary Authors:  Lu Wang (lulu_stat@hotmail.com) and Brianna Lauren 
#Model Summary: Programmed in R-4.1.0, the DOC-M model is a probabilistic and dynamic microsimulation model that projects obesity, diabetes, CVD, 
#and their associated complications for guiding population health and policy decisions. Using US population-based transition probabilities, the model
#tracks a person's annual likelihood of experiencing health events (e.g., developing diabetes and CVD) and death based on individual factors: 
#age, sex, race, blood pressure, total cholesterol, smoking status, and others. Each individual in the DOC-M model can transit through multiple health states
#each year: no CVD or diabetes, diabetes without CVD, CVD without diabetes, both CVD and diabetes, and death, plus four CVD-related events 
#(first or recurrent stroke or coronary heart disease, with an option for revascularization. The model also captures the incidence and prevalence of overweight 
#and obesity based on each individual's BMI. Our modeled population is US adults aged 40-79, corresponding to the sample populations from which diabetes and 
#CVD risk predictions were derived, and our model provides US population estimates by aggregating individual trajectories with appropriate survey weights. 
#Analyses were conducted using probability distributions for all input parameters to capture parameter uncertainty. 
#The comprehensive development and validation of the DOC-M model have been described elsewhere. (Kim et al. Med. Decis. Mak. 2023) 

#################
# 0.Preparation #
#################

# Start the clock
ptm <- proc.time()

#memory.size(max=T) # this is Windows-specific and will not work on the HPC running Linux

# 0.1 Install/Read R Packages

library(survey)
library(svMisc)
library(psych)
library(gdata)
library(dplyr)
library(data.table)   # needed for fread/fwrite
library(foreach)
library(doParallel)
library(abind)
#library(ggplot2)
library(dampack)
library(writexl)

# remove all data from memory
rm(list = ls())

# 0.2 Create Functions 

# 0.2-1 To estimate gamma/beta parameters based on Mean and SD
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

# 0.2-2 Converting multi-year risk (e.g., 10-year ASCVD risk) to annual probabilities
Multi_yr_Risk_to_annual_prob <- function(time, risk) {
  annual_rate <- -(1/time)*log(1-risk)
  annual_prob <- 1 - exp(-annual_rate)
}

# 0.3 Creating Working Directory 


setwd("C:\\Users\\lwang18\\Box\\Model Development - ODC-M\\DOC-M_Validation_Git\\Final")
#setwd("/cluster/tufts/kimlab/lwang18/DOC-M_Validation_Git")

# 0.4 Source other key scripts
source("02_Programs/1a - Diabetes_risk_prediction_FHS (categorical points).R")
source("02_Programs/1b - ASCVD_risk_calculator (calibrated_v2).R")
source("02_Programs/1c - FHS Subsequent_CVD_risk_calculator.R")
source("02_Programs/2a - HrQOL estimator for US general population.R")
source("02_Programs/2b - HCE estimator for US general population.R")

# 0.5 Model settings (set manually or read from command line when submitted through cluster)
args <- commandArgs(trailingOnly = TRUE)  # get all arguments after script.R
# If no arguments are read from command line, set modeling choices manually. 
# Otherwise, read modeling choices from command line.
if (length(args) == 0) {
  seed <- 1234
  n.sim <-  10 #Number of probablistic samplings
  n.cycle <- 5 #Number of Cycle Length: How long does the model run (i.e., analytic time horizon)
  n.sample <- "ALL" #Number of individuals to be selected from the full sample; if full sample, enter "ALL"
  n.loop=1 ##set to 1 for model validation
  #n.loop is the Number of replicates for each individual, each person is replicated for n.loop times and the outcomes are averaged across these replicates 
  intervention="No Policy"
} else {
  # expecting 4 arguments
  if (length(args) != 4) {
    stop("ERROR: Incorrect number of command line arguments", call. = FALSE)
  }
  seed <- as.numeric(args[1]) # extracting first arg as seed and attempt to cast to number
  n.sim <- as.numeric(args[2])
  n.cycle <- as.numeric(args[3])
  # Assume entire sample
  n.sample = "ALL"
  n.loop=as.numeric(args[4])
  ##set to 1 for model validation
  #n.loop is the Number of replicates for each individual, each person is replicated for n.loop times and the outcomes are averaged across these replicates 
  intervention="No Policy"
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
}

output_pathRDS =  paste("03_Output/sim_out", n.sim, intervention, "SEED", seed,  ".rda", sep = "_")
output_pathCVS = "03_Output/Processed/final_output.csv"

set.seed(seed)

################################################################################################
# 1 Defining and Importing Necessary Impute parameters, creating n.sim random draws
################################################################################################

##$$Policy specific settings

# 1.0 Policy-effect size, costs, and discounting rate

# policy effect not exist for model validation 
# Policy costs not exist for model validation

beta_cost <- beta_QALY <- 0.03 #Annual discounting rate of costs and QALYs


# 1.1 Read in master input file

print('Importing data')

data_for_analysis<- fread("01_Input/NHANES/NHANES0102_input.csv", stringsAsFactors = TRUE, data.table = FALSE)

Agesex<-ifelse(data_for_analysis$Age<65 & data_for_analysis$Female==1,1,
               ifelse(data_for_analysis$Age<65 & data_for_analysis$Female==0,2,
                      ifelse(data_for_analysis$Age>=65 & data_for_analysis$Female==1,3, 
                             ifelse(data_for_analysis$Age>=65 & data_for_analysis$Female==0,4,0))))

# 1.1.1 DM risk adjustment for non-whites
data_for_analysis$risk_adjustment.DM <- ifelse(data_for_analysis$DEMO %in% c("NHWM", "NHWF", "Female", "Male"), 1.0,
                                               ifelse(data_for_analysis$DEMO %in% c("NHBM", "NHBM"), 1.5, 2.4))
# 1.2 health states

##Initial health states.

name.health.state <- c("No CVD, No Diabetes", "No CVD, With Diabetes", "First Stroke", "First CHD w/o RVSC", "First CHD with RVSC",
                       "CVD History, No Diabetes", "CVD History, With Diabetes", "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC", 
                       "DM_Death", "Stroke_Death", "CHD_Death", "Non_DM_Non_CVD_Death") 
n.health.state <- length(name.health.state) # number of health state that individuals can transition over time
name.death.states <- c("DM_Death", "Stroke_Death", "CHD_Death", "Non_DM_Non_CVD_Death")


# 1.3 Health-state/Event specific mortality data  

#1.3-3a: Re-calibration
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

Non_CVD_DM_mortality <- fread("01_Input/Non_DM_IHD_Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ corrected).csv", stringsAsFactors = TRUE, data.table = FALSE)
Non_CVD_DM_mortality_Male <- Adj_Non_CVD_DM_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Male, Non_CVD_DM_mortality$Male_SE))
Non_CVD_DM_mortality_Female <- Adj_Non_CVD_DM_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$Female, Non_CVD_DM_mortality$Female_SE))
Non_CVD_DM_mortality_NHWM <- Adj_Non_CVD_DM_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWM, Non_CVD_DM_mortality$NHWM_SE))
Non_CVD_DM_mortality_NHWF <- Adj_Non_CVD_DM_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHWF, Non_CVD_DM_mortality$NHWF_SE))
Non_CVD_DM_mortality_NHBM <- Adj_Non_CVD_DM_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBM, Non_CVD_DM_mortality$NHBM_SE))
Non_CVD_DM_mortality_NHBF <- Adj_Non_CVD_DM_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$NHBF, Non_CVD_DM_mortality$NHBF_SE))
Non_CVD_DM_mortality_HM <- Adj_Non_CVD_DM_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HM, Non_CVD_DM_mortality$HM_SE))
Non_CVD_DM_mortality_HF <- Adj_Non_CVD_DM_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, Non_CVD_DM_mortality$HF, Non_CVD_DM_mortality$HF_SE))

stroke_mortality <- fread("01_Input/Stroke_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
stroke_mortality_Male <- Adj_CVD_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Male, stroke_mortality$Male_SE))
stroke_mortality_Female <- Adj_CVD_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$Female, stroke_mortality$Female_SE))
stroke_mortality_NHWM <- Adj_CVD_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWM, stroke_mortality$NHWM_SE))
stroke_mortality_NHWF <- Adj_CVD_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHWF, stroke_mortality$NHWF_SE))
stroke_mortality_NHBM <- Adj_CVD_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBM, stroke_mortality$NHBM_SE))
stroke_mortality_NHBF <- Adj_CVD_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$NHBF, stroke_mortality$NHBF_SE))
stroke_mortality_HM <- Adj_CVD_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HM, stroke_mortality$HM_SE))
stroke_mortality_HF <- Adj_CVD_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, stroke_mortality$HF, stroke_mortality$HF_SE))

CHD_mortality <- fread("01_Input/IHD_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
CHD_mortality_Male <- Adj_CVD_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Male, CHD_mortality$Male_SE))
CHD_mortality_Female <- Adj_CVD_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$Female, CHD_mortality$Female_SE))
CHD_mortality_NHWM <- Adj_CVD_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWM, CHD_mortality$NHWM_SE))
CHD_mortality_NHWF <- Adj_CVD_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHWF, CHD_mortality$NHWF_SE))
CHD_mortality_NHBM <- Adj_CVD_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBM, CHD_mortality$NHBM_SE))
CHD_mortality_NHBF <- Adj_CVD_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$NHBF, CHD_mortality$NHBF_SE))
CHD_mortality_HM <- Adj_CVD_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HM, CHD_mortality$HM_SE))
CHD_mortality_HF <- Adj_CVD_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, CHD_mortality$HF, CHD_mortality$HF_SE))

DM_mortality <- fread("01_Input/DM_Cause_Mortality (Annual probability) (2001-2016 & 85+ adjusted).csv", stringsAsFactors = TRUE, data.table = FALSE)
DM_mortality_Male <- Adj_DM_mortality_Male*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Male, DM_mortality$Male_SE))
DM_mortality_Female <- Adj_DM_mortality_Female*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$Female, DM_mortality$Female_SE))
DM_mortality_NHWM <- Adj_DM_mortality_NHWM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWM, DM_mortality$NHWM_SE))
DM_mortality_NHWF <- Adj_DM_mortality_NHWF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHWF, DM_mortality$NHWF_SE))
DM_mortality_NHBM <- Adj_DM_mortality_NHBM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBM, DM_mortality$NHBM_SE))
DM_mortality_NHBF <- Adj_DM_mortality_NHBF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$NHBF, DM_mortality$NHBF_SE))
DM_mortality_HM <- Adj_DM_mortality_HM*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HM, DM_mortality$HM_SE))
DM_mortality_HF <- Adj_DM_mortality_HF*t(mapply(calc_nsims_rbeta, n.sim, DM_mortality$HF, DM_mortality$HF_SE))


# 1.4 Secular trends in major risk factors: estimated as average annual percent change over 1999-2016 

risk_factor_trend <- fread("01_Input/Risk_Factors_Trends_Summary (Final)_calibrated.csv", stringsAsFactors = TRUE, data.table = FALSE)
risk_factor_trend_sim <- matrix(NA, nrow=nrow(risk_factor_trend), ncol=n.sim)
colnames(risk_factor_trend_sim) <- paste("simulation_", 1:n.sim, sep = " ")

risk_factor_trend_sim <- t(mapply(rnorm, n=n.sim, risk_factor_trend$APC, risk_factor_trend$APC_SE))

for (i in 1:nrow(risk_factor_trend)) {
  risk_factor_trend_sim[i,] <- as.matrix(rnorm(n.sim, risk_factor_trend$APC[i], risk_factor_trend$APC_SE[i]))
}

risk_factor_trend_sim <- cbind(risk_factor_trend, risk_factor_trend_sim)


# 1.5 Gender- and race-specific proportiona of CHD cases Among ASCVD cases (CHD + Stroke)
Prop_CHD <- t(fread("01_Input/Prop_CHD.csv", stringsAsFactors = FALSE, data.table = FALSE)[,-1])

data_for_analysis$p_CHD_first<-Prop_CHD[data_for_analysis[,"DEMO"],1]
data_for_analysis$p_CHD_recurrent<-Prop_CHD[data_for_analysis[,"DEMO"],2]

#1.6 Health Care Expenditure Model
HCE_parameter_sim <- matrix(NA, nrow=nrow(HCE_parameters), ncol=n.sim)
rownames(HCE_parameter_sim) <- rownames(HCE_parameters)
colnames(HCE_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

HCE_parameter_sim_diabe <- HCE_parameter_sim

for (i in 1:nrow(HCE_parameters)) {
  HCE_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HCE_parameters$Beta[i], HCE_parameters$SE[i]))
}

#1.7 HrQOL Prediction
HRQOL_parameter_sim <- matrix(NA, nrow=nrow(HRQOL_parameters), ncol=n.sim)
rownames(HRQOL_parameter_sim) <- rownames(HRQOL_parameters)
colnames(HRQOL_parameter_sim) <- paste("simulation_", 1:n.sim, sep = " ")

for (i in 1:nrow(HRQOL_parameters)) {
  HRQOL_parameter_sim[i,] <- as.matrix(rnorm(n.sim, HRQOL_parameters$Beta[i], HRQOL_parameters$SE[i]))
}

#1.8 additional health care cost parameters

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


# 1.9 Productivity costs associated with CHD, stroke
#ppp Add productivity costs associated with stroke, and CHD respectively 
#average annual productivity costs for stroke and CHD cases was calculated based on total annual productivity costs and population counts of stroke and CHD in the US in 2019 respectively.
#Source: https://www.heart.org/idc/groups/heart-public/@wcm/@adv/documents/downloadable/ucm_491513.pdf
#$ were inflated from 2015 dollar to 2021 dollar

c_prod_stroke =4661
c_prod_stroke_sim = calc_nsims_rgamma(n.sim, c_prod_stroke, c_prod_stroke*0.2)

c_prod_CHD =6833
c_prod_CHD_sim = calc_nsims_rgamma(n.sim, c_prod_CHD, c_prod_CHD*0.2)

c_prod_CVD =6163
c_prod_CVD_sim = calc_nsims_rgamma(n.sim, c_prod_CVD, c_prod_CVD*0.2)

# 1.10 Import risk factor-disease etiologic effects data inputs : Age-specific relative risk estimates between risk factor and disease outcomes
## this can be updated when updated evidence available 
## not relevant for model validation



#################################################################################################################################
# 2 Estimating disease-specific risk, health-related quality of life (HrQOL), and healthcare expenditures (HCE) at the baseline #
#################################################################################################################################

# 2.1 FHS 8-year Diabetes Risk Prediction (With BMI)
variable_for_raw.input <- names(data_for_analysis)
raw.input.data <- data_for_analysis
data_for_analysis$DM_risk_8yr <- calc_DM_risk(raw.input.data)
data_for_analysis$DM_prob <- Multi_yr_Risk_to_annual_prob(time=8, risk=data_for_analysis$DM_risk_8yr)

# 2.2 ACC/AHA ASCVD 10-year Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
data_for_analysis$CVD_prob <- Multi_yr_Risk_to_annual_prob(time=10, risk=data_for_analysis$ASCVD_Risk_10yr)

# 2.3 FHS 2-year CVD recurrent Risk Prediction
raw.input.data <- data_for_analysis
data_for_analysis$CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
data_for_analysis$CVD_recurrent_prob <- Multi_yr_Risk_to_annual_prob(time=2, risk=data_for_analysis$CVD_Recurrent_risk_2yr)

# 2.4 Individual HrQOL prediction
raw.input.data <- data_for_analysis
data_for_analysis$HRQOL_scores <- calc_HRQOL(raw.input.data, HRQOL_parameters[,1])

# 2.5 Individual HCE prediction
raw.input.data <- data_for_analysis
data_for_analysis$HCE_predict <- calc_HCE(raw.input.data, HCE_parameters[,1])


n.individual <- nrow(data_for_analysis)*n.loop # number of individuals in the model

#########################################################################
# 3 Main Model - PLEASE DO NOT MODIFY UNLESS YOU KNOW WHAT YOU'RE DOING #
#########################################################################
# Define output variables, this can be modified depending on project needs. 
##fullset of output variables that is available for output
##Out_variables<-c( "SEQN","Subject_ID","Obesity","Diabetes", "Incident First CVD", "Incident Recurrent CVD","Incident CVD","CVD_history", "Non_DM_Non_CVD_Death", "CHD_Death", "Stroke_Death", "DM_Death", "Life Years","effect_disc", "HCE_disc", "Prod_cost", "Food_cost", "Food_cost_disc",  "Prod_cost_disc","Admin_costs", "Total_cost_health","Total_cost_societ" )
Out_variables<-c( "Subject_ID","Obesity","Diabetes","CVD_history", "All_Deaths","DM_Death", "CVD_Death", "Non_CMD_Death" )

# 3.1 source the scrit that defined the simulation function for running n.sim simulations 
source("02_Programs/3_sim_function_V.R")

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
  run_sim(s, "No Policy", 1)
}

stopCluster(cl)



#########################################################################
# 6 Summarize and save output                                           #
#########################################################################
print("Summarizing and saving output")
summ_start = proc.time()

##Save model output, a four dimension array

  saveRDS(sim_out, file = output_pathRDS )
  

# Process output for model validation
  
  model_out_means = readRDS(output_pathRDS)[,,"mean",]
  model_out_vars = readRDS(output_pathRDS)[,,"var",]
  
  vars <- c("Obesity", "Diabetes", "CVD_history", "All_Deaths", "DM_Death", "CVD_Death", "Non_CMD_Death" )
  races <- c("all", "NHW", "NHB", "Hisp")
  year1 = 2001
  
  
  # Initialize final output (means and 95% CIs for each variable, race, and cycle - long format)
  final_out = data.frame()
  for (race in races) {
    # Combine output from several runs
    r=switch(race, "all"=1, "NHW" = 2, "NHB" = 3, "Hisp" = 4)
    for (i_var in vars) {
      temp.means =matrix(as.numeric(model_out_means[i_var,((n.cycle+1)*(r-1)+1):((n.cycle+1)*r),]), nrow=n.cycle+1)
      temp.within.var = matrix( as.numeric(model_out_vars[i_var,((n.cycle+1)*(r-1)+1):((n.cycle+1)*r),]),nrow=n.cycle+1)
      temp_out <- data.frame(Year = year1:(year1+n.cycle))
      temp_out$Outcome <- i_var
      temp_out$RE <- race
      # Mean of means across simulations (final mean per cycle)
      temp_out$mean = apply(temp.means, 1, mean)
      # Mean of variances across simulations (final within-simulation variance per cycle)
      within.var <- apply(temp.within.var, 1, mean)
      # Variance (SE^2) of means across simulations (final between-simulation variance per cycle)
      between.var <- apply(temp.means, 1, function(x) sd(x)^2/n.sim)
      # Final variance calculation per cycle (Rubin's rule formula)
      # Source: Dakin et al. Accurately Reflecting Uncertainty When Using Patient-Level Simulation Models to Extrapolate Clinical Trial Data. MDM. 2020.
      final.var <- within.var + (1 + (1/n.sim)) * between.var
      final.se <- sqrt(final.var)
      # Calculate lower limit of 95% CI
      temp_out$LL = temp_out$mean - 1.96 * final.se
      # Calculate upper limit of 95% CI
      temp_out$UL = temp_out$mean + 1.96 * final.se
      final_out = rbind(final_out, temp_out)
    }
  }
  
  # Save output
  write.csv(final_out, file = output_pathCVS, row.names = F)
  
print("Time to summarize and save output:")
proc.time() - summ_start

