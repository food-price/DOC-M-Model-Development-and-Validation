##10-year ASCVD Risk Calculator##
##Primary Source: Goff et al. (Circulation, 2014) Appendix 7, Table A
##Author: David D. Kim, Tufts Medical Center, DKim3@tuftsmedicalcenter.org

##Risk Equation Parameters
ASCVD_risk_params <- read.csv("01_Input/ASCVD Risk Equation Parameters.csv", row.names = 1)

calc_ASCVD_risk <- function(raw.input.data) {
  ##Input Data
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(ASCVD_risk_params))
  colnames(risk.input.data) <- rownames(ASCVD_risk_params)
  
  risk_score <- matrix(nrow=nrow(raw.input.data), ncol=4)
  colnames(risk_score) <- c("individual Sum", "Sex_Race_Overall_Mean", "Baseline_Survival", "ASCVD_Risk_10yr")
  
  risk.input.data[,1] = log(raw.input.data[,"Age"])
  risk.input.data[,2] = (log(raw.input.data[,"Age"]))^2
  risk.input.data[,3] = log(raw.input.data[,"Total_Chol"])
  risk.input.data[,4] = log(raw.input.data[,"Age"])*log(raw.input.data[,"Total_Chol"])
  risk.input.data[,5] = log(raw.input.data[,"HDL"])
  risk.input.data[,6] = log(raw.input.data[,"Age"])*log(raw.input.data[,"HDL"])
  risk.input.data[,7] = log(raw.input.data[,"SBP"])*raw.input.data[,"HPT_Txt"]
  risk.input.data[,8] = log(raw.input.data[,"Age"])*log(raw.input.data[,"SBP"])*raw.input.data[,"HPT_Txt"]
  risk.input.data[,9] = log(raw.input.data[,"SBP"])*(1-raw.input.data[,"HPT_Txt"])
  risk.input.data[,10] = log(raw.input.data[,"Age"])*log(raw.input.data[,"SBP"])*(1-raw.input.data[,"HPT_Txt"])
  risk.input.data[,11] = raw.input.data[,"Smoking"]
  risk.input.data[,12] = log(raw.input.data[,"Age"])*raw.input.data[,"Smoking"]
  risk.input.data[,13] = raw.input.data[,"Diabetes"]
  
  #Risk_calculation
  
  #Part 1: Individual Sum (=coefficients*value)
  individual_sum <- ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]!=2, risk.input.data[,1:13] %*% ASCVD_risk_params[1:13,"F_White"],
                           ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]==2, risk.input.data[,1:13] %*% ASCVD_risk_params[1:13,"F_Black"],
                                  ifelse(raw.input.data[,"Female"]==0 & raw.input.data[,"Race"]!=2, risk.input.data[,1:13] %*% ASCVD_risk_params[1:13,"M_White"],
                                         risk.input.data[,1:13] %*% ASCVD_risk_params[1:13,"M_Black"] )))
  
  overall_mean <- ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]!=2, ASCVD_risk_params[14,"F_White"],
                         ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]==2, ASCVD_risk_params[14,"F_Black"],
                                ifelse(raw.input.data[,"Female"]==0 & raw.input.data[,"Race"]!=2, ASCVD_risk_params[14,"M_White"],
                                       ASCVD_risk_params[14,"M_Black"] )))
  
  baseline_survival <- ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]!=2, ASCVD_risk_params[15,"F_White"], 
                              ifelse(raw.input.data[,"Female"]==1 & raw.input.data[,"Race"]==2, ASCVD_risk_params[15,"F_Black"],
                                     ifelse(raw.input.data[,"Female"]==0 & raw.input.data[,"Race"]!=2, ASCVD_risk_params[15,"M_White"],
                                            ASCVD_risk_params[15,"M_Black"] )))
  
  raw.input.data$ASCVD_Risk_10yr <- 1 - baseline_survival^exp(individual_sum-overall_mean)
  
  raw.input.data$ASCVD_Risk_10yr <- ifelse(raw.input.data[,"Race"]==2, raw.input.data$ASCVD_Risk_10yr*0.9, raw.input.data$ASCVD_Risk_10yr)
  
  return(raw.input.data$ASCVD_Risk_10yr)
}

