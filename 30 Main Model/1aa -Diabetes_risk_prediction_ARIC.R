##9-year Diabetes Risk Calculator: ARIC .##
##Primary Source: by Lacy et al and Mann et al

##Risk Equation Parameters

dm_risk_params2 <- read.csv("00 Input Data/ARIC_Diabetes_prediction_model.csv", row.names = 1)

calc_DM_risk_ARIC1 <- function(raw.input.data) {
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(dm_risk_params2))
  colnames(risk.input.data) <- rownames(dm_risk_params2)
  
  risk.input.data[,1] = 1.0
  risk.input.data[,2] = raw.input.data[,"Age"]
  risk.input.data[,3] = raw.input.data[,"Black"]
  risk.input.data[,4] = raw.input.data[,"DM_family"]
  risk.input.data[,5] = raw.input.data[,"Glucose"]
  risk.input.data[,6] = raw.input.data[,"SBP"]
  risk.input.data[,7] = raw.input.data[,"WC"]
  risk.input.data[,8] = raw.input.data[,"Height"]
  risk.input.data[,9] = raw.input.data[,"HDL"]
  risk.input.data[,10] = raw.input.data[,"Trig"]
  #Treating all missing information as zero
  risk.input.data[is.na(risk.input.data)] <- 0
  
  ARIC_Lacy <-  risk.input.data[,1:10] %*% as.matrix(dm_risk_params2[1:10,"Lacy"])

  raw.input.data$DM_ARIC1 = (exp(ARIC_Lacy)/(1+exp(ARIC_Lacy)))
  return(raw.input.data$DM_ARIC1)
}

calc_DM_risk_ARIC2 <- function(raw.input.data) {
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(dm_risk_params2))
  colnames(risk.input.data) <- rownames(dm_risk_params2)
  risk.input.data[,1] = 1.0
  risk.input.data[,2] = raw.input.data[,"Age"]
  risk.input.data[,3] = raw.input.data[,"Black"]
  risk.input.data[,4] = raw.input.data[,"DM_parent"]
  risk.input.data[,5] = raw.input.data[,"Glucose"]
  risk.input.data[,6] = raw.input.data[,"SBP"]
  risk.input.data[,7] = raw.input.data[,"WC"]
  risk.input.data[,8] = raw.input.data[,"Height"]
  risk.input.data[,9] = raw.input.data[,"HDL"]
  risk.input.data[,10] = raw.input.data[,"Trig"]
  #Treating all missing information as zero
  risk.input.data[is.na(risk.input.data)] <- 0
  
  ARIC_Mann <-  risk.input.data[,1:10] %*% dm_risk_params2[1:10,"Mann"]
  raw.input.data$DM_ARIC2 = (exp(ARIC_Mann)/(1+exp(ARIC_Mann)))
  
  return(raw.input.data$DM_ARIC2)
}

dm_risk_params4 <- read.csv("00 Input Data/ARIC_Diabetes_prediction_model_1.csv", row.names = 1)

calc_DM_risk_ARIC3 <- function(raw.input.data) {
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(dm_risk_params4))
  colnames(risk.input.data) <- rownames(dm_risk_params4)
  
  risk.input.data[,1] = 1.0
  risk.input.data[,2] = raw.input.data[,"Age"]
  risk.input.data[,3] = raw.input.data[,"Black"]
  risk.input.data[,4] = raw.input.data[,"DM_parent"]
  risk.input.data[,5] = raw.input.data[,"SBP"]
  risk.input.data[,6] = raw.input.data[,"WC"]
  risk.input.data[,7] = raw.input.data[,"Height"]

  #Treating all missing information as zero
  risk.input.data[is.na(risk.input.data)] <- 0
  
  ARIC_Lacy <-  risk.input.data[,1:7] %*% as.matrix(dm_risk_params4[1:7,"Parm"])
  
  raw.input.data$DM_ARIC1 = (exp(ARIC_Lacy)/(1+exp(ARIC_Lacy)))
  return(raw.input.data$DM_ARIC1)
}

dm_risk_params5 <- read.csv("00 Input Data/ARIC_Diabetes_prediction_model_2.csv", row.names = 1)

calc_DM_risk_ARIC4 <- function(raw.input.data) {
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(dm_risk_params5))
  colnames(risk.input.data) <- rownames(dm_risk_params5)
  
  risk.input.data[,1] = 1.0
  risk.input.data[,2] = raw.input.data[,"Age"]
  risk.input.data[,3] = raw.input.data[,"Black"]
  risk.input.data[,4] = raw.input.data[,"DM_parent"]
  risk.input.data[,5] = raw.input.data[,"Glucose"]
  risk.input.data[,6] = raw.input.data[,"SBP"]
  risk.input.data[,7] = raw.input.data[,"WC"]
  risk.input.data[,8] = raw.input.data[,"Height"]
  
  #Treating all missing information as zero
  risk.input.data[is.na(risk.input.data)] <- 0
  
  ARIC_Lacy <-  risk.input.data[,1:8] %*% as.matrix(dm_risk_params5[1:8,"Parm"])
  
  raw.input.data$DM_ARIC1 = (exp(ARIC_Lacy)/(1+exp(ARIC_Lacy)))
  return(raw.input.data$DM_ARIC1)
}


