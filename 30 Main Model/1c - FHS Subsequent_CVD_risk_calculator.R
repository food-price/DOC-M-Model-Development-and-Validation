##2-year Subsequent CHD Risk Calculator##
##Primary Source: D'Agostino et al. (Am Heart J, 2000) TABLE IV

##Risk Equation Parameters
recur_CVD_risk_params <- read.csv("00 Input Data/FHS Recurrent CVD Risk Equation Parameters.csv", row.names = 1)

calc_recur_CVD_risk <- function(raw.input.data) {
  risk.input.data <- matrix(nrow=nrow(raw.input.data), ncol=nrow(recur_CVD_risk_params))
  colnames(risk.input.data) <- rownames(recur_CVD_risk_params)
  
  risk_score <- matrix(nrow=nrow(raw.input.data), ncol=4)
  colnames(risk_score) <- c("m_parameter", "extreme_value", "u_parameter","CVD_Recurrent_risk_2yr")
  
  risk.input.data[,1] = 1.0
  risk.input.data[,2] = raw.input.data[,"Age"]
  #risk.input.data[,2] = (log(raw.input.data[,"Age"]))^2
  risk.input.data[,3] = log((raw.input.data[,"Total_Chol"])/(raw.input.data[,"HDL"]))
  risk.input.data[,4] = log(raw.input.data[,"SBP"])
  risk.input.data[,5] = raw.input.data[,"Diabetes"]
  risk.input.data[,6] = raw.input.data[,"Smoking"]
  
  #Treating all missing information as zero
  risk.input.data[is.na(risk.input.data)] <- 0
  
  m_parameter <- ifelse(raw.input.data[,"Female"]==1, risk.input.data[,1:6] %*% recur_CVD_risk_params[1:6,"Women"], risk.input.data[,1:6] %*% recur_CVD_risk_params[1:6,"Men"])
  extreme_value <- ifelse(raw.input.data[,"Female"]==1, recur_CVD_risk_params[7,"Women"], recur_CVD_risk_params[7,"Men"])
  u_parameter = (log(2)-m_parameter)/extreme_value
  raw.input.data$CVD_Recurrent_risk_2yr = 1-exp(-exp(u_parameter))
  
  return(raw.input.data$CVD_Recurrent_risk_2yr)
}

