##8-year Diabetes Incidence Score##
##Primary Source: Wilson et al (2007, Arch Int Med) Table 6
##Author: David D. Kim, Tufts Medical Center, DKim3@tuftsmedicalcenter.org


dm_risk_params <- read.csv("00 Input Data/FHS Diabetes Risk Model Parameters.csv", row.names = 1)

calc_DM_risk <- function(raw.input.data) {
  names <- c(rownames(dm_risk_params))
  
  ###Diabetes
  DM_risk_score <- matrix(nrow=nrow(raw.input.data), ncol=1)
  colnames(DM_risk_score) <- c("Diabete risk score")
  
  #Predictor 1: glucose level
  raw.input.data$Glucose_cat<-0
  raw.input.data$Glucose_cat[raw.input.data$Glucose>=100 & raw.input.data$Glucose<126] <- 1
  
  #Predictor 2: BMI
  raw.input.data$BMI_cat1<-0
  raw.input.data$BMI_cat1[raw.input.data$BMI>=25 & raw.input.data$BMI<30]<-1
  
  raw.input.data$BMI_cat2<-0
  raw.input.data$BMI_cat2[raw.input.data$BMI>=30]<-1
  
  #Predictor 3: HDL
  raw.input.data$HDL_cat<-0
  raw.input.data$HDL_cat[raw.input.data$HDL<40 & raw.input.data$Female==0]<-1
  raw.input.data$HDL_cat[raw.input.data$HDL<50 & raw.input.data$Female==1]<-1
  
  #Predictor 4: DM Parent
  # Already NHANES variable
  
  #Predictor 5: Trig_level
  raw.input.data$Trig_cat<-0
  raw.input.data$Trig_cat[raw.input.data$Trig>=150]<-1
  
  #Predictor 6: HBP
  raw.input.data$HBP_cat<-0
  raw.input.data$HBP_cat[raw.input.data$SBP>=130 | raw.input.data$HPT_Txt == 1 | raw.input.data$DBP>=85] <- 1
  
  DM_risk_input<-raw.input.data[c("seqn",names)]
  
  ############################################
  # If value is missing, place 0
  ############################################
  DM_risk_input[is.na(DM_risk_input)] <- 0
  
  ############################################
  # Diabetes Risk Calculation
  ############################################
  DM_risk_input$DM_scores<- as.matrix(DM_risk_input[2:8]) %*% as.matrix(dm_risk_params[1:7,1])
  
  DM_risk_input$DM_risk <- 0
  DM_risk_input$DM_risk[DM_risk_input$DM_scores<=10] <- 3
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==11 | DM_risk_input$DM_scores==12] <- 4
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==13] <- 5
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==14] <- 6
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==15] <- 7
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==16] <- 9
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==17] <- 11
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==18] <- 13
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==19] <- 15
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==20] <- 18
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==21] <- 21
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==22] <- 25
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==23] <- 29
  DM_risk_input$DM_risk[DM_risk_input$DM_scores==24] <- 33
  DM_risk_input$DM_risk[DM_risk_input$DM_scores>=25] <- 35
  
  DM_risk_input$DM_risk_8yr <- DM_risk_input$DM_risk/100
  return (DM_risk_input$DM_risk_8yr)
}

