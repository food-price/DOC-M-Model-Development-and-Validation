##Predicting HRQOL Among the non-institutionalized US population##
##Primary Source: Lubetkin et al (2005, Quality of Life Research) Table 2
##Author: David D. Kim, Tufts Medical Center, DKim3@tuftsmedicalcenter.org

###HRQOL
HRQOL_parameters <- read.csv("00 Input Data/HrQOL Estimation Parameters.csv", row.names = 1)

calc_HRQOL <- function(raw.input.data, params) {
  # Input: NHANES data for analysis and HRQOL parameters
  names <- c(rownames(HRQOL_parameters))
  
  #Age Categories
  raw.input.data$Age_18_39 <- 0
  raw.input.data$Age_18_39[raw.input.data$Age>=18 & raw.input.data$Age<40] <- 1
  
  raw.input.data$Age_40_59 <- 0
  raw.input.data$Age_40_59[raw.input.data$Age>=40 & raw.input.data$Age<60] <- 1
  
  raw.input.data$Age_60_69 <- 0
  raw.input.data$Age_60_69[raw.input.data$Age>=60 & raw.input.data$Age<70] <- 1
  
  raw.input.data$Age_over70 <- 0
  raw.input.data$Age_over70[raw.input.data$Age>=70] <- 1
  
  #Gender
  raw.input.data$Male <- 0
  raw.input.data$Male[raw.input.data$Female==0] <- 1
  
  #Race
  raw.input.data$White<-0
  raw.input.data$White[raw.input.data$Race==1] <- 1
  
  raw.input.data$Black<-0
  raw.input.data$Black[raw.input.data$Race==2] <- 1
  
  raw.input.data$Hispanic<-0
  raw.input.data$Hispanic[raw.input.data$Race==3] <- 1
  
  raw.input.data$HBP<-0
  raw.input.data$HBP[raw.input.data$SBP>=130 | raw.input.data$HPT_Txt == 1 | raw.input.data$DBP>=85] <- 1
  
  HRQOL_input<-raw.input.data[c("seqn",names)]
  
  ############################################
  # If value is missing, place 0
  ############################################
  HRQOL_input[is.na(HRQOL_input)] <- 0
  
  HRQOL_baseline_mean <- 0.901 # Lubetkin et al (2005, Quality of Life Research) Table 1 for Age 18-39
  HRQOL_baseline_SE <- 0.026
  #Baseline HrQOL of 0.9 is based on Average EQ-5D index for age 18-39 in the same study
  #Theoretically, the baseline HRQOL should be based on age 18-39, female, white without any clinical conditions. 
  #without such detail baseline information, we use the average HRQOL among age 18-39 as a baseline HRQOL. 
  
  ############################################
  # Diabetes Risk Calculation
  ############################################
  HRQOL_input$HRQOL_scores<- HRQOL_baseline_mean + (as.matrix(HRQOL_input[2:13]) %*% as.matrix(params)) 
  
  return(HRQOL_input$HRQOL_scores)
}



