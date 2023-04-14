##Predicting HCE among the non-institutionalized US Population##
##Source: Author's analysis of the 2014-2016 MEPS data
##Author: David D. Kim, Tufts Medical Center, DKim3@tuftsmedicalcenter.org

HCE_parameters <- read.csv("01_Input/HCE Estimation Parameters.csv", row.names = 1)

calc_HCE <- function(raw.input.data, params) {
  # Input: NHANES data for analysis and HCE parameters
  names <- c(rownames(HCE_parameters))
  
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
  
  HCE_input<-raw.input.data[c("SEQN",names)]
  HCE_input$Age <- HCE_input$Age-40
  
  ############################################
  # If value is missing, place 0
  ############################################
  HCE_input[is.na(HCE_input)] <- 0
  
  HCE_baseline_mean <- 2878.84 # Baseline HCE were estimated for individuals of age 40, BMI 28, male, non-Hispanic White, without any clinical condition
  HCE_baseline_SE <- 103.5559
  
  ############################################
  # Diabetes Risk Calculation
  ############################################
  HCE_input$HCE_predict<- HCE_baseline_mean + (as.matrix(HCE_input[2:11]) %*% as.matrix(params)) 
}

