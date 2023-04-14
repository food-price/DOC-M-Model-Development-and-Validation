#Updating Disease history status
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



run_sim <- function(s, intervention, timehz) {
 

###################################
####1. set up initial values########
###################################
  # Returns array with variables for each individual after n.cycle
  
  # Define inpute variables
  Input_variables <- c("SEQN","Subject_ID","DEMO", "Age", "Total_Chol","HDL", "SBP", "DBP", "BMI",  "Obesity", "weight", "height","Trig", "Glucose", 
                       "Diabetes", "CVD_history","DM_prob", "CVD_prob", "CVD_recurrent_prob","HRQOL_scores", "HCE_predict")

  #Set up initial values
  #replicate each individual in the data for n.loop times, to minimize stochastic error
  sim_out_t<-do.call("rbind", replicate(n.loop, data_for_analysis[,Input_variables], simplify = FALSE))
  sim_out_t$state=rep(data_for_analysis$initial_H, n.loop)
  sim_out_t[,"Age_cycle"]=sim_out_t[,"Age"]
  sim_out_t[,"BaseBMI"]=sim_out_t[,"BMI"]
  sim_out_t[,"HCE_disc"] <- sim_out_t$HCE_predict
  sim_out_t[,"effect_disc"] <- sim_out_t$HRQOL_scores
  sim_out_t[,"Life Years"]=0
  sim_out_t[,"Incident First CVD"]=0
  sim_out_t[,"Incident Recurrent CVD"]=0
  sim_out_t[,"Incident CVD"]=0
  sim_out_t[,"All_Deaths"] <- 0
  sim_out_t[,"CVD_Death"]<- 0
  sim_out_t[,"DM_Death"]<- 0
  sim_out_t[,"Non_CMD_Death"]<- 0
  
  if (timehz==1){
  # Initialize array to hold output from single simulation
  sim_out<-array(NA, dim=c(nrow(data_for_analysis), length(Out_variables), n.cycle+1),
                 dimnames = list(data_for_analysis$Subject_ID,
                                Out_variables, paste("cycle", 0:n.cycle, sep = " ")))
  sim_out[,,1]<-as.matrix(sim_out_t[Out_variables])
  } 
########################################
######2. run model for n.cycle###############
########################################  
    
  for (t in 1:n.cycle) {
    print(t)
    ###################################################################################  
    # 2.1 Updating Risk Factors over time: Applying secular trends based on age, gender, R/E
    ###################################################################################  
    
    #Time-varying data inputs
    
    sim_out_t[,"Age_cycle"] <- sim_out_t[,"Age"] + t
    ##update BMI##

  sim_out_t[,"BMI"] <- update_risk_factor(sim_out_t, "BMI", "BMI", s)

    # Truncate BMI to the interval 12-70 kg/m2. Source: ???	https://doi.org/10.2105/AJPH.2008.137364  
   sim_out_t[,"BMI"] <- ifelse(sim_out_t[,"BMI"] < 12, 12,
                                ifelse(sim_out_t[,"BMI"]>70, 70, sim_out_t[,"BMI"]))
   sim_out_t[,"Total_Chol"] <- update_risk_factor(sim_out_t, "Total_Chol", "Total_Chol", s)
   sim_out_t[,"HDL"] <- update_risk_factor(sim_out_t, "HDL", "HDL", s)
   sim_out_t[,"SBP"] <- update_risk_factor(sim_out_t, "SBP", "SBP", s)
   sim_out_t[,"DBP"] <- update_risk_factor(sim_out_t, "DBP", "SBP", s)
   sim_out_t[,"Trig"] <- update_risk_factor(sim_out_t, "Trig", "Trig", s)
   sim_out_t[,"Glucose"] <- update_risk_factor(sim_out_t, "Glucose", "Glucose", s)
    
   
   ###################################################################################  
   # 2.2 Updating RRs of changing risk factors
   ###################################################################################  
   
    # this is study specific, changed to one for model validation 
   
    RR_diff_total_CHD <- 1
    RR_diff_total_DM <- 1
    RR_diff_total_Stroke <- 1

    ###################################################################################  
    # 2.3 Updating disease risks###########
    ###################################################################################  
    
    vars<-c("SEQN", "Age_cycle","Glucose","BMI","Total_Chol","HDL","Trig","SBP","DBP","Diabetes","CVD_history","Age")
    raw.input.data<-cbind(do.call("rbind", replicate(n.loop, data_for_analysis[,c("Female", "Race", "DM_parent",  "HPT_Txt",  "Smoking")], simplify = FALSE)),sim_out_t[,vars])
    raw.input.data$Age<-raw.input.data$Age_cycle
    
    if (t%/%2 > 0 & t%%2 == 0)  {
      CVD_Recurrent_risk_2yr <- calc_recur_CVD_risk(raw.input.data)
     sim_out_t[,"CVD_recurrent_prob"] <- Multi_yr_Risk_to_annual_prob(time=2, risk=CVD_Recurrent_risk_2yr)
    } else {
     sim_out_t[,"CVD_recurrent_prob"] <-sim_out_t[,"CVD_recurrent_prob"]
    }
    
    if (t%/%8 > 0 & t%%8 == 0)  {
      DM_risk_8yr <- calc_DM_risk(raw.input.data)
     sim_out_t[,"DM_prob"] <- Multi_yr_Risk_to_annual_prob(time=8, risk=DM_risk_8yr)
    } else {
     sim_out_t[,"DM_prob"] <-sim_out_t[,"DM_prob"]
    }
    
    if (t%/%10 > 0 & t%%10 == 0) {
      ASCVD_Risk_10yr <- calc_ASCVD_risk(raw.input.data)
     sim_out_t[,"CVD_prob"] <- Multi_yr_Risk_to_annual_prob(time=10, risk=ASCVD_Risk_10yr)
    } else {
     sim_out_t[,"CVD_prob"] <-sim_out_t[,"CVD_prob"]
    }
    
    ###################################################################################  
    # 2.4 Updating transition probabilities###########
    ###################################################################################    
    #Disaggregating ASCVD
    
    p.CHD_first <- rep(data_for_analysis$p_CHD_first, n.loop) #Prob. of CHD event (Angina, MI, Fatal MI, Fatal CHD) Source: Benjamin et al, Circulation, 2018 Table 13-1 & 18-1, and 18-2
    p.CHD_recurrent <-  rep(data_for_analysis$p_CHD_recurrent, n.loop)
    
    #Defining transition probabilities
    # Subtract 19 from age because indexing starts at age 20
    p.death.DM <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ DM_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ DM_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ DM_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ DM_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ DM_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ DM_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death.CHD <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ CHD_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ CHD_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ CHD_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ CHD_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ CHD_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ CHD_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ CHD_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ CHD_mortality_HM[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death.Stroke <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ stroke_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ stroke_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ stroke_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ stroke_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ stroke_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ stroke_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ stroke_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ stroke_mortality_HF[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    p.death <- case_when(
      sim_out_t[,"DEMO"] == "Male" ~ Non_CVD_DM_mortality_Male[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "Female" ~ Non_CVD_DM_mortality_Female[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWM" ~ Non_CVD_DM_mortality_NHWM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHWF" ~ Non_CVD_DM_mortality_NHWF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBM" ~ Non_CVD_DM_mortality_NHBM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "NHBF" ~ Non_CVD_DM_mortality_NHBF[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HM" ~ Non_CVD_DM_mortality_HM[sim_out_t[,"Age_cycle"]-19,s],
      sim_out_t[,"DEMO"] == "HF" ~ Non_CVD_DM_mortality_HF[sim_out_t[,"Age_cycle"]-19,s]
    )
    
    #Markov State #1: "No CVD, No Diabetes"
    p.H.2.DM <- as.numeric(sim_out_t[,"DM_prob"])*rep(data_for_analysis$risk_adjustment.DM, n.loop)
    p.H.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out_t[,"CVD_prob"])
    p.H.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out_t[,"CVD_prob"])
    p.H.2.H <- ifelse((p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death) > 1, 0,
                      1 - (p.H.2.DM + p.H.2.initial_CHD + p.H.2.initial_Stroke + p.death))
    
    #Markov State #2: "No CVD, With Diabetes
    p.DM.2.death <- (1-exp(-(p.death+p.death.DM)))
    p.DM.2.initial_CHD <- RR_diff_total_CHD*p.CHD_first*as.numeric(sim_out_t[,"CVD_prob"])
    p.DM.2.initial_Stroke <- RR_diff_total_Stroke*(1-p.CHD_first)*as.numeric(sim_out_t[,"CVD_prob"])
    p.DM.2.DM <- as.numeric(ifelse((p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death) > 1, 0,
                                   1 - (p.DM.2.initial_CHD + p.DM.2.initial_Stroke + p.DM.2.death)))
    
    #Markov State #3: "First Stroke"
    p.initial_Stroke.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke)))))
    p.initial_Stroke.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_Stroke.2.death))
    p.initial_Stroke.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_Stroke.2.death, 0))
    
    #Markov State #4: "First CHD w/o RVSC" 
    p.initial_CHD_No_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD)))))
    p.initial_CHD_No_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_No_RVSC.2.death))
    p.initial_CHD_No_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_No_RVSC.2.death, 0))
    
    #Markov State #5: "First CHD with RVSC"
    p.initial_CHD_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s])))))
    p.initial_CHD_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.initial_CHD_RVSC.2.death))
    p.initial_CHD_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.initial_CHD_RVSC.2.death, 0))
    
    #Markov State #6: "CVD History, No Diabetes"
    p.CVD.2.DM <- as.numeric(sim_out_t[,"DM_prob"])*rep(data_for_analysis$risk_adjustment.DM, n.loop)
    p.CVD.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke)))
    p.CVD.2.CVD <- as.numeric(ifelse((p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death) > 1, 0,
                                     1 - (p.CVD.2.DM + p.CVD.2.Sub_CHD + p.CVD.2.Sub_Stroke + p.CVD.2.death)))
    
    #Markov State #7: "CVD History, With Diabetes"
    p.CVD_DM.2.Sub_CHD <- RR_diff_total_CHD*p.CHD_recurrent*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD_DM.2.Sub_Stroke <- RR_diff_total_Stroke*(1-p.CHD_recurrent)*as.numeric(sim_out_t[,"CVD_recurrent_prob"])
    p.CVD_DM.2.death <- (1-exp(-(p.death+p.death.CHD+p.death.Stroke+p.death.DM)))
    p.CVD_DM.2.CVD_DM <- as.numeric(ifelse((p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke) > 1, 0,
                                           1 - (p.CVD_DM.2.death + p.CVD_DM.2.Sub_CHD + p.CVD_DM.2.Sub_Stroke)))
    
    #Markov State #8: "Subsequent Stroke"
    p.sub_Stroke.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.Stroke+p.death.DM))), (1-exp(-(p.death+p.death.Stroke)))))
    p.sub_Stroke.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_Stroke.2.death))
    p.sub_Stroke.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_Stroke.2.death, 0))
    
    #Markov State #9: "Subsequent CHD w/o RVSC" 
    p.sub_CHD_No_RVSC.2.death <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM))), (1-exp(-(p.death+p.death.CHD)))))
    p.sub_CHD_No_RVSC.2.CVD_No_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_No_RVSC.2.death))
    p.sub_CHD_No_RVSC.2.CVD_DM <- as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_No_RVSC.2.death, 0))
    
    #Markov State #10: "Subsequent CHD with RVSC"
    p.sub_CHD_RVSC.2.death <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, (1-exp(-(p.death+p.death.CHD+p.death.DM+p.death.RVSC_sim[s]))), (1-exp(-(p.death+p.death.CHD+p.death.RVSC_sim[s])))))
    p.sub_CHD_RVSC.2.CVD_No_DM <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 0, 1 - p.sub_CHD_RVSC.2.death))
    p.sub_CHD_RVSC.2.CVD_DM <-  as.numeric(ifelse(as.numeric(sim_out_t[,"Diabetes"]) ==1, 1 - p.sub_CHD_RVSC.2.death, 0))
    
    ###################################################################################  
    ####2.5 Assign transition probablities to the transition matrix ##########
    ###################################################################################
    p.transition <- array(NA, dim=c(n.individual, n.health.state),
                          dimnames =  list(c(1:n.individual),name.health.state))
    p.transition[,"No CVD, No Diabetes"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.H, rep(0, n.individual))
    p.transition[,"No CVD, With Diabetes"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.DM, 
                                                     ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.DM.2.DM, rep(0, n.individual)))
    p.transition[,"First Stroke"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.H.2.initial_Stroke, 
                                            ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.DM.2.initial_Stroke, rep(0, n.individual)))
    p.transition[,"First CHD w/o RVSC"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", (1-p.RVSC)*p.H.2.initial_CHD, 
                                                  ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", (1-p.RVSC)*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"First CHD with RVSC"] <- ifelse(sim_out_t[,"state"] == "No CVD, No Diabetes", p.RVSC*p.H.2.initial_CHD, 
                                                   ifelse(sim_out_t[,"state"] == "No CVD, With Diabetes", p.RVSC*p.DM.2.initial_CHD, rep(0, n.individual)))
    p.transition[,"CVD History, No Diabetes"] <- case_when(
      sim_out_t[,"state"] == "First Stroke" ~ p.initial_Stroke.2.CVD_No_DM, 
      sim_out_t[,"state"] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "CVD History, No Diabetes" ~ p.CVD.2.CVD, 
      sim_out_t[,"state"] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_No_DM, 
      sim_out_t[,"state"] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_No_DM, 
      sim_out_t[,"state"] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_No_DM, 
      TRUE ~ rep(0, n.individual)
    )
    p.transition[,"CVD History, With Diabetes"] <- case_when(
      sim_out_t[,"state"] == "First Stroke" ~ p.initial_Stroke.2.CVD_DM, 
      sim_out_t[,"state"] == "First CHD w/o RVSC" ~ p.initial_CHD_No_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "First CHD with RVSC" ~ p.initial_CHD_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "CVD History, No Diabetes" ~ p.CVD.2.DM, 
      sim_out_t[,"state"] == "CVD History, With Diabetes" ~ p.CVD_DM.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent Stroke" ~ p.sub_Stroke.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent CHD w/o RVSC" ~ p.sub_CHD_No_RVSC.2.CVD_DM, 
      sim_out_t[,"state"] == "Subsequent CHD with RVSC" ~ p.sub_CHD_RVSC.2.CVD_DM, 
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Subsequent Stroke"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", p.CVD.2.Sub_Stroke,
                                                 ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", p.CVD_DM.2.Sub_Stroke, rep(0, n.individual)))
    p.transition[,"Subsequent CHD w/o RVSC"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", (1-p.RVSC)*p.CVD.2.Sub_CHD, 
                                                       ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", (1-p.RVSC)*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    p.transition[,"Subsequent CHD with RVSC"] <- ifelse(sim_out_t[,"state"] == "CVD History, No Diabetes", p.RVSC*p.CVD.2.Sub_CHD, 
                                                        ifelse(sim_out_t[,"state"] == "CVD History, With Diabetes", p.RVSC*p.CVD_DM.2.Sub_CHD, rep(0, n.individual)))
    
    p.transition[,"DM_Death"] <- case_when(
     sim_out_t[,"state"] == "DM_Death" ~ rep(1, n.individual),
     sim_out_t[,"state"] %in% c("Non_DM_Non_CVD_Death", "Stroke_Death", "CHD_Death") ~ rep(0, n.individual),
      sim_out[,"Diabetes",t] == 1 ~ p.death.DM, 
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Stroke_Death"] <- case_when(
     sim_out_t[,"state"] %in% 
        c("First Stroke", "Subsequent Stroke", "CVD History, No Diabetes", "CVD History, With Diabetes") ~ p.death.Stroke, 
     sim_out_t[,"state"] == "Stroke_Death" ~ rep(1, n.individual),
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"CHD_Death"] <- case_when(
     sim_out_t[,"state"] %in%
        c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes") ~ p.death.CHD,
     sim_out_t[,"state"] %in% c("First CHD with RVSC", "Subsequent CHD with RVSC") ~ p.death.CHD+p.death.RVSC_sim[s],
     sim_out_t[,"state"] == "CHD_Death" ~ rep(1, n.individual),
      TRUE ~ rep(0, n.individual)
    )
    
    p.transition[,"Non_DM_Non_CVD_Death"] <- case_when(
     sim_out_t[,"state"] %in% c("DM_Death", "Stroke_Death", "CHD_Death") ~ rep(0, n.individual),
     sim_out_t[,"state"] == "Non_DM_Non_CVD_Death" ~ rep(1, n.individual),
      TRUE ~ p.death
    )
    
    
    # Check that all transition probabilities add to 1 (rounded to 3 digits)
    # if (sum(round(rowSums(p.transition),3)!=1) != 0) {
    #   p_sums = round(rowSums(p.transition),3)
    #   error_out = sim_out[p_sums!=1,"state",t] # Output state of person(s) with error
    #   stop("Transition probabilities do not add to 1. ", paste("Simulation", s, ", Time", t, ". "), "State(s) with error: ", error_out)
    # }
    
 
    ###################################################################################  
    ####2.6 # Transition to the next health state ##########
    ################################################################################### 
    
   sim_out_t[,"state"]<- apply(p.transition, 1, function(x) sample(name.health.state, 1, prob = x))
    
   
    ###################################################################################  
    ####2.7 Update health status, disease outcome, QALY, and economic outcomes##########
    ################################################################################### 
    
    
   sim_out_t[,"Diabetes"] <- ifelse(sim_out_t[,"state"]%in% c("No CVD, With Diabetes", "CVD History, With Diabetes"), 1, 
                                     ifelse(sim_out_t[,"state"] %in% name.death.states, NA, 
                                            ifelse(rep(data_for_analysis$Glucose, n.loop) > 126, 1,sim_out_t[,"Diabetes"])))
    
   sim_out_t[,"CVD_history"] <- ifelse( sim_out_t[,"state"] %in% c("First Stroke", "First CHD w/o RVSC", "First CHD with RVSC", "CVD History, No Diabetes", "CVD History, With Diabetes",
                                                                   "Subsequent Stroke", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), 1, 
                                                ifelse(sim_out_t[,"state"] %in% name.death.states, NA, sim_out_t[,"CVD_history"]))
    
   
   
   sim_out_t[,"BMI"] <- ifelse(sim_out_t[,"state"] %in% name.death.states, NA ,sim_out_t[,"BMI"])
   sim_out_t[,"Obesity"] <- ifelse(sim_out_t[,"BMI"] < 30, 0, 
                                    ifelse(sim_out_t[,"BMI"] >= 30, 1, sim_out_t[,"Obesity"]))
    
    
    raw.input.data[,vars]<-sim_out_t[,vars]
       
    # Update QALYs and costs
    HRQOL_scores_t1 <- calc_HRQOL(raw.input.data, HRQOL_parameter_sim[,s])
    HRQOL_scores_t1 <- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"), HRQOL_scores_t1+u_stroke_sim[s], 
                                           ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "First CHD with RVSC", "Subsequent CHD w/o RVSC", "Subsequent CHD with RVSC"), HRQOL_scores_t1+u_CHD_sim[s],
                                                  ifelse(sim_out_t[,"state"] %in% name.death.states, 0, HRQOL_scores_t1)))
    
    sim_out_t[,"HRQOL_scores"] =sim_out_t[,"HRQOL_scores"]+HRQOL_scores_t1
    sim_out_t[,"effect_disc"] <- sim_out_t[,"effect_disc"]+HRQOL_scores_t1/((1+beta_QALY)^(t-1))
    

    HCE_predict_t1 <- calc_HCE(raw.input.data, HCE_parameter_sim[,s])
    HCE_predict_t1 <- ifelse(sim_out_t[,"state"]%in% c("First Stroke", "Subsequent Stroke"), HCE_predict_t1 +c_stroke_sim[s] , 
                                          ifelse(sim_out_t[,"state"]%in% c("First CHD w/o RVSC", "Subsequent CHD w/o RVSC"), HCE_predict_t1+c_CHD_sim[s], 
                                                 ifelse(sim_out_t[,"state"]%in% c("First CHD with RVSC", "Subsequent CHD with RVSC"), HCE_predict_t1+c_CHD_sim[s]+c_RVSC_sim[s],
                                                        ifelse(sim_out_t[,"state"] %in% name.death.states, 0, HCE_predict_t1))))
    sim_out_t[,"HCE_predict"]=sim_out_t[,"HCE_predict"]+HCE_predict_t1
    sim_out_t[,"HCE_disc"] <- sim_out_t[,"HCE_disc"]+HCE_predict_t1/((1+beta_cost)^(t-1))
    
    sim_out_t[,"Life Years"]<-sim_out_t[,"Life Years"]+(!sim_out_t[,"state"] %in% name.death.states)
    
    sim_out_t[,"Incident First CVD"] = sim_out_t[,"Incident First CVD"]+ ifelse(sim_out_t[,"state"]=="First Stroke"| sim_out_t[,"state"]=="First CHD with RVSC" |sim_out_t[,"state"]=="First CHD w/o RVSC",1, 0)
    
    sim_out_t[,"Incident Recurrent CVD"] = sim_out_t[,"Incident Recurrent CVD"]+ ifelse(sim_out_t[,"state"]=="Subsequent Stroke"|sim_out_t[,"state"]=="Subsequent CHD with RVSC"|sim_out_t[,"state"]=="Subsequent CHD w/o RVSC",1,0)
    
    sim_out_t[,"Incident CVD"] = sim_out_t[,"Incident Recurrent CVD"]+sim_out_t[,"Incident First CVD"]
    
    sim_out_t[,"All_Deaths"] <-  ifelse(sim_out_t[,"state"] %in% c("Non_DM_Non_CVD_Death", "CHD_Death", "Stroke_Death", "DM_Death"), 1,0)
    
    sim_out_t[,"CVD_Death"] <-  ifelse(sim_out_t[,"state"]  %in% c("CHD_Death", "Stroke_Death" ), 1,0)
  
    sim_out_t[,"DM_Death"] <-  ifelse(sim_out_t[,"state"] =="DM_Death", 1,0)
    
    sim_out_t[,"Non_CMD_Death"] <-  ifelse(sim_out_t[,"state"]=="Non_DM_Non_CVD_Death", 1,0)
    
    
    sim_out_mean<-sim_out_t[ ,c(Out_variables)]  %>% 
      group_by(Subject_ID)%>% 
      summarize_at(Out_variables[-1], mean, na.omit=TRUE)
    if (timehz==1) {
    sim_out[,,t+1]<-as.matrix(sim_out_mean)
    } else sim_out=sim_out_mean
   
  }

###output specific for model validation      
          
  options(survey.lonely.psu = "adjust")
  design <- svydesign(id=~sdmvpsu, strata=~sdmvstra, weights=~wtmec2yr, nest=TRUE, data=data_for_analysis)
  design_NHW <- subset(design, (DEMO=="NHWM" | DEMO=="NHWF"))
  design_NHB <- subset(design, (DEMO=="NHBM" | DEMO=="NHBF"))
  design_Hisp <- subset(design, (DEMO=="HM" | DEMO=="HF"))
  
  vars <- c("Obesity", "Diabetes", "CVD_history", "All_Deaths","DM_Death", "CVD_Death", "Non_CMD_Death")
  races <- c("all", "NHW", "NHB", "Hisp")
  stat<-c("mean", "var")
  mean_out <- array(NA, dim = c(length(vars)+1, (n.cycle+1)*4, 2), 
                    dimnames = list(c(vars,"race"), rep(paste("cycle", 0:n.cycle, sep = " "), 4), stat))
  
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
        
        for (var in vars) {
          mean_out[var,t2, "mean"] <- svymean(as.numeric(sim_out[,var,t]), design, deff=F, na.rm=T)
          mean_out [var,t2, "var"] <- SE(svymean(as.numeric(sim_out[,var,t]), design, deff=F, na.rm=T))^2
        }
        mean_out["race",t2, "mean"] <- race
        mean_out["race",t2, "var"] <- race
 
      } else {
        
        for (var in vars) {
          mean_out[var,t2, "mean"] <- svymean(as.numeric(sim_out[which(data_for_analysis$Race==race_num),var,t]), race_design, deff=F, na.rm=T)
          mean_out[var,t2, "var"] <- SE(svymean(as.numeric(sim_out[which(data_for_analysis$Race==race_num),var,t]), race_design, deff=F, na.rm=T))^2
        }
        mean_out["race",t2, "mean"] <- race
        mean_out["race",t2, "var"] <- race
     }
   }
  }
  return(mean_out)
}

  
    