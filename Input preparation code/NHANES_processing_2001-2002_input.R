# Script to process NHANES data as model input file
# Combines separate NHANES data sets
# Cleans data
# Imputes data using multiple imputation
# Subsets population of interest

# Adapted by Lu Wang (brianna.lauren@tufts.edu) from NHANES processng codes initially created by Brianna Lauren (brianna.lauren@tufts.edu)
# Date: Jan 26, 2021
# Last revision: June 24, 2021


##### LIBRARIES ####################################################################################
library(dplyr)
library(purrr)
library(foreign)
library(Hmisc) # multiple imputation
library(survey)

##### GENERAL SETTINGS #############################################################################

setwd("C:/Users/lwang18/Documents/Github/ODC-M_Validation/00 Input Data/NHANES")

# File paths for saving intermediate and final data sets
raw_filepath = "./NHANES_0118_raw.csv"
preImp_filepath = "./NHANES_0102_preImp.csv"
Imp_filepath = "./NHANES_0102_Imp.csv"

##### MERGE AND APPEND NHANES DATA #################################################################



# Define NHANES files and variables that are needed
# Excludes drug (RXQ_RX) and dietary data because handled separately (multiple rows per person)

# Cycles of NHANES to include
cycles <- c("2001-2002")
cycle_letters <- c("B")
files <- c("DEMO", "DIQ", "BMX","BPQ","BPX", "SMQ", "MCQ", "CDQ", "L10_2", "L13_2", "FSQ")

nhanes_vars <- c("SEQN", "SDDSRVYR", "WTINT2YR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "RIDAGEYR", 
                 "RIAGENDR", "RIDRETH1", "RIDRETH3", "DMDEDUC2", "INDFMPIR", "LBXTC", "LBXTR", "LB2TC", "LB2TR", "LB2HDL",
                 "LBDHDD", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4", "BPXDI1", "BPXDI2", "BPXDI3", 
                 "BPXDI4", "BPQ020", "BPQ040A", "SMQ020", "SMQ040", "MCQ300C", "BMXBMI", "BMXHT", 
                 "BMXWT", "DIQ010", "LB2GH","LB2GLU", "LBXGH", "LBXGLU", "LBXGLT",  
                 "MCQ160B", "MCQ160C", "MCQ160E", "MCQ160F", "CDQ001", "CDQ002", "CDQ003", "CDQ004", 
                 "CDQ005", "CDQ006", "CDQ009D", "CDQ009E", "CDQ009F", "CDQ009G", 
                 "MCQ260AA", "MCQ260AB", "MCQ260AC", "MCQ260AD", "MCQ260AE", "MCQ260AF", "MCQ260AG", "MCQ260AH", "MCQ260AI")

read_data <- function(cycle_i, file_base) {
  cycle = cycles[cycle_i]
  cycle_letter = cycle_letters[cycle_i]
  file = paste(file_base, "_", cycle_letter, ".XPT", sep = "")
  data = read.xport(paste('.', cycle, file, sep = '/'))
  return(data)
}

data_cycle = lapply(files, read_data, cycle_i = 1)

# Merge each file for a single cycle (left merge)
cycle_df = reduce(data_cycle, merge, by = "SEQN", all.x = T)
# Select variables
nhanes_raw = select(cycle_df, any_of(nhanes_vars))

nhanes_raw <- nhanes_raw  %>%  rename(LBXGH=LB2GH, LBXGLU=LB2GLU, LBXTC=LB2TC, LBXTR=LB2TR, LBDHDD=LB2HDL) 
nhanes_raw$LBXGLT= NA




##### PREPARING FOR IMPUTATION: RENAMING, RECODING, AND CREATING VARIABLES #########################


# Prepare data for imputation
nhanes <- nhanes_raw %>% 
  rename( 
    Age = RIDAGEYR,
    edu = DMDEDUC2,  # Education level
    pir = INDFMPIR,  # Income to poverty ratio
    Total_Chol = LBXTC,   # Total cholesterol
    Trig = LBXTR,    # Triglycerides
    HDL = LBDHDD,    # HDL cholesterol
    BMI = BMXBMI,    # Body mass index
    height = BMXHT,  # Height (cm)
    weight = BMXWT,  # Weight (kg)
    dm_self = DIQ010, # Self-reported diabetes
    HbA1c = LBXGH,   # Hemoglobin A1c
    Glucose = LBXGLU,    # Fasting plasma glucose
    ogtt = LBXGLT,   # Oral glucose tolerance test
    chd = MCQ160C,   # Coronary heart disease
    mi = MCQ160E,    # Heart attack
    Stroke = MCQ160F, # Stroke
    # DM_family = MCQ300C # not for 2001-2002 cycle
    ) %>%
  mutate(
    # Recode "Refused" or "Don't Know" to missing
    edu = ifelse(edu == 7 | edu == 9, NA, edu),
    dm_self = ifelse(dm_self == 7 | dm_self == 9, NA, ifelse(dm_self == 3, 2, dm_self)), # assume borderline is negative
    DM_family= ifelse(MCQ260AA==1 |MCQ260AB==2|MCQ260AC==3|MCQ260AD==4|MCQ260AE==5|MCQ260AF==6|MCQ260AG==7|MCQ260AH==8, 1, 0),
    chd = ifelse(chd == 7 | chd == 9, NA, chd),
    mi = ifelse(mi == 7 | mi == 9, NA, mi),
    Stroke = ifelse(Stroke == 7 | Stroke == 9, NA, Stroke),
    CDQ001 = ifelse(CDQ001 == 7 | CDQ001 == 9, NA, CDQ001),
    CDQ002 = ifelse(CDQ002 == 7 | CDQ002 == 9, NA, CDQ002),
    CDQ003 = ifelse(CDQ003 == 7 | CDQ003 == 9, NA, CDQ003),
    CDQ004 = ifelse(CDQ004 == 7 | CDQ004 == 9, NA, CDQ004),
    CDQ005 = ifelse(CDQ005 == 7 | CDQ005 == 9, NA, CDQ005),
    CDQ006 = ifelse(CDQ006 == 7 | CDQ006 == 9, NA, CDQ006),
    BPQ040A = ifelse(BPQ040A == 7 | BPQ040A == 9, NA, BPQ040A),
    DM_family = ifelse(is.na(DM_family), 0, DM_family),
    SMQ020 = ifelse(SMQ020 == 7 | SMQ020 == 9, NA, SMQ020),
    SMQ040 = ifelse(SMQ040 == 7 | SMQ040 == 9, NA, SMQ040),
    # Define additional variables
    Female = case_when(RIAGENDR == 1 ~ 0, RIAGENDR == 2 ~ 1),
    Race = case_when(RIDRETH1 == 1 | RIDRETH1 == 2 ~ 3, # Hispanic
                     RIDRETH1 == 3 ~ 1,                 # White
                     RIDRETH1 == 4 ~ 2,                 # Black
                     RIDRETH1 == 5 ~ 4),                # Other
    SBP = rowMeans(select(nhanes_raw, starts_with("BPXSY")), na.rm = T),
    DBP = rowMeans(select(nhanes_raw, starts_with("BPXDI")), na.rm = T),
    HPT_Txt = case_when(BPQ020 == 2 | BPQ040A == 2 ~ 0, 
                        BPQ020 == 1 & BPQ040A == 1 ~ 1),
    Smoking = case_when(SMQ020 == 2 | (SMQ020 == 1 & SMQ040 == 3) ~ 0,
                        SMQ020 == 1 & (SMQ040 == 1 | SMQ040 == 2) ~ 1),
    # NOTE: Rose Questionnaire was only asked to participants 40 years and older
    # Participants less than 40 years assumed to be negative
    roseQ = ifelse(Age < 40, 0, ifelse(CDQ001==1 & CDQ002==1 & CDQ004==1 &  CDQ003!=1 & CDQ005==1 & CDQ006==1 & 
                                         ((CDQ009D==4 | CDQ009E==5) | (CDQ009F==6 & CDQ009G==7)), 1, 0))
)
# Read drug information file (can be used for all NHANES cycles)
rx_info = read.xport('./RXQ_DRUG.XPT')
rx_info = rx_info %>%
  select(RXDDRGID, RXDDRUG, RXDDCI1A, RXDDCI1B, RXDDCI1C) %>%
  arrange(RXDDRGID)

# Append individual medication files from all cycles
rx_df = data.frame()
for (cycle_i in 1:length(cycles)) {
  cycle = cycles[cycle_i]
  cycle_rx = read_data(cycle_i, "RXQ_RX")
  cycle_rx = select(cycle_rx, SEQN, RXD030, RXDDRGID)
  rx_df = bind_rows(rx_df, cycle_rx)
}

# Sort by RXDDRGID
rx_df = arrange(rx_df, RXDDRGID, SEQN)
# Merge with rx_info
rx_df_w_info = merge(rx_df, rx_info, all.x = T)
rx_df_w_info = arrange(rx_df_w_info, SEQN)

# Create medication boolean variables
rx_df_w_info = rx_df_w_info %>%
  mutate(
    dm_rx_temp = case_when(RXDDCI1A == 358 & RXDDCI1B == 99 ~ 1, 
                           RXD030 == 2 | (RXD030 == 1 & !is.na(RXDDRUG)) ~ 0),
    angina_rx_temp = case_when(RXDDCI1A == 40 & RXDDCI1B == 45 ~ 1,
                               RXD030 == 2 | (RXD030 == 1 & !is.na(RXDDRUG)) ~ 0)
  )

rx_vars = rx_df_w_info %>%
  group_by(SEQN) %>%
  summarise(dm_rx = any(dm_rx_temp==1),
            angina_rx = any(angina_rx_temp==1))

# Merge processed medication variable
nhanes = merge(nhanes, rx_vars, all.x = T)

# Angina
nhanes = mutate(nhanes, angina = ifelse(roseQ == 1 | angina_rx==1, 1, 0))

# Define additional variables (these variables will not be imputed, but will be used for pre-imputation population stats)
# Note the way ifelse works: if some criteria are negative and some are NA, it will return NA.
# If at least one criterium is positive, it will return positive (even if others are NA)
# Change to: if all are NA, return NA. If some are negative and some are NA, return negative.
# If some are positive and some are NA, return positive.
# To implement this change, I will sum across the criteria, removing NA values.
# First, create indicator variables
nhanes$dm_self[nhanes$dm_self==7|nhanes$dm_self==9] =NA
nhanes$dm_self=nhanes$dm_self==1
nhanes$glucose_dm = nhanes$Glucose > 125
nhanes$ogtt_dm = nhanes$ogtt >= 200
nhanes$hba1c_dm = nhanes$HbA1c >= 6.5
# Calculate the sum across the criteria, removing NA values
nhanes$dm_score = rowSums(select(nhanes, dm_self, dm_rx, glucose_dm, ogtt_dm, hba1c_dm), na.rm = T)
# If all criteria are NA, Diabetes status is NA. Otherwise, if at least one criterion is met, Diabetes status is 1.
nhanes$Diabetes = ifelse(is.na(nhanes$dm_self)&nhanes$glucose_dm & nhanes$ogtt_dm & nhanes$hba1c_dm , NA,
                         #& nhanes$glucose_dm!=1 & nhanes$ogtt_dm!=1 & nhanes$hba1c_dm!=1 &  nhanes$dm_rx!=1
                         ifelse(nhanes$dm_score > 0, 1, 0))

nhanes = nhanes %>%
  mutate(
    HBP = ifelse(HPT_Txt == 1 | SBP >= 130 | DBP >= 85, 1, 0),
    CVD_history = ifelse(angina == 1 | chd == 1 | mi == 1 | Stroke == 1, 1, 0),
    # Define survey weights over 3 cycles
    # Before imputation, weights should be based on the subset of the population with examination data
    WT_TOTAL = 1/n.cycles * WTMEC2YR
  )

nhanes = nhanes %>%
  filter(Age >= 18) %>%
  select(SEQN, SDDSRVYR, WTINT2YR, WTMEC2YR, SDMVPSU, SDMVSTRA, Age, edu, pir, Total_Chol, Trig, HDL, DM_family, BMI,
         height, weight, dm_self, HbA1c, Glucose, ogtt,  chd, mi, Stroke, Female, Race, SBP, DBP, HPT_Txt, Smoking, 
         dm_rx, angina, Diabetes, HBP, CVD_history, WT_TOTAL)

# Dietary data
diet_df = data.frame()
for (cycle_i in 1:length(cycles)) {
  cycle = cycles[cycle_i]
  cycle_letter = cycle_letters[cycle_i]
  
  # Day 1
  file1 = paste("fped_dr1iff_", cycle_letter, ".csv", sep = "")
  diet1 = read.csv(paste('./', cycle, file1, sep = '/'))
  # Calculate totals
  diet1_sums = diet1 %>%
    group_by(SEQN) %>%
    summarise(fruit = sum(F_TOTAL, na.rm = T),
              veg = sum(V_TOTAL, na.rm = T),
              wh_grain = sum(G_WHL, na.rm = T),
              seafd = sum(M_FISH_HI, na.rm = T) + sum(M_FISH_LO, na.rm = T),
              nutsds = sum(M_NUTSD, na.rm = T),
              oils = sum(DISCFAT_OIL, na.rm = T))
  
  # # Day 2
  # file2 = paste("fped_dr2iff_", cycle_letter, ".csv", sep = "")
  # diet2 = read.csv(paste('./01_Input/NHANES/Raw', cycle, file2, sep = '/'))
  # diet2_sums = diet2 %>%
  #   group_by(SEQN, WTDRD1, WTDR2D) %>%
  #   summarise(fruit_d2 = sum(DR2I_F_TOTAL, na.rm = T),
  #             veg_d2 = sum(DR2I_V_TOTAL, na.rm = T),
  #             wh_grain_d2 = sum(DR2I_G_WHOLE, na.rm = T),
  #             seafd_d2 = sum(DR2I_PF_SEAFD_HI, na.rm = T) + sum(DR2I_PF_SEAFD_LOW, na.rm = T),
  #             nutsds_d2 = sum(DR2I_PF_NUTSDS, na.rm = T),
  #             oils_d2 = sum(DR2I_OILS, na.rm = T))
  
  # # Merge Day 1 and Day 2
  # diet_d1_d2 = merge(diet1_sums, diet2_sums, all = T)
  # # Average Day 1 and Day 2
  # diet_avg = diet_d1_d2 %>%
  #   mutate(fruit = rowMeans(select(diet_d1_d2, fruit_d1, fruit_d2), na.rm = T),
  #          veg = rowMeans(select(diet_d1_d2, veg_d1, veg_d2), na.rm = T),
  #          wh_grain = rowMeans(select(diet_d1_d2, wh_grain_d1, wh_grain_d2), na.rm = T),
  #          seafd = rowMeans(select(diet_d1_d2, seafd_d1, seafd_d2), na.rm = T),
  #          nutsds = rowMeans(select(diet_d1_d2, nutsds_d1, nutsds_d2), na.rm = T),
  #          oils = rowMeans(select(diet_d1_d2, oils_d1, oils_d2), na.rm = T)) %>%
  #   select(SEQN, WTDRD1, WTDR2D, fruit, veg, wh_grain, seafd, nutsds, oils)
  # Append to final df
  diet_df = rbind(diet_df, diet1_sums)
}

# Merge dietary data with other data
nhanes_diet = merge(nhanes, diet_df, all.x = T)


##### MULTIPLE IMPUTATION ##########################################################################

set.seed(1899)
orig.data = nhanes
# orig.data = dm_fi_sample
orig.data$Subject_ID <- c(1:nrow(orig.data))
n.participants <- nrow(orig.data)

num.imp <-10 # Set number of imputations
num.knots <- 3 # set 0 for linear predictions of continous variables
orig.data$Race = factor(orig.data$Race)
#orig.data$hhfs = factor(orig.data$hhfs)
orig.data$edu = factor(orig.data$edu)
imputed.data <- aregImpute(~Age+Female+Race+BMI+height+weight+Total_Chol+Trig+HDL+SBP+DBP+HPT_Txt+Smoking+Glucose+HbA1c+DM_family+dm_self+
                             pir+edu+angina+dm_rx+chd+mi+Stroke, orig.data, n.impute=num.imp, nk=num.knots, x=T)
imputed.dataset <- data.frame()
for(j in 1:num.imp) {  
  print(j)
  completed <- orig.data
  imputed <- impute.transcan(imputed.data, imputation = j, data = orig.data, list.out=TRUE, pr=FALSE, check = FALSE)
  completed[names(imputed)] <- imputed # Replace variables from the original dataset with imputed values
  imputed.dataset <- rbind(imputed.dataset, as.data.frame(completed))
}  

# Define remaining variables
imputed.dataset = imputed.dataset %>%
  mutate(
    HBP = ifelse(HPT_Txt == 1 | SBP >= 130 | DBP >= 85, 1, 0),
    CVD_history = ifelse(angina == 1 | chd == 1 | mi == 1 | Stroke == 1, 1, 0),
    Diabetes = ifelse(dm_self == 1 | dm_rx == 1 | Glucose > 125 | ogtt >= 200 | HbA1c >= 6.5, 1, 0),
    Prediabetes = ifelse(dm_self == 1 | dm_rx == 1 | Glucose >= 100 | ogtt >= 140 | HbA1c >= 5.7, 1, 0),
    # After imputation, we've filled in missing examination data, so we should use interview weights
    WT_TOTAL = 1/n.cycles * WTINT2YR
  )

# Save data set
write.csv(imputed.dataset, Imp_filepath, row.names = F)



