# Process output from DOCM model
# Current output: means and variances across individuals for each cycle and simulation

setwd("")

library(abind)

# Settings
vars <- c("Obesity", "Diabetes", "CVD_history", "All-cause mortality","DM mortality","CVD mortality", "Non CMD mortality" )
races <- c("all", "NHW", "NHB", "Hisp")
# seeds = c("1234", "2345", "23456", "234567")
seed= c("1234")
date = "2022-06-22"
year1 = 2001
n.cycle = 18
n.sim = 100
output_path = "final_output_2022_06_22b.csv"

# Initialize final output (means and 95% CIs for each variable, race, and cycle - long format)
final_out = data.frame()
model_out_means = readRDS(paste("sim_out_1000_No_Policy", "SEED", seed, date, ".rda", sep = "_"))[,,"mean",]
model_out_vars = readRDS(paste("sim_out_1000_No_Policy", "SEED", seed, date, ".rda", sep = "_"))[,,"var",]

for (race in races) {
  # Combine output from several runs
  r=switch(race, "all"=1, "NHW" = 2, "NHB" = 3, "Hisp" = 4)
  for (i_var in vars) {
    temp.means =matrix(as.numeric(model_out_means[i_var,(n.cycle*(r-1)+1):(n.cycle*r),]), nrow=n.cycle)
    temp.within.var = matrix( as.numeric(model_out_vars[i_var,(n.cycle*(r-1)+1):(n.cycle*r),]),nrow=n.cycle)
    temp_out <- data.frame(Year = year1:(year1+n.cycle-1))
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
write.csv(final_out, file = output_path, row.names = F)

