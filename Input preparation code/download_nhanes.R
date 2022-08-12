# Code to download NHANES data
# Created by Brianna Lauren (brianna.lauren@tufts.edu)
# Date: Nov 3 2020
# Last Updated by Lu Wang June24 2022

library(data.table)
library(foreign)
library(dplyr)

setwd('....')

# Download NHANES datasets
cycles <- c("2013-2014", "2015-2016")
cycle_letters <- c("H", "I")
files <- c("DEMO", "TCHOL", "HDL", "BMX", "BPQ", "BPX", "SMQ", "DIQ", "GHB", "GLU", "OGTT", "RXQ_RX", "MCQ", "CDQ", "FSQ", "DR1IFF", "DR2IFF")
base_url = "https://wwwn.cdc.gov/nchs/nhanes"
#Important, need to predefine the folder for each cycle to store the data, and the folder name ned to be consistent with that specified in cycles : e.g. "2013-2014", "2015-2016"
for (cycle_i in 1:length(cycles)) {
  for (file_base in files) {
    cycle = cycles[cycle_i]
    cycle_letter = cycle_letters[cycle_i]
    file = paste(file_base, "_", cycle_letter, ".XPT", sep = "")
    if (!file.exists(paste(".", cycle, file, sep = "/"))) {
      download.file(paste(base_url, cycle, file, sep = "/"), paste(".", cycle, file, sep = "/"), mode = "wb")
    }
  }
}
