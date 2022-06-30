#------------------------------------------------------------------------------#
# Extreme values -- ADAPTA 
# 
# Alliance Bioversity-CIAT, 2022
#------------------------------------------------------------------------------#
# R options
g <- gc(reset = T); rm(list = ls()) # Emptying the garbage collector
# .rs.restartR()                    # Restart R session
options(warn = -1, scipen = 999)    # Remove warning alerts and scientific notation
suppressMessages(library(pacman))   # Loading R-packages
suppressMessages(pacman::p_load(tidyverse,terra,extRemes)) # rtsVis
# reticulate::use_miniconda(condaenv = 'Miniconda3')
# suppressMessages(pacman::p_load(tensorflow, keras))
#------------------------------------------------------------------------------#
# main function
source("//CATALOGUE/Workspace14/WFP_ClimateRiskPr/0.Project_Documents/ext_val_main.R", encoding = 'UTF-8')
#------------------------------------------------------------------------------#
# root
root <- paste0('//catalogue/Workspace14/WFP_ClimateRiskPr/0.Project_Documents/')
# inputs
country <- 'kenya'
season <- c('S1', 'S2')
#------------------------------------------------------------------------------#
# indices
indx  <- c("AT","CDD","CSDI","HSI","IRR","NDWS","NTx35",
           "NWLD","NWLD50","NWLD90","P5D","P95","SHI","THI","TR")

# calc extreme values per indice
for(i in 1:length(indx)){
  for(s in 1:length(season)){
      infile <- terra::rast(paste0(root,country,"/results/Seasonal/",season[s],"/",indx[i],".tif"))
      outfile <- paste0(root,country,"/","results/Extreme_values/",indx[i],"_",season[s],"_extreme_values.tif")
      if (!file.exists(outfile)) {
      ExtVal(cdd1 = infile,
             outfile = outfile) 
    }
    else {
      cat(paste0("extreme value already calculated: ",indx[i],"_",season[s],"\n"))
    }
  }
} 
#------------------------------------------------------------------------------#