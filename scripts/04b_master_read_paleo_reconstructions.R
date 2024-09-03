

# ------------------------------------------------------------------------------------
# Read Proxy data and proxy-based products:
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023
library(matrixStats)
# requires to run 04a_master_load_reconstructions.R beforehand!

# setwd to project folder:
setwd(".")




# 00. Load functions & code:
# ------------------------------------------------------------------------------------
source("code/_post-processing_reconstructions.R")


# 01. Read Crowley et al. 2014 reconstruction
# ------------------------------------------------------------------------------------
{
  proxy.Crowley2014 = read.table("data/00_DATASET/paleoclimate/Crowley_etal_2014/eft234-sup-0007-ts05.txt", header = F, skip=1, nrows = 203)
  names(proxy.Crowley2014) <- c("Year",	"30N-90N comp",	"30S-30N comp",	"30S-30N enso comp",	"30S-30N ext comp",	"30S-30N ext comp enso",	"30S-90S comp",	"Global proxy (1782)",
                                "Global proxy/instrument (appendix, 1894)",	"Global instrument composite (HadCruNoaaNasa)",	"HadCRUt3 30N-90N",	"HadCRUt3 30S-30N",	"HadCRUt3 30S-90S",	"GHG [temp scaled]",
                                "Volcanic [temp scaled]",	"(Global Proxy 1782)-GHG",	"(Global Proxy 1782)-(GHG+Volc)",	"(Global proxy/instrum 1894)-(GHG+Volc)",	"(Global instrument composite)-(GHG+Volc)")
}


# 02. Read PAGES2k / Neukom et al. (2019) temperature proxies:
# ------------------------------------------------------------------------------------

  neukom2019_BHM = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/BHM.txt", header = T)
  neukom2019_CPS_new = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/CPS_new.txt", header = F, skip = 1)
  colnames(neukom2019_CPS_new) = colnames(neukom2019_BHM)
  neukom2019_DA = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/DA.txt", header = F, skip = 1)
  colnames(neukom2019_DA) = colnames(neukom2019_BHM)
  neukom2019_M08 = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/M08.txt", header = F, skip = 1) # str(neukom2019_M08)
  colnames(neukom2019_M08) = colnames(neukom2019_BHM)
  neukom2019_OIE = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/OIE.txt", header = F, skip = 1) # str(neukom2019_OIE)
  colnames(neukom2019_OIE) = colnames(neukom2019_BHM)
  neukom2019_PAI = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/PAI.txt", header = F, skip = 1) # str(neukom2019_PAI)
  colnames(neukom2019_PAI) = colnames(neukom2019_BHM)
  neukom2019_PCR = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/PCR.txt", header = F, skip = 1) # str(neukom2019_PCR)
  colnames(neukom2019_PCR) = colnames(neukom2019_BHM)
  
  # get full ensemble range:
  neukom2019_full_ensemble = read.table(file = "data/00_DATASET/paleoclimate/neukom2019temp/Full_ensemble_median_and 95pct_range.txt", header = F, skip = 5) # str(neukom2019_PCR)
  colnames(neukom2019_full_ensemble) = c("Year", "CowtanWay_instrumental_target", "Full_ensemble_median", "Full_ensemble_2.5th_percentile",	
                                         "Full_ensemble_97.5th_percentile",	"CowtanWay_instrumental_target_31year_filtered",	"31year_filtered_full_ensemble_median", 
                                         "31year_filtered_full_ensemble_2.5th_percentile",	"31year_filtered_full_ensemble_97.5th_percentile")


# 03. Read Ocean2k / Tierney et al. (2015) temperature proxies:
# ------------------------------------------------------------------------------------
source("scripts/03d_read_ocean2k.R")
  
  
  
