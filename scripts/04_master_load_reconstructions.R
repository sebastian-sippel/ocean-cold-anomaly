
# ------------------------------------------------------------------------------------
# Load all reconstructions, including filtering:
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023
library(matrixStats)

# 00. Load functions & code:
# ------------------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_attribution_hildreth-lu.R")


# 01. Load global observations:
# ------------------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations_SST.R")
# load additional SST datasets:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/SST_datasets.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/TLand_datasets.RData")


# 02. Load NEW reconstructions:
# ------------------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03b_load_global_obs_reconstruction.R")


# 03. Load CMIP reconstructions:
# ------------------------------------------------------------------------------------
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tas_land_all.df_v4.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tos_all.df_v4.RData")
CMIP6.tos_all.df$ann$Yhat$GSAT[which(CMIP6.tos_all.df$ann$Yhat$GSAT > 20)] = NA



# 04. Load CMIP6-piControl reconstructions:
# ------------------------------------------------------------------------------------
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
trends_piControl = extract.nyear.trend_(XAX = CMIP6.tas_ann_piControl_ct, nyears = 50, trend.sep = 20, 
                                        var.names = c("GMST_FM", "GMSST", "GMLSAT_NI", "TMMSAT_40S_40N", "TMLSAT_40S_40N"))


# 05. Load forced response estimates from DAMIP simulations:
# ------------------------------------------------------------------------------------
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_ALL_ct_f.RData")
CMIP6.MMM = rowMeans((CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist))
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HISTssp245_ct_f.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HIST_ct_f.RData")
names = names(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist)

CMIP6.GMSST.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GMMSAT.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GMLSAT.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GMST.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GSAT.tos_yhat.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GSAT.tas_land_yhat.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GMST.tos_yhat.f = matrix(data = NA, nrow=165, ncol = length(names))
CMIP6.GMST.tas_land_yhat.f = matrix(data = NA, nrow=165, ncol = length(names))

for (i in 1:length(names)) {
  print(i)
  cur.mod = names[i]
  mod.ix = which(cur.mod == CMIP6.tos_all.df$M$mod)
  # CMIP6.tos_all.df$ann$Y$GMSST
  
  for (cur.year in 1850:2014) {
    CMIP6.GMSST.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Y$GMSST[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GMLSAT.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Y$GMLSAT_NI[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GMMSAT.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Y$GMMSAT[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GMST.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Y$GMST_FM[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GSAT.tos_yhat.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Yhat$GSAT[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GSAT.tas_land_yhat.f[cur.year-1849,i] = mean(CMIP6.tas_land_all.df$ann$Yhat$GSAT[mod.ix][which(CMIP6.tas_land_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GMST.tos_yhat.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Yhat$GMST_FM[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
    CMIP6.GMST.tas_land_yhat.f[cur.year-1849,i] = mean(CMIP6.tas_land_all.df$ann$Yhat$GMST_FM[mod.ix][which(CMIP6.tas_land_all.df$M$year[mod.ix] == cur.year)])
  }
}


# 07. Read Proxy data and proxy-based products:
# ------------------------------------------------------------------------------------

# Read Proxy data:
{
  setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/proxy/Crowley_etal_2014/")
  
  proxy.Crowley2014 = read.table("eft234-sup-0007-ts05.txt", header = F, skip=1, nrows = 203)
  names(proxy.Crowley2014) <- c("Year",	"30N-90N comp",	"30S-30N comp",	"30S-30N enso comp",	"30S-30N ext comp",	"30S-30N ext comp enso",	"30S-90S comp",	"Global proxy (1782)",
                                "Global proxy/instrument (appendix, 1894)",	"Global instrument composite (HadCruNoaaNasa)",	"HadCRUt3 30N-90N",	"HadCRUt3 30S-30N",	"HadCRUt3 30S-90S",	"GHG [temp scaled]",
                                "Volcanic [temp scaled]",	"(Global Proxy 1782)-GHG",	"(Global Proxy 1782)-(GHG+Volc)",	"(Global proxy/instrum 1894)-(GHG+Volc)",	"(Global instrument composite)-(GHG+Volc)")
  
  
  # Read PAGES2k / Neukom et al. (2019) temperature proxies:
  setwd("/net/h2o/climphys1/sippels/_DATASET/neukom2019temp/recons/")
  neukom2019_BHM = read.table(file = "BHM.txt", header = T)
  neukom2019_CPS_new = read.table(file = "CPS_new.txt", header = F, skip = 1)
  colnames(neukom2019_CPS_new) = colnames(neukom2019_BHM)
  neukom2019_DA = read.table(file = "DA.txt", header = F, skip = 1)
  colnames(neukom2019_DA) = colnames(neukom2019_BHM)
  neukom2019_M08 = read.table(file = "M08.txt", header = F, skip = 1) # str(neukom2019_M08)
  colnames(neukom2019_M08) = colnames(neukom2019_BHM)
  neukom2019_OIE = read.table(file = "OIE.txt", header = F, skip = 1) # str(neukom2019_OIE)
  colnames(neukom2019_OIE) = colnames(neukom2019_BHM)
  neukom2019_PAI = read.table(file = "PAI.txt", header = F, skip = 1) # str(neukom2019_PAI)
  colnames(neukom2019_PAI) = colnames(neukom2019_BHM)
  neukom2019_PCR = read.table(file = "PCR.txt", header = F, skip = 1) # str(neukom2019_PCR)
  colnames(neukom2019_PCR) = colnames(neukom2019_BHM)
  
  # get full ensemble range:
  neukom2019_full_ensemble = read.table(file = "Full_ensemble_median_and 95pct_range.txt", header = F, skip = 5) # str(neukom2019_PCR)
  colnames(neukom2019_full_ensemble) = c("Year", "CowtanWay_instrumental_target", "Full_ensemble_median", "Full_ensemble_2.5th_percentile",	
                                         "Full_ensemble_97.5th_percentile",	"CowtanWay_instrumental_target_31year_filtered",	"31year_filtered_full_ensemble_median", 
                                         "31year_filtered_full_ensemble_2.5th_percentile",	"31year_filtered_full_ensemble_97.5th_percentile")
}


# 00.(g) Get additional proxy datasets:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/03d_read_ocean2k.R")

get.RE_score <- function(cur.region ) {
  cur.year = c(cur.region [[5]])
  year.ix = which(cur.year %in% 1870:2000)
  RE = sapply(X = 1:dim(cur.region[[2]])[2], FUN=function(i) {
    sum(cur.region[[2]][year.ix,i], na.rm = T)
  })
  return(RE) # which(RE > 0)
}




# 08. Attribution and Filtering of time series:
# ------------------------------------------------------------------------------------
{
library("dplR") # package for band-pass filtering
    
    CMIP6.MMM = rowMeans(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist)[1:171]
    # f = rowMeans(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist)[1:171]
    # f = cbind(rowMeans(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histGHG), rowMeans(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histAER), rowMeans(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histNAT))
    
    # GSAT Prepare data for plotting:
    GSAT.tas_land = get.df(Y = OBS.tas_land_$GSAT$ann$mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2020)
    GSAT.tos = get.df(Y = OBS.tos_$GSAT$ann$mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2020)
    GSAT.HadISST = get.df(Y = HadISST.annual$GSAT_mod_p1_min[1:151], f = CMIP6.MMM, years = HadISST.annual$Year[1:151], center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.COBE_SST2 = get.df(Y = COBE_SST2.annual$GSAT_mod_p1_min, f = CMIP6.MMM, years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.ERSSTv5 = get.df(Y = ERSSTv5.annual$GSAT_mod_p1_min[1:167], f = CMIP6.MMM, years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.ERSSTv4 = get.df(Y = ERSSTv4.annual$GSAT_mod_p1_min, f = CMIP6.MMM, years = ERSSTv4.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.BEST_Land = get.df(Y = BEST_Land.annual$GSAT_mod_p1_min[101:271], f = CMIP6.MMM, years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GSAT$ann$mod_p1_min), f = CMIP6.MMM, years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GSAT$ann$mod_p1_min), f = CMIP6.MMM, years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GSAT.tos_ua = get.df(Y = HadSST_ua.annual$GSAT_mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2020)
    
    
    # GMSST Prepare data for plotting:
    GMSST.HadSST4 = get.df(Y = HadSST4.global.annual$Anomaly[1:171], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.tas_land = get.df(Y = OBS.tas_land_$GMSST$ann$mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMSST.tos = get.df(Y = OBS.tos_$GMSST$ann$mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMSST.HadISST = get.df(Y = HadISST.annual$GMSST_mod_p1_min[1:151], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = HadISST.annual$Year[1:151], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMSST_mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMSST_mod_p1_min[1:167], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.ERSSTv4 = get.df(Y = ERSSTv4.annual$GMSST_mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = ERSSTv4.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.BEST_Land = get.df(Y = BEST_Land.annual$GMSST_mod_p1_min[101:271], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMSST$ann$mod_p1_min), f = CMIP6.MMM, years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMSST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMSST$ann$mod_p1_min), f = CMIP6.MMM, years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMSST.tos_ua = get.df(Y = HadSST_ua.annual$GMSST_mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2020)
    
    
    # GMLSAT Prepare data for plotting:
    GMLSAT_NI.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.tas_land = get.df(Y = OBS.tas_land_$GMLSAT_NI$ann$mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMLSAT_NI.tos = get.df(Y = OBS.tos_$GMLSAT_NI$ann$mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMLSAT_NI.HadISST = get.df(Y = HadISST.annual$GMLSAT_NI_mod_p1_min[1:151], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = HadISST.annual$Year[1:151], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMLSAT_NI_mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMLSAT_NI_mod_p1_min[1:167], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.ERSSTv4 = get.df(Y = ERSSTv4.annual$GMLSAT_NI_mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = ERSSTv4.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.BEST_Land = get.df(Y = BEST_Land.annual$GMLSAT_NI_mod_p1_min[101:271], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMLSAT_NI$ann$mod_p1_min), f = CMIP6.MMM, years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMLSAT_NI.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMLSAT_NI$ann$mod_p1_min), f = CMIP6.MMM, years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMLSAT_NI.tos_ua = get.df(Y = HadSST_ua.annual$GMLSAT_NI_mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2020)
    
    
    # GMST Prepare data for plotting:
    GMST.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.tas_land = get.df(Y = OBS.tas_land_$GMST_FM$ann$mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMST.tos = get.df(Y = OBS.tos_$GMST$ann$mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2014)
    GMST.HadISST = get.df(Y = HadISST.annual$GMST_FM_mod_p1_min[1:151], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = HadISST.annual$Year[1:151], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMST_FM_mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMST_FM_mod_p1_min[1:167], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.ERSSTv4 = get.df(Y = ERSSTv4.annual$GMST_FM_mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = ERSSTv4.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.BEST_Land = get.df(Y = BEST_Land.annual$GMST_FM_mod_p1_min[101:271], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMST_FM$ann$mod_p1_min), f = CMIP6.MMM, years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMST_FM$ann$mod_p1_min), f = CMIP6.MMM, years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2020)
    GMST.tos_ua = get.df(Y = HadSST_ua.annual$GMST_FM_mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
}




# 09. Define function for window correlation and trend analysis:
# ------------------------------------------------------------------------------------

get.diff.window.cor <- function(x, y, years = 1850:2020, w.width = 51, W = 20, center = T, ens.means = NA) {
  
  x.low = pass.filt(y = x, W = W, type = "low", method = "Butterworth")
  x.high = pass.filt(y = x, W = W, type = "high", method = "Butterworth")
  y.low = pass.filt(y = y, W = W, type = "low", method = "Butterworth")
  y.high = pass.filt(y = y, W = W, type = "high", method = "Butterworth")
  if (center == T & !is.numeric(ens.means)) {
    x.low = x.low - mean(x.low[match(x = 1961:1990, table = years)])
    x.high = x.high - mean(x.high[match(x = 1961:1990, table = years)])
    y.low = y.low - mean(y.low[match(x = 1961:1990, table = years)])
    y.high = y.high - mean(y.high[match(x = 1961:1990, table = years)])
  } else if (center == T & is.numeric(ens.means)) {
    x.low = x.low - mean(ens.means$x.low[match(x = 1961:1990, table = years)])
    x.high = x.high - mean(ens.means$x.high[match(x = 1961:1990, table = years)])
    y.low = y.low - mean(ens.means$y.low[match(x = 1961:1990, table = years)])
    y.high = y.high - mean(ens.means$y.high[match(x = 1961:1990, table = years)])
  }
  # Get Differences & Correlations:
  # plot(y); lines(x-y); lines(x.low - y.low)
  ret.df = data.frame(x, x.high, x.low, y, y.high, y.low, 
                      x-y, x.high-y.high, x.low-y.low,  
                      rollapply(data = cbind(x, y), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA),
                      rollapply(data = cbind(x.high, y.high), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA), 
                      rollapply(data = cbind(x.low, y.low), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA))
  names(ret.df) = c("x", "x.high", "x.low", "y", "y.high", "y.low", 
                    "diff_x_y", "diff_x.high_y.high", "diff_x.low_y.low", 
                    "cor_x_y", "cor_x.high_y.high", "cor_x.low_y.low")
  return(ret.df)
}

get.diff.window.cor_ens <- function(x, y, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix) {
  
  if (center.ens.means == T) {
    ens.means = get.diff.window.cor(colMedians(x[ens.ix,]), colMedians(y[ens.ix,]), years = years, w.width = w.width, W = W, center = T)
  } else {
    ens.means = NA
  }
  
  dat = list()
  for (en in 1:length(ens.ix)) {
    if (center.ens.means == T) {
      dat[[en]] = get.diff.window.cor(x[en,], y[en,], years = years, w.width = w.width, W = W, center = center, ens.means = ens.means) 
    } else {
      dat[[en]] = get.diff.window.cor(x[en,], y[en,], years = years, w.width = w.width, W = W, center = center, ens.means = ens.means) 
    }
  }
  
  # Get quantile estimates:
  diff_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x_y)), probs = c(0.025, 0.5, 0.975))
  diff_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.high_y.high)), probs = c(0.025, 0.5, 0.975))
  diff_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.low_y.low)), probs = c(0.025, 0.5, 0.975))
  
  cor_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x_y)), probs = c(0.025, 0.5, 0.975))
  cor_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.high_y.high)), probs = c(0.025, 0.5, 0.975))
  cor_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.low_y.low)), probs = c(0.025, 0.5, 0.975))
  
  ret.df = data.frame(diff_x_y, diff_x.high_y.high, diff_x.low_y.low, cor_x_y, cor_x.high_y.high, cor_x.low_y.low)
  names(ret.df) = c("diff_x_y_2.5", "diff_x_y_50", "diff_x_y_97.5", "diff_x.high_y.high_2.5", "diff_x.high_y.high_50", "diff_x.high_y.high_97.5", "diff_x.low_y.low_2.5", "diff_x.low_y.low_50", "diff_x.low_y.low_97.5", 
                    "cor_x_y_2.5", "cor_x_y_50", "cor_x_y_97.5", "cor_x.high_y.high_2.5", "cor_x.high_y.high_50", "cor_x.high_y.high_97.5", "cor_x.low_y.low_2.5", "cor_x.low_y.low_50", "cor_x.low_y.low_97.5")
  return(ret.df)
}

get.trend <- function(x, trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014) {
  
  if (all(is.na(x))) return(rep(NA, length(trend.years)))
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", 1:length(trend.years), sep="")
  return(trends.out)
}


# get.trend(x = OBS.tos_$GMSST$ann$mod_p1_min[1,], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2020)
# -> further changes to do: Good name for trends, or meta-file...  
get.trend_ <- function(x, trend.length = 50, years = 1850:2014) {
  
  trend.years = lapply(X = years[1]:(tail(years, 1)-trend.length+1), FUN = function(cur.year) cur.year:(cur.year+trend.length-1))
    
  if (all(is.na(x))) return(rep(NA, length(trend.years)))
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", seq(1850, 1850+length(trend.years)-1), sep="")
  return(trends.out)
}






# 10. Derive window correlation:
# ------------------------------------------------------------------------------------
OBS_mod_p1_min = get.diff.window.cor_ens(x = OBS.tos_$GSAT$ann$mod_p1_min, y = OBS.tas_land_$GSAT$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 94:200)

# implement CMIP trends here...


# 11. Derive OBS trends & low-freq. correction:
# ------------------------------------------------------------------------------------

# correction with hybrid36 time series:
# plot(GMSST.tos$mod_p1_min_50, type='l', ylim = c(-0.8, 0.8))
# lines(GMSST.hybrid36$mod_p1_min_50, col = "blue")
# lines(GMSST.tos$lp_50, col = "black", lty = 2)
# lines(GMSST.hybrid36$lp_50, col = "blue", lty = 2)

# lines(x=c(40,40), y = c(-1,1)); lines(x=c(90,90), y = c(-1,1))
hybrid.cor = rep(0, 171)
hybrid.cor[41:100] = -GMSST.tos$res_lp_50[41:100] + GMSST.hybrid36$res_lp_50[41:100]
# plot(hybrid.cor)
OBS.tos_GMSST_cor = OBS.tos_$GMSST$ann$mod_p1_min + rep.row(x = hybrid.cor, n = 200)

# get trends:
probs = c(0.025, 0.5, 0.975)

tos_trends_GMSST = apply(X = OBS.tos_$GMSST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_, trend.length = 50, years = 1850:2020)
tos_trends_GMSST_q = apply(X = tos_trends_GMSST, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))
tas_land_trends_GMSST = apply(X = OBS.tas_land_$GMSST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_, trend.length = 50, years = 1850:2020)
tas_land_trends_GMSST_q = apply(X = tas_land_trends_GMSST, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

tas_land_trends_GMLSAT = apply(X = OBS.tas_land_$GMLSAT_NI$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_, trend.length = 50, years = 1850:2020)
tas_land_trends_GMLSAT_q = apply(X = tas_land_trends_GMLSAT, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

tos_hybrid36_trends_GMSST = get.trend_(x = hybrid36.annual$Anomaly, trend.length = 50, years = 1850:2020)
tos_hybrid36.cor_trends_GMSST = apply(X = OBS.tos_GMSST_cor, MARGIN = 1, FUN=get.trend_, trend.length = 50, years = 1850:2020)
tos_hybrid36.cor_trends_GMSST_q = apply(X = tos_hybrid36.cor_trends_GMSST, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))




## 12. Derive CMIP6 trends (and CMIP window correlation):
## ----------------------------------------
CMIP6.Tropics = cbind(CMIP6.tos_all.df$ann$Y$IndianOcean, CMIP6.tos_all.df$ann$Y$WPacific, CMIP6.tos_all.df$ann$Y$WAtlantic)
CMIP6.Tropics_ = apply(X = CMIP6.Tropics, MARGIN = 1, FUN=function(x) weighted.mean(x = x, w = c(w_IOcean, w_WPacific, w_WAtlantic)))

{
  
  # select CMIP6 historical members:
  CMIP6.tas_land.all = data.frame(cbind(all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  CMIP6.tos.all = data.frame(cbind(all = paste(CMIP6.tos_all.df$M$mod, "_", CMIP6.tos_all.df$M$scen, "_", CMIP6.tos_all.df$M$ens.mem, sep="")))
  ens.mem = data.frame(cbind(mod=CMIP6.tas_land_all.df$M$mod, scen = CMIP6.tas_land_all.df$M$scen, ens.mem = CMIP6.tas_land_all.df$M$ens.mem, 
                             all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  ens.mem.un = unique(ens.mem)
  ens.mem.un = ens.mem.un[which(ens.mem.un$scen == "historical"),]
  # remove all ensemble members that contain NA's:
  na.mems = unique(CMIP6.tos.all$all[which(is.na(CMIP6.tos_all.df$ann$Yhat$GSAT))])
  omit.ix=na.omit(match(x = na.mems, table = ens.mem.un$all))
  if (length(omit.ix) > 0)  ens.mem.un = ens.mem.un[-omit.ix,]
  
  CMIP6.tas_land_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tos_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tas_land_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tos_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  
  # Get trends for each ensemble member:
  CMIP6.TMLSAT_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.TMMSAT_tos.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMLSAT_NI_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMSST_tos.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMSST_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMSST_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMLSAT_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.Tropics_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  CMIP6.GMST_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 116))
  
  
  for (en in 1:dim(ens.mem.un)[1]) {
    print(en)
    
    ix_tas_land = which(CMIP6.tas_land.all$all == ens.mem.un$all[en] & CMIP6.tas_land_all.df$M$year %in% 1850:2014)
    ix_tos = which(CMIP6.tos.all$all == ens.mem.un$all[en] & CMIP6.tos_all.df$M$year %in% 1850:2014)
    
    CMIP6.tas_land_hist_mod_p1[en,] = CMIP6.tas_land_all.df$ann$Yhat$GSAT[ix_tas_land]
    CMIP6.tos_hist_mod_p1[en,] = CMIP6.tos_all.df$ann$Yhat$GSAT[ix_tos]
    
    CMIP6.tas_land_hist_mod_p1_pt1[en,] = CMIP6.tas_land_all.df$ann$Yhat_pt1$GSAT[ix_tas_land]
    CMIP6.tos_hist_mod_p1_pt1[en,] = CMIP6.tos_all.df$ann$Yhat_pt1$GSAT[ix_tos]
    
    # Get trends: 
    CMIP6.TMLSAT_tas_land.trends[en,] = get.trend_(x = CMIP6.tas_land_all.df$ann$Yhat$TMLSAT_40S_40N[ix_tas_land], trend.length = 50, years = 1850:2014)
    CMIP6.TMMSAT_tos.trends[en,] = get.trend_(x = CMIP6.tos_all.df$ann$Yhat$TMMSAT_40S_40N[ix_tos], trend.length = 50, years = 1850:2014)
    CMIP6.GMLSAT_NI_tas_land.trends[en,] = get.trend_(x = CMIP6.tas_land_all.df$ann$Yhat$GMLSAT_NI[ix_tas_land], trend.length = 50, years = 1850:2014)
    CMIP6.GMSST_tos.trends[en,] = get.trend_(x = CMIP6.tos_all.df$ann$Yhat$GMSST[ix_tos], trend.length = 50, years = 1850:2014)
    CMIP6.GMSST_tas_land.trends[en,] = get.trend_(x = CMIP6.tas_land_all.df$ann$Yhat$GMSST[ix_tas_land], trend.length = 50, years = 1850:2014)
    CMIP6.GMSST_true.trends[en,] = get.trend_(x = CMIP6.tos_all.df$ann$Y$GMSST[ix_tos], trend.length = 50, years = 1850:2014)
    CMIP6.GMLSAT_true.trends[en,] = get.trend_(x = CMIP6.tos_all.df$ann$Y$GMLSAT_NI[ix_tas_land], trend.length = 50, years = 1850:2014)
    CMIP6.Tropics_true.trends[en,] = get.trend_(x = CMIP6.Tropics_[ix_tas_land], trend.length = 50, years = 1850:2014)
    # CMIP6.GMST_true.trends[en,] = get.trend(x = CMIP6.tos_all.df$ann$Y$GMST_FM[ix_tos], trend.length = 50, years = 1850:2014)
  }
  
  CMIP6_mod_p1_min = get.diff.window.cor_ens(x = CMIP6.tos_hist_mod_p1, y = CMIP6.tas_land_hist_mod_p1, years = 1850:2014, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix = 1:599)
  CMIP6_mod_p1_min_pt1 = get.diff.window.cor_ens(x = CMIP6.tos_hist_mod_p1_pt1, y = CMIP6.tas_land_hist_mod_p1_pt1, years = 1850:2014, w.width = 51, W = 20, center = F, center.ens.means = F, ens.ix = 1:599)
}



## 13. Define linear constraints based on ensemble or individual point estimate:
## ----------------------------------------

# Function to get linear model + constraint based on model relationship:
# y = CMIP6.trends$GMSST3_true
# x = CMIP6.trends$Tropics3_true
# x_new = ocean2k_trends[3]



get.linear.model.constraint <- function(y, x, x_new, plot.constraint = T) {
  dat = data.frame(y = y, x = x)
  new.dat = data.frame(x = x_new)
  reg.mod = lm(y ~ x, data = dat)
  
  n=length(reg.mod$residuals)
  RSS = sum((reg.mod$residuals)^2)
  s2 = 1/(n - 2) * RSS  # MSE
  sigma_x = sqrt( var(dat$x, na.rm=T) )
  new.dat.pred = predict(reg.mod, newdata = new.dat)
  sigma_f_x = sqrt(s2) * sqrt(1 + 1/n + (median(new.dat.pred) - mean(dat$x, na.rm = T)) / (n * sigma_x^2) )
  
  # combine variances in quadrature:
  mean_out =  median(new.dat.pred)
  sd_out =  sqrt(sigma_f_x^2) # + var(new.dat.pred))
  
  
  # plot emergent constraint for an ensemble:
  if(plot.constraint == T) {
    plot(x, y, col="red")
    abline(reg.mod, col ="red")
    
    plotCI(x = median(x_new), y = 0, li = quantile(x_new, 0.025), ui = quantile(x_new, 0.975), 
           err = "x", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
    plotCI(x = median(x_new), y = mean_out, li = mean_out - 2 * sd_out, ui = mean_out + 2 * sd_out, 
           err = "y", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
  }
  
  return(data.frame(mean_out=mean_out, sd_out = sd_out))
}


get.linear.model.constraint_ens <- function(y, x, x_new, plot.constraint = T) {
  
  dat = data.frame(y = y, x = x)
  new.dat = data.frame(x = x_new)
  reg.mod = lm(y ~ x, data = dat)
  
  n=length(reg.mod$residuals)
  RSS = sum((reg.mod$residuals)^2)
  s2 = 1/(n - 2) * RSS  # MSE
  sigma_x = sqrt( var(dat$x) )
  new.dat.pred = predict(reg.mod, newdata = new.dat)
  sigma_f_x = sqrt(s2) * sqrt(1 + 1/n + (median(new.dat.pred) - mean(dat$x)) / (n * sigma_x^2) )
  
  # combine variances in quadrature:
  mean_out =  median(new.dat.pred)
  sd_out =  sqrt(sigma_f_x^2 + var(new.dat.pred))
  
  
  # plot emergent constraint for an ensemble:
  if(plot.constraint == T) {
    plot(x, y, col="red")
    abline(reg.mod, col ="red")
    
    plotCI(x = median(x_new), y = 0, li = quantile(x_new, 0.025), ui = quantile(x_new, 0.975), 
           err = "x", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
    plotCI(x = median(x_new), y = mean_out, li = mean_out - 2 * sd_out, ui = mean_out + 2 * sd_out, 
           err = "y", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
  }
  
  return(data.frame(mean_out=mean_out, sd_out = sd_out))
}



