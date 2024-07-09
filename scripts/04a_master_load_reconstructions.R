
# ------------------------------------------------------------------------------------
# Load all reconstructions, including filtering:
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023
library(matrixStats)

# setwd to project folder:
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/")

# 00. Load functions & code:
# ------------------------------------------------------------------------------------
source("code/_convenience/frenchcolormap.R")
source("code/_functions_CMIP6.R")
source("code/_attribution_hildreth-lu.R")
source("code/_post-processing_reconstructions.R")


# 01. Load global observations:
# ------------------------------------------------------------------------------------

# load standard GMST datasets (+CRUTEM5 raw / HadSST4 raw):
load("data/03_processedOBS_reconstr/GMST_datasets.RData")

# load SST and Land Datasets, including their reconstructions for target variables:
load("data/03_processedOBS_reconstr/SST_datasets.RData")
load("data/03_processedOBS_reconstr/TLand_datasets.RData")


# 02. Load NEW OBS reconstructions (CRUTEM5, HadSST4, CLASSNMAT, ...):
# ------------------------------------------------------------------------------------
source("scripts/03b_load_global_obs_reconstruction.R")


# 03. Load long CMIP reconstructions (for >100.000 CMIP years):
# ------------------------------------------------------------------------------------
load("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/03_processedCMIP6_reconstr/CMIP6.tas_land_all.df_v5.RData")
load("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/03_processedCMIP6_reconstr/CMIP6.tos_all.df.RData")
CMIP6.tos_all.df$ann$Yhat$GSAT[which(CMIP6.tos_all.df$ann$Yhat$GSAT > 20)] = NA



# 04. Load CMIP6-piControl reconstructions:
# ------------------------------------------------------------------------------------
# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
# trends_piControl = extract.nyear.trend_(XAX = CMIP6.tas_ann_piControl_ct, nyears = 50, trend.sep = 20, 
#                                        var.names = c("GMST_FM", "GMSST", "GMLSAT_NI", "TMMSAT_40S_40N", "TMLSAT_40S_40N"))


# 05. Estimate forced response:
# ------------------------------------------------------------------------------------

# select models with at least 5 ensemble members to calculate the forced response:
names = unique(CMIP6.tos_all.df$M$mod)[which(sapply(unique(CMIP6.tos_all.df$M$mod), 
                                                    FUN=function(cur.mod) length(unique(CMIP6.tos_all.df$M$ens.mem[which(cur.mod == CMIP6.tos_all.df$M$mod & CMIP6.tos_all.df$M$scen == "historical")]))) >= 3)]

CMIP6.GSAT.f = matrix(data = NA, nrow=165, ncol = length(names))
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
  
  for (cur.year in 1850:2014) {
    CMIP6.GSAT.f[cur.year-1849,i] = mean(CMIP6.tos_all.df$ann$Y$GSAT[mod.ix][which(CMIP6.tos_all.df$M$year[mod.ix] == cur.year)])
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




# 08. Attribution and Filtering of time series:
# ------------------------------------------------------------------------------------
{
library("dplR") # package for band-pass filtering

    # GSAT Prepare data for plotting:
    GSAT.tas_land = get.df(Y = OBS.tas_land_$GSAT$ann$mod_p1_min, f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GSAT.tos = get.df(Y = OBS.tos_$GSAT$ann$mod_p1_min, f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GSAT.COBE_SST2 = get.df(Y = COBE_SST2.annual$GSAT_mod_p1_min, f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GSAT.ERSSTv5 = get.df(Y = ERSSTv5.annual$GSAT_mod_p1_min[1:167], f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GSAT.BEST_Land = get.df(Y = BEST_Land.annual$GSAT_mod_p1_min[101:271], f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GSAT.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GSAT$ann$mod_p1_min), f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GSAT.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GSAT$ann$mod_p1_min), f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GSAT.tos_ua = get.df(Y = HadSST_ua.annual$GSAT_mod_p1_min, f = rowMeans(CMIP6.GSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
    
    # GMSST Prepare data for plotting:
    GMSST.HadSST4 = get.df(Y = HadSST4.global.annual$Anomaly[1:171], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.tas_land = get.df(Y = OBS.tas_land_$GMSST$ann$mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMSST.tos = get.df(Y = OBS.tos_$GMSST$ann$mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMSST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMSST_mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMSST_mod_p1_min[1:167], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.BEST_Land = get.df(Y = BEST_Land.annual$GMSST_mod_p1_min[101:271], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMSST$ann$mod_p1_min), f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMSST$ann$mod_p1_min), f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMSST.tos_ua = get.df(Y = HadSST_ua.annual$GMSST_mod_p1_min, f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
    
    
    # GMLSAT Prepare data for plotting:
    GMLSAT_NI.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.tas_land = get.df(Y = OBS.tas_land_$GMLSAT_NI$ann$mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMLSAT_NI.tos = get.df(Y = OBS.tos_$GMLSAT_NI$ann$mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMLSAT_NI.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMLSAT_NI_mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMLSAT_NI_mod_p1_min[1:167], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.BEST_Land = get.df(Y = BEST_Land.annual$GMLSAT_NI_mod_p1_min[101:271], f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMLSAT_NI$ann$mod_p1_min), f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMLSAT_NI$ann$mod_p1_min), f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMLSAT_NI.tos_ua = get.df(Y = HadSST_ua.annual$GMLSAT_NI_mod_p1_min, f = rowMeans(CMIP6.GMLSAT.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
    
    
    # GMST Prepare data for plotting:
    GMST.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.tas_land = get.df(Y = OBS.tas_land_$GMST_FM$ann$mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMST.tos = get.df(Y = OBS.tos_$GMST$ann$mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
    GMST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMST_FM_mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMST_FM_mod_p1_min[1:167], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.BEST_Land = get.df(Y = BEST_Land.annual$GMST_FM_mod_p1_min[101:271], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMST_FM$ann$mod_p1_min), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMST_FM$ann$mod_p1_min), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
    GMST.tos_ua = get.df(Y = HadSST_ua.annual$GMST_FM_mod_p1_min, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
}



# 09. Derive window correlation for OBS and CMIP (Figure 2):
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/")

# OBS Window correlation:
OBS_mod_p1_min = get.diff.window.cor_ens(x = OBS.tos_$GMST_FM$ann$mod_p1_min, y = OBS.tas_land_$GMST_FM$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 1:200)
OBS_mod_p0 = get.diff.window.cor_ens(x = OBS.tos_$GMST_FM$ann$mod_p0, y = OBS.tas_land_$GMST_FM$ann$mod_p0, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 1:200)

CMIP6_matrix = get.CMIP6.recon.matrix(CMIP6.tas_land_all.df, CMIP6.tos_all.df)
CMIP6_mod_p1_min = get.diff.window.cor_ens(y = CMIP6_matrix$tas_land_mod_p1, x = CMIP6_matrix$tos_mod_p1, years = 1850:2014, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix = 1:602)
CMIP6_mod_p1_min_pt1 = get.diff.window.cor_ens(y = CMIP6_matrix$tas_land_mod_p1_pt1, x = CMIP6_matrix$tos_mod_p1_pt1, years = 1850:2014, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix = 1:602)
CMIP6_mod_p0 = get.diff.window.cor_ens(y = CMIP6_matrix$tas_land_mod_p0, x = CMIP6_matrix$tos_mod_p0, years = 1850:2014, w.width = 51, W = 20, center = F, center.ens.means = F, ens.ix = 1:602)
CMIP6_mod_p0_pt1 = get.diff.window.cor_ens(y = CMIP6_matrix$tas_land_mod_p0_pt1, x = CMIP6_matrix$tos_mod_p0_pt1, years = 1850:2014, w.width = 51, W = 20, center = F, center.ens.means = F, ens.ix = 1:602)



