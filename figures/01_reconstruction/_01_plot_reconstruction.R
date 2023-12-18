
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021
library(matrixStats)

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//code/_attribution_hildreth-lu.R")

# Load global observations:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations_SST.R")
# load additional SST datasets:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/SST_datasets.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/TLand_datasets.RData")

# 00.(c) Load *new* reconstructions:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03b_load_global_obs_reconstruction.R")

# 00.(d) Load forced response estimates from DAMIP simulations:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_ALL_ct_f.RData")
CMIP6.MMM = rowMeans((CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist))
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HISTssp245_ct_f.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HIST_ct_f.RData")


# 00.(e) Load CMIP reconstructions for forced responses:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tas_land_all.df_v5.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tos_all.df_v5.RData")

# plot(CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist$`ACCESS-ESM1-5`)
# plot(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist$`ACCESS-ESM1-5`[1:180])

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


i = 7
plot(x = 1850:2014, y = CMIP6.GSAT.tas_land_yhat.f[,i], type='n')
lines(x = 1850:2014, y = rowMeans((CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist)), lty = 1)
lines(x = 1850:2014, y = rowMeans(CMIP6.GSAT.tos_yhat.f), type='l', col = "blue", lty = 2)
lines(x = 1850:2014, y = rowMeans(CMIP6.GSAT.tas_land_yhat.f), col = "brown",  lty = 1)
lines(x = 1850:2014, y = rowMeans(CMIP6.GMST.f, na.rm = T), col = "red",  lty = 1)




# ------------------------------------------------------------------------------------
# Plot reconstruction(s):
# ------------------------------------------------------------------------------------


# 01_main_reconstruction: AGMT anomalies | the !main! reconstructions (no uncertainties/sensitivities!):
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/01_reconstruction/")

library(RColorBrewer)
col = brewer.pal(n = 8, name = "Dark2")

library("dplR") # package for band-pass filtering




## ----------------------------------------------------------------------------------------
# 00. Run filtering on all time series here:
## ----------------------------------------------------------------------------------------
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




## ----------------------------------------------------------------------------------------
## Get Numbers for Abstract:
## ----------------------------------------------------------------------------------------


get.range.trends <- function(obs.data, cur.name = "GSAT", mod_type = "mod_p1", late.ix = 162:171, early.ix = 1:51) {
  
  trend.est = sapply(X = 1:200, FUN=function(ix) mean(obs.data[[cur.name]]$ann[[mod_type]][ix,late.ix]) - mean(obs.data[[cur.name]]$ann[[mod_type]][ix,early.ix]))
  trend.range = c(mean(trend.est), quantile(trend.est, probs = c(0.025, 0.5, 0.975)))
  return(trend.range)
}

get.OLS.trends <- function(obs.data, cur.name = "GSAT", mod_type = "mod_p1", years = 1880:2012) {
  
  year.ix = which(c(1850:2020) %in% years)
  trend.est = sapply(X = 1:200, FUN=function(ix) lm(c(obs.data[[cur.name]]$ann[[mod_type]][ix,year.ix]) ~ years)$coefficients[2] * length(years))
  trend.range = c(mean(trend.est), quantile(trend.est, probs = c(0.025, 0.5, 0.975)))
  return(trend.range)
}


## Warming trends:
## ----------------------------------------------------------------------------------------
# GMST reconstructions, 2011-2020 vs. 1850-1900:
get.range.trends(obs.data = OBS.tas_land_, cur.name = "GMST_FM", mod_type = "mod_p1_min", late.ix = 162:171, early.ix = 1:51)  # 1.06°C
get.range.trends(obs.data = OBS.tos_, cur.name = "GMST_FM", mod_type = "mod_p1_min", late.ix = 162:171, early.ix = 1:51)

## Correlations:
## ----------------------------------------------------------------------------------------
#### GMST:
### ANNUAL:
# 1850-1900
cor(CRUTEM5.global.annual$Anomaly[1:51], HadSST4.global.annual$Anomaly[1:51], use = "complete.obs")  # 0.47
cor(GMST.tas_land$mod_p1_min_50[1:51], GMST.tos$mod_p1_min_50[1:51]) # 0.70 ... hooray!

cor(CRUTEM5.global.annual$Anomaly[101:171], HadSST4.global.annual$Anomaly[101:171], use = "complete.obs") 

cor(GSAT.tas_land$hp_50[1:51], GSAT.tos$hp_50[1:51])
OBS_mod_p1_min = get.diff.window.cor_ens(x = OBS.tos_$GSAT$ann$mod_p1_min, y = OBS.tas_land_$GSAT$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 94:200)

OBS_mod_p1_min$cor_x.high_y.high_50


raw.data.scale$cor_x.high_y.high

cor(GMST.tas_land$mod_p1_min_50[101:171], GMST.tos$mod_p1_min_50[101:171])


mean(GMST.tas_land$mod_p1_min_50[51:(51+35)]) - mean(GMST.tos$mod_p1_min_50[51:(51+35)])



### MONTHLY:
cor(CRUTEM5.global.monthly$Anomaly[1:(51*12)], HadSST4.global.monthly$Anomaly[1:(51*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[(0*12):(51*12)], OBS.tos$GMST_FM$mon$mod_p1[(0*12):(51*12)])  # 0.37 ... hooray!
cor(OBS.tas_land$GMST_FM$mon$mod_p0[1:(51*12)], OBS.tos$GMST_FM$mon$mod_p0[1:(51*12)])  # 0.29 ... hooray!
cor(OBS.tas_land$GMST_FM$mon$mod_gta[1:(51*12)], OBS.tos$GMST_FM$mon$mod_gta[1:(51*12)])  # 0.30 ... hooray!

# 1850-2020
cor(CRUTEM5.global.monthly$Anomaly[1:(170*12)], HadSST4.global.monthly$Anomaly[1:(170*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[1:(170*12)], OBS.tos$GMST_FM$mon$mod_p1_min[1:(170*12)])












## ----------------------------------------------------------------------------------------
## 01a. Plot GSAT filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GSAT_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)

  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")

  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }

  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$mod_p1_min_2.5, 
                                               rev(GSAT.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$mod_p1_min_2.5, 
                                               rev(GSAT.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GSAT.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$lp_2.5, 
                                               rev(GSAT.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$lp_2.5, 
                                               rev(GSAT.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GSAT.tos$lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
      c = 1+0.5*2; sc = 1
      polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$hp_2.5, 
                                                 rev(GSAT.tas_land$hp_97.5)) * sc + c, 
              col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
      
      polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$hp_2.5, 
                                                 rev(GSAT.tos$hp_97.5)) * sc + c, 
              col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$hp_50 * sc + c, col = "blue", lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
    }

  # Forced response:
  {
      c = -1 + 0.5*3; sc = 1

      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$forced * sc + c, col = "darkorange", lwd = 3)
      lines(x = 1850:2020, y = GSAT.tos$forced * sc + c, col = "blue", lwd = 2)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$forced * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$forced * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$forced * sc + c, col = "darkblue", lty = 4)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$forced * sc + c, col = "darkgoldenrod4", lty = 2)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$forced * sc + c, col = "darkgoldenrod4", lty = 4)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$forced * sc + c, col = "blueviolet", lty = 2)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$forced * sc + c, col = "magenta4", lty = 1)
  }
    
  # Unforced, LP:
  {
      c = -3 + 0.5*4; sc = 1
    
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tas_land$res_lp_2.5[1:165], 
                                                 rev(GSAT.tas_land$res_lp_97.5[1:165])) * sc + c, 
              col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
      
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tos$res_lp_2.5[1:165], 
                                                 rev(GSAT.tos$res_lp_97.5[1:165])) * sc + c, 
              col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$res_lp_50 * sc + c, col = "darkorange", lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$res_lp_50 * sc + c, col = "blue", lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$res_lp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$res_lp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$res_lp_50 * sc + c, col = "darkblue", lty = 4)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$res_lp_50 * sc + c, col = "blueviolet", lty = 2)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$res_lp_50 * sc + c, col = "magenta4", lty = 1)
    }
    
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tas_land$res_hp_2.5[1:165], 
                                                 rev(GSAT.tas_land$res_hp_97.5[1:165])) * sc + c, 
              col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
      
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tos$res_hp_2.5[1:165], 
                                                 rev(GSAT.tos$res_hp_97.5[1:165])) * sc + c, 
              col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$res_hp_50 * sc + c, col = "darkorange", lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$res_hp_50 * sc + c, col = "blue", lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$res_hp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$res_hp_50 * sc + c, col = "darkblue", lty = 2)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$res_hp_50 * sc + c, col = "darkblue", lty = 4)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$res_hp_50 * sc + c, col = "blueviolet", lty = 2)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$res_hp_50 * sc + c, col = "magenta4", lty = 2)
    }
    
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
           lty = c(1, 1, 1, 2, 4, 2, 2, 4), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
           cex = 0.75, ncol = 3, bg = "white")
    
    ### 
    # cor(scale(pass.filt(y = OBS.tas_land$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
    #    scale(pass.filt(y = OBS.tos$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
    
    # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
    #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
    
  

}
dev.off()



## ----------------------------------------------------------------------------------------
## 01b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GSAT_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$mod_p1_min_2.5, 
                                               rev(GSAT.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$mod_p1_min_2.5, 
                                               rev(GSAT.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2, lwd = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4, lwd = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2, lwd = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4, lwd = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2, lwd = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$lp_2.5, 
                                               rev(GSAT.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$lp_2.5, 
                                               rev(GSAT.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$lp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$hp_2.5, 
                                               rev(GSAT.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$hp_2.5, 
                                               rev(GSAT.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$hp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()





## ----------------------------------------------------------------------------------------
## 02a. Plot GMST filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "d. Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "e. Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "f. Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$mod_p1_min_2.5, 
                                               rev(GMST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$mod_p1_min_2.5, 
                                               rev(GMST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$lp_2.5, 
                                               rev(GMST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$lp_2.5, 
                                               rev(GMST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$hp_2.5, 
                                               rev(GMST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$hp_2.5, 
                                               rev(GMST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$forced * sc + c, col = "darkorange", lwd = 3)
    lines(x = 1850:2020, y = GMST.tos$forced * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$forced * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$forced * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$forced * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$forced * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$forced * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tas_land$res_lp_2.5[1:165], 
                                               rev(GMST.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tos$res_lp_2.5[1:165], 
                                               rev(GMST.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$res_lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tas_land$res_hp_2.5[1:165], 
                                               rev(GMST.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tos$res_hp_2.5[1:165], 
                                               rev(GMST.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = "magenta4", lty = 2)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 02b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$mod_p1_min_2.5, 
                                               rev(GMST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$mod_p1_min_2.5, 
                                               rev(GMST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2, lwd = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4, lwd = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2, lwd = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4, lwd = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2, lwd = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$lp_2.5, 
                                               rev(GMST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$lp_2.5, 
                                               rev(GMST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$hp_2.5, 
                                               rev(GMST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$hp_2.5, 
                                               rev(GMST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()





## ----------------------------------------------------------------------------------------
## 03a. Plot GMSST filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMSST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean sea surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$mod_p1_min_2.5, 
                                               rev(GMSST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$mod_p1_min_2.5, 
                                               rev(GMSST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$lp_2.5, 
                                               rev(GMSST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$lp_2.5, 
                                               rev(GMSST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$hp_2.5, 
                                               rev(GMSST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$hp_2.5, 
                                               rev(GMSST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$forced * sc + c, col = "darkorange", lwd = 3)
    lines(x = 1850:2020, y = GMSST.tos$forced * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$forced * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$forced * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$forced * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$forced * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$forced * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tas_land$res_lp_2.5[1:165], 
                                               rev(GMSST.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tos$res_lp_2.5[1:165], 
                                               rev(GMSST.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$res_lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$res_lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$res_lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$res_lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$res_lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tas_land$res_hp_2.5[1:165], 
                                               rev(GMSST.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tos$res_hp_2.5[1:165], 
                                               rev(GMSST.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$res_hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$res_hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$res_hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$res_hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$res_hp_50 * sc + c, col = "magenta4", lty = 2)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 03b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMSST_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean sea surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$mod_p1_min_2.5, 
                                               rev(GMSST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$mod_p1_min_2.5, 
                                               rev(GMSST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2, lwd = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4, lwd = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2, lwd = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4, lwd = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2, lwd = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$lp_2.5, 
                                               rev(GMSST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$lp_2.5, 
                                               rev(GMSST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$lp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$hp_2.5, 
                                               rev(GMSST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$hp_2.5, 
                                               rev(GMSST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$hp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()




  



## ----------------------------------------------------------------------------------------
## 03a. Plot GMLSAT_NI filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMLSAT_NI_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean land surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$lp_2.5, 
                                               rev(GMLSAT_NI.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$lp_2.5, 
                                               rev(GMLSAT_NI.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$hp_2.5, 
                                               rev(GMLSAT_NI.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$hp_2.5, 
                                               rev(GMLSAT_NI.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$forced * sc + c, col = "darkorange", lwd = 3)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$forced * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$forced * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$forced * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$forced * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$forced * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$forced * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$forced * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tas_land$res_lp_2.5[1:165], 
                                               rev(GMLSAT_NI.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tos$res_lp_2.5[1:165], 
                                               rev(GMLSAT_NI.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$res_lp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$res_lp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$res_lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$res_lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$res_lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$res_lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$res_lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tas_land$res_hp_2.5[1:165], 
                                               rev(GMLSAT_NI.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tos$res_hp_2.5[1:165], 
                                               rev(GMLSAT_NI.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$res_hp_50 * sc + c, col = "darkorange", lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$res_hp_50 * sc + c, col = "blue", lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$res_hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$res_hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$res_hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$res_hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$res_hp_50 * sc + c, col = "magenta4", lty = 2)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 03b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMLSAT_NI_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean land surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$mod_p1_min_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$mod_p1_min_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$mod_p1_min_50 * sc + c, col = "darkblue", lty = 2, lwd = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$mod_p1_min_50 * sc + c, col = "darkblue", lty = 4, lwd = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 2, lwd = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$mod_p1_min_50 * sc + c, col = "darkgoldenrod4", lty = 4, lwd = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$mod_p1_min_50 * sc + c, col = "blueviolet", lty = 2, lwd = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$mod_p1_min_50 * sc + c, col = "magenta4", lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$lp_2.5, 
                                               rev(GMLSAT_NI.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$lp_2.5, 
                                               rev(GMLSAT_NI.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$lp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$lp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$lp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$lp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$lp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$lp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$lp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$lp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$hp_2.5, 
                                               rev(GMLSAT_NI.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color("darkorange", alpha = 30), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$hp_2.5, 
                                               rev(GMLSAT_NI.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color("blue", alpha = 30), border = make.transparent.color("blue", 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$hp_50 * sc + c, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$hp_50 * sc + c, col = "blue", lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$hp_50 * sc + c, col = "darkblue", lty = 2)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$hp_50 * sc + c, col = "darkblue", lty = 4)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$hp_50 * sc + c, col = "darkgoldenrod4", lty = 2)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$hp_50 * sc + c, col = "darkgoldenrod4", lty = 4)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$hp_50 * sc + c, col = "blueviolet", lty = 2)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$hp_50 * sc + c, col = "magenta4", lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c("darkorange", "blue", "magenta4", "darkblue", "darkblue", "blueviolet", "darkgoldenrod4", "darkgoldenrod4"),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()




