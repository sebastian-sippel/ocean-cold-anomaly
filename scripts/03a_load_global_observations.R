
# ------------------------------------------------------------------------------------
# Load CRU global observational estimates for comparison.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021


# 00.(b) load CRUTEM5 global(!!) data:
{
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/CRUTEM5/5d00_monthly/")
  
  CRUTEM5.global.annual = read.table("CRUTEM.5.0.1.0.summary_series.global.annual.csv", sep = ",", header = T)
  names(CRUTEM5.global.annual) <- c("Time", "Anomaly", "lower_CI2.5", "upper_CI97.5")
  CRUTEM5.global.annual_comp = read.table("CRUTEM.5.0.1.0.component_series.global.annual.csv", sep = ",", header = T)
  
  CRUTEM5.global.monthly = read.table("CRUTEM.5.0.1.0.summary_series.global.monthly.csv", sep = ",", header = T)
  names(CRUTEM5.global.monthly) <- c("Time", "Anomaly", "lower_CI2.5", "upper_CI97.5")
  CRUTEM5.global.monthly_comp = read.table("CRUTEM.5.0.1.0.component_series.global.monthly.csv", sep = ",", header = T)
}


# 00.(d) load HadSST4 global(!!) data:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")
{
  HadSST4.global.annual = read.table("HadSST.4.0.1.0_annual_GLOBE.csv", sep = ",", header = T)
  names(HadSST4.global.annual) <- c("Year", "Anomaly", "total_uncertainty", "uncorrelated_uncertainty", 
                                    "correlated_uncertainty", "bias_uncertainty", "coverage_uncertainty", "lower_CI2.5", "upper_CI97.5")
  HadSST4.global.monthly = read.table("HadSST.4.0.1.0_monthly_GLOBE.csv", sep = ",", header = T)
  names(HadSST4.global.monthly) <- c("Year", "Month", "Anomaly", "total_uncertainty", "uncorrelated_uncertainty", 
                                     "correlated_uncertainty", "bias_uncertainty", "coverage_uncertainty", "lower_CI2.5", "upper_CI97.5")
}


# 00.(e) load HadCRUT5:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT5/time_series/")
{
  HadCRUT5.global.annual = read.table("HadCRUT.5.0.1.0.analysis.summary_series.global.annual.csv", sep = ",", header = T)
  names(HadCRUT5.global.annual) <- c("Year", "Anomaly", "lower_CI2.5", "upper_CI97.5") 
  
  HadCRUT5.global.monthly = read.table("HadCRUT.5.0.1.0.analysis.summary_series.global.monthly.csv", sep = ",", header = T)
  names(HadCRUT5.global.monthly) <- c("Year", "Anomaly", "lower_CI2.5", "upper_CI97.5") 
}

# plot(HadSST4.global.annual$anomaly, type='l')
# lines(HadSST4.global.annual$upper_bound_95pct_bias_uncertainty_range, type='l', col = "blue")
# lines(HadSST4.global.annual$lower_bound_95pct_bias_uncertainty_range, type='l', col = "blue")






# ------------------------------------------------------------------------------------
# Load *ALL OTHER* Observational Estimates:
# ------------------------------------------------------------------------------------

get_annual_average <- function(cur.ts) {
  return(colMeans(matrix(cur.ts, nrow = 12), na.rm=T))
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}



# 02. Cowtan+Way Updated to 08/2022:
# ---------------------------------------------
# https://www-users.york.ac.uk/~kdc3/papers/coverage2013/series.html
setwd("/net/h2o/climphys1/sippels/_DATASET/Cowtan_Way2014/_orig/v202208/")

# CW2014-HadCRUT4-kriging:
CW2014.global.monthly_had4_krig = read.table("had4_krig_v2_0_0.csv", sep=",", header = F)
names(CW2014.global.monthly_had4_krig) = c("Year", "Anomaly", "V3", "V4", "V5")

CW2014.global.annual_had4_krig = read.table("had4_krig_annual_v2_0_0.csv", sep=",", header = F)
names(CW2014.global.annual_had4_krig) = c("Year", "Anomaly", "V3", "V4", "V5")

# CW2014-cobe2cru-kriging:
CW2014.global.monthly_cobe2cru_krig = read.table("cobe2cru_krig_v2_0_0.csv", sep=",", header = F)
names(CW2014.global.monthly_cobe2cru_krig) = c("Year", "Anomaly")
# plot(CW2014.global.annual_had4_krig$Year, CW2014.global.annual_had4_krig$Anomaly)
CW2014.global.annual_cobe2cru_krig = data.frame(matrix(data = NA, nrow = length(CW2014.global.monthly_cobe2cru_krig$Anomaly)/12, ncol = 2))
names(CW2014.global.annual_cobe2cru_krig) = c("Year", "Anomaly")
CW2014.global.annual_cobe2cru_krig$Year = get_annual_average(as.numeric(CW2014.global.monthly_cobe2cru_krig$Year)) - 0.5
CW2014.global.annual_cobe2cru_krig$Year[1] = 1850
CW2014.global.annual_cobe2cru_krig$Anomaly = get_annual_average(CW2014.global.monthly_cobe2cru_krig$Anomaly)

# CW2014-had4sst4-kriging:
CW2014.global.monthly_had4sst4_krig = read.table("had4sst4_krig_v2_0_0.csv", sep=",", header = F)
names(CW2014.global.monthly_had4sst4_krig) = c("Year", "Anomaly")
CW2014.global.annual_had4sst4_krig = data.frame(matrix(data = NA, nrow = length(CW2014.global.monthly_had4sst4_krig$Anomaly)/12, ncol = 2))
names(CW2014.global.annual_had4sst4_krig) = c("Year", "Anomaly")
CW2014.global.annual_had4sst4_krig$Year = get_annual_average(as.numeric(CW2014.global.monthly_had4sst4_krig$Year)) - 0.5
CW2014.global.annual_had4sst4_krig$Year[1] = 1850
CW2014.global.annual_had4sst4_krig$Anomaly = get_annual_average(CW2014.global.monthly_had4sst4_krig$Anomaly)

## CW2014 ENSEMBLE:
setwd("/net/h2o/climphys1/sippels/_DATASET/CW2014/")
CW2014.global.monthly_had4_krig_ENS = read.table("had4_krig_ensemble_v2_0_0.txt", header = F)
CW2014.global.annual_had4_krig_ENS = apply(X = CW2014.global.monthly_had4_krig_ENS[1:2052,], MARGIN = 2, FUN=get_annual_average)


# 03. NASA GISSTEMP:
# ---------------------------------------------
# https://data.giss.nasa.gov/gistemp/
setwd("/net/h2o/climphys1/sippels/_DATASET/GISS/_orig/v202208/")

GISS.global.monthly = read.table("GLB.Ts+dSST.csv", sep=",", header = T, skip = 1, nrows=142)
GISS.global.annual = data.frame(matrix(data = NA, nrow = length(GISS.global.monthly$Year), ncol = 2))
names(GISS.global.annual) = c("Year", "Anomaly")
GISS.global.annual$Year = GISS.global.monthly$Year
GISS.global.annual$Anomaly = GISS.global.monthly$J.D

# AIRS Temperature data:
AIRS_v6.global.monthly = read.table("GISSTEMP_AIRS_GLB.Ts+dSST.csv", sep=",", header = T, skip = 1, nrows=21)
AIRS_v7.global.monthly = read.table("GISSTEMP_AIRS_GLB.Ts+dSST.csv", sep=",", header = T, skip = 24, nrows=21)

# AIRS_v7.global.monthly = data.frame(matrix(data = NA, nrow = length(GISS.global.monthly$Year), ncol = 2))
# get_annual_average(cur.ts = as.numeric(t(as.matrix(AIRS_v7.global.monthly[,2:13]))))
GISS.global.monthly.2002_2022 = read.table("GISSTEMP_AIRS_GLB.Ts+dSST.csv", sep=",", header = T, skip = 47, nrows=21)


# 04. BerkeleyEarth:
# ---------------------------------------------
# http://berkeleyearth.org/data/
setwd("/net/h2o/climphys1/sippels/_DATASET/BEST/_orig/v202208/")
BEST.global.monthly = read.table("Land_and_Ocean_complete.csv", sep=",", header = F, skip = 86, nrows = 2064)
names(BEST.global.monthly) = c("Year", "Month", "Anomaly", "Unc.", "1y-Anomaly", "1y-Unc.", "5y-Anomaly", "5y-Unc.", "10y-Anomaly", "10y-Unc.", "20y-Anomaly", "20y-Unc.")
BEST.global.annual = data.frame(matrix(data = NA, nrow = length(BEST.global.monthly$Year)/12, ncol = 2))
names(BEST.global.annual) = c("Year", "Anomaly")
BEST.global.annual$Year = get_annual_average(cur.ts = BEST.global.monthly$Year)
BEST.global.annual$Anomaly = get_annual_average(cur.ts = BEST.global.monthly$Anomaly)

BEST.land.monthly = read.table("Complete_TAVG_complete.csv", sep=",", header = F, skip = 35, nrows = 3264)
names(BEST.land.monthly) = c("Year", "Month", "Anomaly", "Unc.", "1y-Anomaly", "1y-Unc.", "5y-Anomaly", "5y-Unc.", "10y-Anomaly", "10y-Unc.", "20y-Anomaly", "20y-Unc.")
BEST.land.annual = data.frame(matrix(data = NA, nrow = length(BEST.land.monthly$Year)/12, ncol = 2))
names(BEST.land.annual) = c("Year", "Anomaly")
BEST.land.annual$Year = get_annual_average(cur.ts = BEST.land.monthly$Year)
BEST.land.annual$Anomaly = get_annual_average(cur.ts = BEST.land.monthly$Anomaly)


# 05. NOAA _ GobalTemperatureData
# ------------------------------------------
# Version 4: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00934
# Version 5: https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v5/access/timeseries/
# Version 5 is the correct one 

setwd("/net/h2o/climphys1/sippels/_DATASET/NOAA_MLOST/v5/_orig/")
NOAA.global.annual = read.table("aravg.ann.land_ocean.90S.90N.v5.0.0.202207.csv", sep=",", header = F)
names(NOAA.global.annual) = c("Year", "Anomaly", "V3", "V4", "V5", "V6")
NOAA.global.monthly = read.table("aravg.mon.land_ocean.90S.90N.v5.0.0.202207.csv", sep=",", header = F)
names(NOAA.global.monthly) = c("Year", "Month", "Anomaly", "V3", "V4", "V5", "V6", "V7", "V8", "V9")

# read (separate) land- and ocean data:
# setwd("/net/h2o/climphys1/sippels/_DATASET/NOAA_MLOST/v5/_orig/")

# 06. ERA5.
# setwd("/net/h2o/")

# -> only for late period...


# 07. Japanese Temperature Average
# ------------------------------------------
# https://ds.data.jma.go.jp/tcc/tcc/products/gwp/temp/list/year_wld.html
setwd("/net/h2o/climphys1/sippels/_DATASET/JMA_GTA/_orig/v202208/")
JMA.global.annual = read.table("year_wld.csv", sep=",", header = T)

# 08. Kadow et al., 2020
# ------------------------------------------
# https://www.nature.com/articles/s41561-020-0582-5#Sec14
setwd("/net/h2o/climphys1/sippels/_DATASET/Kadow/FREVA-CLINT-climatereconstructionAI-75e7231/reconstructions/")

library(raster)
AI20cr = brick("20crAI_HadCRUT4_4.6.0.0_tas_mon_185001-201812.nc") + 0 
AIcmip = brick("cmipAI_HadCRUT4_4.6.0.0_tas_mon_185001-201812.nc") + 0 
Kadow.global.monthly <- data.frame(matrix(data = NA, nrow = 2028, ncol = 4))
names(Kadow.global.monthly) <- c("Year", "Month", "Anomaly_20crAI", "Anomaly_cmipAI")
Kadow.global.monthly$Year <- c(rep.row(x = 1850:2018, n = 12))
Kadow.global.monthly$Month <- c(rep(1:12, 2028/12))
Kadow.global.monthly$Anomaly_20crAI = t(values(AI20cr)) %*% values(raster::area(AI20cr)) / sum(values(raster::area(AI20cr)))
Kadow.global.monthly$Anomaly_cmipAI = t(values(AIcmip)) %*% values(raster::area(AIcmip)) / sum(values(raster::area(AIcmip)))
# Annual Anomaly: 
Kadow.global.annual = data.frame(matrix(data = NA, nrow = length(Kadow.global.monthly$Month)/12, ncol = 3))
names(Kadow.global.annual) = c("Year", "Anomaly_20crAI", "Anomaly_cmipAI")
Kadow.global.annual$Year = get_annual_average(cur.ts = Kadow.global.monthly$Year)
Kadow.global.annual$Anomaly_20crAI = get_annual_average(cur.ts = Kadow.global.monthly$Anomaly_20crAI)
Kadow.global.annual$Anomaly_cmipAI = get_annual_average(cur.ts = Kadow.global.monthly$Anomaly_cmipAI)



# 09. CR20 for comparison:
# ------------------------------------------
# setwd("/net/h2o/climphys1/sippels/_DATASET/20CR/V3/tas/2mMO/")
# test = nc_open("air.2m.mon.mean_gl.nc")
# ncvar_get(nc = test, varid = "air") - 273.15


# 10. Load Cowtan hybrid36 dataset:
## Read Cowtan-hybrid36:
hybrid36 = read.table(file = "/net/h2o/climphys1/sippels/_DATASET/Cowtan_etal_2017_SSTBIAS/results/coastal_wt_y0/hybrid_36m.temp")
names(hybrid36) = c("Year", "Anomaly")
hybrid36.annual = data.frame(cbind(1850:2016, sapply(X = 1850:2016, FUN=function(year) mean(hybrid36$Anomaly[which(floor(hybrid36$Year) == year)]))))
names(hybrid36.annual) = c("Year", "Anomaly")



setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/figures/")


