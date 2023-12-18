

# ------------------------------------------------------------------------------------
# Load global observational SST data for comparison.
# ------------------------------------------------------------------------------------

# +project onto reconstruction coefficients to compare on equal coverage...

# Sebastian Sippel
# 25.10.2022

library(raster)
library(ncdf4)



get_annual_average <- function(cur.ts) {
  return(colMeans(matrix(cur.ts, nrow = 12), na.rm=T))
}

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# read reconstruction coefficients \beta:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tos.df_v4.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tas_land.df_v4.RData")
# length(CMIP6.tos.df$mon[[1]]$beta$GMSST)


# Get names for coefficients:
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N", 
              "TMSST_40S_40N", "TMSST_25S_25N_", "IndianOcean", "WPacific", "EPacific", "WAtlantic")


# 00 Berkeley Earth Land :
# ---------------------------------------------
# http://berkeleyearth.lbl.gov/regions/global-land

setwd("/net/h2o/climphys1/sippels/_DATASET/BEST/monthly_1d00/")
# system("cdo -ymonsub Complete_TAVG_LatLong1_20211208.nc -ymonmean -seldate,1961-01-01,1990-12-31 Complete_TAVG_LatLong1_20211208.nc Complete_TAVG_LatLong1_20211208_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt Complete_TAVG_LatLong1_20211208_anom.nc Complete_TAVG_LatLong1_20211208_anom_5d00.nc")
BEST_Land_anom = brick("Complete_TAVG_LatLong1_20211208_anom.nc") + 0
BEST_Land_anom_5d00 = brick("Complete_TAVG_LatLong1_20211208_anom_5d00.nc") + 0
BEST_Land_anom_5d00_ = t(values(BEST_Land_anom_5d00))

BEST_Land.monthly = data.frame(matrix(NA, nrow = 3264, ncol = 59))
names(BEST_Land.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
BEST_Land.monthly$Year = c(rep.row(1750:2021, 12))
BEST_Land.monthly$Month = c(rep.col(1:12, 272))

aw = values(raster::area(subset(BEST_Land_anom, 1)))


# get weighted LAND average:
for (i in 1:3261) {
  print(i)
  w = aw
  BEST_Land.monthly$Anomaly[i] = weighted.mean(x = values(subset(BEST_Land_anom, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(BEST_Land.monthly$Year[i], 1850:2020); mon.ix = BEST_Land.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  tas = c(matrix(BEST_Land_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tas_land.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix]    # image.plot(matrix(HadISST_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tas))) {
    print("NA")
    tas[which(is.na(tas))] = mean(tas, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    BEST_Land.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    BEST_Land.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    BEST_Land.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    BEST_Land.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}


BEST_Land.annual = data.frame(matrix(NA, nrow = 272, ncol = 59))
names(BEST_Land.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
BEST_Land.annual$Year = c(1750:2021)

for (i in 3:59) {
  BEST_Land.annual[[i]] = get_annual_average(cur.ts = BEST_Land.monthly[[i]])
}




# 01 HadISST:
# ---------------------------------------------
# https://www.metoffice.gov.uk/hadobs/hadisst/data/download.html

# load.files("/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt")

setwd("/net/h2o/climphys1/sippels/_DATASET/HadISST/")

# system("cd /net/h2o/climphys1/sippels/_DATASET/HadISST/")
# system("cdo setrtomiss,-10000,-3 HadISST_sst.nc HadISST_sst_.nc")
# system("cdo -ymonsub HadISST_sst_.nc -ymonmean -seldate,1961-01-01,1990-12-31 HadISST_sst_.nc HadISST_sst_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt HadISST_sst_anom.nc HadISST_sst_anom_5d00.nc")

HadISST_ice = brick("HadISST_ice.nc") + 0
HadISST_sst_anom = brick("HadISST_sst_anom.nc") + 0
HadISST_sst_anom_5d00 = brick("HadISST_sst_anom_5d00.nc") + 0
HadISST_sst_anom_5d00_ = t(values(HadISST_sst_anom_5d00))

HadISST.monthly = data.frame(matrix(NA, nrow = 1836, ncol = 59))
names(HadISST.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
HadISST.monthly$Year = c(rep.row(1870:2022, 12))
HadISST.monthly$Month = c(rep.col(1:12, 153))

aw = values(raster::area(subset(HadISST_sst_anom, 1)))

# get weighted SST average:
for (i in 1:1833) {
  print(i)
  w = aw * values((1 - subset(HadISST_ice, i)))
  HadISST.monthly$Anomaly[i] = weighted.mean(x = values(subset(HadISST_sst_anom, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(HadISST.monthly$Year[i], 1850:2020); mon.ix = HadISST.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  tos = c(matrix(HadISST_sst_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(HadISST_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tos))) {
    print("NA")
    tos[which(is.na(tos))] = mean(tos, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    HadISST.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    HadISST.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    HadISST.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    HadISST.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}

HadISST.annual = data.frame(matrix(NA, nrow = 153, ncol = 59))
names(HadISST.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
HadISST.annual$Year = c(1870:2022)

for (i in 3:59) {
  HadISST.annual[[i]] = get_annual_average(cur.ts = HadISST.monthly[[i]])
}


# plot(x = HadISST.annual$Year, y = HadISST.annual$Anomaly, type='l')
# lines(x = HadISST.annual$Year, y = HadISST.annual$mod_p1_min, col = "red")




# 02. COBE-SST
# ---------------------------------------------
# https://psl.noaa.gov/data/gridded/data.cobe2.html

setwd("/net/h2o/climphys1/sippels/_DATASET/COBE-SST2/")

#system("cd /net/h2o/climphys1/sippels/_DATASET/COBE-SST2/")
# system("cdo -ymonsub sst.mon.mean.nc -ymonmean -seldate,1961-01-01,1990-12-31 sst.mon.mean.nc sst.mon.mean_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt sst.mon.mean_anom.nc sst.mon.mean_anom_5d00.nc")

COBE_SST2_ice = brick("icec.mon.mean.nc") + 0
COBE_SST2_sst_anom = brick("sst.mon.mean_anom.nc") + 0
COBE_SST2_sst_anom_5d00 = brick("sst.mon.mean_anom_5d00.nc") + 0
COBE_SST2_sst_anom_5d00_ = t(values(COBE_SST2_sst_anom_5d00))

COBE_SST2.monthly = data.frame(matrix(NA, nrow = 2040, ncol = 59))
names(COBE_SST2.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
COBE_SST2.monthly$Year = c(rep.row(1850:2019, 12))
COBE_SST2.monthly$Month = c(rep.col(1:12, 170))

aw = values(raster::area(subset(COBE_SST2_sst_anom, 1)))

# get weighted SST average:
for (i in 1:2040) {
  print(i)
  w = aw * values((1 - subset(COBE_SST2_ice, i)))
  COBE_SST2.monthly$Anomaly[i] = weighted.mean(x = values(subset(COBE_SST2_sst_anom, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(COBE_SST2.monthly$Year[i], 1850:2020); mon.ix = COBE_SST2.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  tos = c(matrix(COBE_SST2_sst_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(COBE_SST2_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tos))) {
    print("NA")
    tos[which(is.na(tos))] = mean(tos, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    COBE_SST2.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    COBE_SST2.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    COBE_SST2.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    COBE_SST2.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}

COBE_SST2.annual = data.frame(matrix(NA, nrow = 170, ncol = 59))
names(COBE_SST2.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
COBE_SST2.annual$Year = c(1850:2019)

for (i in 3:59) {
  COBE_SST2.annual[[i]] = get_annual_average(cur.ts = COBE_SST2.monthly[[i]])
}



# 03. ERSST
# ---------------------------------------------
# Version 5: https://psl.noaa.gov/data/gridded/data.noaa.ersst.v5.html
setwd("/net/h2o/climphys1/sippels/_DATASET/ERSST/v5/")
# system("cd /net/h2o/climphys1/sippels/_DATASET/ERSST/v5/")
# system("cdo -ymonsub /net/h2o/climphys1/sippels/_DATASET/ERSST/v5/sst.mnmean.nc -ymonmean -seldate,1961-01-01,1990-12-31 /net/h2o/climphys1/sippels/_DATASET/ERSST/v5/sst.mnmean.nc /net/h2o/climphys1/sippels/_DATASET/ERSST/v5/sst.mnmean_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt sst.mnmean_anom.nc sst.mnmean_anom_5d00.nc")


ERSSTv5_sst_anom = brick("sst.mnmean_anom.nc") + 0
ERSSTv5_sst_anom_5d00 = brick("sst.mnmean_anom_5d00.nc") + 0
ERSSTv5_sst_anom_5d00_ = t(values(ERSSTv5_sst_anom_5d00))

ERSSTv5.monthly = data.frame(matrix(NA, nrow = 2028, ncol = 59))
names(ERSSTv5.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
ERSSTv5.monthly$Year = c(rep.row(1854:2022, 12))
ERSSTv5.monthly$Month = c(rep.col(1:12, 169))

aw = values(raster::area(subset(ERSSTv5_sst_anom, 1)))

# get weighted SST average:
for (i in 1:2025) {
  print(i)
  w = aw # * values((1 - subset(ERSSTv5_ice, i)))
  ERSSTv5.monthly$Anomaly[i] = weighted.mean(x = values(subset(ERSSTv5_sst_anom, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(ERSSTv5.monthly$Year[i], 1850:2020); mon.ix = ERSSTv5.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  tos = c(matrix(ERSSTv5_sst_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(ERSSTv5_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tos))) {
    print("NA")
    tos[which(is.na(tos))] = mean(tos, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    ERSSTv5.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    ERSSTv5.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    ERSSTv5.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    ERSSTv5.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}

ERSSTv5.annual = data.frame(matrix(NA, nrow = 169, ncol = 59))
names(ERSSTv5.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
ERSSTv5.annual$Year = c(1854:2022)

for (i in 3:59) {
  ERSSTv5.annual[[i]] = get_annual_average(cur.ts = ERSSTv5.monthly[[i]])
}


# 04. ERSST
# ---------------------------------------------
# Version 4: https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v4/
setwd("/net/h2o/climphys1/sippels/_DATASET/ERSST/v4/")
#system("cd /net/h2o/climphys1/sippels/_DATASET/ERSST/v4/")
# system("cdo -ymonsub /net/h2o/climphys1/sippels/_DATASET/ERSST/v4/sst.mnmean.nc -ymonmean -seldate,1961-01-01,1990-12-31 /net/h2o/climphys1/sippels/_DATASET/ERSST/v4/sst.mnmean.nc /net/h2o/climphys1/sippels/_DATASET/ERSST/v4/sst.mnmean_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt sst.mnmean_anom.nc sst.mnmean_anom_5d00.nc")

ERSSTv4_sst_anom = brick("sst.mnmean_anom.nc") + 0
ERSSTv4_sst_anom_5d00 = brick("sst.mnmean_anom_5d00.nc") + 0
ERSSTv4_sst_anom_5d00_ = t(values(ERSSTv4_sst_anom_5d00))

ERSSTv4.monthly = data.frame(matrix(NA, nrow = 2004, ncol = 59))
names(ERSSTv4.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
ERSSTv4.monthly$Year = c(rep.row(1854:2020, 12))
ERSSTv4.monthly$Month = c(rep.col(1:12, 167))

aw = values(raster::area(subset(ERSSTv4_sst_anom, 1)))

# get weighted SST average:
for (i in 1:1994) {
  print(i)
  w = aw # * values((1 - subset(ERSSTv5_ice, i)))
  ERSSTv4.monthly$Anomaly[i] = weighted.mean(x = values(subset(ERSSTv4_sst_anom, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(ERSSTv4.monthly$Year[i], 1850:2020); mon.ix = ERSSTv4.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  tos = c(matrix(ERSSTv4_sst_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(ERSSTv5_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tos))) {
    print("NA")
    tos[which(is.na(tos))] = mean(tos, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    ERSSTv4.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    ERSSTv4.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    ERSSTv4.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    ERSSTv4.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}


ERSSTv4.annual = data.frame(matrix(NA, nrow = 167, ncol = 59))
names(ERSSTv4.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
ERSSTv4.annual$Year = c(1854:2020)

for (i in 3:59) {
  ERSSTv4.annual[[i]] = get_annual_average(cur.ts = ERSSTv4.monthly[[i]])
}






# 05 HadSST4 - UNADJUSTED:
# ---------------------------------------------
# https://www.metoffice.gov.uk/hadobs/hadsst4/data/download.html

## Get conversion to CMIP6 grid and distance matrix for rasterbrick format:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_code/tools/convert.to.eurocentric.R")

setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT5/5d00_monthly_noninfilled/")
raster.template_CRU = raster(subset(brick("HadCRUT.5.0.1.0.weights.nc"), 1000))
raster.template.coords = coordinates(raster.template_CRU)

# match coordinates... and transfer:
x.coord.CRU = c(matrix(coordinates(raster.template_CRU)[,1], 72, 36)[,36:1])
y.coord.CRU = c(matrix(coordinates(raster.template_CRU)[,2], 72, 36)[,36:1])
x.coord = c(matrix(coordinates(raster.template)[,1], 72, 36)[,36:1])
y.coord = c(matrix(coordinates(raster.template)[,2], 72, 36)[,36:1])

transfer.CRU.to.CMIP6.grid = rep(NA, 72*36)
for (i in 1:(72*36)) {
  # ix = which( (x.coord[i] == x.coord.CRU | x.coord[i] == x.coord.CRU - 360 ) & y.coord[i] == y.coord.CRU )
  transfer.CRU.to.CMIP6.grid[i] = which( (x.coord == x.coord.CRU[i] | x.coord == x.coord.CRU[i] + 360 ) & y.coord == y.coord.CRU[i] )
}

# load.files("/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt")

setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/other/")

# system("cd /net/h2o/climphys1/sippels/_DATASET/HadISST/")
# system("cdo setrtomiss,-10000,-3 HadISST_sst.nc HadISST_sst_.nc")
# system("cdo -ymonsub HadISST_sst_.nc -ymonmean -seldate,1961-01-01,1990-12-31 HadISST_sst_.nc HadISST_sst_anom.nc")
# system("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt HadISST_sst_anom.nc HadISST_sst_anom_5d00.nc")

# HadSST_ua = brick("HadSST.4.0.1.0_unadjusted.nc") + 0
HadSST_ua_sst_anom = brick("HadSST.4.0.1.0_unadjusted.nc") + 0
HadSST_ua_sst_anom_5d00 = convert.to.pacificcentric(brick("HadSST.4.0.1.0_unadjusted.nc") + 0)
HadSST_ua_sst_anom_5d00_ = t(apply(get_CMIP5_array(file = paste("HadSST.4.0.1.0_unadjusted.nc", sep=""), var="tos", time.count = -1), MARGIN = 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
# image.plot(matrix(test[2000,], 72, 36))

HadSST_ua.monthly = data.frame(matrix(NA, nrow = 2052, ncol = 59))
names(HadSST_ua.monthly) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
HadSST_ua.monthly$Year = c(rep.row(1850:2020, 12))
HadSST_ua.monthly$Month = c(rep.col(1:12, 171))

aw = values(raster::area(subset(HadSST_ua_sst_anom_5d00, 1)))

# get weighted SST average:
for (i in 1:2052) {
  print(i)
  w = aw
  HadSST_ua.monthly$Anomaly[i] = weighted.mean(x = values(subset(HadSST_ua_sst_anom_5d00, i)), w = w, na.rm=T)
  
  # PROJECT ON DIFFERENT FINGERPRINTS:
  year.ix = match(HadSST_ua.monthly$Year[i], 1850:2020); mon.ix = HadSST_ua.monthly$Month[i]
  if(is.na(year.ix)) next;
  
  # tos = c(matrix(HadSST_ua_sst_anom_5d00_[i,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(HadISST_sst_anom_5d00_[i,], 72, 36)[,36:1])
  tos = c(matrix(HadSST_ua_sst_anom_5d00_[i,], 72, 36))[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(HadISST_sst_anom_5d00_[i,], 72, 36)[,36:1])
  if (any(is.na(tos))) {
    print("NA")
    tos[which(is.na(tos))] = mean(tos, na.rm=T)
  }
  # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
  
  for (cur.name.ix in 1:length(names.vec)) {
    HadSST_ua.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
    HadSST_ua.monthly[[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
    HadSST_ua.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
    HadSST_ua.monthly[[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
  }
}

HadSST_ua.annual = data.frame(matrix(NA, nrow = 171, ncol = 59))
names(HadSST_ua.annual) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
HadSST_ua.annual$Year = c(1850:2020)

for (i in 3:59) {
  HadSST_ua.annual[[i]] = get_annual_average(cur.ts = HadSST_ua.monthly[[i]])
}

# plot(HadSST_ua.annual$GSAT_mod_p1_min)




# plot different SST analysis:
plot(HadISST.annual$Year, HadISST.annual$Anomaly, xlim = c(1850, 2022), col = "darkgray", type='l', ylim = c(-0.6, 0.6))
lines(x = HadISST.annual$Year, y = HadISST.annual$GMSST_mod_p1_min, col = "darkgray", lty = 2)
lines(COBE_SST2.annual$Year, COBE_SST2.annual$Anomaly, col = "darkblue")
lines(COBE_SST2.annual$Year, COBE_SST2.annual$GMSST_mod_p1_min, col = "darkblue", lty = 2)
lines(ERSSTv5.annual$Year, ERSSTv5.annual$Anomaly, col = "darkred")
lines(ERSSTv5.annual$Year, ERSSTv5.annual$GMSST_mod_p1_min, col = "darkred", lty = 2)
lines(ERSSTv4.annual$Year, ERSSTv4.annual$Anomaly, col = "red")
lines(ERSSTv4.annual$Year, ERSSTv4.annual$GMSST_mod_p1_min, col = "red", lty = 2)


save(list = c("HadISST.annual", "HadISST.monthly", "COBE_SST2.annual", "COBE_SST2.monthly",
              "ERSSTv5.annual", "ERSSTv5.monthly", "ERSSTv4.annual", "ERSSTv4.monthly",
              "HadSST_ua.annual", "HadSST_ua.monthly"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/SST_datasets.RData")

save(list = c("BEST_Land.annual", "BEST_Land.monthly"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/TLand_datasets.RData")





### ------------------------------------------------------------
### Load PAGES2k global field reconstructions (???)
### ------------------------------------------------------------






