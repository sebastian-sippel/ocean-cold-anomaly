
# ------------------------------------------------------------------------------------
# Load global observational SST data for comparison.
# ------------------------------------------------------------------------------------

# +project onto reconstruction coefficients to compare on equal coverage...

# Sebastian Sippel
# 25.10.2022

library(raster)
library(ncdf4)

source("/net/h2o/climphys1/sippels/_code/tools/convert.to.eurocentric.R")

get_annual_average <- function(cur.ts) {
  return(colMeans(matrix(cur.ts, nrow = 12), na.rm=T))
}

get_annual_average_Apr_Mar <- function(cur.ts) {
  cur.ts = c(cur.ts[-c(1:3)], NA, NA, NA)
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





# 01 Load SST-based Millenium Reanalysis:
# ---------------------------------------------
# https://www.ncei.noaa.gov/access/paleo-search/study/27850
# https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/tardif2019lmr/v2_0/

library(foreach)
library(doParallel)
registerDoParallel(cores=20)

setwd("/net/h2o/climphys1/sippels/_DATASET/LMR_v2/")
# system("ncks --mk_rec_dmn time sst_MCruns_ensemble_mean_LMRv2.0.nc sst_MCruns_ensemble_mean_LMRv2.0_.nc")

# system("cdo -sellevidx,1 sst_MCruns_ensemble_mean_LMRv2.0_.nc test1.nc")
# system("ncwa -a MCrun sst_MCruns_ensemble_mean_LMRv2.0_.nc test2.nc")
# "ncks -F -d MCrun,1 sst_MCruns_ensemble_mean_LMRv2.0_.nc test3.nc"

# Get ENSEMBLE MEAN DATASET:
# system(paste("ncwa -a MCrun sst_MCruns_ensemble_mean_LMRv2.0_.nc sst_MCruns_ensemble_mean_LMRv2.0_ensmean.nc", sep=""))
# system(paste("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt sst_MCruns_ensemble_mean_LMRv2.0_ensmean.nc sst_MCruns_ensemble_mean_LMRv2.0_ensmean1.nc", sep=""))


# Process each file to 5d00:
foreach(mc=1:20) %dopar% {
  print(mc)
  system(paste("ncks -F -d MCrun,", mc, " sst_MCruns_ensemble_mean_LMRv2.0_.nc test", mc,".nc", sep=""))
  system(paste("ncwa -a MCrun test", mc,".nc test", mc,"_.nc", sep=""))
  system(paste("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt -seldate,1850-01-01,2000-12-31 test", mc,"_.nc sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,".nc -ymonmean -seldate,1961-01-01,1990-12-31 sst_MCruns_ensemble_mean_LMRv2.0_ensmean1.nc ", 
               "sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep=""))
}
# Ensemble Mean:
# system("cdo ensmean test1_.nc test2_.nc test3_.nc test4_.nc test5_.nc test6_.nc test7_.nc test8_.nc test9_.nc test10_.nc test11_.nc test12_.nc test13_.nc test14_.nc test15_.nc test16_.nc test17_.nc test18_.nc test19_.nc test20_.nc sst_MCruns_ensemble_mean_LMRv2.0_ensmean.nc")
system("rm test*")

### READ AIR/LAND MASK:

sftlf = convert.to.pacificcentric(raster("/net/h2o/climphys1/sippels/_DATA/grid/5d00_static/cmip5_masks/sftlf_g025.nc"))
# plot(sftlf)


LMRv2_GMSST = matrix(data = NA, nrow = 20, ncol = 151)

for (mc in 1:20) {
  print(mc)  
  test = brick(paste("sst_MCruns_ensemble_mean_LMRv2.0_MC", mc, "_anom.nc", sep="")) + 0
  # test = brick("sst_MCruns_ensemble_mean_LMRv2.0_.nc", var = "sst", level = mc) + 0
  test.values = values(test)  # str(test.values)
  areaw = values(area(test)) * (1- values(sftlf))
  LMRv2_GMSST[mc,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw, na.rm = T))
}
LMRv2_GMSST_years = c(1850:2000)


# Project onto different reconstructions:
LMRv2.tos.monthly = list()

for (mc in 1:20) {
  print(mc)
  
  LMRv2_sst_anom_5d00 = brick(paste("sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep="")) + 0
  LMRv2_sst_anom_5d00_ = t(values(LMRv2_sst_anom_5d00))  
  
  LMRv2.tos.monthly[[mc]] = data.frame(matrix(NA, nrow = 1812, ncol = 59))
  names(LMRv2.tos.monthly[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  LMRv2.tos.monthly[[mc]]$Year = c(rep.row(1850:2000, 12))
  LMRv2.tos.monthly[[mc]]$Month = c(rep.col(1:12, 151))
  
  aw = values(raster::area(subset(LMRv2_sst_anom_5d00, 1)))
  
  # get weighted SST average:
  for (i in 1:1812) {
    print(i)
    w = aw # * values((1 - subset(ERSSTv5_ice, i)))
    LMRv2.tos.monthly[[mc]]$Anomaly[i] = weighted.mean(x = values(subset(LMRv2_sst_anom_5d00, floor((i-1)/12)+1)), w = w, na.rm=T)
    
    # PROJECT ON DIFFERENT FINGERPRINTS:
    year.ix = match(LMRv2.tos.monthly[[mc]]$Year[i], 1850:2020); mon.ix = LMRv2.tos.monthly[[mc]]$Month[i]
    if(is.na(year.ix)) next;
    
    tos = c(matrix(LMRv2_sst_anom_5d00_[floor((i-1)/12)+1,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(LMRv2_sst_anom_5d00_[i,], 72, 36)[,36:1])
    if (any(is.na(tos))) {
      print("NA")
      tos[which(is.na(tos))] = mean(tos, na.rm=T)
    }
    # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
    
    for (cur.name.ix in 1:length(names.vec)) {
      LMRv2.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
      LMRv2.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
      LMRv2.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
      LMRv2.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
    }
  }
}


LMRv2.tos.annual = list()
for (mc in 1:20) {
  LMRv2.tos.annual[[mc]] = data.frame(matrix(NA, nrow = 151, ncol = 59))
  names(LMRv2.tos.annual[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  LMRv2.tos.annual[[mc]]$Year = c(1850:2000)
  
  for (i in 3:59) {
    LMRv2.tos.annual[[mc]][[i]] = get_annual_average(cur.ts = LMRv2.tos.monthly[[mc]][[i]])
  }
}

save(list = c("LMRv2_GMSST", "LMRv2_GMSST_years", "LMRv2.tos.monthly", "LMRv2.tos.annual"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/LMRv2.tos.RData")

# plot(LMRv2.annual[[1]]$Year, LMRv2.annual[[1]]$GSAT_mod_p1_min, type="l")
# plot(LMRv2.monthly[[1]]$Year, LMRv2.monthly[[1]]$GSAT_mod_p1_min)
# LMRv2_sstmean = LMR_sstmean
# LMRv2_sstmean_years = LMR_sstmean_years






# 02 Load Tair-based Millenium Reanalysis:
# ---------------------------------------------
# https://www.ncei.noaa.gov/access/paleo-search/study/27850
# https://www.ncei.noaa.gov/pub/data/paleo/reconstructions/tardif2019lmr/v2_0/

library(foreach)
library(doParallel)
registerDoParallel(cores=20)

setwd("/net/h2o/climphys1/sippels/_DATASET/LMR_v2/")
# system("ncks --mk_rec_dmn time air_MCruns_ensemble_mean_LMRv2.0.nc air_MCruns_ensemble_mean_LMRv2.0_.nc")

# Get ENSEMBLE MEAN DATASET:
system(paste("ncwa -a MCrun air_MCruns_ensemble_mean_LMRv2.0_.nc air_MCruns_ensemble_mean_LMRv2.0_ensmean.nc", sep=""))
system(paste("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt air_MCruns_ensemble_mean_LMRv2.0_ensmean.nc air_MCruns_ensemble_mean_LMRv2.0_ensmean1.nc", sep=""))


# Process each file to 5d00:
foreach(mc=1:20) %dopar% {
  print(mc)
  system(paste("ncks -F -d MCrun,", mc, " air_MCruns_ensemble_mean_LMRv2.0_.nc test", mc,".nc", sep=""))
  system(paste("ncwa -a MCrun test", mc,".nc test", mc,"_.nc", sep=""))
  system(paste("cdo -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -seldate,1850-01-01,2000-12-31 test", mc,"_.nc air_MCruns_ensemble_mean_LMRv2.0_MC", mc,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 air_MCruns_ensemble_mean_LMRv2.0_MC", mc,".nc -ymonmean -seldate,1961-01-01,1990-12-31 air_MCruns_ensemble_mean_LMRv2.0_ensmean1.nc ", 
               "air_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep=""))
}
# Ensemble Mean:
# system("cdo ensmean test1_.nc test2_.nc test3_.nc test4_.nc test5_.nc test6_.nc test7_.nc test8_.nc test9_.nc test10_.nc test11_.nc test12_.nc test13_.nc test14_.nc test15_.nc test16_.nc test17_.nc test18_.nc test19_.nc test20_.nc sst_MCruns_ensemble_mean_LMRv2.0_ensmean.nc")
system("rm test*")


LMRv2_GSAT = matrix(data = NA, nrow = 20, ncol = 151)
LMRv2_GMLSAT = matrix(data = NA, nrow = 20, ncol = 151)
LMRv2_GMST = matrix(data = NA, nrow = 20, ncol = 151)

for (mc in 1:20) {
  print(mc)
  
  test = brick(paste("air_MCruns_ensemble_mean_LMRv2.0_MC", mc, "_anom.nc", sep="")) + 0
  # test = brick("sst_MCruns_ensemble_mean_LMRv2.0_.nc", var = "sst", level = mc) + 0
  test.values = values(test)  # str(test.values)
  areaw = values(area(test)) * (values(sftlf))
  LMRv2_GMLSAT[mc,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw, na.rm = T))
  areaw_all = values(area(test))
  LMRv2_GSAT[mc,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw_all, na.rm = T))
}
LMRv2_GSAT_years = c(1850:2000)
LMRv2_GMLSAT_years = c(1850:2000)

# Project onto different reconstructions:
LMRv2.tas_land.monthly = list()

for (mc in 1:20) {
  print(mc)
  
  LMRv2_air_anom_5d00 = brick(paste("air_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep="")) + 0
  LMRv2_air_anom_5d00_ = t(values(LMRv2_air_anom_5d00))
  
  LMRv2.tas_land.monthly[[mc]] = data.frame(matrix(NA, nrow = 1812, ncol = 59))
  names(LMRv2.tas_land.monthly[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  LMRv2.tas_land.monthly[[mc]]$Year = c(rep.row(1850:2000, 12))
  LMRv2.tas_land.monthly[[mc]]$Month = c(rep.col(1:12, 151))
  
  aw = values(raster::area(subset(LMRv2_air_anom_5d00, 1)))
  
  # get weighted SST average:
  for (i in 1:1812) {
    print(i)
    w = aw # * values((1 - subset(ERSSTv5_ice, i)))
    LMRv2.tas_land.monthly[[mc]]$Anomaly[i] = weighted.mean(x = values(subset(LMRv2_air_anom_5d00, floor((i-1)/12)+1)), w = w, na.rm=T)
    
    # PROJECT ON DIFFERENT FINGERPRINTS:
    year.ix = match(LMRv2.tas_land.monthly[[mc]]$Year[i], 1850:2020); mon.ix = LMRv2.tas_land.monthly[[mc]]$Month[i]
    if(is.na(year.ix)) next;
    
    tas = c(matrix(LMRv2_air_anom_5d00_[floor((i-1)/12)+1,], 72, 36)[,36:1])[CMIP6.tas_land.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(LMRv2_sst_anom_5d00_[i,], 72, 36)[,36:1])
    if (any(is.na(tas))) {
      print("NA")
      tas[which(is.na(tas))] = mean(tas, na.rm=T)
    }
    # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
    
    for (cur.name.ix in 1:length(names.vec)) {
      LMRv2.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
      LMRv2.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
      LMRv2.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
      LMRv2.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tas %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
    }
  }
}


LMRv2.tas_land.annual = list()
for (mc in 1:20) {
  LMRv2.tas_land.annual[[mc]] = data.frame(matrix(NA, nrow = 151, ncol = 59))
  names(LMRv2.tas_land.annual[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  LMRv2.tas_land.annual[[mc]]$Year = c(1850:2000)
  
  for (i in 3:59) {
    LMRv2.tas_land.annual[[mc]][[i]] = get_annual_average(cur.ts = LMRv2.tas_land.monthly[[mc]][[i]])
  }
}

## aggregate GMLSAT and LMRv2_GMSST to GMST:
LMRv2_GMST = 0.67 * LMRv2_GMSST + 0.33 * LMRv2_GMLSAT
LMRv2_GMST_years = LMRv2_GMLSAT_years
# plot(LMRv2.annual[[1]]$Year, LMRv2.annual[[1]]$GSAT_mod_p1_min, type="l")
# plot(LMRv2.monthly[[1]]$Year, LMRv2.monthly[[1]]$GSAT_mod_p1_min)

save(list = c("LMRv2_GSAT", "LMRv2_GSAT_years", "LMRv2_GMLSAT", "LMRv2_GMLSAT_years", 
              "LMRv2_GMST", "LMRv2_GMST_years",
              "LMRv2.tas_land.monthly", "LMRv2.tas_land.annual"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/LMRv2.tas_land.RData")





# 03 Load Temperature Reconstructions based on Neukom etal:
# ---------------------------------------------


# 03 Load Neukom field reconstructions:
# ---------------------------------------------
setwd("/net/h2o/climphys1/sippels/_DATASET/Neukom2019_Nature_field_reconstruction/")
## Temperature reconstructions are in exact same CMIP5 grid!!


## subtract 1961-1990 average from reconstruction:
system("cdo setgrid,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt DA.nc DA_.nc")
system("cdo setgrid,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt GraphEM.nc GraphEM_.nc")

system(paste("ncwa -a lev AM.nc AM_ensmean.nc", sep=""))
system(paste("ncwa -a lev CCA.nc CCA_ensmean.nc", sep=""))
system(paste("ncwa -a lev CPS.nc CPS_ensmean.nc", sep=""))
system(paste("ncwa -O -a ens DA_.nc DA_ensmean.nc", sep=""))
system(paste("cdo -O -settaxis,0001-01-01,12:00:00,1year -settunits,days DA_ensmean.nc DA_ensmean.nc", sep=""))

# system(paste("ncatted -O -a units,time,o,c,'years since 0000-01-01 12:00:00' -a long_name,time,o,c,time DA_ensmean.nc DA_ensmean.nc", sep=""))
#system(paste("ncks -O --mk_rec_dmn time DA_ensmean.nc DA_ensmean.nc", sep=""))
#system(paste("cdo -O -settunits,days -settaxis,0001-01-01,12:00:00,1year DA_ensmean.nc DA_ensmean_test.nc", sep=""))
system(paste("ncwa -a ens GraphEM_.nc GraphEM_ensmean.nc", sep=""))
#system(paste("cdo -O -settunits,days -settaxis,0001-01-01,12:00:00,1year GraphEM_ensmean.nc GraphEM_ensmean.nc", sep=""))
system(paste("ncwa -a lev PCR.nc PCR_ensmean.nc", sep=""))



### CONTINUE HERE !!!

for (i in 1:100) {
  print(i)
  # AM.nc
  system(paste("ncks -O -F -d lev,", i, " AM.nc test", i,".nc", sep=""))
  system(paste("ncwa -O -a lev test", i,".nc test", i,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 AM_ensmean.nc AM", i, "_anom.nc", sep=""))
  system(paste("rm test", i,".nc", sep=""))
  
  # CCA.nc
  system(paste("ncks -O -F -d lev,", i, " CCA.nc test", i,".nc", sep=""))
  system(paste("ncwa -O -a lev test", i,".nc test", i,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 CCA_ensmean.nc CCA", i, "_anom.nc", sep=""))
  system(paste("rm test", i,".nc", sep=""))
  
  # CPS.nc
  system(paste("ncks -O -F -d lev,", i, " CPS.nc test", i,".nc", sep=""))
  system(paste("ncwa -O -a lev test", i,".nc test", i,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 CPS_ensmean.nc CPS", i, "_anom.nc", sep=""))
  system(paste("rm test", i,".nc", sep=""))
  
  # DA.nc
  ## DA.nc has problems!!
  
  # system(paste("ncks -O -F -d ens,", i, " DA_.nc test", i,".nc", sep=""))
  # system(paste("ncwa -O -a ens test", i,".nc test", i,".nc", sep=""))
  # system(paste("cdo -O -settaxis,0001-01-01,12:00:00,1year -settunits,days test", i,".nc test", i,".nc", sep=""))
  # system(paste("ncatted -O -a units,time,o,c,'years since 0000-01-01 12:00:00' -a long_name,time,o,c,time test", i,".nc test", i,".nc", sep=""))
  # system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 DA_ensmean.nc DA", i, "_anom.nc", sep=""))
  # system(paste("rm test", i,".nc", sep=""))

  # GraphEM.nc
  ## GraphEM has problems!!
  
  # system(paste("ncks -O -F -d ens,", i, " GraphEM.nc test", i,".nc", sep=""))
  # system(paste("ncwa -O -a ens test", i,".nc test", i,".nc", sep=""))
  # system(paste("cdo -O -settunits,days -settaxis,0001-01-01,12:00:00,1year test", i,".nc test", i,".nc", sep=""))
  # system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 GraphEM_ensmean.nc GraphEM", i, "_anom.nc", sep=""))
  # system(paste("rm test", i,".nc", sep=""))

  # PCR.nc
  system(paste("ncks -O -F -d lev,", i, " PCR.nc test", i,".nc", sep=""))
  system(paste("ncwa -O -a lev test", i,".nc test", i,".nc", sep=""))
  system(paste("cdo -ymonsub -seldate,1850-01-01,2000-12-31 test", i, ".nc -ymonmean -seldate,1961-01-01,1990-12-31 PCR_ensmean.nc PCR", i, "_anom.nc", sep=""))
  system(paste("rm test", i,".nc", sep=""))
}


sftlf = convert.to.pacificcentric(raster("/net/h2o/climphys1/sippels/_DATA/grid/5d00_static/cmip5_masks/sftlf_g025.nc"))

file.list = c(paste("AM", 1:100, "_anom.nc", sep=""), paste("CCA", 1:100, "_anom.nc", sep=""),
  paste("CPS", 1:100, "_anom.nc", sep=""), paste("PCR", 1:100, "_anom.nc", sep=""))

Neukom_GMSST = matrix(data = NA, nrow = 400, ncol = 151)
Neukom_GMLSAT = matrix(data = NA, nrow = 400, ncol = 151)
Neukom_GMST = matrix(data = NA, nrow = 400, ncol = 151)
Neukom_years = 1850:2000

for (i in 1:400) {
  print(i)  
  # test = brick(paste("AM", i, "_anom.nc", sep="")) + 0
  test = brick(file.list[i]) + 0  # plot(test, 1)
  test.values = values(test)  # str(test.values)
  areaw_SST = values(area(test)) * (1- values(sftlf))
  areaw_LAND = values(area(test)) * (values(sftlf))
  areaw = values(area(test))
  Neukom_GMSST[i,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw_SST, na.rm = T))
  Neukom_GMLSAT[i,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw_LAND, na.rm = T))
  Neukom_GMST[i,] = sapply(X = 1:151, FUN=function(i) weighted.mean(x = test.values[,i], w = areaw, na.rm = T))
}

save(list = c("Neukom_GMSST", "Neukom_GMLSAT", "Neukom_GMST"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/Neukom_fieldrec.RData")




# Project onto tos reconstruction:
Neukom.tos.monthly = list()

for (mc in 1:600) {
  print(mc)
  
  Neukom_anom_5d00 = brick(file.list[mc]) + 0   #brick(paste("sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep="")) + 0
  Neukom_anom_5d00_ = t(values(Neukom_anom_5d00))  
  
  Neukom.tos.monthly[[mc]] = data.frame(matrix(NA, nrow = 1812, ncol = 59))
  names(Neukom.tos.monthly[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  Neukom.tos.monthly[[mc]]$Year = c(rep.row(1850:2000, 12))
  Neukom.tos.monthly[[mc]]$Month = c(rep.col(1:12, 151))
  
  aw = values(raster::area(subset(Neukom_anom_5d00, 1)))
  
  # get weighted SST average:
  for (i in 1:1812) {
    # print(i)
    w = aw # * values((1 - subset(ERSSTv5_ice, i)))
    Neukom.tos.monthly[[mc]]$Anomaly[i] = weighted.mean(x = values(subset(Neukom_anom_5d00, floor((i-1)/12)+1)), w = w, na.rm=T)
    
    # PROJECT ON DIFFERENT FINGERPRINTS:
    year.ix = match(Neukom.tos.monthly[[mc]]$Year[i], 1850:2020); mon.ix = Neukom.tos.monthly[[mc]]$Month[i]
    if(is.na(year.ix)) next;
    
    tos = c(matrix(Neukom_anom_5d00_[floor((i-1)/12)+1,], 72, 36)[,36:1])[CMIP6.tos.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(LMRv2_sst_anom_5d00_[i,], 72, 36)[,36:1])
    if (any(is.na(tos))) {
      print("NA")
      tos[which(is.na(tos))] = mean(tos, na.rm=T)
    }
    # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tos_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
    
    for (cur.name.ix in 1:length(names.vec)) {
      Neukom.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
      Neukom.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
      Neukom.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
      Neukom.tos.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tos %*% CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tos.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
    }
  }
}


Neukom.tos.annual = list()
for (mc in 1:600) {
  Neukom.tos.annual[[mc]] = data.frame(matrix(NA, nrow = 151, ncol = 59))
  names(Neukom.tos.annual[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  Neukom.tos.annual[[mc]]$Year = c(1850:2000)
  
  for (i in 3:59) {
    Neukom.tos.annual[[mc]][[i]] = get_annual_average_Apr_Mar(cur.ts = Neukom.tos.monthly[[mc]][[i]])
  }
}


### save stuff...
save(list = c("Neukom.tos.monthly", "Neukom.tos.annual"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/Neukom.tos.RData")


load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/Neukom.tos.RData")

# Project onto tas_land reconstruction:
Neukom.tas_land.monthly = list()

for (mc in 1:600) {
  print(mc)
  
  Neukom_anom_5d00 = brick(file.list[mc]) + 0   #brick(paste("sst_MCruns_ensemble_mean_LMRv2.0_MC", mc,"_anom.nc", sep="")) + 0
  Neukom_anom_5d00_ = t(values(Neukom_anom_5d00))  
  
  Neukom.tas_land.monthly[[mc]] = data.frame(matrix(NA, nrow = 1812, ncol = 59))
  names(Neukom.tas_land.monthly[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  Neukom.tas_land.monthly[[mc]]$Year = c(rep.row(1850:2000, 12))
  Neukom.tas_land.monthly[[mc]]$Month = c(rep.col(1:12, 151))
  
  aw = values(raster::area(subset(Neukom_anom_5d00, 1)))
  
  # get weighted SST average:
  for (i in 1:1812) {
    # print(i)
    w = aw # * values((1 - subset(ERSSTv5_ice, i)))
    Neukom.tas_land.monthly[[mc]]$Anomaly[i] = weighted.mean(x = values(subset(Neukom_anom_5d00, floor((i-1)/12)+1)), w = w, na.rm=T)
    
    # PROJECT ON DIFFERENT FINGERPRINTS:
    year.ix = match(Neukom.tas_land.monthly[[mc]]$Year[i], 1850:2020); mon.ix = Neukom.tas_land.monthly[[mc]]$Month[i]
    if(is.na(year.ix)) next;
    
    tas_land = c(matrix(Neukom_anom_5d00_[floor((i-1)/12)+1,], 72, 36)[,36:1])[CMIP6.tas_land.df$mon[[mon.ix]]$beta$GMSST[[year.ix]]$grid.ix] # image.plot(matrix(LMRv2_sst_anom_5d00_[i,], 72, 36)[,36:1])
    if (any(is.na(tas_land))) {
      print("NA")
      tas_land[which(is.na(tas_land))] = mean(tas_land, na.rm=T)
    }
    # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMSST_v3/", HadISST.monthly$Year[i], "-", formatC(HadISST.monthly$Month[i], width = 2, format="d", flag = "0"), ".RData", sep=""))
    
    for (cur.name.ix in 1:length(names.vec)) {
      Neukom.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_min", sep="")]][i] = tas_land %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[1]
      Neukom.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p0_1se", sep="")]][i] = tas_land %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p0_a0[2]
      Neukom.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_min", sep="")]][i] = tas_land %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,1] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[1]
      Neukom.tas_land.monthly[[mc]][[paste(names.vec[cur.name.ix], "_mod_p1_1se", sep="")]][i] = tas_land %*% CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1[,2] + CMIP6.tas_land.df$mon[[mon.ix]]$beta[[names.vec[cur.name.ix]]][[year.ix]]$mod_p1_a0[2]
    }
  }
}


Neukom.tas_land.annual = list()
for (mc in 1:600) {
  Neukom.tas_land.annual[[mc]] = data.frame(matrix(NA, nrow = 151, ncol = 59))
  names(Neukom.tas_land.annual[[mc]]) = c("Year", "Month", "Anomaly", c(sapply(X = names.vec, FUN=function(x) paste(x, "_", c("mod_p0_min", "mod_p0_1se", "mod_p1_min", "mod_p1_1se"), sep=""))))
  Neukom.tas_land.annual[[mc]]$Year = c(1850:2000)
  
  for (i in 3:59) {
    Neukom.tas_land.annual[[mc]][[i]] = get_annual_average_Apr_Mar(cur.ts = Neukom.tas_land.monthly[[mc]][[i]])
  }
}


### save stuff...
save(list = c("Neukom.tas_land.monthly", "Neukom.tas_land.annual"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/Neukom.tas_land.RData")


