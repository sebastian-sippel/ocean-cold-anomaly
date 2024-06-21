

# ------------------------------------------------------------
### 1. GET ALL CRU data + land masks: 
# ------------------------------------------------------------

# Sebastian Sippel
# 10.08.2019
# 09 / 2021

library(raster)
library(ncdf4)
library(MASS)
library(geosphere)

# Read functions: 
source("/net/h2o/climphys1/sippels/_code/tools/convert.to.eurocentric.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr/code/_functions_CMIP5_extr.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr/code/_functions_CRU.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_functions_CMIP5_extr.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")

## Get conversion to CMIP6 grid and distance matrix for rasterbrick format:
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

# cur.land.weights.CMIP6 = cur.land.weights[transfer.CRU.to.CMIP6.grid]
# image.plot(matrix(cur.land.weights[transfer.CRU.to.CMIP6.grid], 72, 36))


# get distance matrix:
# distGeo(p1 = raster.template.coords[1,], p2 = raster.template.coords[2000,])/1000
# raster.template.distm = distm(x = raster.template.coords, y = raster.template.coords)
# ?distm
# image.plot(raster.template.distm, main="Distance Matrix - raster")

## Get distance matrix for array format:
# x.coord = c(matrix(raster.template.coords[,1], 72, 36)[,36:1])
# y.coord = c(matrix(raster.template.coords[,2], 72, 36)[,36:1])
# xy.coord = data.frame(cbind(x = x.coord, y = y.coord))
# distm_ar = distm(x = xy.coord, y = xy.coord) / 1000 ## distance in km
# image.plot(distm_ar, main="Distance Matrix - array")

## GET RANDOM POINT FROM DISTANCE MATRIX (raster):
# dist.random.point = raster.template
# point.ix = 1500
# values(dist.random.point) = raster.template.distm[point.ix,]
# plot(dist.random.point)




# CRUTEM5 
# ------------------------------
{
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/CRUTEM5/5d00_monthly/")
  # CRUTEM5 = t(apply(X = get_CMIP5_array("CRUTEM.5.0.1.0.anomalies.nc", var="tas", time.count = -1), 3, as.vector))  # image.plot(matrix(CRUTEM5[1,], nrow=72, ncol = 36))
  # CRUTEM5_cr = t(apply(X = get_CMIP5_array("CRUTEM.5.0.1.0.anomalies.nc", var="tas", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]  # image.plot(matrix(CRUTEM5[1,], nrow=72, ncol = 36))
  # image.plot(matrix(CRUTEM5[2000,], 72, 36))
  # image.plot(matrix(CRUTEM5_cr[2000,], 72, 36))
  CRUTEM5 = t(apply(X = get_CMIP5_array("CRUTEM.5.0.1.0.anomalies.nc", var="tas", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CRUTEM5_brick = convert.to.pacificcentric(brick("CRUTEM.5.0.1.0.anomalies.nc", varname="tas"))
  CRUTEM5_nobs = t(apply(get_CMIP5_array(file = "CRUTEM.5.0.1.0.station_counts.nc", var = "tas_nobs", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CRUTEM5_nobs_brick = convert.to.pacificcentric(brick("CRUTEM.5.0.1.0.station_counts.nc", varname="tas_nobs") + 0)
  CRUTEM5_calendar = as.Date(substring(names(CRUTEM5_nobs_brick), first = 2), format = "%Y.%m.%d")
  
  # Get CRUTEM5 uncertainties:
  CRUTEM5_sampling_unc = t(apply(X = get_CMIP5_array("CRUTEM.5.0.1.0.measurement_sampling.nc", var = "tas_unc", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CRUTEM5_station_unc = t(apply(X = get_CMIP5_array("CRUTEM.5.0.1.0.station_uncertainty.nc", var = "tas_corr", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]   # not needed as this appears to be encoded in the Ensemble.
  ## CRUTEM4_unc = sqrt(CRUTEM4_sampling_unc ^ 2 + CRUTEM4_station_unc ^ 2)
  
  
  # these are measurement uncertainties, uncorrelated across the land dataset:
  # plot(HadCRUT4_landunc[cur.date.ix,cur.mask.ix], CRUTEM4_sampling_unc[cur.date.ix,cur.mask.ix])
  # lines(c(0,1), c(0,1))
  # plot(HadCRUT4_landunc[cur.date.ix,cur.mask.ix],   sqrt(CRUTEM4_sampling_unc[cur.date.ix,cur.mask.ix]^2+CRUTEM4_station_unc[cur.date.ix,cur.mask.ix]^2))
  # lines(c(0,1), c(0,1))
}  
  

# HadCRUT5 data:
# ------------------------------
{
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT5/5d00_monthly_noninfilled/")
  HadCRUT5 = t(apply(get_CMIP5_array("HadCRUT.5.0.1.0.anomalies.ensemble_mean.nc", var="tas_mean", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadCRUT5_brick = convert.to.pacificcentric(brick("HadCRUT.5.0.1.0.anomalies.ensemble_mean.nc", varname="tas_mean") + 0)
  HadCRUT5_calendar = as.Date(substring(names(HadCRUT5_brick), first = 2), format = "%Y.%m.%d")
  
  
  # HadCRUT5 land mask: 
  HadCRUT5_land_weights = t(apply(get_CMIP5_array(file = "HadCRUT.5.0.1.0.weights.nc", var = "weights", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadCRUT5_land_weights_brick = convert.to.pacificcentric(brick("HadCRUT.5.0.1.0.weights.nc") + 0)
  # image.plot(matrix(HadCRUT5_land_weights[1000,], 72, 36))
  # plot(HadCRUT5_land_weights_brick, 2000)    # weights over land. land = 1; water = 0
  # which(values(subset(HadCRUT5_land_weights_brick, 2000)) == 0.25)  # small islands receive w = 0.25
  
  
  # HadCRUT4 ensembles and uncertainties:
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT5/5d00_monthly_noninfilled/")
  HadCRUT5_unc = t(apply(get_CMIP5_array("_uncertainties/HadCRUT.5.0.1.0.uncorrelated.nc", var="tas_unc", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadCRUT5_unc_brick = convert.to.pacificcentric(brick("_uncertainties/HadCRUT.5.0.1.0.uncorrelated.nc") + 0)  # Uncorrelated Measurement and sampling uncertainties; simply to add onto ensembles: https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
  
  # HadCRUT4_SSTunc_uncorr = t(apply(get_CMIP5_array("HadCRUT.4.6.0.0.uncorrelated_supplementary.nc", var="standard_error", time.count = -1), 3, as.vector))
  # HadCRUT4_SSTunc_uncorr_brick = brick("HadCRUT.4.6.0.0.uncorrelated_supplementary.nc") + 0  # Uncorrelated Measurement and sampling uncertainties for SST DATA
  
  # Write function to get land+ocean cov. matrix for any month / or uncertainty realization:
  # setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT/5d00_monthly/")
  # for (i in 1:length(HadCRUT4_calendar)) {
  #  cur.year = format(HadCRUT4_calendar[i], "%Y")
  #  cur.month = format(HadCRUT4_calendar[i], "%m")
  #  print(paste(cur.year, cur.month))
  #  test.sample = sample.HadCRUT4.uncertainties(year = cur.year, month = cur.month, n.sim = 100, land_unc = HadCRUT4_landunc, ret.cov = F, cal = HadCRUT4_calendar)
  # }
  
  # HadCRUT5 - Ensemble:  
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT5/5d00_monthly_noninfilled/_HadCRUT.5.0.1.0.anomalies_1_to_200_netcdf/")
  HadCRUT5_ENS = list()
  for (i in 1:200) {
    print(i)
    HadCRUT5_ENS[[i]] = t(apply(get_CMIP5_array(file = paste("HadCRUT.5.0.1.0.anomalies.", i,".nc", sep=""), var="tas", time.count = -1), MARGIN = 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  }
  
  # image.plot(matrix(HadCRUT4_ENS[[99]][1500,], 72, 36), zlim = c(-15, 15))  # ensemble shows very little spread over the Ocean...
  # image.plot(matrix(HadCRUT4[1,], 72, 36))
  
  # Generate ensemble that includes uncertainties for *LAND* data:
  # image.plot(matrix(HadCRUT4_ENS_anom[[1]][801,], 72, 36))
}

  
# HadSST4:
# ------------------------------
{
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")
  HadSST4 = t(apply(get_CMIP5_array("HadSST.4.0.1.0_median.nc", var="tos", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadSST4[which(HadSST4 > 10^36)] = NA
  HadSST4_brick = convert.to.pacificcentric(brick("HadSST.4.0.1.0_median.nc", varname="tos") + 0)
  values(HadSST4_brick)[which(values(HadSST4_brick) > 10^36)] = NA
  HadSST4_nobs = t(apply(get_CMIP5_array(file = "HadSST.4.0.1.0_number_of_observations.nc", var = "numobs", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadSST4_nobs[which(HadSST4_nobs > 10^36)] = NA
  HadSST4_nobs_brick = convert.to.pacificcentric(brick("HadSST.4.0.1.0_number_of_observations.nc", varname="numobs") + 0)
  values(HadSST4_nobs_brick)[which(values(HadSST4_nobs_brick) > 10^36)] = NA
  HadSST4_calendar = as.Date(substring(names(HadSST4_nobs_brick), first = 2), format = "%Y.%m.%d")
  
  # Get HadSST4 uncertainties:
  HadSST4_sampling_unc = t(apply(X = get_CMIP5_array("HadSST.4.0.1.0_sampling_uncertainty.nc", var = "tos_unc", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  HadSST4_measurement_unc = t(apply(X = get_CMIP5_array("HadSST.4.0.1.0_uncorrelated_measurement_uncertainty.nc", var = "tos_unc", time.count = -1), 3, as.vector))[,transfer.CRU.to.CMIP6.grid]   # not needed as this appears to be encoded in the Ensemble.
  # HadSST4_sampling_unc = t(apply(X = get_CMIP5_array("HadSST.4.0.1.0_sampling_uncertainty.nc", var = "tos_unc", time.count = -1), 3, as.vector))
  # HadSST4_measurement_unc = t(apply(X = get_CMIP5_array("HadSST.4.0.1.0_uncorrelated_measurement_uncertainty.nc", var = "tos_unc", time.count = -1), 3, as.vector))
  
  # image.plot(matrix(HadSST4_sampling_unc[1,], 72, 36))
  # image.plot(matrix(HadSST4_sampling_unc[1,transfer.CRU.to.CMIP6.grid][transfer.CRU.to.CMIP6.grid], 72, 36))
  
  # HadSST4 - Ensemble:
  setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/HadSST.4.0.1.0_ensemble/")
  HadSST4_ENS = list()
  for (i in 1:200) {
    print(i)
    HadSST4_ENS[[i]] = t(apply(get_CMIP5_array(file = paste("HadSST.4.0.1.0_ensemble_member_", i,".nc", sep=""), var="tos", time.count = -1), MARGIN = 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
    HadSST4_ENS[[i]][which(HadSST4_ENS[[i]] > 10^36)] = NA
    # plot(brick(paste("HadSST.3.1.1.0.anomalies.", i,".nc", sep="")) + 0, 1000)
  }
}
  

# CLASSnmat:
# ------------------------------
{
  setwd("/net/h2o/climphys1/sippels/_DATASET/CLASSnmat/")
  CLASSnmat = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_anomaly", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  # image.plot(matrix(CLASSnmat[1,], 72, 36))
  CLASSnmat_brick = convert.to.pacificcentric(brick("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", varname="t2m_anomaly") + 0)
  CLASSnmat_calendar = as.Date(substring(names(CLASSnmat_brick), first = 2), format = "%Y.%m.%d")
  
  # Get CLASSnmat uncertainties:
  CLASSnmat_total_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_total_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CLASSnmat_uncorrelated_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_uncorrelated_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CLASSnmat_sampling_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_sampling_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CLASSnmat_correlated_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_correlated_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CLASSnmat_climatology_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_climatology_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
  CLASSnmat_anomaly_unc = t(apply(get_CMIP5_array("gridded/CLASSnmat_1.0.0.0_t2m_base_1961-1990.nc", var="t2m_anomaly_uncertainty", time.count = -1)[,36:1,], 3, as.vector))[,transfer.CRU.to.CMIP6.grid]
}



# Reconstruct CRUTEM5 ENSEMBLE from HadSST4 and HadCRUT5 ensemble: 
# ------------------------------
{
  # find all grid cells that have land+ocean values for each time step:
  CRUTEM5_ENS = list()
  for (i in 1:200) CRUTEM5_ENS[[i]] = matrix(data=NA, nrow = 2058, ncol = 72*36)
  # CRUTEM5_ENS_NH = matrix(data = NA, nrow = 2058, ncol = 200)
  # CRUTEM5_ENS_SH = matrix(data = NA, nrow = 2058, ncol = 200)
  # CRUTEM5_ENS_GL = matrix(data = NA, nrow = 2058, ncol = 200)
  
  CRUTEM5_ENS_anom_NH = matrix(data = NA, nrow = 2058, ncol = 200)
  CRUTEM5_ENS_anom_SH = matrix(data = NA, nrow = 2058, ncol = 200)
  CRUTEM5_ENS_anom_GL = matrix(data = NA, nrow = 2058, ncol = 200)
  
  NH.w = c(matrix(values(raster::area(subset(CRUTEM5_brick, 1))), 72, 36)[,36:1])[which(y.coord > 0)]
  SH.w = c(matrix(values(raster::area(subset(CRUTEM5_brick, 1))), 72, 36)[,36:1])[which(y.coord < 0)]
  NH.ix = which(y.coord > 0)
  SH.ix = which(y.coord < 0)
  
  for (date.ix in 1:2055) {
    print(date.ix)
    
    ## fur current time step, estimate land vs. water fraction:
    # ens.ix = 1
    # plot(HadCRUT4_ENS[[ens.ix]][date.ix, ], HadSST3_ENS[[ens.ix]][date.ix,], xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))
    # cor(HadCRUT4_ENS[[ens.ix]][date.ix, ], HadSST3_ENS[[ens.ix]][date.ix, ], use="complete.obs")
    # plot(HadCRUT4_ENS[[ens.ix]][date.ix, coast.ix], HadSST3_ENS[[ens.ix]][date.ix,coast.ix], xlim = c(-0.2, 0.2), ylim = c(-0.2, 0.2))
    # cor(HadCRUT4_ENS[[ens.ix]][date.ix, coast.ix], HadSST3_ENS[[ens.ix]][date.ix, coast.ix], use="complete.obs")
    # image.plot(matrix(HadCRUT4_ENS[[ens.ix]][date.ix, ], 72, 36), zlim = c(-8, 12))
    # image.plot(matrix(HadSST3_ENS[[ens.ix]][date.ix, ], 72, 36), zlim = c(-8, 12))
  
    ## Derive land-only CRUTEM5 ensemble:
    for (ens.ix in 1:200) {
      # print(ens.ix)
      
      cur.land.weights = HadCRUT5_land_weights[date.ix,] 
      
      # fill land grid cells:
      land.ix = which(cur.land.weights == 1)  # grid cells that contain ONLY land. ALTERNATIVE: HadCRUT5_land_weights[2000, which(is.na(HadSST4_ENS[[ens.ix]][2000,]) & !is.na(HadCRUT5_land_weights[2000,]))]
      CRUTEM5_ENS[[ens.ix]][date.ix,land.ix] = HadCRUT5_ENS[[ens.ix]][date.ix,land.ix]
      # image.plot(matrix(CRUTEM5_ENS[[i]][date.ix,], 72, 36))
      
      # fill combined grid cells:
      coast.ix = which(cur.land.weights < 1 & cur.land.weights > 0)
      CRUTEM5_ENS[[ens.ix]][date.ix,coast.ix] = (HadCRUT5_ENS[[ens.ix]][date.ix,coast.ix] - (1 - cur.land.weights[coast.ix]) * HadSST4_ENS[[ens.ix]][date.ix,coast.ix]) / cur.land.weights[coast.ix]
      # image.plot(matrix(CRUTEM5_ENS[[i]][date.ix,], 72, 36))  # these are only the coast pixels. 
      
      # CRUTEM5_ENS[[ens.ix]][date.ix,land.ix] = HadCRUT4_ENS[[ens.ix]][date.ix,land.ix]
      # CRUTEM5_ENS[[ens.ix]][date.ix,coast.ix] = (HadCRUT4_ENS[[ens.ix]][date.ix,coast.ix] - (1 - cur.land.fraction) * HadSST3_ENS[[ens.ix]][date.ix,coast.ix]) / cur.land.fraction
      # CRUTEM5_ENS_NH[date.ix,ens.ix] = stats::weighted.mean(x = CRUTEM5_ENS[[ens.ix]][date.ix,NH.ix],
      #                                                           w = NH.w, na.rm = T)
      
      CRUTEM5_ENS_anom_NH[date.ix,ens.ix] = stats::weighted.mean(x = CRUTEM5_ENS[[ens.ix]][date.ix,NH.ix] - CRUTEM5[date.ix,NH.ix],
                           w = NH.w, na.rm = T)
      if (any(!is.na(CRUTEM5_ENS[[ens.ix]][date.ix,SH.ix]))) {
        # CRUTEM5_ENS_SH[date.ix,ens.ix] = stats::weighted.mean(x = CRUTEM5_ENS[[ens.ix]][date.ix,SH.ix],
        #                                                           w = SH.w, na.rm = T) 
        CRUTEM5_ENS_anom_SH[date.ix,ens.ix] = stats::weighted.mean(x = CRUTEM5_ENS[[ens.ix]][date.ix,SH.ix] - CRUTEM5[date.ix,SH.ix],
                                                                 w = SH.w, na.rm = T) 
      }
      # CRUTEM5_ENS_GL[date.ix,ens.ix] = mean(c(CRUTEM5_ENS_NH[date.ix,ens.ix], CRUTEM5_ENS_SH[date.ix,ens.ix]), na.rm = T)
      CRUTEM5_ENS_anom_GL[date.ix,ens.ix] = mean(c(CRUTEM5_ENS_anom_NH[date.ix,ens.ix], CRUTEM5_ENS_anom_SH[date.ix,ens.ix]), na.rm = T)
    }
      
    #
    # plot(apply(sapply(X = CRUTEM4_ENS, FUN = function(x) x[date.ix,coast.ix]), 1, sd), apply(sapply(X = HadCRUT4_ENS, FUN = function(x) x[date.ix,coast.ix]), 1, sd))
    # lines(c(0,1), c(0,1), col = "red")
    # CRUTEM4 ensemble has much larger uncertainties than HadCRUT4-ensemble...
  }
  
  # GET ENSEMBLE ANOMALIES TO SUBTRACT:
  HadCRUT5_ENS_anom = list()
  for (i in 1:200) HadCRUT5_ENS_anom[[i]] = HadCRUT5_ENS[[i]] - HadCRUT5
  
  CRUTEM5_ENS_anom = list()
  for (i in 1:200) CRUTEM5_ENS_anom[[i]] = CRUTEM5_ENS[[i]][1:2055,] - CRUTEM5
  
  HadSST4_ENS_anom = list()
  for (i in 1:200) HadSST4_ENS_anom[[i]] = HadSST4_ENS[[i]][1:2059,] - HadSST4[1:2059,]
  
  ## check if ensemble mean corresponds to CRUTEM5 dataset:
  # test = sapply(X = CRUTEM5_ENS, FUN=function(x) x[2005,which(HadCRUT5_land_weights[2005,] == 1)])
  # plot(rowMeans(test), CRUTEM5[2005,which(HadCRUT5_land_weights[2005,] == 1)])
  # abline(0,1, col = "red")
  # test = sapply(X = CRUTEM5_ENS, FUN=function(x) x[2005, which(HadCRUT5_land_weights[2005,] < 1)])
  # plot(rowMeans(test), CRUTEM5[2005, which(HadCRUT5_land_weights[2005,] < 1)])
  # abline(0,1, col = "red")
  
  # image.plot(matrix(CRUTEM5_ENS_anom[[i]][1000,], 72, 36))
  
  ## plot uncertainty in CRUTEM5:
  
  
  date.ix = 2000
  test = sapply(X = CRUTEM5_ENS_anom, FUN=function(x) x[date.ix,])
  
  # SD across CRUTEM5 ensemble:
  library(RColorBrewer)
  col = (brewer.pal(n = 9, name = "Oranges"))
  image.plot(matrix(rowSds(test), 72, 36), zlim = c(0, 1), col = col)
  
  # Uncertainty: 
  image.plot(matrix(CRUTEM5_sampling_unc[date.ix,], 72, 36), zlim = c(0, 1), col = col)
  
  # plot ensemble:
  plot(apply(X = CRUTEM5_ENS_anom_NH, MARGIN = 1, FUN = function(x) diff(range(x))))
  
  plot(CRUTEM5_ENS_anom_NH[,1], type='l', ylim = c(-0.5, 0.5))
  for (i in 2:200) {
    lines(CRUTEM5_ENS_anom_NH[,i], type='l')
  }
  
  # from 1929 to 1930: massive reduction in ensemble spread ??
  
  plot(rowSds(CRUTEM5_ENS_anom_SH))
}


# save small junks of data to keep learning task "data-efficient":
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU")

for (mon in 1:12) {
  print(mon)
  time.ix = seq(mon, 2052, by = 12)
  
  
  # CRUTEM5
    CRUTEM5_ = CRUTEM5[time.ix,]
    CRUTEM5_nobs_ = CRUTEM5_nobs[time.ix,]
    CRUTEM5_ENS_anom_ = lapply(X = 1:200, FUN=function(i) CRUTEM5_ENS_anom[[i]][time.ix,])
    CRUTEM5_ENS_anom_GL_ = CRUTEM5_ENS_anom_GL[time.ix,]
    CRUTEM5_sampling_unc_ = CRUTEM5_sampling_unc[time.ix,]
    CRUTEM5_station_unc_ = CRUTEM5_station_unc[time.ix,]
    CRUTEM5_calendar_ = CRUTEM5_calendar[time.ix]
    land.weights = HadCRUT5_land_weights[time.ix,]
    
    save(list = c("CRUTEM5_", "CRUTEM5_nobs_", "CRUTEM5_ENS_anom_", "CRUTEM5_ENS_anom_GL_", 
                  "CRUTEM5_sampling_unc_", "CRUTEM5_station_unc_", 
                  "CRUTEM5_calendar_", "land.weights", "transfer.CRU.to.CMIP6.grid"),
         file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/", "CRUTEM5_mon", mon,".RData", sep=""))
    
    
    # HadSST4
    time.ix = seq(mon, 2052, by = 12)
    HadSST4_ = HadSST4[time.ix,]
    HadSST4_nobs_ = HadSST4_nobs[time.ix,]
    HadSST4_ENS_anom_ = lapply(X = 1:200, FUN=function(i) HadSST4_ENS_anom[[i]][time.ix,])
    HadSST4_sampling_unc_ = HadSST4_sampling_unc[time.ix,]
    HadSST4_measurement_unc_ = HadSST4_measurement_unc[time.ix,]
    # HadSST4_ENS_anom_GL_ = HadSST4_ENS_anom_GL[time.ix,]
    # HadSST4_sampling_unc_ = HadSST4_sampling_unc[time.ix,]
    # HadSST4_station_unc_ = HadSST4_station_unc[time.ix,]
    HadSST4_calendar_ = HadSST4_calendar[time.ix]
    land.weights = HadCRUT5_land_weights[time.ix,]
    
    save(list = c("HadSST4_", "HadSST4_nobs_", "HadSST4_ENS_anom_", 
                  "HadSST4_sampling_unc_", "HadSST4_measurement_unc_", 
                  "HadSST4_calendar_", "land.weights", "transfer.CRU.to.CMIP6.grid"),
         file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/", "HadSST4_mon", mon,".RData", sep=""))
    

    # HadCRUT5
    time.ix = seq(mon, 2052, by = 12)
    HadCRUT5_ = HadCRUT5[time.ix,]
    HadCRUT5_ENS_anom_ = lapply(X = 1:200, FUN=function(i) HadCRUT5_ENS_anom[[i]][time.ix,])
    # HadCRUT5_ENS_anom_GL_ = HadCRUT5_ENS_anom_GL[time.ix,]
    HadCRUT5_unc_ = HadCRUT5_unc[time.ix,]
    # HadCRUT5_station_unc_ = HadCRUT5_station_unc[time.ix,]
    HadCRUT5_calendar_ = HadCRUT5_calendar[time.ix]
    land.weights = HadCRUT5_land_weights[time.ix,]
    
    save(list = c("HadCRUT5_", "HadCRUT5_ENS_anom_", "HadCRUT5_calendar_", 
                  "land.weights", "transfer.CRU.to.CMIP6.grid"),
         file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/", "HadCRUT5_mon", mon,".RData", sep=""))
}


# save small junks of data to keep learning task "data-efficient":
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU")

for (mon in 1:12) {
  print(mon)
  time.ix = seq(mon, dim(CLASSnmat)[1], by = 12)
  
  # CLASSnmat
  CLASSnmat_ = CLASSnmat[time.ix,]
  # CLASSnmat_nobs_ = CLASSnmat_nobs[time.ix,]
  # CLASSnmat_ENS_anom_ = lapply(X = 1:200, FUN=function(i) CLASSnmat_ENS_anom[[i]][time.ix,])
  # CLASSnmat_ENS_anom_GL_ = CLASSnmat_ENS_anom_GL[time.ix,]
  CLASSnmat_total_unc_ = CLASSnmat_total_unc[time.ix,]
  CLASSnmat_calendar_ = CLASSnmat_calendar[time.ix]
  
  save(list = c("CLASSnmat_", 
                "CLASSnmat_total_unc_", 
                "CLASSnmat_calendar_", "transfer.CRU.to.CMIP6.grid"),
       file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/", "CLASSnmat_mon", mon,".RData", sep=""))
}


# ----------------------------------------------------------------------------------------------------------------------

