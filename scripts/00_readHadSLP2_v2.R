
## ------------------------------------
## 00_readHadSLP2_v2
## ------------------------------------

# Sebastian Sippel
# 16.08.2022
# run on n2o-server with R version R-4.1.2

## cmip_split session on n2o, 06.04.2021:
# screen -S _process_HadSLP2
# module load R/4.0.3-openblas 
# R (-> not R-3.6.1)


## read uninterpolated product.
library(ncdf4)
library(raster)
library(fields)
# install.packages("remotes")
# remotes::install_github("SEEG-Oxford/seegSDM")
library("seegSDM")

source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")

# https://www.metoffice.gov.uk/hadobs/hadslp2/data/Read_instructions_HADSLP2_5dg
setwd("/net/h2o/climphys1/sippels/_DATASET/HadSLP2/metoffice/")


# Read Full Dataset:
HadSLP2_interpol_r = brick("/net/h2o/climphys1/sippels/_DATASET/HadSLP2/ncdf_from_NOAA/slp.mnmean.real.nc") + 0
HadSLP2_interpol_test = get_CMIP5_array(file = "/net/h2o/climphys1/sippels/_DATASET/HadSLP2/ncdf_from_NOAA/slp.mnmean.real.nc", var = "slp")
str(HadSLP2_interpol_test)
# image.plot(HadSLP2_interpol_test[,,1860])

HadSLP2_calendar = as.Date(substring(names(HadSLP2_interpol_r), first = 2), format = "%Y.%m.%d")
HadSLP2_calendar = HadSLP2_calendar[1:1860]

## Read HadSLP2_v2 uninterpolated:
library(foreach)
library(doParallel)
registerDoParallel(cores=40)

ret.list = foreach(i=1:1860) %dopar% {
  print(i)
  
  # read HadSLP2 uninterpolated:
  test = t(as.matrix(read.table(file = "hadslp2.0_acts.asc", skip = 37*(i-1) + (i - 1) + 1, nrows = 37)))[,37:1]
  test[which(test == -99990)] = NA
  HadSLP2_uninterpol = c(rbind(test[37:72,], test[1:36,]) / 100) # convert to hPa and to 0-360°
  
  # read HadSLP2 interpolated product:
  test = t(as.matrix(read.table(file = "hadslp2.asc", skip = 37*(i-1) + (i - 1) + 1, nrows = 37)))[,37:1]
  # test[which(test == -99990)] = NA
  HadSLP2_interpol = c(rbind(test[37:72,], test[1:36,]) / 100) # convert to hPa
  
  # read observational errors:
  test = t(as.matrix(read.table(file = "hadslp2_obs-error.asc", skip = 37*(i-1) + (i - 1) + 1, nrows = 37)))[,37:1]
  test[which(test == -99990)] = NA
  HadSLP2_obs_error = c(rbind(test[37:72,], test[1:36,]) / 100) # convert to hPa and to 0-360°
  
  # read number of observations:
  test = t(as.matrix(read.table(file = "hadslp2_nobs.asc", skip = 37*(i-1) + (i - 1) + 1, nrows = 37)))[,37:1]
  test[which(test == -99)] = NA
  HadSLP2_nobs = c(rbind(test[37:72,], test[1:36,])) 
  
  
  ## find HadSLP2_obs_error estimates closest to respective grid cell:
  raster.template = raster(xmn = 0, xmx=360, ymn = -90, ymx=90, nrows = 37, ncols = 72)
  grid.ix = which(!is.na(c(matrix(HadSLP2_uninterpol, 72, 37)[,37:1])))
  
  miss_err_grid = grid.ix[which(is.na(c(matrix(HadSLP2_obs_error, 72, 37)[,37:1])[grid.ix]))]  # missing error estimates on raster grid cells.
  err_raster <- err_raster_f <- raster.template
  values(err_raster) = c(matrix(HadSLP2_obs_error, 72, 37)[,37:1])
  test = data.frame(x = coordinates(err_raster)[,1], y = coordinates(err_raster)[,2])
  test0 = nearestLand(points = test, raster = err_raster, max_distance = 3000 * 10^3)
  values(err_raster_f) <- raster::extract(x = err_raster, y = test0, method = "simple")
  
  HadSLP2_obs_error_f = c(matrix(values(err_raster_f), 72, 37)[,37:1])
  miss_err_grid_f = grid.ix[which(is.na(c(matrix(HadSLP2_obs_error_f, 72, 37)[,37:1])[grid.ix]))]  # missing error estimates on raster grid cells.
  # plot(err_raster); # plot(err_raster_f)
  
  ## set nearest grid cell error and multiply error by x2:
  HadSLP2_obs_error_f2x = HadSLP2_obs_error_f * 2
  HadSLP2_obs_error_f2x[which(!is.na(HadSLP2_obs_error))] = HadSLP2_obs_error[which(!is.na(HadSLP2_obs_error))]
  
  ret.list = (list(HadSLP2_uninterpol = HadSLP2_uninterpol, HadSLP2_interpol = HadSLP2_interpol, HadSLP2_obs_error = HadSLP2_obs_error, 
                HadSLP2_nobs = HadSLP2_nobs, HadSLP2_obs_error_f = HadSLP2_obs_error_f, HadSLP2_obs_error_f2x = HadSLP2_obs_error_f2x,
                grid.ix = grid.ix, miss_err_grid = miss_err_grid, miss_err_grid_f = miss_err_grid_f))
  return(ret.list)
}


## put from ret.list into matrix format:
HadSLP2_uninterpol = matrix(data = NA, nrow = 1860, ncol = 72*37)
HadSLP2_interpol = matrix(data = NA, nrow = 1860, ncol = 72*37)
HadSLP2_nobs = matrix(data = NA, nrow = 1860, ncol = 72*37)
HadSLP2_obs_error = matrix(data = NA, nrow = 1860, ncol = 72*37)
HadSLP2_obs_error_f = matrix(data = NA, nrow = 1860, ncol = 72*37)
HadSLP2_obs_error_f2x = matrix(data = NA, nrow = 1860, ncol = 72*37)

for (i in 1:1860) {
  print(i)
  HadSLP2_uninterpol[i,] = ret.list[[i]]$HadSLP2_uninterpol
  HadSLP2_interpol[i,] = ret.list[[i]]$HadSLP2_interpol
  HadSLP2_nobs[i,] = ret.list[[i]]$HadSLP2_nobs
  HadSLP2_obs_error[i,] = ret.list[[i]]$HadSLP2_obs_error
  HadSLP2_obs_error_f[i,] = ret.list[[i]]$HadSLP2_obs_error_f
  HadSLP2_obs_error_f2x[i,] = ret.list[[i]]$HadSLP2_obs_error_f2x
}


setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_HadSLP2/")
ref.year.ix = match(x = 1961:1990, table = 1850:2004)

for (mon in 1:12) {
  print(mon)
  time.ix = seq(mon, length(HadSLP2_calendar), by = 12)
  
    HadSLP2_uninterpol_mean = apply(X = HadSLP2_uninterpol[time.ix[ref.year.ix],], MARGIN=2, FUN=mean, na.rm=T)
    HadSLP2_uninterpol_mean[which(apply(X = HadSLP2_uninterpol[time.ix[ref.year.ix],], MARGIN = 2, FUN = function(x) length(which(is.na(x))) >= 15))] = NA
    HadSLP2_uninterpol_ = HadSLP2_uninterpol[time.ix,] - rep.row(x = HadSLP2_uninterpol_mean, n = length(time.ix))
    
    HadSLP2_interpol_mean = apply(X = HadSLP2_interpol[time.ix[ref.year.ix],], MARGIN=2, FUN=mean, na.rm=T)
    HadSLP2_interpol_mean[which(apply(X = HadSLP2_interpol[time.ix[ref.year.ix],], MARGIN = 2, FUN = function(x) length(which(is.na(x))) >= 15))] = NA
    HadSLP2_interpol_ = HadSLP2_interpol[time.ix,] - rep.row(x = HadSLP2_interpol_mean, n = length(time.ix))
    
    HadSLP2_obs_error_ = HadSLP2_obs_error[time.ix,]
    HadSLP2_obs_error_f_ = HadSLP2_obs_error_f[time.ix,]
    HadSLP2_obs_error_f2x_ = HadSLP2_obs_error_f2x[time.ix,]
    
    HadSLP2_nobs_ = HadSLP2_nobs[time.ix,]
    HadSLP2_calendar_ = HadSLP2_calendar[time.ix]
  
  
  save(list = c("HadSLP2_uninterpol_", "HadSLP2_interpol_", "HadSLP2_obs_error_", "HadSLP2_obs_error_f_", "HadSLP2_obs_error_f2x_", "HadSLP2_nobs_", "HadSLP2_calendar_",
                "HadSLP2_uninterpol_mean", "HadSLP2_interpol_mean"),
       file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_HadSLP2/", "HadSLP2_mon", mon,".RData", sep=""))
}









## Check maps:
{
i = i
image.plot(matrix(HadSLP2_interpol, 72, 37))
test = raster(xmn = 0, xmx=360, ymn = -90, ymx=90, nrows = 37, ncols = 72)
values(test) = c(matrix(HadSLP2_interpol, 72, 37)[,37:1])
plot(test)
lines(coastsCoarse)

image.plot(matrix(values(subset(HadSLP2_interpol_r, i)), 72, 37)[,37:1])
image.plot(matrix(HadSLP2_uninterpol[i,], 72, 37))
image.plot(matrix(HadSLP2_obs_error[i,], 72, 37))
image.plot(matrix(HadSLP2_nobs[i,], 72, 37))

grid.ix = (which(!is.na(HadSLP2_uninterpol[i,])))
length(grid.ix)
HadSLP2_obs_error[i,]


plot(c(matrix(HadSLP2_uninterpol[i,], 72, 37)), c(matrix(values(subset(HadSLP2_interpol_r, i)), 72, 37)[,37:1]))
cor(c(matrix(HadSLP2_uninterpol[i,], 72, 37)), c(matrix(values(subset(HadSLP2_interpol_r, i)), 72, 37)[,37:1]), use = "complete.obs")
plot(c(matrix(HadSLP2_interpol[i,], 72, 37)), c(matrix(values(subset(HadSLP2_interpol_r, i)), 72, 37)[,37:1]))
}
