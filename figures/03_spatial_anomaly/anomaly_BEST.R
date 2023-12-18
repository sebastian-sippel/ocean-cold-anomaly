

## Create observational maps based on 1870-1900 reference period:
# Sebastian Sippel
# 30.05.2022

library(raster)
library(ncdf4)
library(fields)
library(matrixStats)
library(bigmemory)
library(SpatialEpi)
library(RColorBrewer)
library(rworldmap)


# source("/net/h2o/climphys1/sippels/_projects/low_freq_anchor_v2/code/_convenience/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/low_freq_anchor_v2/code/_convenience/project_raster.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_plot_projected_worldmap_v2.R")
# source("/net/h2o/climphys1/sippels/_projects/low_freq_anchor_v2/code/_convenience/_plot_projected_worldmap_v2.R")


raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw=c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]))

col = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(101))
RdBu = brewer.pal(n = 11, name = "RdBu")

## Define functions:
get.change <- function(cur.raster, cur.nobs, ref.period = 1870:1895, test.period = 1900:1920, frac.data = 0.5) {
  
  year = as.numeric(substring(text = names(cur.raster), 2, 5))
  month = as.numeric(substring(text = names(cur.raster), 7, 8))
  
  ref.period.ix = which(year %in% ref.period)
  test.period.ix = which(year %in% test.period)
  
  ref.raster = raster(cur.raster); values(ref.raster) = NA
  
  # plot(sum(!is.na(subset(cur.nobs, ref.period.ix)))) # number of months with data!
  cur.mask = sum(!is.na(subset(cur.nobs, ref.period.ix))) > (length(ref.period.ix) * frac.data)
  
  ## CONTINUE here with simple averaging with >0.6 data points avail
  # plot(mean(subset(cur.raster, ref.period.ix), na.rm=T))
  cur.raster.ref = mask(x = mean(subset(cur.raster, ref.period.ix), na.rm=T), mask = cur.mask, maskvalue = 0)
  cur.raster.test = mask(x = mean(subset(cur.raster, test.period.ix), na.rm=T), mask = cur.mask, maskvalue = 0)
  
  # plot(cur.raster.test - cur.raster.ref, zlim = c(-1.5, 1.5), col = col)
  raster.diff = cur.raster.test - cur.raster.ref
  return(raster.diff)
}



## ------------------------------------------------------------------------
### 00. Berkeley Earth Comparison Land vs. Ocean: 
## ------------------------------------------------------------------------

# read 1.0 degree land-ocean mask from Koeppen: 
land_fraction = raster("/net/h2o/climphys1/sippels/_DATA/grid/0d50_static/Koeppen/land_sea_mask_0d50.nc")

# read Berkeley Earth combined product with SSTatSeaIce:

setwd("/net/h2o/climphys1/sippels/_DATASET/BEST/_orig/v202302/average_temperature_SSTatSeaIce/")
BEST_anom = brick("Complete_TAVG_LatLong1.nc") + 0
plot(BEST_anom, 5)


# 









