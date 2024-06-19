


## Create observational maps based on 1870-1900 reference period:
# Sebastian Sippel
# 30.05.2022

library(raster)
library(ncdf4)
# library(hydroGOF)
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
get.change <- function(cur.raster, cur.nobs, ref.period = 1871:1890, test.period = 1900:1920, frac.data = 0.5) {
  
  year = as.numeric(substring(text = names(cur.raster), 2, 5))
  month = as.numeric(substring(text = names(cur.raster), 7, 8))
  
  ref.period.ix = which(year %in% ref.period)
  test.period.ix = which(year %in% test.period)
  
  ref.raster = raster(cur.raster); values(ref.raster) = NA
  
  # plot(sum(!is.na(subset(cur.nobs, ref.period.ix)))) # number of months with data!
  cur.mask = sum(!is.na(subset(cur.nobs, ref.period.ix))) > (length(ref.period.ix) * frac.data)
  
  ## CONTINUE here with simple averaging with >0.6 data points avail
  # plot(mean(subset(cur.raster, ref.period.ix), na.rm=T))
  cur.raster.ref = raster::mask(x = mean(subset(cur.raster, ref.period.ix), na.rm=T), mask = cur.mask, maskvalue = 0)
  cur.raster.test = raster::mask(x = mean(subset(cur.raster, test.period.ix), na.rm=T), mask = cur.mask, maskvalue = 0)
  
  # plot(cur.raster.test - cur.raster.ref, zlim = c(-1.5, 1.5), col = col)
  raster.diff = cur.raster.test - cur.raster.ref
  return(raster.diff)
}






## ------------------------------------------------------------------------
### 00. Berkeley Earth Comparison Land vs. Ocean:  -> FIX LATER
## ------------------------------------------------------------------------
{
### Look at Berkeley Earth dataset and 1900-1920 residuals relative to 1880-1900 and 1920-1940:
setwd("/net/h2o/climphys1/sippels/_DATASET/BEST/_orig/v202302/")
system("cdo -yearmean Land_and_Ocean_LatLong1.nc Land_and_Ocean_LatLong1_ann.nc")
BEST_anom = brick("Land_and_Ocean_LatLong1_ann.nc") + 0


## where is no data in HadSST3 ?

# HadSST3 dataset:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")
HadSST4 = brick("HadSST.4.0.1.0_median.nc") + 0
HadSST4_nobs = brick("HadSST.4.0.1.0_number_of_observations.nc") + 0
# get land mask of HadSST4:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/CRUTEM4v/5d00_monthly/")
cur.mask = brick("/net/h2o/climphys1/sippels/_DATA/grid/5d00_static/cmip5_masks/sftlf_g025.nc") + 0
#cur.mask = any(!is.na(HadSST3))
# plot(land.mask, 1)


#### Get BerkeleyEarth comparison for 1901-1920:
#
  cold_anomaly = calc(subset(x = BEST_anom, c(52:71)), fun=mean, na.rm=T)
  ref_period = calc(subset(x = BEST_anom, c(22:41)), fun = mean, na.rm=T)
  
  # covered for at least 50% of the time during 1901-1920:
  test.obs = (subset(x = HadSST4_nobs, (51*12 + 1):(71*12)) > 1)
  noncov = (calc(test.obs, fun = sum, na.rm = T)) < (nlayers(test.obs) * 0.5)
  
  setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/02_spatial_anomaly/")

  points.to.add = SpatialPoints(data.frame(coordinates(noncov)[which(values(noncov) == 1 & values(cur.mask) < 0.5),]))
  plot_projected_worldmap_eucentric(file.name = "_BEST_1901-20.png", beta = cold_anomaly - ref_period, 
                          disagg = 1, col = col, zlim = c(-3, 3), legend.text = "Temperature Difference [°C]", main = "Mean Temperature Anomaly 1901-20 w.r.t. 1871-1890 (BerkeleyEarth)",
                          points.to.add = points.to.add)

  #  png("")
  # par(mfrow=c(1,1), oma = c(0,0,0,0), mar = c(4, 4, 3, 4))
  # plot(cold_anomaly - ref_period, zlim = c(-3, 3), col = col, 
  #     main = "Mean Temperature Anomaly 1901-20 w.r.t. 1871-1900 (BerkeleyEarth)")
  # lines(coastsCoarse, col = "grey")
  # points(SpatialPoints(data.frame(coordinates(noncov)[which(values(noncov) == 1 & values(cur.mask) < 0.5),])), pch = 16, cex = 0.4)
}








## ------------------------------------------------------------------------
### 01. HadSST4 vs. CRUTEM5:
## ------------------------------------------------------------------------

### Look at CRUTEM5 and HadSST4 datasets:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/CRUTEM5/5d00_monthly/")
CRUTEM5 = brick("CRUTEM.5.0.1.0.anomalies.nc") + 0
CRUTEM5_nobs = brick("CRUTEM.5.0.1.0.station_counts.nc") + 0

# HadSST4 dataset:
setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")
HadSST4 = brick("HadSST.4.0.1.0_median.nc") + 0
HadSST4_nobs = brick("HadSST.4.0.1.0_number_of_observations.nc") + 0

## should be at least observations 50% of the time:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)

# plot(HadSST4_diff, zlim = c(-1.5, 1.5), col = col)
# plot(CRUTEM5_diff, zlim = c(-1.5, 1.5), col = col)

# get overlapping grid cells:
ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
ix1 = which(values(!is.na(CRUTEM5_diff)))
ix2 = which(values(!is.na(HadSST4_diff)))

combi.raster = raster(HadSST4_diff)
values(combi.raster) = NA
values(combi.raster)[ix2] = values(HadSST4_diff)[ix2]
values(combi.raster)[ix1] = values(CRUTEM5_diff)[ix1]
# values(combi.raster)[ix] = NA

plot(CRUTEM5_diff, zlim = c(-2, 2), col = col)
lines(coastsCoarse, col = "grey")

plot(HadSST4_diff, zlim = c(-2, 2), col = col)
lines(coastsCoarse, col = "grey")

setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/02_spatial_anomaly/")
plot_projected_worldmap_eucentric_DIFF(file.name = "_HadSST4_vs_CRUTEM5_1901-20_frac50.png", beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = T, 
                                                   legend.text = "Temperature difference [°C] 1901-1920 vs. 1871-1890", main = "", nlab = 5, legend = T, asp = 1, 
                                                   raster.to.add = CRUTEM5_diff)

# with at least 20% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)

plot_projected_worldmap_eucentric_DIFF(file.name = "_HadSST4_vs_CRUTEM5_1901-20_frac20.png", beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = T, 
                                       legend.text = "Temperature difference [°C] 1901-1920 vs. 1871-1890", main = "", nlab = 5, legend = T, asp = 1, 
                                       raster.to.add = CRUTEM5_diff)




##### --------------------------------------------------------------------------------------------- ######
# 02. boxplots with direct comparison of grid cells:
##### --------------------------------------------------------------------------------------------- ######

par(mar=c(5,5,1,1))
plot(values(CRUTEM5_diff)[ix])
lines(values(HadSST4_diff)[ix], col="red")

plot(values(CRUTEM5_diff)[ix] - values(HadSST4_diff)[ix])

## make collection of boxplots in different regions:

# top row: Boxplots of relative differences to 1871-1900 period over land and ocean
# bottom row: Relative differences between land and ocean.

# Regions:
# Global, NH-extra, Tropics., SH-extra
# Pacific, Atlantic, Indian Ocean.

# Domain definition, Division of global oceans: 
# https://www.researchgate.net/publication/297661864_Fast_multidimensional_ensemble_empirical_mode_decomposition_for_the_analysis_of_big_spatio-temporal_datasets

setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/02_spatial_anomaly/")

# with at least 50% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)

pdf("_temperature_differences_frac50.pdf", width = 8, height= 4)
{
  par(mfrow = c(1, 1), mar = c(4,5,1,1))
  plot(c(1,1), type="n", xlim = c(0.5, 7.5), ylim = c(-1.1, 1.1),
       xlab = "", xaxt="n", ylab = "Temperature difference in coastal grid cells [°C], \n 1901-1920 w.r.t. 1871-1890", las = 1, xaxs="i")
  axis(side = 2, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(x = c(0, 10), y = c(0, 0), col = "black", lty = 2)
  
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 0.8, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 1.2, varwidth = T, notch = F, boxwex = 0.3, frame=F, xaxt="n", axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 1.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 1, col = "grey40")
    
    lines(x = c(1.5, 1.5), y = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 1.8, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 2.2, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 2.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 2, col = "grey40")
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 2.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 3.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 3.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 3, col = "grey40")
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 3.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 4.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 4.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 4, col = "grey40")
    
    lines(x = c(4.5, 4.5), y = c(-2, 2), col = "grey40")
    
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 4.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 5.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 5.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 5, col = "grey40")
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 5.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 6.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 6.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 6, col = "grey40")
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 6.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 7.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 7.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 7, col = "grey40")
  }
  mtext(text = "Globe", side = 1, line = 1, at = 1)
  mtext(text = "NH Extra-\nTropics", side = 1, line = 1, at = 2)
  mtext(text = "Tropics", side = 1, line = 1, at = 3)
  mtext(text = "SH Extra-\nTropics", side = 1, line = 1, at = 4)
  mtext(text = "Pacific\n Ocean", side = 1, line = 1, at = 5)
  mtext(text = "Atlantic\n Ocean", side = 1, line = 1, at = 6)
  mtext(text = "Indian\n Ocean", side = 1, line = 1, at = 7)
  
  legend("top", c("Land", "Ocean"), col = c("darkorange", "blue"), lty = 1, lwd = 8, inset = 0.1, bg = "white")
}
dev.off()


## try whether violin plot looks nicer:
install.packages("vioplot")
library(vioplot)
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")


pdf("_temperature_differences_frac50_violin.pdf", height = 9, width = 5)
{
  par(mfrow = c(1, 1), mar = c(4,1,2,1))
  plot(c(1,1), type="n", ylim = c(1, 8), xlim = c(-1.1, 1.1),
       ylab = "", yaxt="n", xlab = "Temperature difference [°C], \n 1901-1920 w.r.t. 1871-1890", las = 1, xaxs="i", main = "Coastal Grid Cells")
  axis(side = 1, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(y = c(0, 10), x = c(0, 0), col = "black", lty = 2)
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 7.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 7.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 7.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 7.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 7.9, labels = c("Global"), adj = 0)
    text(x = 0.9, y = 7.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    lines(y = c(7.1, 7.1), x = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 6.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 6.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 6.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 6.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 6.9, labels = c("Northern Hemisphere \n Extra-Tropics"), adj = 0)
    text(x = 0.9, y = 6.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 5.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 5.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 5.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 5.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 5.9, labels = c("Tropics"), adj = 0)
    text(x = 0.9, y = 5.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    at = 4.5; label = "Southern Hemisphere \n Extra-Tropics"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    lines(y = c(4.1, 4.1), x = c(-2, 2), col = "grey40")
  
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    at = 3.5; label = "Coastal Pacific"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 2.5; label = "Coastal Atlantic"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 1.5; label = "Coastal Indian Ocean"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
  }
  legend("bottomright", c("Land", "Ocean"), col = c("orange", "lightblue"), lty = 1, lwd = 8, inset = 0.02, bg = "white")
}
dev.off()



# with at least 20% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)



pdf("_temperature_differences_frac20_violin.pdf", height = 9, width = 5)
{
  par(mfrow = c(1, 1), mar = c(4,1,2,1))
  plot(c(1,1), type="n", ylim = c(1, 8), xlim = c(-1.1, 1.1),
       ylab = "", yaxt="n", xlab = "Temperature difference [°C], \n 1901-1920 w.r.t. 1871-1890", las = 1, xaxs="i", main = "Coastal Grid Cells")
  axis(side = 1, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(y = c(0, 10), x = c(0, 0), col = "black", lty = 2)
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 7.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 7.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 7.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 7.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 7.9, labels = c("Global"), adj = 0)
    text(x = 0.9, y = 7.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    lines(y = c(7.1, 7.1), x = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 6.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 6.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 6.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 6.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 6.9, labels = c("Northern Hemisphere \n Extra-Tropics"), adj = 0)
    text(x = 0.9, y = 6.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = 5.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 5.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = 5.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 5.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = 5.9, labels = c("Tropics"), adj = 0)
    text(x = 0.9, y = 5.9, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    at = 4.5; label = "Southern Hemisphere \n Extra-Tropics"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    lines(y = c(4.1, 4.1), x = c(-2, 2), col = "grey40")
    
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    at = 3.5; label = "Coastal Pacific"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 2.5; label = "Coastal Atlantic"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 1.5; label = "Coastal Indian Ocean"
    vioplot(x = values(CRUTEM5_diff)[ix], col = "orange", side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = "lightblue", side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -0.9, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.9, y = at + 0.4, labels = paste("n = ", length(ix), sep=""), col = "grey40")
  }
  legend("bottomright", c("Land", "Ocean"), col = c("orange", "lightblue"), lty = 1, lwd = 8, inset = 0.02, bg = "white")
}
dev.off()




pdf("_temperature_differences_frac20.pdf", width = 8, height= 4)
{
  par(mfrow = c(1, 1), mar = c(4,5,1,1))
  plot(c(1,1), type="n", xlim = c(0.5, 7.5), ylim = c(-1.1, 1.1),
       xlab = "", xaxt="n", ylab = "Temperature difference in coastal grid cells [°C], \n 1901-1920 w.r.t. 1871-1900", las = 1, xaxs="i")
  axis(side = 2, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(x = c(0, 10), y = c(0, 0), col = "black", lty = 2)
  
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 0.8, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 1.2, varwidth = T, notch = F, boxwex = 0.3, frame=F, xaxt="n", axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 1.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 1, col = "grey40")
    
    lines(x = c(1.5, 1.5), y = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 1.8, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 2.2, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 2.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 2, col = "grey40")
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 2.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 3.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 3.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 3, col = "grey40")
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 3.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 4.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 4.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 4, col = "grey40")
    
    lines(x = c(4.5, 4.5), y = c(-2, 2), col = "grey40")
    
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 4.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 5.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 5.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 5, col = "grey40")
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 5.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 6.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 6.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 6, col = "grey40")
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    boxplot(values(CRUTEM5_diff)[ix], col = "darkorange", add=T, at = 6.8, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    boxplot(values(HadSST4_diff)[ix], col = "blue", add=T, at = 7.2, varwidth = T, notch = F, boxwex = 0.3, frame = F, axes = F)
    # boxplot(values(HadSST4_diff)[ix] - values(CRUTEM5_diff)[ix], col = "lightgrey", add=T, at = 7.25, varwidth = T, notch = F, boxwex = 0.3, frame=F, axes = F)
    mtext(text = paste("n = ", length(ix), sep=""), side = 3, line = -1, at = 7, col = "grey40")
  }
  mtext(text = "Globe", side = 1, line = 1, at = 1)
  mtext(text = "NH Extra-\nTropics", side = 1, line = 1, at = 2)
  mtext(text = "Tropics", side = 1, line = 1, at = 3)
  mtext(text = "SH Extra-\nTropics", side = 1, line = 1, at = 4)
  mtext(text = "Pacific\n Ocean", side = 1, line = 1, at = 5)
  mtext(text = "Atlantic\n Ocean", side = 1, line = 1, at = 6)
  mtext(text = "Indian\n Ocean", side = 1, line = 1, at = 7)
  
  legend("top", c("Land", "Ocean"), col = c("darkorange", "blue"), lty = 1, lwd = 8, inset = 0.1, bg = "white")
}
dev.off()







