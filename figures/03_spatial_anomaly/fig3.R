
# ------------------------------------------------------------------------
## Fig. 3: 
# ------------------------------------------------------------------------

# Create observational maps based on 1871-1890 reference period
# Create violin plots for coastal grid cells

# Sebastian Sippel
# 04.07.2024

library(raster)
library(ncdf4)
# library(hydroGOF)
library(fields)
library(matrixStats)
library(bigmemory)
library(SpatialEpi)
library(RColorBrewer)
library(rworldmap)


setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/")

source("code/_convenience/project_raster.R")
source("code/_convenience/_plot_projected_worldmap_v2.R")
source("code/_convenience/frenchcolormap.R")

# Define raster and colours:
raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw=c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]))
col = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(101))
RdBu = brewer.pal(n = 11, name = "RdBu")




## Define functions:
## ------------------------------------------------------------------------

# get.change: calculate the mean change between two periods.
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
### 00. Berkeley Earth Comparison Land vs. Ocean: 
## ------------------------------------------------------------------------
{
### Look at Berkeley Earth dataset and 1900-1920 residuals relative to 1880-1900 and 1920-1940:
# setwd("/net/h2o/climphys1/sippels/_DATASET/BEST/_orig/v202302/")
# system("cdo -yearmean Land_and_Ocean_LatLong1.nc Land_and_Ocean_LatLong1_ann.nc")
BEST_anom = brick("data/00_DATASET/obs/BEST/Land_and_Ocean_LatLong1_ann.nc", varname = "temperature") + 0


#### Get BerkeleyEarth comparison for 1901-1920:

# set working directory:
  cold_anomaly = calc(subset(x = BEST_anom, c(52:71)), fun=mean, na.rm=T)
  ref_period = calc(subset(x = BEST_anom, c(22:41)), fun = mean, na.rm=T)
  
  pdf(file = "figures/03_spatial_anomaly/fig3b.pdf", width = 3.5, height = 2.333, pointsize = 5)
  par(mar=c(1,1,1,1), oma=c(0,0,0,0))
  plot_projected_worldmap_eucentric(file.name = NULL, beta = cold_anomaly - ref_period, 
                          disagg = 1, col = col, zlim = c(-2, 2), legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "Berkeley Earth",
                          points.to.add = NULL)
  dev.off()
  
  plot_projected_worldmap_eucentric(file.name = "figures/03_spatial_anomaly/fig3b.png", beta = cold_anomaly - ref_period, 
                                    disagg = 1, col = col, zlim = c(-2, 2), legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "Berkeley Earth",
                                    points.to.add = NULL)
}



## ------------------------------------------------------------------------
### 01. HadSST4 vs. CRUTEM5:
## ------------------------------------------------------------------------

### Look at CRUTEM5 and HadSST4 datasets:
CRUTEM5 = brick("data/00_DATASET/obs/CRUTEM5/CRUTEM.5.0.1.0.anomalies.nc") + 0
CRUTEM5_nobs = brick("data/00_DATASET/obs/CRUTEM5/CRUTEM.5.0.1.0.station_counts.nc") + 0

# HadSST4 dataset:
#setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")
HadSST4 = brick("data/00_DATASET/obs/HadSST.4.0.1.0/HadSST.4.0.1.0_median.nc") + 0
HadSST4_nobs = brick("data/00_DATASET/obs/HadSST.4.0.1.0/HadSST.4.0.1.0_number_of_observations.nc") + 0

## should be at least observations 50% of the time:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)

# get overlapping grid cells:
ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
ix1 = which(values(!is.na(CRUTEM5_diff)))
ix2 = which(values(!is.na(HadSST4_diff)))

combi.raster = raster(HadSST4_diff)
values(combi.raster) = NA
values(combi.raster)[ix2] = values(HadSST4_diff)[ix2]
values(combi.raster)[ix1] = values(CRUTEM5_diff)[ix1]

par(mfrow=c(1,1))
plot(CRUTEM5_diff, zlim = c(-2, 2), col = col)
lines(coastsCoarse, col = "grey")

plot(HadSST4_diff, zlim = c(-2, 2), col = col)
lines(coastsCoarse, col = "grey")


# set working directory:
# setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/03_spatial_anomaly/")


pdf(file = "figures/03_spatial_anomaly/fig3a_frac50.pdf", width = 3.5, height = 2.3333, pointsize = 5)
par(mar=c(1,0.5,1,0.5), oma=c(0,0.5,0,0.5))
plot_projected_worldmap_eucentric_DIFF(file.name = NULL, beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = T, 
                                                   legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "CRUTEM5 and HadSST4", nlab = 5, legend = T, asp = 1, 
                                                   raster.to.add = CRUTEM5_diff)
dev.off()


plot_projected_worldmap_eucentric_DIFF(file.name = "figures/03_spatial_anomaly/fig3a_frac50.png", beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = T, 
                                       legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "CRUTEM5 and HadSST4", nlab = 5, legend = T, asp = 1, 
                                       raster.to.add = CRUTEM5_diff)

# with at least 20% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)

pdf(file = "figures/03_spatial_anomaly/fig3a_frac20.pdf", width = 3.5, height = 2.3333, pointsize = 5)
par(mar=c(1,1,1,1), oma=c(0,0,0,0))
plot_projected_worldmap_eucentric_DIFF(file.name = NULL, beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = F, 
                                       legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "CRUTEM5 and HadSST4", nlab = 5, legend = T, asp = 1, 
                                       raster.to.add = CRUTEM5_diff)
dev.off()

plot_projected_worldmap_eucentric_DIFF(file.name = "figures/03_spatial_anomaly/fig3a_frac20.png", beta = HadSST4_diff, disagg = 1, col = rev(colorRampPalette(RdBu)(101)), zlim = c(-2, 2), useRaster = F, 
                                       legend.text = "Temperature difference [°C], [1901-20] vs. [1871-90]", main = "CRUTEM5 and HadSST4", nlab = 5, legend = T, asp = 1, 
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

# setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/03_spatial_anomaly/")

# with at least 50% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.5)

## try whether violin plot looks nicer:
library(vioplot)

pdf("figures/03_spatial_anomaly/fig3c_frac50.pdf", height = 4.5, width = 2, pointsize = 5)
{
  par(mfrow = c(1, 1), mar = c(4,1,2,1))
  plot(c(1,1), type="n", ylim = c(1, 8), xlim = c(-1.1, 1.1),
       ylab = "", yaxt="n", xlab = "Temperature difference [°C], \n [1901-1920] w.r.t. [1871-1890]", las = 1, xaxs="i", main = "Coastal Grid Cells")
  axis(side = 1, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(y = c(0, 10), x = c(0, 0), col = "black", lty = 2)
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 7.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 7.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 7.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 7.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 7.9, labels = c("Global"), adj = 0)
    text(x = 0.2, y = 7.9, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    lines(y = c(7.1, 7.1), x = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 6.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 6.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 6.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 6.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 6.9, labels = c("Northern Hemisphere \n Extra-Tropics"), adj = 0)
    text(x = 0.2, y = 6.9, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 5.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 5.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 5.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 5.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 5.9, labels = c("Tropics"), adj = 0)
    text(x = 0.2, y = 5.9, labels = paste("n = ", length(ix), "   p<.01 **", sep=""), col = "grey40", adj = 0)
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    at = 4.5; label = "Southern Hemisphere \n Extra-Tropics"
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    lines(y = c(4.1, 4.1), x = c(-2, 2), col = "grey40")
  
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    at = 3.5; label = "Coastal Pacific"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 2.5; label = "Coastal Atlantic"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 1.5; label = "Coastal Indian Ocean"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.05 *", sep=""), col = "grey40", adj = 0)
  }
  legend("bottomright", c("Land", "Ocean"), col = c(make.transparent.color("darkorange", alpha = 75), make.transparent.color("blue", alpha = 75)), lty = 1, lwd = 8, inset = 0.02, bg = "white")
}
dev.off()



# with at least 20% of data:
HadSST4_diff = get.change(cur.raster = HadSST4, cur.nobs = HadSST4_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)
CRUTEM5_diff = get.change(cur.raster = CRUTEM5, cur.nobs = CRUTEM5_nobs, ref.period = 1871:1890, test.period = 1901:1920, frac.data = 0.2)


pdf("figures/03_spatial_anomaly/fig3c_frac20.pdf", height = 4.5, width = 2, pointsize = 5)
{
  par(mfrow = c(1, 1), mar = c(4,1,2,1))
  plot(c(1,1), type="n", ylim = c(1, 8), xlim = c(-1.1, 1.1),
       ylab = "", yaxt="n", xlab = "Temperature difference [°C], \n [1901-1920] w.r.t. [1871-1890]", las = 1, xaxs="i", main = "Coastal Grid Cells")
  axis(side = 1, at = seq(-1, 1, 0.1), labels=F, tcl=0.2)
  
  lines(y = c(0, 10), x = c(0, 0), col = "black", lty = 2)
  {
    # Global:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)))
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 7.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 7.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 7.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 7.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 7.9, labels = c("Global"), adj = 0)
    text(x = 0.2, y = 7.9, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    lines(y = c(7.1, 7.1), x = c(-2, 2), col = "grey40")
    
    # NH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] > 30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 6.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 6.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 6.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 6.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 6.9, labels = c("Northern Hemisphere \n Extra-Tropics"), adj = 0)
    text(x = 0.2, y = 6.9, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    # Tropics < 30°N & > -30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = 5.5, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = 5.5, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = 5.45, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = 5.45, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = 5.9, labels = c("Tropics"), adj = 0)
    text(x = 0.2, y = 5.9, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    # SH Extratropics > 30°N:
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & coordinates(CRUTEM5_diff)[,2] < -30)
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    at = 4.5; label = "Southern Hemisphere \n Extra-Tropics"
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    lines(y = c(4.1, 4.1), x = c(-2, 2), col = "grey40")
    
    # Pacific Ocean:
    # (1) 120°E, 60°N, 180°E, -75°N
    # (2) -180°E, 60°N, -75°E, -75°N
    ix = c(which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > 120 & coordinates(CRUTEM5_diff)[,1] < 180 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75),
           which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                   coordinates(CRUTEM5_diff)[,1] > -180 & coordinates(CRUTEM5_diff)[,1] < -75 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75))
    at = 3.5; label = "Coastal Pacific"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    
    # Atlantic: 
    # -75°E, 60°N, 20°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > -75 & coordinates(CRUTEM5_diff)[,1] < 20 & coordinates(CRUTEM5_diff)[,2] < 60 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 2.5; label = "Coastal Atlantic"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.001 ***", sep=""), col = "grey40", adj = 0)
    
    # Indian Ocean: 
    # 20°E, 30°N, 120°E, -75°N
    ix = which(values(!is.na(HadSST4_diff) & !is.na(CRUTEM5_diff)) & 
                 coordinates(CRUTEM5_diff)[,1] > 20 & coordinates(CRUTEM5_diff)[,1] < 120 & coordinates(CRUTEM5_diff)[,2] < 30 & coordinates(CRUTEM5_diff)[,2] > -75)
    at = 1.5; label = "Coastal Indian Ocean"
    t.test(values(CRUTEM5_diff)[ix], values(HadSST4_diff)[ix])
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("darkorange", alpha = 75), side="right", at = at, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(CRUTEM5_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("blue", alpha = 75), side="left", at = at - 0.05, wex = 0.5, axes = F, add = T, pchMed = 21, horizontal = T)
    vioplot(x = values(HadSST4_diff)[ix], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at - 0.05, wex = 1.5, axes = F, add = T, pchMed = 21, horizontal = T)
    text(x = -1, y = at + 0.4, labels = label, adj = 0)
    text(x = 0.2, y = at + 0.4, labels = paste("n = ", length(ix), "   p<.05 *", sep=""), col = "grey40", adj = 0)
  }
  legend("bottomright", c("Land", "Ocean"), col = c(make.transparent.color("darkorange", alpha = 75), make.transparent.color("blue", alpha = 75)), lty = 1, lwd = 8, inset = 0.02, bg = "white")
}
dev.off()



