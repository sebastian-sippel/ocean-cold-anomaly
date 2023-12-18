
# ------------------------------------------------------------------------------------
# Evaluate tos and tas reconstruction(s) based on CMIP6.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 23.03.2022


# screen -S make_plots
# R

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_functions_CMIP5_extr.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")

# get plotting for fingerprints:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_plot_projected_worldmap_v2.R")


# 00.(c) Load GSAT/GMST reconstructions:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tas_land.df_v5.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tos.df_v5.RData")



### 1a. Evaluation of reconstruction MSE:
### -----------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/")

library(RColorBrewer)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(99)
map.col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(99)

paired.col = brewer.pal(n = 12, name = "Paired")



### Time series of ANNUAL MSE: 

## GMST_FM:
cur.var = "GMST_FM"


pdf("01_annual_evaluation_GMST-reconstr_MSE.pdf", height=7, width=8)
par(mfrow=c(2,2), mar=c(3,4.5,2,0.5))
{
  ylim = c(0, 0.04)
  xlim = c(1850, 2020)
  
  ### Land-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Annual temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "Land-based GMST Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p0_pt0, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt0, col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p1_pt0, col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p0_pt1, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt1, col = "grey50", lwd = 2)
  lines(x = 1850:2020, CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p1_pt1, col = paired.col[6], lwd = 2)
  
  legend("topright", c("Evaluation, without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation, incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
         col = c("black", "grey50", paired.col[5], paired.col[6], "black", "black", "grey50", paired.col[5], paired.col[6]), 
         lty = c(NA, 2, 2, 2, NA, NA, 1, 1, 1), lwd = c(NA, 2,2,2,NA, NA, 2,2,2), inset = 0.02, cex = 0.8)
  
  
  ### SST-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Annual temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "SST-based GMST Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p0_pt0, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt0, col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p1_pt0, col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p0_pt1, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt1, col = "grey50", lwd = 2)
  lines(x = 1850:2020, CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p1_pt1, col = paired.col[6], lwd = 2)
  
  
  
  ### Land-based reconstruction, relative error:
  ylim = c(-80, 80)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = ((CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p0_pt0 / CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt0) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = ((CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p1_pt0 / CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt0) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = ((CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p0_pt1 / CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt1) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = ((CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_p1_pt1 / CMIP6.tas_land.df$ann_MSE[[cur.var]]$mod_gta_pt1) - 1) * 100, col = paired.col[6], lwd = 2)
  
  # ylim = c(-75, 75)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = ((CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p0_pt0 / CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt0) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = ((CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p1_pt0 / CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt0) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = ((CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p0_pt1 / CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt1) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = ((CMIP6.tos.df$ann_MSE[[cur.var]]$mod_p1_pt1 / CMIP6.tos.df$ann_MSE[[cur.var]]$mod_gta_pt1) - 1) * 100, col = paired.col[6], lwd = 2)
}
dev.off()







#### MONTHLY ERROR ESTIMATES:

pdf("01_evaluation_GMST-reconstr_MSE.pdf", height=7, width=8)
par(mfrow=c(2,2), mar=c(3,4.5,2,0.5))
{
  ylim = c(0, 0.06)
  xlim = c(1850, 2020)
  
  ### Land-based reconstruction:
plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
     ylab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "Land-based GMST Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  legend("topright", c("Evaluation, without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation, incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
         col = c("black", "grey50", paired.col[5], paired.col[6], "black", "black", "grey50", paired.col[5], paired.col[6]), 
         lty = c(NA, 2, 2, 2, NA, NA, 1, 1, 1), lwd = c(NA, 2,2,2,NA, NA, 2,2,2), inset = 0.02, cex = 0.8)
  
  
  ### SST-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "SST-based GMST Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  

  ### Land-based reconstruction, relative error:
  ylim = c(-60, 60)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
     ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
     main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 3)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 3)

  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
  
  # ylim = c(-75, 75)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 3)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 3)
  
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
}
dev.off()






## GSAT:
cur.var = "GSAT"


pdf("01_evaluation_GSAT-reconstr_MSE.pdf", height=7, width=8)
par(mfrow=c(2,2), mar=c(3,4.5,2,0.5))
{
  ylim = c(0, 0.06)
  xlim = c(1850, 2020)
  
  ### Land-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "Land-based GSAT Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  legend("topright", c("Evaluation, without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation, incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
         col = c("black", "grey50", paired.col[5], paired.col[6], "black", "black", "grey50", paired.col[5], paired.col[6]), 
         lty = c(NA, 2, 2, 2, NA, NA, 1, 1, 1), lwd = c(NA, 2,2,2,NA, NA, 2,2,2), inset = 0.02, cex = 0.8)
  
  
  ### SST-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "SST-based GSAT Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  
  
  ### Land-based reconstruction, relative error:
  ylim = c(-60, 60)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
  
  # ylim = c(-75, 75)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
}
dev.off()







### Annual GSAT reconstruction:


## GSAT:
cur.var = "GSAT"
Y = rowMeans(sapply(X = 1:12, FUN=function(mon.ix) CMIP6.tas_land.df$mon[[mon.ix]]$Y[[cur.var]]))
# CMIP6.tas_land.df$mon_MSE$GSAT


pdf("01_evaluation_GSAT-reconstr_MSE_annual.pdf", height=7, width=8)
par(mfrow=c(2,2), mar=c(3,4.5,2,0.5))
{
  ylim = c(0, 0.06)
  xlim = c(1850, 2020)
  
  ### Land-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Annual temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "Land-based GSAT Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  CMIP6.tas_land.df$mon[[1]]$Yhat$GMSST$mod_p1$pt0.1se
  
  
  CMIP6.tas_land.df$ann[[cur.var]]$mod_p0$pt1.min
  
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  legend("topright", c("Evaluation, without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation, incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
         col = c("black", "grey50", paired.col[5], paired.col[6], "black", "black", "grey50", paired.col[5], paired.col[6]), 
         lty = c(NA, 2, 2, 2, NA, NA, 1, 1, 1), lwd = c(NA, 2,2,2,NA, NA, 2,2,2), inset = 0.02, cex = 0.8)
  
  
  ### SST-based reconstruction:
  plot(c(1,1), xlim = xlim, ylim = ylim, type='l', las = 1,
       ylab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), xlab = "", main = "SST-based GSAT Reconstruction")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.002), tcl=0.2, labels=F)
  
  # Evaluation without uncertainties+biases:
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x))), col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x))), col = "grey50", lwd = 2, lty = 2)
  lines(x = 1850:2020, (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x))), col = paired.col[6], lwd = 2, lty = 2)
  
  # Evaluation With uncertainties+biases:
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)), col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)), col = "grey50", lwd = 2)
  lines(x = 1850:2020, rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)), col = paired.col[6], lwd = 2)
  
  
  
  ### Land-based reconstruction, relative error:
  ylim = c(-60, 60)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tas_land.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
  
  # ylim = c(-75, 75)
  plot(c(1,1), xlim = c(1850,2020), ylim = ylim, type='l', las = 1,
       ylab = "MSE [in %] relative to GTA-estimator", xlab = "", 
       main = "")
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 5), tcl=0.2, labels=F)
  
  
  lines(x = c(1850,2020), y = c(0,0), col = "grey50", lty = 1)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2, lty = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt0, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt0, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2, lty = 2)
  
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p0_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[5], lwd = 2)
  lines(x = 1850:2020, y = (rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_p1_pt1, FUN=function(x) x)) / rowMeans(sapply(CMIP6.tos.df$mon_MSE[[cur.var]]$mod_gta_pt1, FUN=function(x) x)) - 1) * 100, col = paired.col[6], lwd = 2)
}
dev.off()










