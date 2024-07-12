
# ------------------------------------------------------------------------------------
# Illustrate scatterplot of temperature reconstructions
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 23.11.2023

##


# 00.(a) load  respective functions & code:
source("code/_convenience/frenchcolormap.R")
source("code/_functions_CMIP6.R")

# get plotting for fingerprints:
source("code/_convenience/_plot_projected_worldmap_v2.R")

# setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/")

png(filename = "figures/00_evaluation/04_illustration_scatterplot.png", width = 8, height = 7, units = "in", res = 200)
{
  par(mar=c(4,5,2,1), mfrow=c(2,2))
  
  xlim <- ylim <- c(-1, 1.5)

  ### tas_land: annual 1895 Predicted vs. observed:
  load("data/03_processedCMIP6_reconstr/CMIP6.tas_land_CMIP6_illustration_v5.RData")
  obs = CMIP_illustration_1895$Y_ann
  pred = CMIP_illustration_1895$Yhat_ann_p1_pt1
  # obs1 = obs[which(!is.na(pred))]
  # pred1 = lm(obs ~ pred)$fitted
  plot(y = obs, x = pred, 
       xlim = xlim, ylim = ylim, main = "Land-based 1895 Coverage incl. Unc./Bias",
       ylab = "CMIP6 GMST [°C]", xlab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16, cex = 0.6)   # Monthly predicted vs. observed for 1895-06 time step.
  mtext("a", side = 3, line = 0.5, at = -1.4, adj = 0, cex = 0.9, font = 2)
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs)^2, 2), sep="")), inset = 0.02)
  
  xlim <- ylim <- c(-1, 1.5)

  ### tos: annual 1895 Predicted vs. observed:
  load("data/03_processedCMIP6_reconstr/CMIP6.tos_CMIP6_illustration.RData")
  obs = CMIP_illustration_1895$Y_ann
  pred = CMIP_illustration_1895$Yhat_ann_p1_pt1
  # obs1 = obs[which(!is.na(pred))]
  # pred1 = lm(obs ~ pred)$fitted  
  plot(y = obs, x = pred, 
       xlim = xlim, ylim = ylim, main = "SST-based 1895 Coverage incl. Unc./Bias",
       ylab = "CMIP6 GMST [°C]", xlab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16, cex = 0.6)   # Monthly predicted vs. observed for 1895-06 time step.
  mtext("b", side = 3, line = 0.5, at = -1.4, adj = 0, cex = 0.9, font = 2)
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs, use = "complete.obs")^2, 2), sep="")), inset = 0.02)
  
  ### tas_land: annual 1995 Predicted vs. observed:
  load("data/03_processedCMIP6_reconstr/CMIP6.tas_land_CMIP6_illustration_v5.RData")
  obs = CMIP_illustration_1995$Y_ann
  pred = CMIP_illustration_1995$Yhat_ann_p1_pt1
  # obs1 = obs[which(!is.na(pred))]
  # pred1 = lm(obs ~ pred)$fitted
  plot(y = obs, x = pred, 
       xlim = xlim, ylim = ylim, main = "Land-based 1995 Coverage incl. Unc./Bias",
       ylab = "CMIP6 GMST [°C]", xlab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16, cex = 0.6)   # Monthly predicted vs. observed for 1895-06 time step.
  mtext("c", side = 3, line = 0.5, at = -1.4, adj = 0, cex = 0.9, font = 2)
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs, use = "complete.obs")^2, 2), sep="")), inset = 0.02)
  
  ### tos: annual 1995 Predicted vs. observed:
  load("data/03_processedCMIP6_reconstr/CMIP6.tos_CMIP6_illustration.RData")
  obs = CMIP_illustration_1995$Y_ann
  pred = CMIP_illustration_1995$Yhat_ann_p1_pt1
  # obs1 = obs[which(!is.na(pred))]
  # pred1 = lm(obs ~ pred)$fitted  
  plot(y = obs, x = pred, 
       xlim = xlim, ylim = ylim, main = "SST-based 1995 Coverage incl. Unc./Bias",
       ylab = "CMIP6 GMST [°C]", xlab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16, cex = 0.6)   # Monthly predicted vs. observed for 1895-06 time step.
  mtext("d", side = 3, line = 0.5, at = -1.4, adj = 0, cex = 0.9, font = 2)
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs, use = "complete.obs")^2, 2), sep="")), inset = 0.02)
}
dev.off()


#### -THE END-
