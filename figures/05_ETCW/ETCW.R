
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023

## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04_master_load_reconstructions.R")



# ------------------------------------------------------------------------------------
# 02. Land-ocean warming ratio during ETCW:
# ------------------------------------------------------------------------------------



### INCLUDE GMST CONSTRAINT!!
## try to constrain 
plot(x = CMIP6.trends$GMST3_true, y = CMIP6.trends$GMSST3_true)
## look at Neukom:
test1 = apply(X = neukom2019_BHM[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test2 = apply(X = neukom2019_CPS_new[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test3 = apply(X = neukom2019_DA[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test4 = apply(X = neukom2019_M08[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test5 = apply(X = neukom2019_OIE[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test6 = apply(X = neukom2019_PAI[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))
test7 = apply(X = neukom2019_PCR[1850:2000,2:1001], MARGIN = 2, FUN=function(x) get.trend(x = x, trend.years = trend.years, years = 1850:2000))

lines(x = quantile(test1[3,], probs = c(0.025, 0.975)), y = rep(-0.2, 2), col = "red")
lines(x = quantile(test2[3,], probs = c(0.025, 0.975)), y = rep(-0.15, 2), col = "red")
lines(x = quantile(test3[3,], probs = c(0.025, 0.975)), y = rep(-0.1, 2), col = "red")
lines(x = quantile(test4[3,], probs = c(0.025, 0.975)), y = rep(-0.05, 2), col = "red")
lines(x = quantile(test5[3,], probs = c(0.025, 0.975)), y = rep(0, 2), col = "red")
lines(x = quantile(test6[3,], probs = c(0.025, 0.975)), y = rep(0.05, 2), col = "red")
lines(x = quantile(test7[3,], probs = c(0.025, 0.975)), y = rep(0.1, 2), col = "red")



setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/05_ETCW/")

## ... continue here... !

library(plotrix)
library(RColorBrewer)
library(mvtnorm)
library(ellipse)

red.cols = brewer.pal(n = 9, name = "Reds")

ylim = c(-0.5, 1)
xlim = c(-0.5, 2)


pdf(file = "_01_land_ocean_v2.pdf", width = 6, height = 4.5)
par(mar = c(4,4,1,0), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
layout(matrix(1:2, ncol = 2, byrow = T), widths = c(9,3), heights = c(3,9))

{  
plot(x = 1850:2020, y = 1850:2020, type="n", 
     xlab = "Land warming [°C per n years]", ylab = "Ocean warming [°C per n years]", main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxs="i")

axis(side = 2, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
axis(side = 1, at = seq(-1, 2, 0.05), tcl=0.2, labels=F)
axis(side=4, tick = T, labels = F)
axis(side = 4, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)

# points(y = trends_piControl$Y$GMSST * 50, x = trends_piControl$Y$GMLSAT_NI * 50, col = make.transparent.color(red.cols[2], alpha = 250), pch = 16, cex = 0.6)

# (1) CMIP6 models provide strong constraint between land and ocean:
# points(y = CMIP6.trends$GMSST5_true, x = CMIP6.trends$GMLSAT5_true, col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
points(y = CMIP6.trends$GMSST4_true, x = CMIP6.trends$GMLSAT4_true, col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
points(y = CMIP6.trends$GMSST3_true, x = CMIP6.trends$GMLSAT3_true, col = make.transparent.color(red.cols[7], alpha = 80), pch = 16, cex = 0.6)
points(y = CMIP6.trends$GMSST2_true, x = CMIP6.trends$GMLSAT2_true, col = make.transparent.color(red.cols[5], alpha = 100), pch = 16, cex = 0.6)
points(y = CMIP6.trends$GMSST1_true, x = CMIP6.trends$GMLSAT1_true, col = make.transparent.color(red.cols[3], alpha = 100), pch = 16, cex = 0.6)

# show CMIP6 models in terms of Ellipses:
# data = cbind(CMIP6.trends$GMLSAT5_true, CMIP6.trends$GMSST5_true)
# data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
# lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
data = cbind(CMIP6.trends$GMLSAT4_true, CMIP6.trends$GMSST4_true)
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
data = cbind(CMIP6.trends$GMLSAT3_true, CMIP6.trends$GMSST3_true)
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[7])
data = cbind(CMIP6.trends$GMLSAT2_true, CMIP6.trends$GMSST2_true)
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[5])
data = cbind(CMIP6.trends$GMLSAT1_true, CMIP6.trends$GMSST1_true)
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[3])



## Look at land/ocean warming ratio based on Byrne et al. theory:
dq_O = 0.09 / 1000  # in kg/kg per decade
dT_O = 0.13 # in K pre decade
gamma = 0.72 # ERA-Interim
Lv = 2257  # kJ / kg latent heat of vaporization
cp =  1.0 # kJ/kg/K specific heat capactiy of air at constant pressure
A = 1 + Lv / cp * (1 - gamma) * dq_O / dT_O # A = land-ocean amplification factor
# lines(x = seq(-0.5, 2, 0.5) * A, y = seq(-0.5, 2, 0.5))

## Look at land/ocean warming from raw observations:
# SST_slope_raw = lm(HadSST4.global.annual$Anomaly[which(HadSST4.global.annual$Year %in% 1901:1950)] ~ c(1:50))$coefficients[2] * 50
# Land_slope_raw = lm(CRUTEM5.global.annual$Anomaly[which(CRUTEM5.global.annual$Time %in% 1901:1950)] ~ c(1:50))$coefficients[2] * 50
# points(x = Land_slope_raw, y = SST_slope_raw, pch = 16)

## Look at Cowtan hybrid36:
#hybrid36_slope = lm(hybrid36.annual$Anomaly[which(hybrid36.annual$Year %in% 1901:1950)] ~ c(1:50))$coefficients[2] * 50
#BEST_Land_slope = lm(BEST_Land.annual$Anomaly[which(BEST_Land.annual$Year %in% 1901:1950)] ~ c(1:50))$coefficients[2] * 50
#COBE_SST2_slope = lm(COBE_SST2.annual$Anomaly[which(COBE_SST2.annual$Year %in% 1901:1950)] ~ c(1:50))$coefficients[2] * 50
# points(x = Land_slope_raw, y = hybrid36_slope, col = "blue", pch = 16)


# Plot different land-ocean warming ratios:
for (i in 1:4) {
  plotCI(x = tas_land_trends_GMLSAT_q[2,i], y = tos_trends_GMSST_q[2,i], li = tas_land_trends_GMLSAT_q[1,i], ui = tas_land_trends_GMLSAT_q[3,i], 
         err = "x", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
  plotCI(x = tas_land_trends_GMLSAT_q[2,i], y = tos_trends_GMSST_q[2,i], li = tos_trends_GMSST_q[1,i], ui = tos_trends_GMSST_q[3,i], 
         err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
  
  plotCI(x = tas_land_trends_GMLSAT_q[2,i], y = tos_hybrid36.cor_trends_GMSST_q[2,i], li = tas_land_trends_GMLSAT_q[1,i], ui = tas_land_trends_GMLSAT_q[3,i], 
         err = "x", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
  plotCI(x = tas_land_trends_GMLSAT_q[2,i]+0.01, y = tos_hybrid36.cor_trends_GMSST_q[2,i], li = tos_hybrid36.cor_trends_GMSST_q[1,i], ui = tos_hybrid36.cor_trends_GMSST_q[3,i], 
         err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
}
lines(x = tas_land_trends_GMLSAT_q[2,1:4], y = tos_trends_GMSST_q[2,1:4], col = "blue", lty = 2, lwd = 2)
lines(x = tas_land_trends_GMLSAT_q[2,1:4], y = tos_hybrid36.cor_trends_GMSST_q[2,1:4], col = "darkgoldenrod4", lty = 2, lwd = 2)

text(x = tas_land_trends_GMLSAT_q[2,1], y = tos_trends_GMSST_q[1,1] - 0.08, labels = c("1850-1890"), 
     adj = 0.5, col = "black", cex = 0.8)# red.cols[2])
text(x = tas_land_trends_GMLSAT_q[2,2], y = tos_trends_GMSST_q[1,2] - 0.08, labels = c("1876-1925"), 
     adj = 0.5, col = "black", cex = 0.8)
text(x = tas_land_trends_GMLSAT_q[2,3], y = tos_trends_GMSST_q[3,3] + 0.08, labels = c("1901-1950"), 
     adj = 0.5, col = "black", cex = 0.8)
text(x = tas_land_trends_GMLSAT_q[2,4], y = tos_trends_GMSST_q[1,4] - 0.08, labels = c("1965-2014"), 
     adj = 0.5, col = "black", cex = 0.8)


legend("bottomright", c("Original Reconstruction", "(CRUTEM5+HadSST4)",
                        "CRUTEM5+Coastal Hybrid SST", "", 
                        "CMIP6, 1850-1890", "CMIP6, 1876-1925",
                        "CMIP6, 1901-1950", "CMIP6, 1965-2014"), col = c("blue", NA, "darkgoldenrod4", NA, 
                                                     make.transparent.color(red.cols[3], alpha = 200), make.transparent.color(red.cols[5], alpha = 100),
                                                     make.transparent.color(red.cols[7], alpha = 100), make.transparent.color(red.cols[9], alpha = 100)), 
       pch = c(3, NA, 3, NA, 16, 16, 16, 16), 
       inset = 0.02, cex = 0.7)
}

# second part of plot (constraints):
{
par(mar = c(4,0,2,0))
plot(x = c(0.5, 2.5), y = c(1,1), type="n", 
     xlab = "", ylab = "", main = "", ylim = ylim, xlim = c(0.5, 2.5), las=1, yaxs="i", xaxs="i", xaxt= "n", yaxt="n", bty = "n")

text(x = 1, y = 0.85, labels = c("1876-1925"), srt=45, adj = 0.5, col = "black", cex = 0.8)
text(x = 2, y = 0.85, labels = c("1901-1950"), srt=45, adj = 0.5, col = "black", cex = 0.8)
lines(x = c(1.5,1.5), y = c(-1, 3), lty = 3)
ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)

# 1876-1925
# Recon:
plotCI(x = 0.7, y = tos_trends_GMSST_q[2,2], li = tos_trends_GMSST_q[1,2], ui = tos_trends_GMSST_q[3,2], 
       err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.02)
# Ocean2k-based constraint:
ocean2k_constraint2 = get.linear.model.constraint(y = CMIP6.trends$GMSST2_true, x = CMIP6.trends$Tropics2_true, x_new = rep(ocean2k_trends[2], 2), plot.constraint = F)
plotCI(x = 0.9, y = ocean2k_constraint2$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, 
       liw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.02)
# CRUTEM5 based constraint:
land_constraint2 = get.linear.model.constraint_ens(y = CMIP6.trends$GMSST2_true, x = CMIP6.trends$GMSST_tas_land2, x_new = tas_land_trends_GMSST[2,], plot.constraint = F)
plotCI(x = 1.1, y = land_constraint2$mean_out, uiw = qnorm(p = 0.975) * land_constraint2$sd_out, 
       liw = qnorm(p = 0.975) * land_constraint2$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.02)
# Cowtan corrections:
plotCI(x = 1.3, y = tos_hybrid36.cor_trends_GMSST_q[2,2], li = tos_hybrid36.cor_trends_GMSST_q[1,2], ui = tos_hybrid36.cor_trends_GMSST_q[3,2], 
       err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.02)



# 1901-1950
# Recon:
plotCI(x = 1.7, y = tos_trends_GMSST_q[2,3], li = tos_trends_GMSST_q[1,3], ui = tos_trends_GMSST_q[3,3], 
       err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.02)
# Ocean2k-based constraint:
ocean2k_constraint3 = get.linear.model.constraint(y = CMIP6.trends$GMSST3_true, x = CMIP6.trends$Tropics3_true, x_new = rep(ocean2k_trends[3], 2), plot.constraint = F)
plotCI(x = 1.9, y = ocean2k_constraint3$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, 
       liw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.02)
# CRUTEM5 based constraint:
land_constraint3 = get.linear.model.constraint_ens(y = CMIP6.trends$GMSST3_true, x = CMIP6.trends$GMSST_tas_land3, x_new = tas_land_trends_GMSST[3,], plot.constraint = F)
plotCI(x = 2.1, y = land_constraint3$mean_out, uiw = qnorm(p = 0.975) * land_constraint3$sd_out, 
       liw = qnorm(p = 0.975) * land_constraint3$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.02)
# Cowtan corrections:
plotCI(x = 2.3, y = tos_hybrid36.cor_trends_GMSST_q[2,3], li = tos_hybrid36.cor_trends_GMSST_q[1,3], ui = tos_hybrid36.cor_trends_GMSST_q[3,3], 
       err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.02)

text(x = 1.7, y = -0.49, labels = c("HadSST4 Reconstr."), srt=90, adj = 0, col = "blue", cex = 0.75)
text(x = 1.9, y = -0.49, labels = c("Ocean2k Constraint"), srt=90, adj = 0, col = "darkorchid4", cex = 0.75)
text(x = 2.1, y = -0.49, labels = c("CRUTEM5 Constraint"), srt=90, adj = 0, col = "darkorange", cex = 0.75)
text(x = 2.3, y = -0.49, labels = c("Coastal Hybrid SST"), srt=90, adj = 0, col = "darkgoldenrod4", cex = 0.75)

## 1901-1950 trend: 
# ix = which(COBE_SST2.annual$Year %in% 1901:1950)
# points(x = tas_land_trends_GMLSAT_q[2,3], y = lm(COBE_SST2.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 1, cex = 1)
# ix = which(ERSSTv5.annual$Year %in% 1901:1950)
# points(x = tas_land_trends_GMLSAT_q[2,3], y = lm(ERSSTv5.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 2, cex = 1)
# ix = which(HadSST4.global.annual$Year %in% 1901:1950)
# points(x = tas_land_trends_GMLSAT_q[2,3], y = lm(HadSST4.global.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 4, cex = 1)
# ix = which(hybrid36.annual$Year %in% 1901:1950)
# points(x = tas_land_trends_GMLSAT_q[2,3], y = lm(hybrid36.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "darkgoldenrod4", pch = 4, cex = 2)
}
  dev.off()
  
  




## plot constraint values for periods 1876-1925 and *1901-1950*
## ------------------------------------------------------------
 getwd() 
  
  
  pdf(file = "_02_SST_constraint.pdf", width = 8, height = 4)
  {
  par(mar = c(4,1,1,1), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
  
  plot(c(1,1), type="n", 
     ylab = "", xlab = "Early 20th Century SST Warming [°C], 1901-1950", main = "",  ylim = c(0.8, 2.8), xlim = c(0,0.8), las=1, yaxs="i", xaxs="i", yaxt="n")

# plot original trends:
  plotCI(y = 1, x = tos_trends_GMSST_q[2,3], ui = tos_trends_GMSST_q[3,3], li = tos_trends_GMSST_q[1,3], col = "grey40", add = T, pch = 16, lwd = 2, cex = 2, err = "x")
  text(x = 0.05, y = 1, labels = "HadSST4-reconstr.", adj = 0)
  
# plot other SST datasets:
  ix = which(COBE_SST2.annual$Year %in% 1901:1950)
  points(y = 1.25, x = lm(COBE_SST2.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 1, cex = 2)
  ix = which(ERSSTv5.annual$Year %in% 1901:1950)
  points(y = 1.25, x = lm(ERSSTv5.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 2, cex = 2)
  ix = which(HadSST4.global.annual$Year %in% 1901:1950)
  points(y = 1.25, x = lm(HadSST4.global.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "blue", pch = 4, cex = 2)
  ix = which(hybrid36.annual$Year %in% 1901:1950)
  points(y = 1.25, x = lm(hybrid36.annual$Anomaly[ix] ~ c(1:50))$coefficients[2] * 50, col = "darkgoldenrod4", pch = 4, cex = 2)
  

# plot land-constraint:
  land_constraint = get.linear.model.constraint_ens(y = CMIP6.trends$GMSST3_true, x = CMIP6.trends$GMSST_tas_land3, x_new = tas_land_trends_GMSST[3,], plot.constraint = F)
  plotCI(y = 1.5, x = land_constraint$mean_out, uiw = qnorm(p = 0.975) * land_constraint$sd_out, liw = qnorm(p = 0.975) * land_constraint$sd_out, col = "darkorange", add = T, pch = 16, lwd = 2, cex = 2, err="x")
  text(y = 1.5, x = 0.5, labels = "CRUTEM5-constraint", adj = 0)

# plot Combined hydro:
  plotCI(y = 1.75, x = tos_hybrid36.cor_trends_GMSST_q[2,3], ui = tos_hybrid36.cor_trends_GMSST_q[3,3], li = tos_hybrid36.cor_trends_GMSST_q[1,3], col = "blue", add = T, pch = 16, lwd = 2, cex = 2, err = "x")
  text(y = 1.75, x = 0.5, labels = "SST Reconstruction, \n combined coastal*", adj = 0)
  

# Tropics2k constraint:
  ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)
  ocean2k_constraint = get.linear.model.constraint(y = CMIP6.trends$GMSST3_true, x = CMIP6.trends$Tropics3_true, x_new = rep(ocean2k_trends[3], 2), plot.constraint = F)
  
  plotCI(y = 2, x = ocean2k_constraint$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint$sd_out, liw = qnorm(p = 0.975) * ocean2k_constraint$sd_out, col = "blueviolet", add = T, pch = 16, cex = 2, err = "x")
  text(y = 2, x = 0.5, labels = "Ocean2k-constraint (no unc.)", adj = 0)

# Global mean sea level constraint:
  MSL.constraint = c(-0.7750074, 0.4407971)
  plotCI(y = 2.25, x = 0, ui = MSL.constraint[2], li = MSL.constraint[1], col = "darkgoldenrod4", add = T, pch = NA, cex = 2, err = "x")
  text(y = 2.25, x = 0.5, labels = "Constraint based on \n thermosteric sea level", adj = 0)
  
  legend("top", c("HadSST4", "COBE2-SST", "ERSST5", "Cowtan-hybrid36"), pch = c(4, 1, 2, 4), col = c("blue", "blue", "blue", "darkgoldenrod4"), ncol = 2)
  
  dev.off()
  }
  
  
  
### OK!
  
# Neukom CFRs:
#  Neukom_GMSST_trends = apply(X = Neukom_GMSST, MARGIN = 1, FUN=function(x)  get.trend(x = x, trend.years = trend.years, years = 1850:2000))
#  plotCI(x = 2.25, y = median(Neukom_GMSST_trends[3, 1:100]), ui = quantile(Neukom_GMSST_trends[3, 1:100], 0.975), li = quantile(Neukom_GMSST_trends[3, 1:100], 0.025), err = "y", add = T, col = "blueviolet")
#  plotCI(x = 2.5, y = median(Neukom_GMSST_trends[3, 101:200]), ui = quantile(Neukom_GMSST_trends[3, 101:200], 0.975), li = quantile(Neukom_GMSST_trends[3, 101:200], 0.025), err = "y", add = T, col = "blueviolet")
#  plotCI(x = 2.75, y = median(Neukom_GMSST_trends[3, 201:300]), ui = quantile(Neukom_GMSST_trends[3, 201:300], 0.975), li = quantile(Neukom_GMSST_trends[3, 201:300], 0.025), err = "y", add = T, col = "blueviolet")
#  plotCI(x = 3, y = median(Neukom_GMSST_trends[3, 301:400]), ui = quantile(Neukom_GMSST_trends[3, 301:400], 0.975), li = quantile(Neukom_GMSST_trends[3, 301:400], 0.025), err = "y", add = T, col = "blueviolet")
# 
  


##### THE END ########
  

# Cold ocean anomaly shows up clearly. Reasons why this cold ocean anomaly may be unrealistic:
# (1) literature: biases in ocean are prevalent; comparison to coastal stations appears to show that ocean anomaly is unrealistic.
# (2) "decoupling" of low-frequency variability (poor match between land- and ocean reconstruction) from high-frequency variability (very good match).
# (3) none of the cmip6 models shows a temporal "decoupling pattern" between low- and high-frequency variability.
# (4) none of the cmip6 models supports an ocean-land temperature difference close to that observed.





