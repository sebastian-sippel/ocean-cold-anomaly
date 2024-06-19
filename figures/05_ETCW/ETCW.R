
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023

## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04b_master_read_paleo_reconstructions.R")
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04c_compute_trends_4paleo-comparison.R")




# ------------------------------------------------------------------------------------
# 02. Land-ocean warming ratio during ETCW:
# ------------------------------------------------------------------------------------

setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/05_ETCW/")


library(plotrix)
library(RColorBrewer)
library(mvtnorm)
library(ellipse)

red.cols = brewer.pal(n = 9, name = "Reds")

ylim = c(-0.5, 1)
xlim = c(-0.5, 2)


# ------------------------------------------------------------------------------------
# 40-years figure
# ------------------------------------------------------------------------------------

pdf(file = "05_land_ocean_40_years.pdf", width = 8, height = 5)
par(mar = c(4,4,1,0), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
layout(matrix(1:2, ncol = 2, byrow = T), widths = c(9,4), heights = c(4,9))

{  
plot(x = 1850:2020, y = 1850:2020, type="n", 
     xlab = "Land warming [°C per 40 years]", ylab = "Ocean warming [°C per 40 years]", main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxs="i")

axis(side = 2, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
axis(side = 1, at = seq(-1, 2, 0.05), tcl=0.2, labels=F)
axis(side=4, tick = T, labels = F)
axis(side = 4, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)

# points(y = trends_piControl$Y$GMSST * 50, x = trends_piControl$Y$GMLSAT_NI * 50, col = make.transparent.color(red.cols[2], alpha = 250), pch = 16, cex = 0.6)

# (1) CMIP6 models provide strong constraint between land and ocean:
points(y = CMIP6.GMSST_true.trends[,4], x = CMIP6.GMLSAT_true.trends[,4], col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
points(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMLSAT_true.trends[,3], col = make.transparent.color(red.cols[7], alpha = 80), pch = 16, cex = 0.6)
points(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMLSAT_true.trends[,2], col = make.transparent.color(red.cols[5], alpha = 100), pch = 16, cex = 0.6)
points(y = CMIP6.GMSST_true.trends[,1], x = CMIP6.GMLSAT_true.trends[,1], col = make.transparent.color(red.cols[3], alpha = 100), pch = 16, cex = 0.6)

# show CMIP6 models in terms of Ellipses:
# data = cbind(CMIP6.trends$GMLSAT5_true, CMIP6.trends$GMSST5_true)
# data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
# lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
data = cbind(CMIP6.GMLSAT_true.trends[,4], CMIP6.GMSST_true.trends[,4])
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
data = cbind(CMIP6.GMLSAT_true.trends[,3], CMIP6.GMSST_true.trends[,3])
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[7])
data = cbind(CMIP6.GMLSAT_true.trends[,2], CMIP6.GMSST_true.trends[,2])
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[5])
data = cbind(CMIP6.GMLSAT_true.trends[,1], CMIP6.GMSST_true.trends[,1])
data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
lines(x = data.el[,1], y = data.el[,2], col = red.cols[3])


# Plot different land-ocean warming ratios:
i.seq = c(1:4)
for (i in i.seq) {
  plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_tos.trends_q[2,i], li = OBS.GMLSAT_tas_land.trends_q[1,i], ui = OBS.GMLSAT_tas_land.trends_q[3,i], 
         err = "x", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
  plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_tos.trends_q[2,i], li = OBS.GMSST_tos.trends_q[1,i], ui = OBS.GMSST_tos.trends_q[3,i], 
         err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
  
  plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_hybrid36cor.trends_q[2,i], li = OBS.GMLSAT_tas_land.trends_q[1,i], ui = OBS.GMLSAT_tas_land.trends_q[3,i], 
         err = "x", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
  plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i]+0.01, y = OBS.GMSST_hybrid36cor.trends_q[2,i], li = OBS.GMSST_hybrid36cor.trends_q[1,i], ui = OBS.GMSST_hybrid36cor.trends_q[3,i], 
         err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
}
# get further estimates CRUTEM5+HadSST4 Raw Averages:
points(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], 
       y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], pch = 8, col = "blue")
# points(x = get.trend(x = OBS.GMST.Kadow_tas_land[1:171], trend.years = trend.years, years = 1850:2020)[1:4], 
#       y = get.trend(x = OBS.GMST.Kadow_tos[1:171], trend.years = trend.years, years = 1850:2020)[1:4], pch = 15, col = "grey")
lines(x = OBS.GMLSAT_tas_land.trends_q[2,1:4], y = OBS.GMSST_tos.trends_q[2,1:4], col = "blue", lty = 2, lwd = 1)
lines(x = OBS.GMLSAT_tas_land.trends_q[2,1:4], y = OBS.GMSST_hybrid36cor.trends_q[2,1:4], col = "darkgoldenrod4", lty = 2, lwd = 1)
lines(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], 
      y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], col = "blue", lty = 2, lwd = 1)

text(x = OBS.GMLSAT_tas_land.trends_q[2,1], y = OBS.GMSST_tos.trends_q[1,1] - 0.04, labels = c("1851-1890"), 
     adj = 0.5, col = "black", cex = 0.8)
text(x = OBS.GMLSAT_tas_land.trends_q[2,2], y = OBS.GMSST_tos.trends_q[1,2] - 0.02, labels = c("1871-1910"), 
     adj = 0.5, col = "black", cex = 0.8)
text(x = OBS.GMLSAT_tas_land.trends_q[2,3], y = OBS.GMSST_tos.trends_q[3,3] + 0.08, labels = c("1901-1940"), 
     adj = 0.5, col = "black", cex = 0.8)
text(x = OBS.GMLSAT_tas_land.trends_q[2,4], y = OBS.GMSST_tos.trends_q[1,4] - 0.04, labels = c("1975-2014"), 
     adj = 0.5, col = "black", cex = 0.8)


legend("bottomright", c("Original CRUTEM5+HadSST4",
                        "Reconstruction", "(CRUTEM5+HadSST4)",
                        "Rec. CRUTEM5+Coastal ", "Hybrid SST Corrected", 
                        "CMIP6, 1851-1890", "CMIP6, 1871-1910",
                        "CMIP6, 1901-1940", "CMIP6, 1975-2014"), col = c("blue", "blue", NA, "darkgoldenrod4", NA, 
                                                     make.transparent.color(red.cols[3], alpha = 200), make.transparent.color(red.cols[5], alpha = 100),
                                                     make.transparent.color(red.cols[7], alpha = 100), make.transparent.color(red.cols[9], alpha = 100)), 
       pch = c(8, 3, NA, 3, NA, 16, 16, 16, 16), 
       inset = 0.02, cex = 0.7)
}

# second part of plot (constraints):
{
par(mar = c(4,0,2,0))
plot(x = c(0.5, 2.5), y = c(1,1), type="n", 
     xlab = "", ylab = "", main = "", ylim = ylim, xlim = c(0.5, 2.5), las=1, yaxs="i", xaxs="i", xaxt= "n", yaxt="n", bty = "n")

text(x = 1, y = 0.65, labels = c("1871-1910"), srt=45, adj = 0.5, col = "black", cex = 0.8)
text(x = 2, y = 0.65, labels = c("1901-1940"), srt=45, adj = 0.5, col = "black", cex = 0.8)
lines(x = c(1.5,1.5), y = c(-1, 3), lty = 3)
legend("top", c("HadSST4", "ERSSTv5", "COBE-SST2"), title = "Original Trends", col = "blue", pch = c(8,6,4), cex = 0.7, bg = "white", inset = 0.00, ncol = 2)
ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)

# 1876-1925
# Recon:
plotCI(x = 0.6, y = OBS.GMSST_tos.trends_q[2,2], li = OBS.GMSST_tos.trends_q[1,2], ui = OBS.GMSST_tos.trends_q[3,2], 
       err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
# Other Ocean Datasets:
points(x = 0.7, y = HadSST4.trends[2], col = "blue", pch = 8)
points(x = 0.7, y = ERSSTv5.trends[2], col = "blue", pch = 6)
points(x = 0.7, y = COBE_SST2.trends[2], col = "blue", pch = 4)

# Ocean2k-based constraint:
ocean2k_constraint2 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.Tropics_true.trends[,2], x_new = rep(ocean2k_trends[2], 2), plot.constraint = F)
plotCI(x = 0.8, y = ocean2k_constraint2$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, 
       liw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)

# Neukom Paleo-constraint:
# lm(CMIP6.GMSST_true.trends[,2] ~ CMIP6.GMST_true.trends[,2])
# cor(CMIP6.GMSST_true.trends[,2], CMIP6.GMST_true.trends[,2], use = "complete.obs")
# c(sapply(X = neukom2019_trend_all, FUN=function(x) x[2,]))
for (i in 1:7) {
  neukom_constraint = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMST_true.trends[,2], x_new = neukom2019_trend_all[[i]][2,], plot.constraint = F)
  plotCI(x = 0.85+i/20, y = neukom_constraint$mean_out, uiw = qnorm(p = 0.975) * neukom_constraint$sd_out, 
         liw = qnorm(p = 0.975) * neukom_constraint$sd_out, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
}
# Crowley:
# crowley_constraint = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMST_true.trends[,2], x_new = crowley2014.trend[2], plot.constraint = F)
#plotCI(x = 1.25, y = crowley_constraint$mean_out, uiw = qnorm(p = 0.975) * crowley_constraint$sd_out, 
#       liw = qnorm(p = 0.975) * crowley_constraint$sd_out, col = "darkblue", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
# CRUTEM5 based constraint:
land_constraint2 = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMLSAT_NI_tas_land.trends[,2], x_new = OBS.GMLSAT_tas_land.trends[2,], plot.constraint = F)
plotCI(x = 1.3, y = land_constraint2$mean_out, uiw = qnorm(p = 0.975) * land_constraint2$sd_out, 
       liw = qnorm(p = 0.975) * land_constraint2$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
# Cowtan corrections:
plotCI(x = 1.4, y = OBS.GMSST_hybrid36cor.trends_q[2,2], li = OBS.GMSST_hybrid36cor.trends_q[1,2], ui = OBS.GMSST_hybrid36cor.trends_q[3,2], 
       err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.01)


# 1901-1950
# Recon:
plotCI(x = 1.6, y = OBS.GMSST_tos.trends_q[2,3], li = OBS.GMSST_tos.trends_q[1,3], ui = OBS.GMSST_tos.trends_q[3,3], 
       err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
# Other Ocean Datasets:
points(x = 1.7, y = HadSST4.trends[3], col = "blue", pch = 8)
points(x = 1.7, y = ERSSTv5.trends[3], col = "blue", pch = 6)
points(x = 1.7, y = COBE_SST2.trends[3], col = "blue", pch = 4)

# Ocean2k-based constraint:
ocean2k_constraint3 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.Tropics_true.trends[,3], x_new = rep(ocean2k_trends[3], 2), plot.constraint = F)
plotCI(x = 1.8, y = ocean2k_constraint3$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, 
       liw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
# Neukom Paleo-constraint:
for (i in 1:7) {
  neukom_constraint = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMST_true.trends[,3], x_new = neukom2019_trend_all[[i]][3,], plot.constraint = F)
  plotCI(x = 1.85+i/20, y = neukom_constraint$mean_out, uiw = qnorm(p = 0.975) * neukom_constraint$sd_out, 
         liw = qnorm(p = 0.975) * neukom_constraint$sd_out, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
}
# Crowley:
# crowley_constraint = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMST_true.trends[,3], x_new = crowley2014.trend[3], plot.constraint = F)
#plotCI(x = 2.25, y = crowley_constraint$mean_out, uiw = qnorm(p = 0.975) * crowley_constraint$sd_out, 
#       liw = qnorm(p = 0.975) * crowley_constraint$sd_out, col = "darkblue", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)

# CRUTEM5 based constraint:
land_constraint3 = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMLSAT_NI_tas_land.trends[,3], x_new = OBS.GMLSAT_tas_land.trends[3,], plot.constraint = F)
plotCI(x = 2.3, y = land_constraint3$mean_out, uiw = qnorm(p = 0.975) * land_constraint3$sd_out, 
       liw = qnorm(p = 0.975) * land_constraint3$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
# Cowtan corrections:
plotCI(x = 2.4, y = OBS.GMSST_hybrid36cor.trends_q[2,3], li = OBS.GMSST_hybrid36cor.trends_q[1,3], ui = OBS.GMSST_hybrid36cor.trends_q[3,3], 
       err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.01)


text(x = 1.6, y = -0.49, labels = c("HadSST4 Reconstr."), srt=90, adj = 0, col = "blue", cex = 0.75)
text(x = 1.8, y = -0.49, labels = c("Ocean2k Constraint"), srt=90, adj = 0, col = "darkorchid4", cex = 0.75)
text(x = 2.05, y = -0.49, labels = c("Pages2k GMST Constraint"), srt=90, adj = 0, col = "grey", cex = 0.75)
text(x = 2.3, y = -0.49, labels = c("CRUTEM5 Constraint"), srt=90, adj = 0, col = "darkorange", cex = 0.75)
text(x = 2.4, y = -0.49, labels = c("Coastal Hybrid SST"), srt=90, adj = 0, col = "darkgoldenrod4", cex = 0.75)

}
  dev.off()
  
  


  
  
# ------------------------------------------------------------------------------------
# 50-years figure
# ------------------------------------------------------------------------------------
  
pdf(file = "05_land_ocean_50_years.pdf", width = 8, height = 5)
par(mar = c(4,4,1,0), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
layout(matrix(1:2, ncol = 2, byrow = T), widths = c(9,4), heights = c(4,9))
  
{  
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         xlab = "Land warming [°C per 50 years]", ylab = "Ocean warming [°C per 50 years]", main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxs="i")
    
    axis(side = 2, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    axis(side = 1, at = seq(-1, 2, 0.05), tcl=0.2, labels=F)
    axis(side=4, tick = T, labels = F)
    axis(side = 4, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    
    # points(y = trends_piControl$Y$GMSST * 50, x = trends_piControl$Y$GMLSAT_NI * 50, col = make.transparent.color(red.cols[2], alpha = 250), pch = 16, cex = 0.6)
    
    # (1) CMIP6 models provide strong constraint between land and ocean:
    points(y = CMIP6.GMSST_true.trends[,8], x = CMIP6.GMLSAT_true.trends[,8], col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,7], x = CMIP6.GMLSAT_true.trends[,7], col = make.transparent.color(red.cols[7], alpha = 80), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMLSAT_true.trends[,6], col = make.transparent.color(red.cols[5], alpha = 100), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,5], x = CMIP6.GMLSAT_true.trends[,5], col = make.transparent.color(red.cols[3], alpha = 100), pch = 16, cex = 0.6)
    
    # show CMIP6 models in terms of Ellipses:
    # data = cbind(CMIP6.trends$GMLSAT5_true, CMIP6.trends$GMSST5_true)
    # data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
    # lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
    data = cbind(CMIP6.GMLSAT_true.trends[,8], CMIP6.GMSST_true.trends[,8])
    data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
    lines(x = data.el[,1], y = data.el[,2], col = red.cols[9])
    data = cbind(CMIP6.GMLSAT_true.trends[,7], CMIP6.GMSST_true.trends[,7])
    data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
    lines(x = data.el[,1], y = data.el[,2], col = red.cols[7])
    data = cbind(CMIP6.GMLSAT_true.trends[,6], CMIP6.GMSST_true.trends[,6])
    data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
    lines(x = data.el[,1], y = data.el[,2], col = red.cols[5])
    data = cbind(CMIP6.GMLSAT_true.trends[,5], CMIP6.GMSST_true.trends[,5])
    data.el = ellipse(x = cor(data, use = "complete.obs"), scale = colSds(data, na.rm=T), centre = colMeans(data, na.rm=T), level = 0.95)
    lines(x = data.el[,1], y = data.el[,2], col = red.cols[3])
    
    
    # Plot different land-ocean warming ratios:
    i.seq = c(5:8)
    for (i in i.seq) {
      plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_tos.trends_q[2,i], li = OBS.GMLSAT_tas_land.trends_q[1,i], ui = OBS.GMLSAT_tas_land.trends_q[3,i], 
             err = "x", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
      plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_tos.trends_q[2,i], li = OBS.GMSST_tos.trends_q[1,i], ui = OBS.GMSST_tos.trends_q[3,i], 
             err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
      
      plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i], y = OBS.GMSST_hybrid36cor.trends_q[2,i], li = OBS.GMLSAT_tas_land.trends_q[1,i], ui = OBS.GMLSAT_tas_land.trends_q[3,i], 
             err = "x", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
      plotCI(x = OBS.GMLSAT_tas_land.trends_q[2,i]+0.01, y = OBS.GMSST_hybrid36cor.trends_q[2,i], li = OBS.GMSST_hybrid36cor.trends_q[1,i], ui = OBS.GMSST_hybrid36cor.trends_q[3,i], 
             err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.005)
    }
    # get further estimates CRUTEM5+HadSST4 Raw Averages:
    points(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], 
           y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], pch = 8, col = "blue")
    # points(x = get.trend(x = OBS.GMST.Kadow_tas_land[1:171], trend.years = trend.years, years = 1850:2020)[1:4], 
    #       y = get.trend(x = OBS.GMST.Kadow_tos[1:171], trend.years = trend.years, years = 1850:2020)[1:4], pch = 15, col = "grey")
    lines(x = OBS.GMLSAT_tas_land.trends_q[2,5:8], y = OBS.GMSST_tos.trends_q[2,5:8], col = "blue", lty = 2, lwd = 1)
    lines(x = OBS.GMLSAT_tas_land.trends_q[2,5:8], y = OBS.GMSST_hybrid36cor.trends_q[2,1:4], col = "darkgoldenrod4", lty = 2, lwd = 1)
    lines(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], 
          y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], col = "blue", lty = 2, lwd = 1)
    
    text(x = OBS.GMLSAT_tas_land.trends_q[2,5], y = OBS.GMSST_tos.trends_q[1,5] - 0.04, labels = c("1851-1900"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = OBS.GMLSAT_tas_land.trends_q[2,6], y = OBS.GMSST_tos.trends_q[1,6] - 0.04, labels = c("1871-1920"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = OBS.GMLSAT_tas_land.trends_q[2,7], y = OBS.GMSST_tos.trends_q[3,7] + 0.04, labels = c("1901-1950"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = OBS.GMLSAT_tas_land.trends_q[2,8], y = OBS.GMSST_tos.trends_q[1,8] - 0.04, labels = c("1965-2014"), 
         adj = 0.5, col = "black", cex = 0.8)
    
    
    legend("bottomright", c("Original CRUTEM5+HadSST4",
                            "Reconstruction", "(CRUTEM5+HadSST4)",
                            "Rec. CRUTEM5+Coastal ", "Hybrid SST Corrected", 
                            "CMIP6, 1851-1900", "CMIP6, 1871-1920",
                            "CMIP6, 1901-1950", "CMIP6, 1965-2014"), col = c("blue", "blue", NA, "darkgoldenrod4", NA, 
                                                                             make.transparent.color(red.cols[3], alpha = 200), make.transparent.color(red.cols[5], alpha = 100),
                                                                             make.transparent.color(red.cols[7], alpha = 100), make.transparent.color(red.cols[9], alpha = 100)), 
           pch = c(8, 3, NA, 3, NA, 16, 16, 16, 16), 
           inset = 0.02, cex = 0.7)
  }
  
# second part of plot (constraints):
{
    par(mar = c(4,0,2,0))
    plot(x = c(0.5, 2.5), y = c(1,1), type="n", 
         xlab = "", ylab = "", main = "", ylim = ylim, xlim = c(0.5, 2.5), las=1, yaxs="i", xaxs="i", xaxt= "n", yaxt="n", bty = "n")
    
    text(x = 1, y = 0.65, labels = c("1871-1920"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    text(x = 2, y = 0.65, labels = c("1901-1950"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    lines(x = c(1.5,1.5), y = c(-1, 3), lty = 3)
    legend("top", c("HadSST4", "ERSSTv5", "COBE-SST2"), title = "Original Trends", col = "blue", pch = c(8,6,4), cex = 0.7, bg = "white", inset = 0.00, ncol = 2)
    ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)
    
    # 1876-1925
    # Recon:
    plotCI(x = 0.6, y = OBS.GMSST_tos.trends_q[2,6], li = OBS.GMSST_tos.trends_q[1,6], ui = OBS.GMSST_tos.trends_q[3,6], 
           err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
    # Other Ocean Datasets:
    points(x = 0.7, y = HadSST4.trends[6], col = "blue", pch = 8)
    points(x = 0.7, y = ERSSTv5.trends[6], col = "blue", pch = 6)
    points(x = 0.7, y = COBE_SST2.trends[6], col = "blue", pch = 4)
    
    # Ocean2k-based constraint:
    ocean2k_constraint2 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.Tropics_true.trends[,6], x_new = rep(ocean2k_trends[6], 2), plot.constraint = F)
    plotCI(x = 0.8, y = ocean2k_constraint2$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, 
           liw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    
    # Neukom Paleo-constraint:
    # lm(CMIP6.GMSST_true.trends[,2] ~ CMIP6.GMST_true.trends[,2])
    # cor(CMIP6.GMSST_true.trends[,2], CMIP6.GMST_true.trends[,2], use = "complete.obs")
    # c(sapply(X = neukom2019_trend_all, FUN=function(x) x[2,]))
    for (i in 1:7) {
      neukom_constraint = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMST_true.trends[,6], x_new = neukom2019_trend_all[[i]][6,], plot.constraint = F)
      plotCI(x = 0.85+i/20, y = neukom_constraint$mean_out, uiw = qnorm(p = 0.975) * neukom_constraint$sd_out, 
             liw = qnorm(p = 0.975) * neukom_constraint$sd_out, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }
    # Crowley:
    #crowley_constraint = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMST_true.trends[,6], x_new = crowley2014.trend[6], plot.constraint = F)
    #plotCI(x = 1.25, y = crowley_constraint$mean_out, uiw = qnorm(p = 0.975) * crowley_constraint$sd_out, 
    #       liw = qnorm(p = 0.975) * crowley_constraint$sd_out, col = "darkblue", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    # CRUTEM5 based constraint:
    land_constraint2 = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMLSAT_NI_tas_land.trends[,6], x_new = OBS.GMLSAT_tas_land.trends[6,], plot.constraint = F)
    plotCI(x = 1.3, y = land_constraint2$mean_out, uiw = qnorm(p = 0.975) * land_constraint2$sd_out, 
           liw = qnorm(p = 0.975) * land_constraint2$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
    # Cowtan corrections:
    plotCI(x = 1.4, y = OBS.GMSST_hybrid36cor.trends_q[2,6], li = OBS.GMSST_hybrid36cor.trends_q[1,6], ui = OBS.GMSST_hybrid36cor.trends_q[3,6], 
           err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.01)
    
    
    # 1901-1950
    trend.ix = 7
    # Recon:
    plotCI(x = 1.6, y = OBS.GMSST_tos.trends_q[2,trend.ix], li = OBS.GMSST_tos.trends_q[1,trend.ix], ui = OBS.GMSST_tos.trends_q[3,trend.ix], 
           err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
    # Other Ocean Datasets:
    points(x = 1.7, y = HadSST4.trends[trend.ix], col = "blue", pch = 8)
    points(x = 1.7, y = ERSSTv5.trends[trend.ix], col = "blue", pch = 6)
    points(x = 1.7, y = COBE_SST2.trends[trend.ix], col = "blue", pch = 4)
    
    # Ocean2k-based constraint:
    ocean2k_constraint3 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,trend.ix], x = CMIP6.Tropics_true.trends[,trend.ix], x_new = rep(ocean2k_trends[trend.ix], 2), plot.constraint = F)
    plotCI(x = 1.8, y = ocean2k_constraint3$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, 
           liw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    # Neukom Paleo-constraint:
    for (i in 1:7) {
      neukom_constraint = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,trend.ix], x = CMIP6.GMST_true.trends[,trend.ix], x_new = neukom2019_trend_all[[i]][trend.ix,], plot.constraint = F)
      plotCI(x = 1.85+i/20, y = neukom_constraint$mean_out, uiw = qnorm(p = 0.975) * neukom_constraint$sd_out, 
             liw = qnorm(p = 0.975) * neukom_constraint$sd_out, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }
    # Crowley:
    # crowley_constraint = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,trend.ix], x = CMIP6.GMST_true.trends[,trend.ix], x_new = crowley2014.trend[trend.ix], plot.constraint = F)
    # plotCI(x = 2.25, y = crowley_constraint$mean_out, uiw = qnorm(p = 0.975) * crowley_constraint$sd_out, 
    #       liw = qnorm(p = 0.975) * crowley_constraint$sd_out, col = "darkblue", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    
    
    # CRUTEM5 based constraint:
    land_constraint3 = get.linear.model.constraint_ens(y = CMIP6.GMSST_true.trends[,trend.ix], x = CMIP6.GMLSAT_NI_tas_land.trends[,trend.ix], x_new = OBS.GMLSAT_tas_land.trends[trend.ix,], plot.constraint = F)
    plotCI(x = 2.3, y = land_constraint3$mean_out, uiw = qnorm(p = 0.975) * land_constraint3$sd_out, 
           liw = qnorm(p = 0.975) * land_constraint3$sd_out, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
    # Cowtan corrections:
    plotCI(x = 2.4, y = OBS.GMSST_hybrid36cor.trends_q[2,trend.ix], li = OBS.GMSST_hybrid36cor.trends_q[1,trend.ix], ui = OBS.GMSST_hybrid36cor.trends_q[3,trend.ix], 
           err = "y", add = T, col = "darkgoldenrod4", pch = NA, lwd = 2, sfrac=0.01)
    
    
    text(x = 1.6, y = -0.49, labels = c("HadSST4 Reconstr."), srt=90, adj = 0, col = "blue", cex = 0.75)
    text(x = 1.8, y = -0.49, labels = c("Ocean2k Constraint"), srt=90, adj = 0, col = "darkorchid4", cex = 0.75)
    text(x = 2.05, y = -0.49, labels = c("Pages2k GMST Constraint"), srt=90, adj = 0, col = "grey", cex = 0.75)
    #text(x = 2.25, y = -0.49, labels = c("Crowley Constraint"), srt=90, adj = 0, col = "blue", cex = 0.75)
    text(x = 2.3, y = -0.49, labels = c("CRUTEM5 Constraint"), srt=90, adj = 0, col = "darkorange", cex = 0.75)
    text(x = 2.4, y = -0.49, labels = c("Coastal Hybrid SST"), srt=90, adj = 0, col = "darkgoldenrod4", cex = 0.75)
    
  }
dev.off()
  

  
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





