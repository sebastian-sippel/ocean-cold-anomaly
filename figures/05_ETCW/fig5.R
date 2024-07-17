
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023

## load all data for reconstructions:
source("scripts/04a_master_load_reconstructions.R")
source("scripts/04b_master_read_paleo_reconstructions.R")
source("scripts/04c_compute_trends_4paleo-comparison.R")




# ------------------------------------------------------------------------------------
# 02. Land-ocean warming ratio during ETCW:
# ------------------------------------------------------------------------------------

# setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/05_ETCW/")


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

# old name: 05_land_ocean_40_years.pdf
pdf(file = "figures/05_ETCW/fig5_fl.pdf", width = 8, height = 5)
    par(mar = c(4,4,1,0), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
    layout(matrix(1:2, ncol = 2, byrow = T), widths = c(9,4), heights = c(4,9))
    
    # first part of plot:
    {  
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         xlab = "Land warming [°C per 40 years]", ylab = "Ocean warming [°C per 40 years]", main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxs="i")

    # Add label (a) at the top left outside of the plot
    mtext("a", side = 3, line = 0, at = -0.9, adj = 0, cex = 1.1, font = 2)
    mtext("b", side = 3, line = 0, at = 2.05, adj = 0, cex = 1.1, font = 2)
    mtext("Constraints on \n ocean warming", side = 1, line = 3, at = 2.3, adj = 0, cex = 1, font = 1)    
      
    axis(side = 2, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    axis(side = 1, at = seq(-1, 2, 0.05), tcl=0.2, labels=F)
    axis(side=4, tick = T, labels = F)
    axis(side = 4, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    
    # (1) CMIP6 models provide strong constraint between land and ocean:
    points(y = CMIP6.GMSST_true.trends[,4], x = CMIP6.GMLSAT_true.trends[,4], col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMLSAT_true.trends[,3], col = make.transparent.color(red.cols[7], alpha = 80), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMLSAT_true.trends[,2], col = make.transparent.color(red.cols[5], alpha = 100), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,1], x = CMIP6.GMLSAT_true.trends[,1], col = make.transparent.color(red.cols[3], alpha = 100), pch = 16, cex = 0.6)
    
    # show CMIP6 models in terms of Ellipses:
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
    
    # get different constraints:
    CRUTEM5_constraint = lapply(X = 1:4, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMLSAT_true.trends[,i], x = CMIP6.GMLSAT_NI_tas_land.trends[,i], 
                                                             x_new = OBS.GMLSAT_tas_land.trends[i,], plot.constraint = F))
    HadSST4_constraint = lapply(X = 1:4, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,i], x = CMIP6.GMSST_tos.trends[,i], 
                                                             x_new = OBS.GMSST_tos.trends[i,], plot.constraint = F))
    hybridSST_constraint = lapply(X = 1:4, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,i], x = CMIP6.GMSST_tos.trends[,i], 
                                                             x_new = OBS.GMSST_hybrid36.trends[i,], plot.constraint = F))

    for (i in 1:4) {
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = HadSST4_constraint[[i]]$mean_out, 
             li = CRUTEM5_constraint[[i]]$q025, ui = CRUTEM5_constraint[[i]]$q975, 
             err = "x", add = T, col = "grey40", pch = NA, lwd = 2, lty = 2, sfrac=0.005)      
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = HadSST4_constraint[[i]]$mean_out, 
             li = HadSST4_constraint[[i]]$q025, ui = HadSST4_constraint[[i]]$q975, 
             err = "y", add = T, col = "grey40", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
      
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = hybridSST_constraint[[i]]$mean_out, li = CRUTEM5_constraint[[i]]$q025, ui = CRUTEM5_constraint[[i]]$q975, 
             err = "x", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.005)
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out + 0.01, y = hybridSST_constraint[[i]]$mean_out, li = hybridSST_constraint[[i]]$q025, ui = hybridSST_constraint[[i]]$q975, 
             err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.005)
    }
      lines(x = sapply(CRUTEM5_constraint, FUN=function(x) x$mean_out), y = sapply(HadSST4_constraint, FUN=function(x) x$mean_out), col = "grey40", lty = 5, lwd = 1.5)
      lines(x = sapply(CRUTEM5_constraint, FUN=function(x) x$mean_out), y = sapply(hybridSST_constraint, FUN=function(x) x$mean_out), col = "brown4", lty = 5, lwd = 1.5)
    
    # get further estimates original CRUTEM5+HadSST4 Averages:
    points(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], 
           y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], pch = 8, col = "black")
    lines(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], 
          y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[1:4], col = "black", lty = 5, lwd = 1.5)
    
    text(x = CRUTEM5_constraint[[1]]$mean_out, y = HadSST4_constraint[[1]]$q025 - 0.04, labels = c("1851-1890"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[2]]$mean_out, y = HadSST4_constraint[[2]]$q025 - 0.03, labels = c("1871-1910"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[3]]$mean_out, y = HadSST4_constraint[[3]]$q975 + 0.04, labels = c("1901-1940"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[4]]$mean_out + 0.06, y = HadSST4_constraint[[4]]$q025 - 0.08, labels = c("1975-2014"), 
         adj = 0.5, col = "black", cex = 0.8)
    
    legend("bottomright", c("Original CRUTEM5 / HadSST4",
                            #"CRUTEM5 / HadSST4", "Reconstruction constraint",
                            #"CRUTEM5 / CoastalHybridSST", "Reconstruction Constraint",
                            expression(hat(T)[CRUTEM5]^GMLSAT * " / " * hat(T)[HadSST4]^GMSST * " constraint"),
                            expression(hat(T)[CRUTEM5]^GMLSAT * " / " * hat(T)[CoastalHybridSST]^GMSST * " constraint"),
                            "CMIP6, 1851-1890", "CMIP6, 1871-1910",
                            "CMIP6, 1901-1940", "CMIP6, 1975-2014"), col = c("black", "grey40", NA, "brown4", NA,
                                                         make.transparent.color(red.cols[3], alpha = 200), make.transparent.color(red.cols[5], alpha = 100),
                                                         make.transparent.color(red.cols[7], alpha = 100), make.transparent.color(red.cols[9], alpha = 100)), 
           pch = c(8, 3, NA, 3, NA, 16, 16, 16, 16), pt.cex = c(1, 1.5, NA, 1.5, NA, 1, 1, 1, 1),
           inset = 0.02, cex = 0.65)
    }
    
    # second part of plot (constraints):
    {
    par(mar = c(4,0,1,0))
    plot(x = c(0.5, 2.6), y = c(1,1), type="n", 
         xlab = "", ylab = "", main = "", ylim = ylim, xlim = c(0.5, 2.6), las=1, yaxs="i", xaxs="i", xaxt= "n", yaxt="n", bty = "n")
    
    text(x = 1, y = 0.65, labels = c("1871-1910"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    text(x = 2, y = 0.65, labels = c("1901-1940"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    lines(x = c(1.5,1.5), y = c(-1, 3), lty = 3)
    legend("top", c("HadSST4", "ERSSTv5", "COBE-SST2"), title = "Original SST datasets", col = "blue", pch = c(8,6,4), cex = 0.7, bg = "white", inset = 0.00, ncol = 2)
    ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)
    
    # 1876-1925
    # Recon:
    plotCI(x = 0.6, y = HadSST4_constraint[[2]]$mean_out, li = HadSST4_constraint[[2]]$q025, ui = HadSST4_constraint[[2]]$q975, 
           err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
    
    # Other Ocean Datasets:
    points(x = 0.7, y = HadSST4.trends[2], col = "blue", pch = 8)
    points(x = 0.7, y = ERSSTv5.trends[2], col = "blue", pch = 6)
    points(x = 0.7, y = COBE_SST2.trends[2], col = "blue", pch = 4)
    
    # Ocean2k-based constraint:
    ocean2k_constraint2 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.Tropics_true.trends[,2], x_new = ocean2k_trends[2], plot.constraint = F)
    plotCI(x = 0.8, y = ocean2k_constraint2$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, 
           liw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    
    # Check points for plotting:
    #points(x = 0.55, y = -0.4)
    #points(x = 0.55, y = 0.8)
    
    # Neukom Paleo-constraint:
    for (i in 1:7) {
      neukom_constraint = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMST_true.trends[,2], x_new = neukom2019_trend_all[[i]][2,], plot.constraint = F)
      plotCI(x = 0.85+i/20, y = neukom_constraint$mean_out, ui = neukom_constraint$q975, 
             li = neukom_constraint$q025, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }

    # CRUTEM5 based constraint:
    land_constraint2 = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,2], x = CMIP6.GMLSAT_NI_tas_land.trends[,2], x_new = OBS.GMLSAT_tas_land.trends[2,], plot.constraint = F)
    plotCI(x = 1.3, y = land_constraint2$mean_out, ui = land_constraint2$q975, 
           li = land_constraint2$q025, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)

    # Cowtan corrections:
    plotCI(x = 1.4, y = hybridSST_constraint[[2]]$mean_out, li = hybridSST_constraint[[2]]$q025, ui = hybridSST_constraint[[2]]$q975, 
           err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.01)
    
    
    # 1901-1940
    # Recon:
    plotCI(x = 1.6, y = HadSST4_constraint[[3]]$mean_out, li = HadSST4_constraint[[3]]$q025, ui = HadSST4_constraint[[3]]$q975, 
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
      neukom_constraint = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMST_true.trends[,3], x_new = neukom2019_trend_all[[i]][3,], plot.constraint = F)
      plotCI(x = 1.85+i/20, y = neukom_constraint$mean_out, ui = neukom_constraint$q975, 
             li = neukom_constraint$q025, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }
    
    # CRUTEM5 based constraint:
    land_constraint3 = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,3], x = CMIP6.GMLSAT_NI_tas_land.trends[,3], x_new = OBS.GMLSAT_tas_land.trends[3,], plot.constraint = F)
    plotCI(x = 2.3, y = land_constraint3$mean_out, ui = land_constraint3$q975, 
           li = land_constraint3$q025, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
    # Cowtan corrections:
    plotCI(x = 2.4, y = hybridSST_constraint[[3]]$mean_out, li = hybridSST_constraint[[3]]$q025, ui = hybridSST_constraint[[3]]$q975, 
           err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.01)
    
    text(x = 1.6, y = -0.49, labels = expression(hat(T)[HadSST4] * ""), srt=90, adj = 0, col = "blue", cex = 0.75)
    # text(x = 1.6, y = -0.49, labels = "HadSST4 Reconstruction", srt=90, adj = 0, col = "blue", cex = 0.75)
    text(x = 1.8, y = -0.49, labels = c("Ocean2k ('Best Recons.')"), srt=90, adj = 0, col = "darkorchid4", cex = 0.75)
    text(x = 2.05, y = -0.49, labels = c("Pages2k GMST Recons."), srt=90, adj = 0, col = "grey", cex = 0.75)
    # text(x = 2.28, y = -0.49, labels = "CRUTEM5 Reconstruction", srt=90, adj = 0, col = "darkorange", cex = 0.72)
    # text(x = 2.42, y = -0.49, labels = "CoastalHybridSST Reconstruction", srt=90, adj = 0, col = "brown4", cex = 0.75)
    text(x = 2.26, y = -0.49, labels = expression(hat(T)[CRUTEM5] * ""), srt=90, adj = 0, col = "darkorange", cex = 0.75)
    text(x = 2.47, y = -0.49, labels = expression(hat(T)[CoastalHybridSST] * ""), srt=90, adj = 0, col = "brown4", cex = 0.75)
}
dev.off()
  
  


  
  
# ------------------------------------------------------------------------------------
# 50-years figure
# ------------------------------------------------------------------------------------
  
  
# Supplement Figure ??  
  pdf(file = "figures/05_ETCW/05_land_ocean_50_years_fl.pdf", width = 8, height = 5)
  par(mar = c(4,4,1,0), cex.lab = 1, cex.axis=1, mfrow=c(1,1))
  layout(matrix(1:2, ncol = 2, byrow = T), widths = c(9,4), heights = c(4,9))
  
  # first part of plot:
  {  
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         xlab = "Land warming [°C per 50 years]", ylab = "Ocean warming [°C per 50 years]", main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxs="i")
    
    # Add label (a) at the top left outside of the plot
    mtext("a", side = 3, line = 0, at = -0.9, adj = 0, cex = 1.1, font = 2)
    mtext("b", side = 3, line = 0, at = 2.05, adj = 0, cex = 1.1, font = 2)
    mtext("Constraints on \n ocean warming", side = 1, line = 3, at = 2.3, adj = 0, cex = 1, font = 1)    
    
    axis(side = 2, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    axis(side = 1, at = seq(-1, 2, 0.05), tcl=0.2, labels=F)
    axis(side=4, tick = T, labels = F)
    axis(side = 4, at = seq(-1, 1.5, 0.05), tcl=0.2, labels=F)
    
    # (1) CMIP6 models provide strong constraint between land and ocean:
    points(y = CMIP6.GMSST_true.trends[,8], x = CMIP6.GMLSAT_true.trends[,8], col = make.transparent.color(red.cols[9], alpha = 40), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,7], x = CMIP6.GMLSAT_true.trends[,7], col = make.transparent.color(red.cols[7], alpha = 80), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMLSAT_true.trends[,6], col = make.transparent.color(red.cols[5], alpha = 100), pch = 16, cex = 0.6)
    points(y = CMIP6.GMSST_true.trends[,5], x = CMIP6.GMLSAT_true.trends[,5], col = make.transparent.color(red.cols[3], alpha = 100), pch = 16, cex = 0.6)
    
    # show CMIP6 models in terms of Ellipses:
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
    
    # get different constraints:
    CRUTEM5_constraint = lapply(X = 5:8, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMLSAT_true.trends[,i], x = CMIP6.GMLSAT_NI_tas_land.trends[,i], 
                                                                                             x_new = OBS.GMLSAT_tas_land.trends[i,], plot.constraint = F))
    HadSST4_constraint = lapply(X = 5:8, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,i], x = CMIP6.GMSST_tos.trends[,i], 
                                                                                             x_new = OBS.GMSST_tos.trends[i,], plot.constraint = F))
    hybridSST_constraint = lapply(X = 5:8, FUN=function(i) get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,i], x = CMIP6.GMSST_tos.trends[,i], 
                                                                                               x_new = OBS.GMSST_hybrid36.trends[i,], plot.constraint = F))
    
    for (i in 1:4) {
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = HadSST4_constraint[[i]]$mean_out, 
             li = CRUTEM5_constraint[[i]]$q025, ui = CRUTEM5_constraint[[i]]$q975, 
             err = "x", add = T, col = "grey40", pch = NA, lwd = 2, lty = 2, sfrac=0.005)      
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = HadSST4_constraint[[i]]$mean_out, 
             li = HadSST4_constraint[[i]]$q025, ui = HadSST4_constraint[[i]]$q975, 
             err = "y", add = T, col = "grey40", pch = NA, lwd = 2, lty = 2, sfrac=0.005)
      
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out, y = hybridSST_constraint[[i]]$mean_out, li = CRUTEM5_constraint[[i]]$q025, ui = CRUTEM5_constraint[[i]]$q975, 
             err = "x", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.005)
      plotCI(x = CRUTEM5_constraint[[i]]$mean_out + 0.01, y = hybridSST_constraint[[i]]$mean_out, li = hybridSST_constraint[[i]]$q025, ui = hybridSST_constraint[[i]]$q975, 
             err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.005)
    }
    lines(x = sapply(CRUTEM5_constraint, FUN=function(x) x$mean_out), y = sapply(HadSST4_constraint, FUN=function(x) x$mean_out), col = "grey40", lty = 5, lwd = 1.5)
    lines(x = sapply(CRUTEM5_constraint, FUN=function(x) x$mean_out), y = sapply(hybridSST_constraint, FUN=function(x) x$mean_out), col = "brown4", lty = 5, lwd = 1.5)
    
    # get further estimates CRUTEM5+HadSST4 Raw Averages:
    points(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], 
           y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], pch = 8, col = "black")
    lines(x = get.trend(x = CRUTEM5.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], 
          y = get.trend(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)[5:8], col = "black", lty = 5, lwd = 1.5)
    
    text(x = CRUTEM5_constraint[[1]]$mean_out, y = HadSST4_constraint[[1]]$q025 - 0.04, labels = c("1851-1900"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[2]]$mean_out, y = HadSST4_constraint[[2]]$q025 - 0.03, labels = c("1871-1920"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[3]]$mean_out, y = HadSST4_constraint[[3]]$q975 + 0.08, labels = c("1901-1950"), 
         adj = 0.5, col = "black", cex = 0.8)
    text(x = CRUTEM5_constraint[[4]]$mean_out + 0.06, y = HadSST4_constraint[[4]]$q025 - 0.08, labels = c("1965-2014"), 
         adj = 0.5, col = "black", cex = 0.8)
    
    
    legend("bottomright", c("Original CRUTEM5 / HadSST4",
                            #"CRUTEM5 / HadSST4", "Reconstruction constraint",
                            #"CRUTEM5 / CoastalHybridSST", "Reconstruction Constraint",
                            expression(hat(T)[CRUTEM5]^GMLSAT * " / " * hat(T)[HadSST4]^GMSST * " constraint"),
                            expression(hat(T)[CRUTEM5]^GMLSAT * " / " * hat(T)[CoastalHybridSST]^GMSST * " constraint"),
                            "CMIP6, 1851-1890", "CMIP6, 1871-1910",
                            "CMIP6, 1901-1940", "CMIP6, 1975-2014"), col = c("black", "grey40", NA, "brown4", NA,
                                                                             make.transparent.color(red.cols[3], alpha = 200), make.transparent.color(red.cols[5], alpha = 100),
                                                                             make.transparent.color(red.cols[7], alpha = 100), make.transparent.color(red.cols[9], alpha = 100)), 
           pch = c(8, 3, NA, 3, NA, 16, 16, 16, 16), pt.cex = c(1, 1.5, NA, 1.5, NA, 1, 1, 1, 1),
           inset = 0.02, cex = 0.65)
  }
  
  # second part of plot (constraints):
  {
    par(mar = c(4,0,1,0))
    plot(x = c(0.5, 2.6), y = c(1,1), type="n", 
         xlab = "", ylab = "", main = "", ylim = ylim, xlim = c(0.5, 2.6), las=1, yaxs="i", xaxs="i", xaxt= "n", yaxt="n", bty = "n")

    text(x = 1, y = 0.65, labels = c("1871-1920"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    text(x = 2, y = 0.65, labels = c("1901-1950"), srt=45, adj = 0.5, col = "black", cex = 0.8)
    lines(x = c(1.5,1.5), y = c(-1, 3), lty = 3)
    legend("top", c("HadSST4", "ERSSTv5", "COBE-SST2"), title = "Original SST datasets", col = "blue", pch = c(8,6,4), cex = 0.7, bg = "white", inset = 0.00, ncol = 2)
    ocean2k_trends = get.trend(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = ocean2k_$Year)
    
    # 1876-1925
    # Recon:
    plotCI(x = 0.6, y = HadSST4_constraint[[2]]$mean_out, li = HadSST4_constraint[[2]]$q025, ui = HadSST4_constraint[[2]]$q975, 
           err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
    
    # Other Ocean Datasets:
    points(x = 0.7, y = HadSST4.trends[6], col = "blue", pch = 8)
    points(x = 0.7, y = ERSSTv5.trends[6], col = "blue", pch = 6)
    points(x = 0.7, y = COBE_SST2.trends[6], col = "blue", pch = 4)
    
    # Ocean2k-based constraint:
    ocean2k_constraint2 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.Tropics_true.trends[,6], x_new = ocean2k_trends[6], plot.constraint = F)
    plotCI(x = 0.8, y = ocean2k_constraint2$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, 
           liw = qnorm(p = 0.975) * ocean2k_constraint2$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    
    # Check points for plotting:
    #points(x = 0.55, y = -0.4)
    #points(x = 0.55, y = 0.8)
    
    # Neukom Paleo-constraint:
    for (i in 1:7) {
      neukom_constraint = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMST_true.trends[,6], x_new = neukom2019_trend_all[[i]][6,], plot.constraint = F)
      plotCI(x = 0.85+i/20, y = neukom_constraint$mean_out, ui = neukom_constraint$q975, 
             li = neukom_constraint$q025, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }
    
    # CRUTEM5 based constraint:
    land_constraint2 = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,6], x = CMIP6.GMLSAT_NI_tas_land.trends[,6], x_new = OBS.GMLSAT_tas_land.trends[6,], plot.constraint = F)
    plotCI(x = 1.3, y = land_constraint2$mean_out, ui = land_constraint2$q975, 
           li = land_constraint2$q025, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
    
    # Cowtan corrections:
    plotCI(x = 1.4, y = hybridSST_constraint[[2]]$mean_out, li = hybridSST_constraint[[2]]$q025, ui = hybridSST_constraint[[2]]$q975, 
           err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.01)
    
    # 1901-1950
    # Recon:
    plotCI(x = 1.6, y = HadSST4_constraint[[3]]$mean_out, li = HadSST4_constraint[[3]]$q025, ui = HadSST4_constraint[[3]]$q975, 
           err = "y", add = T, col = "blue", pch = NA, lwd = 2, lty = 2, sfrac=0.01)
    
    # Other Ocean Datasets:
    points(x = 1.7, y = HadSST4.trends[7], col = "blue", pch = 8)
    points(x = 1.7, y = ERSSTv5.trends[7], col = "blue", pch = 6)
    points(x = 1.7, y = COBE_SST2.trends[7], col = "blue", pch = 4)
    
    # Ocean2k-based constraint:
    ocean2k_constraint3 = get.linear.model.constraint(y = CMIP6.GMSST_true.trends[,7], x = CMIP6.Tropics_true.trends[,7], x_new = ocean2k_trends[7], plot.constraint = F)
    plotCI(x = 1.8, y = ocean2k_constraint3$mean_out, uiw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, 
           liw = qnorm(p = 0.975) * ocean2k_constraint3$sd_out, col = "darkorchid4", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    # Neukom Paleo-constraint:
    for (i in 1:7) {
      neukom_constraint = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,7], x = CMIP6.GMST_true.trends[,7], x_new = neukom2019_trend_all[[i]][7,], plot.constraint = F)
      plotCI(x = 1.85+i/20, y = neukom_constraint$mean_out, ui = neukom_constraint$q975, 
             li = neukom_constraint$q025, col = "grey", add = T, pch = NA, cex = 2, err = "y", lwd = 2, sfrac=0.01)
    }
    
    # CRUTEM5 based constraint:
    land_constraint3 = get.linear.model.constraint_ens_emp(y = CMIP6.GMSST_true.trends[,7], x = CMIP6.GMLSAT_NI_tas_land.trends[,7], x_new = OBS.GMLSAT_tas_land.trends[7,], plot.constraint = F)
    plotCI(x = 2.3, y = land_constraint3$mean_out, ui = land_constraint3$q975, 
           li = land_constraint3$q025, col = "darkorange", add = T, pch = NA, lwd = 2, cex = 2, err="y", sfrac=0.01)
    # Cowtan corrections:
    plotCI(x = 2.4, y = hybridSST_constraint[[3]]$mean_out, li = hybridSST_constraint[[3]]$q025, ui = hybridSST_constraint[[3]]$q975, 
           err = "y", add = T, col = "brown4", pch = NA, lwd = 2, sfrac=0.01)
    
    text(x = 1.6, y = -0.49, labels = expression(hat(T)[HadSST4] * ""), srt=90, adj = 0, col = "blue", cex = 0.75)
    # text(x = 1.6, y = -0.49, labels = "HadSST4 Reconstruction", srt=90, adj = 0, col = "blue", cex = 0.75)
    text(x = 1.8, y = -0.49, labels = c("Ocean2k ('Best Recons.')"), srt=90, adj = 0, col = "darkorchid4", cex = 0.75)
    text(x = 2.05, y = -0.49, labels = c("Pages2k GMST Recons."), srt=90, adj = 0, col = "grey", cex = 0.75)
    # text(x = 2.28, y = -0.49, labels = "CRUTEM5 Reconstruction", srt=90, adj = 0, col = "darkorange", cex = 0.72)
    # text(x = 2.42, y = -0.49, labels = "CoastalHybridSST Reconstruction", srt=90, adj = 0, col = "brown4", cex = 0.75)
    text(x = 2.26, y = -0.49, labels = expression(hat(T)[CRUTEM5] * ""), srt=90, adj = 0, col = "darkorange", cex = 0.75)
    text(x = 2.47, y = -0.49, labels = expression(hat(T)[CoastalHybridSST] * ""), srt=90, adj = 0, col = "brown4", cex = 0.75)
  }
  dev.off()
  

  
### OK!




