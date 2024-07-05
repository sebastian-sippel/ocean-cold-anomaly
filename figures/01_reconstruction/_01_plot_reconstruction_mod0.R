
# ------------------------------------------------------------------------------------
# Evaluate and plot reconstructions.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 03.01.2024


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")




# 08. Attribution and Filtering of time series for mod_p0:
# ------------------------------------------------------------------------------------
{
  library("dplR") # package for band-pass filtering
  
  # GMST Prepare data for plotting:
  GMST.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tas_land = get.df(Y = OBS.tas_land_$GMST_FM$ann$mod_p0, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
  GMST.tos = get.df(Y = OBS.tos_$GMST$ann$mod_p0, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
  # GMST.HadISST = get.df(Y = HadISST.annual$GMST_FM_mod_p0[1:151], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = HadISST.annual$Year[1:151], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMSST_mod_p0_1se, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMST_FM_mod_p0_1se[1:167], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.BEST_Land = get.df(Y = BEST_Land.annual$GMST_FM_mod_p0_1se[101:271], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMST_FM$ann$mod_p0), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMST_FM$ann$mod_p0), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tos_ua = get.df(Y = HadSST_ua.annual$GMST_FM_mod_p0_1se, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
}







## 01b. Get correction vector(s):
## ----------------------------------------
{
# get correction vector from Cowtan:
cor_vec_hybrid36 = GMST.hybrid36$mod_p1_min_50[1:167] - GMST.tos$mod_p1_min_50[1:167]
cor_vec_hadsst4 = GMST.tos$mod_p1_min_50 - GMST.tos_ua$mod_p1_min_50
}




# ------------------------------------------------------------------------------------
# Plot reconstruction(s):
# ------------------------------------------------------------------------------------


# 01_main_reconstruction: AGMT anomalies | the !main! reconstructions (no uncertainties/sensitivities!):
# ------------------------------------------------------------------------------------

library(RColorBrewer)
col = brewer.pal(n = 8, name = "Dark2")

library("dplR") # package for band-pass filtering





# Define colours for plotting: 
col.CRUTEM5 = "darkorange"
col.HadSST4 = "blue"
col.COBE2 = "slateblue4"
col.ERSST5 = "deepskyblue4"

col.BEST = "darkgoldenrod4"
col.hybrid36 = "brown4"
col.CLASSNMAT = "lightskyblue"
col.HadSST4_ua = "magenta4"



## ----------------------------------------------------------------------------------------
## 02b. Plot GMST filtered time series with hybrid correction vector:
## ----------------------------------------------------------------------------------------


# mod_p1:
pdf(file = "figures/01_reconstruction/fig1_mod_p0.pdf", width = 8, height=11)
{
  par(mfrow=c(1, 1), mar=c(2,5,1,4))
  ylim = c(-4.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    # Add the rectangle polygon to the plot
    polygon(x = c(1900, 1930, 1930, 1900), y = c(-5, -5, 6.3, 6.3), border = "darkgray", col = make.transparent.color("lightgrey", alpha = 60))
    
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 1)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 1)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 1)
    #lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 1)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 1)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 1)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 1)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 1)
    text(x = 1850, y = 5.5, labels = "a. Original Reconstruction", col = "black", pos = 4, font = 2)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 1)
    text(x = 1850, y = 3.3+0.5, labels = "b. Low-pass filtered \n (>20yr)", col = "black", pos = 4, font = 2)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 1)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "black", pos = 4, font = 2)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 1)
    text(x = 1850, y = -0.5 + 0.5*3, labels = "d. Forced response", col = "black", pos = 4, font = 2)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 1)
    text(x = 1850, y = -2.5 + 0.5*4, labels = "e. Unforced, low-pass filtered (>20yr)", col = "black", pos = 4, font = 2)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 1)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "f. Unforced, high-pass filtered (<20yr)", col = "black", pos = 4, font = 2)
    
    axis(side = 2, at = seq(-7.5, -6, 0.5) + 0.5*6, labels=seq(-0.5, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-7.5, -6, 0.1) + 0.5*6, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-7, -7) + 0.5*6, col = "darkgray", lty = 1)
    text(x = 1850, y = -6.1 + 0.5*6, labels = "g. Implied SST adjustments w.r.t. HadSST4-unadj.", col = "black", pos = 4, font = 2)
    
    mtext(text = "Global implied \n SST Adjustment [°C]", side = 2, line = 3, at = -3.8, col = "black", cex = 0.9)
  }
  
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$mod_p1_min_2.5, 
                                               rev(GMST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$mod_p1_min_2.5, 
                                               rev(GMST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1.5)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1.5)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$lp_2.5, 
                                               rev(GMST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$lp_2.5, 
                                               rev(GMST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1.5)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1.5)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tas_land$hp_2.5, 
                                               rev(GMST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMST.tos$hp_2.5, 
                                               rev(GMST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1.5)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1.5)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMST.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$forced * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tas_land$res_lp_2.5[1:165], 
                                               rev(GMST.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tos$res_lp_2.5[1:165], 
                                               rev(GMST.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1.5)
    lines(x = 1850:2020, y = GMST.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1.5)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tas_land$res_hp_2.5[1:165], 
                                               rev(GMST.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tos$res_hp_2.5[1:165], 
                                               rev(GMST.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1.5)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1.5)
    
    # Lines for other observational datasets:
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 5, lwd = 1.5)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 5, lwd = 1.5)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 5, lwd = 1.5)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1.5)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1.5)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1.5)
  }
  
  # SST Adjustment:
  {
    c = -7+0.5*6; sc = 1
    
    lines(x = 1850:2020, y = cor_vec_hadsst4 * sc + c, type='l', ylim = c(-1, 1), col = col.HadSST4, lwd = 1.5)
    #lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6) * sc + c, col = "blue", lty = 2)
    lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36) * sc + c, col = col.hybrid36, lwd = 1.5)
  }
  
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4-unadj", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "CoastalHybridSST"),
         lty = c(1, 1, 1, 5, 5, 1, 5, 1), lwd = 1.5, col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  
  legend("bottomright", c("HadSST4", "CoastalHybridSST"),
         col = c(col.HadSST4, col.hybrid36),  
         cex = 0.7, ncol = 2, bg = "white", inset = 0.01, lwd = 2)
  
}
dev.off()



