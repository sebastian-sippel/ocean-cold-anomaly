
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023

library(matrixStats)
library(zoo)
library(dplR)


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")


## 01a. Scale observations into GMST:
## ----------------------------------------
wr = c(mean(CRUTEM5.global.annual$Anomaly[142:171]) - mean(CRUTEM5.global.annual$Anomaly[8:51])) / c(mean(HadSST4.global.annual$Anomaly[142:171]) - mean(HadSST4.global.annual$Anomaly[8:51]))
alpha = 0.67
scale.land = alpha / wr + 1 - alpha   # see calculation!!
scale.sst = alpha + (1-alpha) * wr

# GSAT_blended_ann_mod_p1 = colMedians(OBS.tos_$GSAT$ann$mod_p1_min) * sea_fraction + colMedians(OBS.tas_land_$GSAT$ann$mod_p1_min) * land_ice_fraction
ocean.raw.scale = HadSST4.global.annual$Anomaly[1:171] * scale.sst
land.raw.scale = CRUTEM5.global.annual$Anomaly[1:171] * scale.land

# regress CRUTEM5 and HadSST4:
raw.data.scale = get.diff.window.cor(x = ocean.raw.scale[8:171], y = land.raw.scale[8:171], years = 1857:2020, w.width = 51, W = 20, center = T)




# ------------------------------------------------------------------------------------
# 02. Plot cold anomaly plausibility in CMIP6
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/02_cmip6/")


pdf(file = "02_plausibility_CMIP6.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(3,5,1,5))
  ylim = c(-3.5, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  mtext(text = "Temperature difference, \n GMST(ocean) minus GMST(land) [°C]", side = 2, line = 3, at = 3.5, adj = 0.5)
  mtext(text = "Pearson Correlation, \n ocean/land (in 51-year window)", side = 2, line = 3, at = -1.5, adj = 0.5)
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 6, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 4.1, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.7 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    

    lines(x = c(1850, 2020), y = c(1, 1), col = "black", lty = 1)
    
    text(x = 1855, y = 0.7, labels = "d. Original Reconstruction, Correlation", col = "grey40", pos = 4)
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*1, labels=F, tcl=0.2)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*1, labels=F, tcl=0.2)
    text(x = 1855, y = -1.3, labels = "e. High-pass filtered (<20yr), Correlation", col = "grey40", pos = 4)
    
    # axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    # axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    # lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    # text(x = 1855, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x_y_2.5, rev(CMIP6_mod_p1_min$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x_y_2.5, rev(OBS_mod_p1_min$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x.low_y.low_2.5, rev(CMIP6_mod_p1_min$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x.low_y.low_2.5, rev(OBS_mod_p1_min$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.low_y.low * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x.high_y.high_2.5, rev(CMIP6_mod_p1_min$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x.high_y.high_2.5, rev(OBS_mod_p1_min$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  
  # Original Reconstruction, Correlation:
  {
    c = -2 + 0.5*1; sc = 2
    
    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p1_min_pt1$cor_x_y_2.5[26:140], rev(CMIP6_mod_p1_min_pt1$cor_x_y_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x_y_2.5[26:146], rev(OBS_mod_p1_min$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    # Original data:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x_y_2.5[26:146], rev(OBS_mod_p1_min$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered Reconstruction, Correlation:
  {
    # c = -4 + 0.5 - 1.25; sc = 2
    c = -4 + 0.5*1; sc = 2
    
    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p1_min_pt1$cor_x.high_y.high_2.5[26:140], rev(CMIP6_mod_p1_min_pt1$cor_x.high_y.high_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x.high_y.high_2.5[26:146], rev(OBS_mod_p1_min$cor_x.high_y.high_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  legend("bottom", c("Observational reconstruction", "Raw observations, scaled", "CMIP6 reconstruction"), 
         col = c("grey40", "black", "red"),  
         cex = 0.8, ncol = 2, inset = 0.02, lwd = 2, lty = c(1, 2, 1))
}
dev.off()



pdf(file = "02_plausibility_CMIP6_mod_p0.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(3,5,1,5))
  ylim = c(-3.5, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  mtext(text = "Temperature difference, \n GMST(ocean) minus GMST(land) [°C]", side = 2, line = 3, at = 3.5, adj = 0.5)
  mtext(text = "Pearson Correlation, \n ocean/land (in 51-year window)", side = 2, line = 3, at = -1.5, adj = 0.5)
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 6, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 4.1, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.7 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    
    lines(x = c(1850, 2020), y = c(1, 1), col = "black", lty = 1)
    
    text(x = 1855, y = 0.7, labels = "d. Original Reconstruction, Correlation", col = "grey40", pos = 4)
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*1, labels=F, tcl=0.2)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*1, labels=F, tcl=0.2)
    text(x = 1855, y = -1.3, labels = "e. High-pass filtered (<20yr), Correlation", col = "grey40", pos = 4)
    
    # axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    # axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    # lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    # text(x = 1855, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p0$diff_x_y_2.5, rev(CMIP6_mod_p0$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p0$diff_x_y_2.5, rev(OBS_mod_p0$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p0$diff_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p0$diff_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p0$diff_x.low_y.low_2.5, rev(CMIP6_mod_p0$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p0$diff_x.low_y.low_2.5, rev(OBS_mod_p0$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p0$diff_x.low_y.low_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.low_y.low * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p0$diff_x.low_y.low_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p0$diff_x.high_y.high_2.5, rev(CMIP6_mod_p0$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p0$diff_x.high_y.high_2.5, rev(OBS_mod_p0$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p0$diff_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p0$diff_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  
  # Original Reconstruction, Correlation:
  {
    c = -2 + 0.5*1; sc = 2
    
    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p0_pt1$cor_x_y_2.5[26:140], rev(CMIP6_mod_p0_pt1$cor_x_y_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p0$cor_x_y_2.5[26:146], rev(OBS_mod_p0$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    # Original data:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p0$cor_x_y_2.5[26:146], rev(OBS_mod_p0$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    lines(x = 1850:2014, y = CMIP6_mod_p0_pt1$cor_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p0$cor_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered Reconstruction, Correlation:
  {
    c = -4 + 0.5*1; sc = 2
    
    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p0_pt1$cor_x.high_y.high_2.5[26:140], rev(CMIP6_mod_p0_pt1$cor_x.high_y.high_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p0$cor_x.high_y.high_2.5[26:146], rev(OBS_mod_p0$cor_x.high_y.high_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p0_pt1$cor_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p0$cor_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  legend("bottom", c("Observational reconstruction", "Raw observations, scaled", "CMIP6 reconstruction"), 
         col = c("grey40", "black", "red"),  
         cex = 0.8, ncol = 2, inset = 0.02, lwd = 2, lty = c(1, 2, 1))
}
dev.off()




