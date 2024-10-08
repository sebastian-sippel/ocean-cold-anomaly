
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.06.2024

library(matrixStats)
library(zoo)
library(dplR)


## load all data for reconstructions:
source("scripts/04a_master_load_reconstructions.R")


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

pdf(file = "figures/02_cmip6/fig2.pdf", width = 3.5, height=5.2, pointsize = 5)
{
  par(mfrow=c(1, 1), mar=c(3,6,1,5))
  ylim = c(-3.7, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  # mtext(text = "Temperature difference \n GMST(ocean) minus GMST(land) [°C]", side = 2, line = 3, at = 3.5, adj = 0.5)
  #mtext(text = "Temperature difference between ocean- and land \n based GMST reconstruction [°C], ", side = 2, line = 3, at = 3.5, adj = 0.5)
  #mtext(paste("Temperature difference between ocean- and land- \n based GMST reconstruction [°C]"), side = 2, line = 3, at = 3.5, adj = 0.5)
  
  mtext("Temperature difference between ocean- and land- ", side = 2, line = 5, at = 3.5, adj = 0.5)
  mtext(expression("based GMST reconstruction ["*hat(T)[Ocean]-hat(T)[Land]*", in °C]"),
        side = 2, line = 3.6, at = 3.5, adj = 0.5)
  
  mtext("Pearson Correlation between ocean- and land-", side = 2, line = 5, at = -1.5, adj = 0.5)
  mtext(expression("based GMST reconstruction ["*rho*"("*hat(T)[Ocean]*", "*hat(T)[Land]*")]"),
        side = 2, line = 3.6, at = -1.5, adj = 0.5)
  
  
  # mtext(text = "Pearson Correlation between ocean- \n and land- based GMST reconstruction", side = 2, line = 3, at = -1.5, adj = 0.5)
  
  
  {
    # Add the rectangle polygon to the plot
    polygon(x = c(1900, 1930, 1930, 1900), y = c(-5, -5, 6.3, 6.3), border = "darkgray", col = make.transparent.color("lightgrey", alpha = 60), lwd = 0.5)
    
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    # lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 1, lwd = 0.5)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 1, lwd = 0.5)
    text(x = 1850, y = 6, labels = "a. Original Reconstruction", col = "black", pos = 4, font = 2)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 1, lwd = 0.5)
    text(x = 1850, y = 4.1, labels = "b. Low-pass \n filtered (>20yr)", col = "black", pos = 4, font = 2)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 1, lwd = 0.5)
    text(x = 1850, y = 1.7 + 0.5*2, labels = "c. High-pass \n filtered (<20yr)", col = "black", pos = 4, font = 2)
    

    lines(x = c(1850, 2020), y = c(1, 1), col = "black", lty = 1, lwd = 0.5)
    
    text(x = 1855, y = 0.7, labels = "d. Original Reconstruction, Correlation", col = "black", pos = 4, font = 2)
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*1, labels=F, tcl=0.2)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*1, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*1, labels=F, tcl=0.2)
    text(x = 1855, y = -1.3, labels = "e. High-pass filtered (<20yr), Correlation", col = "black", pos = 4, font = 2)
    
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
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x_y_50 * sc + c, col = "red", lwd = 1.5)
    lines(x = 1857:2020, y = raw.data.scale$diff_x_y * sc + c, col = "black", lwd = 1.5, lty = 5)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x_y_50 * sc + c, col = "grey40", lwd = 1.5)
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
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "red", lwd = 1.5)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.low_y.low * sc + c, col = "black", lwd = 1.5, lty = 5)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "grey40", lwd = 1.5)
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
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "red", lwd = 1.5)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.high_y.high * sc + c, col = "black", lwd = 1.5, lty = 5)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "grey40", lwd = 1.5)
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
    
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x_y_50 * sc + c, col = "red", lwd = 1.5)
    lines(x = 1857:2020, y = raw.data.scale$cor_x_y * sc + c, col = "black", lwd = 1.5, lty = 5)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x_y_50 * sc + c, col = "grey40", lwd = 1.5)
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
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x.high_y.high_50 * sc + c, col = "red", lwd = 1.5)
    lines(x = 1857:2020, y = raw.data.scale$cor_x.high_y.high * sc + c, col = "black", lwd = 1.5, lty = 5)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x.high_y.high_50 * sc + c, col = "grey40", lwd = 1.5)
  }
  
  legend("bottomleft", c(expression("Observational GMST reconstructions ("*hat(T)[HadSST4]*" and "*hat(T)[CRUTEM5]*")"), 
                          "Original HadSST4 vs. original CRUTEM5 (scaled into GMST)", 
                          expression("CMIP6 GMST reconstructions ("*hat(T)[CMIP6-Ocean]*" and "*hat(T)[CMIP6-Land]*")")),
                          # "CMIP6 Ocean vs. Land  (GMST reconstructions)"), 
         col = c("grey40", "black", "red"),  
         cex = 1, ncol = 1, inset = 0.01, lwd = 1.5, lty = c(1, 5, 1), bg = "white", xpd = T)
}
dev.off()






