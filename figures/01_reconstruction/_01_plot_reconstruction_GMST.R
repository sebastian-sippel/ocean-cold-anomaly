
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//code/_functions_CMIP6.R")

# Load global observations:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations.R")

# 00.(c) Load *new* reconstructions:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03b_load_global_obs_reconstruction.R")



# ------------------------------------------------------------------------------------
# Plot reconstruction(s):
# ------------------------------------------------------------------------------------


# 01_main_reconstruction: AGMT anomalies | the !main! reconstructions (no uncertainties/sensitivities!):
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/01_reconstruction/")

library(RColorBrewer)
col = brewer.pal(n = 8, name = "Dark2")


## GMST-annual reconstruction:
{
  
  # GMST reconstruction-mod_p0:
  pdf(file = "SI05_GMST_land_vs_ocean_reconstruction_mod_p0.pdf", width = 8, height=4)
  {
    par(mfrow=c(1, 1), mar=c(2,4,1,1))
    ylim = c(-1, 1); xlim = c(1850,2020)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Temperature Anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1)
    axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
    
    ## polygons for uncertainties:
    polygon(x = c(1850:2021, 2021:1850), y = c(HadCRUT5.global.annual$lower_CI2.5, 
                                               rev(HadCRUT5.global.annual$upper_CI97.5)), 
            col = make.transparent.color("wheat4", alpha = 50), border = make.transparent.color("wheat4", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tas_land$GMST_FM$ann$mod_p0_2.5, 
                                               rev(OBS.tas_land$GMST_FM$ann$mod_p0_97.5)), 
            col = make.transparent.color("darkorange", alpha = 50), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tos$GMST_FM$ann$mod_p0_2.5, 
                                               rev(OBS.tos$GMST_FM$ann$mod_p0_97.5)), 
            col = make.transparent.color("blue", alpha = 50), border = make.transparent.color("darkblue", 100))
    
    lines(x = CW2014.global.annual_had4sst4_krig$Year, y = CW2014.global.annual_had4sst4_krig$Anomaly, col = "grey40", lty = 2)
    lines(x = CW2014.global.annual_cobe2cru_krig$Year - 0.5, y = CW2014.global.annual_cobe2cru_krig$Anomaly, col = "grey40", lty = 3)
    
    lines(x = JMA.global.annual$Year, y = JMA.global.annual$Global - mean(JMA.global.annual$Global[match(x = 1961:1990, table = JMA.global.annual$Year)]), col = "black", lty = 3)
    lines(x = BEST.global.annual$Year, y = BEST.global.annual$Anomaly - mean(BEST.global.annual$Anomaly[match(x = 1961:1990, table = BEST.global.annual$Year)]), col = "black", lty = 2)
    lines(x = NOAA.global.annual$Year, y = NOAA.global.annual$Anomaly - mean(NOAA.global.annual$Anomaly[match(x = 1961:1990, table = NOAA.global.annual$Year)]), col = "black", lty = 4)
    lines(x = GISS.global.annual$Year, y = GISS.global.annual$Anomaly - mean(GISS.global.annual$Anomaly[match(x = 1961:1990, table = GISS.global.annual$Year)]), col = "black", lty = 5)
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = OBS.tas_land$GMST_FM$ann$mod_p0, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = OBS.tos$GMST_FM$ann$mod_p0, col = "darkblue", lwd = 2)
    lines(x = 1850:2021, y = HadCRUT5.global.annual$Anomaly, col = "grey40", lwd = 1)
    
    legend("topleft", c("HadCRUT5", "CRUTEM5-reconstr.", "HadSST4-reconstr.", "CW2014", "CW2014-COBE2", "JMA-GMST", "BEST", "NOAA-GlobalTemp-v5", "NASA-GISS"),
           lty = c(1, 1, 1, 2, 3, 3, 2, 4, 5), col = c("wheat4", "darkorange", "darkblue", "grey40", "grey40", "black", "black", "black", "black"), lwd = c(1, 2, 2, 1, 1, 1, 1, 1, 1), inset = 0.02, 
           cex = 0.75, ncol = 2, title = "GMST reconstruction (no unc./bias in training) and GMST datasets")
    dev.off()
  }
  
  # GMST reconstruction-mod_p1:
  pdf(file = "SI05_GMST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=4)
  {
    par(mfrow=c(1, 1), mar=c(2,4,1,1))
    ylim = c(-1, 1); xlim = c(1850,2020)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Temperature Anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1)
    axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
    
    ## polygons for uncertainties:
    polygon(x = c(1850:2021, 2021:1850), y = c(HadCRUT5.global.annual$lower_CI2.5, 
                                               rev(HadCRUT5.global.annual$upper_CI97.5)), 
            col = make.transparent.color("wheat4", alpha = 50), border = make.transparent.color("wheat4", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tas_land$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tas_land$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("darkorange", alpha = 50), border = make.transparent.color("darkorange", 100))

    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tos$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tos$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("blue", alpha = 50), border = make.transparent.color("darkblue", 100))

    lines(x = CW2014.global.annual_had4sst4_krig$Year, y = CW2014.global.annual_had4sst4_krig$Anomaly, col = "grey40", lty = 2)
    lines(x = CW2014.global.annual_cobe2cru_krig$Year - 0.5, y = CW2014.global.annual_cobe2cru_krig$Anomaly, col = "grey40", lty = 3)
    
    lines(x = JMA.global.annual$Year, y = JMA.global.annual$Global - mean(JMA.global.annual$Global[match(x = 1961:1990, table = JMA.global.annual$Year)]), col = "black", lty = 3)
    lines(x = BEST.global.annual$Year, y = BEST.global.annual$Anomaly - mean(BEST.global.annual$Anomaly[match(x = 1961:1990, table = BEST.global.annual$Year)]), col = "black", lty = 2)
    lines(x = NOAA.global.annual$Year, y = NOAA.global.annual$Anomaly - mean(NOAA.global.annual$Anomaly[match(x = 1961:1990, table = NOAA.global.annual$Year)]), col = "black", lty = 4)
    lines(x = GISS.global.annual$Year, y = GISS.global.annual$Anomaly - mean(GISS.global.annual$Anomaly[match(x = 1961:1990, table = GISS.global.annual$Year)]), col = "black", lty = 5)
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = OBS.tas_land$GMST_FM$ann$mod_p1_min, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = OBS.tos$GMST_FM$ann$mod_p1_min, col = "darkblue", lwd = 2)
    lines(x = 1850:2021, y = HadCRUT5.global.annual$Anomaly, col = "grey40", lwd = 1)
    
    
    legend("topleft", c("HadCRUT5", "CRUTEM5-reconstr.", "HadSST4-reconstr.", "CW2014", "CW2014-COBE2", "JMA-GMST", "BEST", "NOAA-GlobalTemp-v5", "NASA-GISS"),
           lty = c(1, 1, 1, 2, 3, 3, 2, 4, 5), col = c("wheat4", "darkorange", "darkblue", "grey40", "grey40", "black", "black", "black", "black"), lwd = c(1, 2, 2, 1, 1, 1, 1, 1, 1), inset = 0.02, 
           cex = 0.75, ncol = 2, title = "GMST reconstruction (train unc./bias incl.) and GMST datasets")
    dev.off()
  }
  
  
  # GMST reconstruction-mod_p1-simplified:
  pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1_simple.pdf", width = 8, height=4)
  {
    par(mfrow=c(1, 1), mar=c(2,4,1,1))
    ylim = c(-1, 1); xlim = c(1850,2020)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Temperature Anomaly [°C]", xlab = "", main = "Annual Global Mean Surface Temperature", ylim = ylim, xlim = xlim, las=1)
    axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
    
    ## polygons for uncertainties:
    polygon(x = c(1850:2021, 2021:1850), y = c(HadCRUT5.global.annual$lower_CI2.5, 
                                               rev(HadCRUT5.global.annual$upper_CI97.5)), 
            col = make.transparent.color("wheat4", alpha = 50), border = make.transparent.color("wheat4", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tas_land$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tas_land$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("darkorange", alpha = 50), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tos$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tos$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("blue", alpha = 50), border = make.transparent.color("darkblue", 100))
    
    # lines(x = CW2014.global.annual_had4sst4_krig$Year, y = CW2014.global.annual_had4sst4_krig$Anomaly, col = "grey40", lty = 2)
    # lines(x = CW2014.global.annual_cobe2cru_krig$Year - 0.5, y = CW2014.global.annual_cobe2cru_krig$Anomaly, col = "grey40", lty = 3)
    
    # lines(x = JMA.global.annual$Year, y = JMA.global.annual$Global - mean(JMA.global.annual$Global[match(x = 1961:1990, table = JMA.global.annual$Year)]), col = "black", lty = 3)
    # lines(x = BEST.global.annual$Year, y = BEST.global.annual$Anomaly - mean(BEST.global.annual$Anomaly[match(x = 1961:1990, table = BEST.global.annual$Year)]), col = "black", lty = 2)
    # lines(x = NOAA.global.annual$Year, y = NOAA.global.annual$Anomaly - mean(NOAA.global.annual$Anomaly[match(x = 1961:1990, table = NOAA.global.annual$Year)]), col = "black", lty = 4)
    # lines(x = GISS.global.annual$Year, y = GISS.global.annual$Anomaly - mean(GISS.global.annual$Anomaly[match(x = 1961:1990, table = GISS.global.annual$Year)]), col = "black", lty = 5)
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = OBS.tas_land$GMST_FM$ann$mod_p1_min, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = OBS.tos$GMST_FM$ann$mod_p1_min, col = "darkblue", lwd = 2)
    lines(x = 1850:2021, y = HadCRUT5.global.annual$Anomaly, col = "grey40", lwd = 2)
    
    
    legend("topleft", c("HadCRUT5", "CRUTEM5-reconstr.", "HadSST4-reconstr."),
           lty = c(1, 1, 1), col = c("wheat4", "darkorange", "darkblue"), lwd = c(2, 2, 2), inset = 0.02, cex = 1)
    dev.off()
  }
}






## CRUTEM5 vs. HadSST4:

pdf(file = "SI06_CRUTEM5_HadSST4.pdf", width = 8, height=4)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,1))
  ylim = c(-1.2, 1.2); xlim = c(1850,2020)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Temperature Anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1)
  axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  
  ## polygons for uncertainties:
  polygon(x = c(1850:2021, 2021:1850), y = c(CRUTEM5.global.annual$lower_CI2.5, 
                                             rev(CRUTEM5.global.annual$upper_CI97.5)), 
          col = make.transparent.color("wheat4", alpha = 50), border = make.transparent.color("wheat4", 100))
  
  polygon(x = c(1850:2021, 2021:1850), y = c(HadSST4.global.annual$lower_CI2.5, 
                                             rev(HadSST4.global.annual$upper_CI97.5)), 
          col = make.transparent.color("blue", alpha = 50), border = make.transparent.color("blue", 100))
  
  # Lines for main reconstructions:
  lines(x = 1850:2021, y = CRUTEM5.global.annual$Anomaly, col = "grey40", lwd = 2)
  lines(x = 1850:2021, y = HadSST4.global.annual$Anomaly, col = "blue", lwd = 2)
  
  lines(x = 1970:2020, y = lm(c(CRUTEM5.global.annual$Anomaly[121:171]) ~ c(1970:2020))$fitted, col = "grey40", lwd = 2)
  lines(x = 1970:2020, y = lm(c(HadSST4.global.annual$Anomaly[121:171]) ~ c(1970:2020))$fitted, col = "blue", lwd = 2)

  lines(x = 1901:1939, y = lm(c(CRUTEM5.global.annual$Anomaly[52:90]) ~ c(1901:1939))$fitted, col = "grey40", lwd = 2)
  lines(x = 1901:1939, y = lm(c(HadSST4.global.annual$Anomaly[52:90]) ~ c(1901:1939))$fitted, col = "blue", lwd = 2)
  
  legend("topleft", c("CRUTEM5", "HadSST4"),
         lty = c(1, 1), col = c("wheat4", "blue"), lwd = c(2, 2), inset = 0.02, 
         cex = 0.75, ncol = 2, title = "Original Data")
  dev.off()
}


# Ratio's of early and late warming trend slopes:
# ------------------------------------------------
lm(c(CRUTEM5.global.annual$Anomaly[121:171]) ~ c(1970:2020))$coefficients[2] /
lm(c(HadSST4.global.annual$Anomaly[121:171]) ~ c(1970:2020))$coefficients[2]

lm(c(CRUTEM5.global.annual$Anomaly[52:90]) ~ c(1901:1939))$coefficients[2] /
lm(c(HadSST4.global.annual$Anomaly[52:90]) ~ c(1901:1939))$coefficients[2]






