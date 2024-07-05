
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2023

## Load MR observations:
load("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/03_processedOBS_reconstr/OBS.tos_MR.RData")
OBS.tos_MR_ = OBS.tos_

## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")



# 09. Read-In of Kadow et al. AI tool for ocean->land and land->ocean extrapolation:
# ------------------------------------------------------------------------------------
# setwd("/net/h2o/climphys1/sippels/_DATASET/Kadow_land-vs-ocean/")
kadow_crutem5_gmst = ncvar_get(nc_open("data/00_DATASET/obs/Kadow-etal-2020/cmip6ctaspadzens-2_tas_mon-gl-72x36_crutem5.0.2_observation_all-avg_1850-2...FM_SUBCLIM"))
kadow_hadsst4_gmst = ncvar_get(nc_open("data/00_DATASET/obs/Kadow-etal-2020/cmip6ctaspadzens-2_tas_mon-gl-72x36_hadsst4_observation_all-avg_1850-2023_...FM_SUBCLIM"))




# ------------------------------------------------------------------------------------
# Plot reconstruction(s):
# ------------------------------------------------------------------------------------


# 01_main_reconstruction: AGMT anomalies | the !main! reconstructions (no uncertainties/sensitivities!):
# ------------------------------------------------------------------------------------
# setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/01_reconstruction/")

library(RColorBrewer)
col = brewer.pal(n = 8, name = "Dark2")




## GMST-annual reconstruction:
{
  
  # GMST reconstruction-mod_p0:
  pdf(file = "figures/01_reconstruction/SI01_GMST_land_vs_ocean_reconstruction.pdf", width = 8, height=8)
  {
    par(mfrow=c(2, 1), mar=c(2,4,1,1))
    ylim = c(-1, 1); xlim = c(1850,2020)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1)
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
    
    # Lines for Kadow reconstructions:
    lines(x = 1850:2023, y = kadow_crutem5_gmst, col = "gold1", type="l", lwd = 2)
    lines(x = 1850:2023, y = kadow_hadsst4_gmst, col = "lightblue", type="l", lwd = 2)
    
    legend("topleft", c("HadCRUT5", "CRUTEM5-reconstr.", "CRUTEM5-reconstr. (Kadow et al.)", "HadSST4-reconstr.", "HadSST4-reconstr. (Kadow et al.)", "CW2014", "CW2014-COBE2", "JMA-GMST", "BEST", "NOAA-GlobalTemp-v5", "NASA-GISS"),
           lty = c(1, 1, 1, 1, 1, 2, 3, 3, 2, 4, 5), col = c("wheat4", "darkorange", "gold1", "darkblue", "lightblue", "grey40", "grey40", "black", "black", "black", "black"), lwd = c(1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1), inset = 0.02, 
           cex = 0.7, ncol = 2, title = "GMST reconstruction (no unc./bias in training), Kadow et al. reconstruction and GMST datasets")

    
    # ------------------------------------------------------------
    # Now replot for mod_p1:
    # ------------------------------------------------------------
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
  
}



## CRUTEM5 vs. HadSST4:
pdf(file = "figures/01_reconstruction/02_CRUTEM5_vs_HadSST4.pdf", width = 8, height=4)
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





## SI05 GMST-annual vs. Mean-removed!
{

  # GMST reconstruction-mod_p1:
  pdf(file = "figures/01_reconstruction/SI05_GMST_land_vs_ocean_correlation_MR_mod_p1.pdf", width = 8, height=8)
  {
    par(mfrow=c(2, 1), mar=c(2,4,1,1))
    ylim = c(-1, 1); xlim = c(1850,2020)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Temperature Anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1)
    axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
    
    ## polygons for uncertainties:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tas_land$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tas_land$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("darkorange", alpha = 50), border = make.transparent.color("darkorange", 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS.tos$GMST_FM$ann$mod_p1_min_2.5, 
                                               rev(OBS.tos$GMST_FM$ann$mod_p1_min_97.5)), 
            col = make.transparent.color("blue", alpha = 50), border = make.transparent.color("darkblue", 100))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = OBS.tas_land$GMST_FM$ann$mod_p1_min, col = "darkorange", lwd = 2)
    lines(x = 1850:2020, y = OBS.tos$GMST_FM$ann$mod_p1_min, col = "darkblue", lwd = 2)

    # Include Mean-Removed reconstruction:
    lines(x = 1850:2020, y = colMedians(OBS.tos_MR_$GMST_FM$ann$mod_p1_min), col = "darkred", lwd = 2)
    
    legend("topleft", c("CRUTEM5-reconstr.", "HadSST4-reconstr.", "HadSST4-reconstr., m.r."),
           lty = c(1, 1, 1), col = c("darkorange", "darkblue", "darkred"), lwd = c(2, 2, 2), inset = 0.02, 
           cex = 0.75, ncol = 2, title = "GMST reconstruction (train unc./bias incl.)")
    
    
    w.width = 51
    dat.frame = data.frame(cbind(OBS.tas_land$GMST_FM$ann$mod_p1_min, OBS.tos$GMST_FM$ann$mod_p1_min, colMedians(OBS.tos_MR_$GMST_FM$ann$mod_p1_min)))
    CRUTEM5.HadSST4.recons.cor = rollapply(data = dat.frame, width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA)
    HadSST4.HadSST4_MR.recons.cor = rollapply(data = dat.frame, width = w.width, FUN=function(x) { cor(x[,2], x[,3], use="complete.obs") }, by.column = F, fill = NA)
    CRUTEM5.HadSST4_MR.recons.cor = rollapply(data = dat.frame, width = w.width, FUN=function(x) { cor(x[,1], x[,3], use="complete.obs") }, by.column = F, fill = NA)
    
    col = c("wheat4", "darkorange", "darkblue")
    col.paired = brewer.pal(n = 8, name = "Paired")
    
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = "Pearson Correlation", xlab = "", main = "", ylim = c(-0.2,1), xlim = xlim, las=1)
    axis(side = 1, at = seq(min(xlim), max(xlim), 5), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
    
    lines(x = 1850:2020, y = CRUTEM5.HadSST4.recons.cor, col = "darkgray", lwd = 2)
    lines(x = 1850:2020, y = HadSST4.HadSST4_MR.recons.cor, col = col.paired[2], lwd = 2)
    lines(x = 1850:2020, y = CRUTEM5.HadSST4_MR.recons.cor, col = col.paired[8], lwd = 2)
    
    legend("bottomright", c("Land- vs. SST reconstr. (OBS)", 
                        "Land- vs. SST(m.r.) reconstr. (OBS)", 
                        "SST- vs. SST(m.r.) reconstr. (OBS)"), 
           col = c("darkgray", col.paired[8], col.paired[2]), bg = "white", lty = c(1, 1, 1), lwd = 2, 
           title = "Pearson Correlation \n in 50-year window", inset=0.02, cex = 0.7)
    
  }
  dev.off()
  

}



