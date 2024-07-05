
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04b_master_read_paleo_reconstructions.R")
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04c_compute_trends_4paleo-comparison.R")



# ------------------------------------------------------------------------------------
# Period Difference Figures
# ------------------------------------------------------------------------------------
# setwd("figures/04_paleo_data/")
library(vioplot)


# 1. Ocean 2k:
# ------------------------------------------------------------------------------------

file.name = "figures/04_paleo_data/04b_ocean2k_period_diff.pdf"
ylim = c(-0.6, 0.4)
xlim = c(0,8)
trend.ix = 9


pdf(file = file.name, width = 5, height = 3)
{
    par(mfrow=c(1, 1), mar=c(1,5,1,1))
    # ylim = c(-0.9, 0.5); xlim = c(0,8)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Regional sea surface temperature \n  difference [°C], 1901-20 vs. 1871-90", sep = ""), xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    
    {
      # Tropics, full reconstruction:
      at = 0.7; wex = 0.6
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.Tropics.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.Tropics.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      points(x = at + 1, y = ocean2k.trends[trend.ix], pch = 25, bg = make.transparent.color("darkorchid4", 100), cex = 1.5)
      lines(x = c(2,2), y = c(-1, 1))
      
      # WPacific:
      at = 2.7
      vioplot(x = OBS.WPacific_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WPacific_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WPacific[which(get.RE_score(cur.region = Tierney_regions$wpacific) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$wpacific[[2]], trend.years = trend.years, years = Tierney_best$wpacific[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = wex, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1.5)
      lines(x = c(4,4), y = c(-1, 1))
      
      # Indian Ocean:  
      at = 4.7
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_IOcean[which(get.RE_score(cur.region = Tierney_regions$indian) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$indian[[2]], trend.years = trend.years, years = Tierney_best$indian[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1.5)
      lines(x = c(6,6), y = c(-1, 1))
      
      # WAtlantic:
      at = 6.7
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WAtlantic[which(get.RE_score(cur.region = Tierney_regions$atlantic) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$atlantic[[2]], trend.years = trend.years, years = Tierney_best$atlantic[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1.5)
      lines(x = c(8,8), y = c(-1, 1))
    }
    
    ### Legend:
    text(x = 1.4, y = ylim[2] - 0.1, labels = "Tropics", pos=2, cex = 0.7)
    text(x = 3.4, y = ylim[2] - 0.1, labels = "West \n Pacific", pos=2, cex = 0.7)
    text(x = 5.4, y = ylim[2] - 0.1, labels = "Indian \n Ocean", pos=2, cex = 0.7)
    text(x = 7.4, y = ylim[2] - 0.1, labels = "West \n Atlantic", pos=2, cex = 0.7)
    # text(x = 10, y = 0.4, labels = "East \n Pacific", pos=2)
    
    legend("bottomleft", c("Land-based reconstruction \n of ocean region", "HadSST4, region average", "Tierney et al., 2015, \n 'Best reconstruction'"), 
           col = c(make.transparent.color("darkorange", 75), 
                                                                   make.transparent.color("blue", 75),
                                                                   "black"), lwd = c(6,6, NA), pch = c(NA, NA, 25), pt.bg = c(NA, NA, make.transparent.color("darkorchid4", 100)), cex = 0.6,
           inset = 0.01, ncol = 2, bg = "white")
    
    # legend("bottomleft", c("Land-based reconstruction \n of ocean region", "HadSST4, region average", "Tierney et al., 2015, \n 'Best reconstruction'",
    #                       "Tierney et al., 2015, Range"), col = c(make.transparent.color("orange", 75), 
    #                                                               make.transparent.color("blue", 75),
    #                                                               rep(make.transparent.color("darkorchid4", 100), 2)), lwd = c(6,6, NA, 6), pch = c(NA, NA, 25, NA), pt.bg = c(NA, NA, "blueviolet", NA), cex = 0.6,
    #       inset = 0.01, ncol = 2, bg = "white")
  }
dev.off()







# 2. Neukom global-scale climate field reconstructions + GMST estimates:
# ------------------------------------------------------------------------------------
file.name = "04a_Pages2k_period_diff.pdf"
ylim = c(-0.5, 0.5)
xlim = c(0,6.3)
trend.ix = 9


  pdf(file = file.name, width = 5, height = 3)
  {
    par(mfrow=c(1, 1), mar=c(0.2,5,0.2,0.2))
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Global mean surface temperature \n difference [°C], 1901-20 vs. 1871-90", sep = ""), xlab = "", 
         main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    
    at = 0.5; wex = 0.6; cex = 0.6
    
    # Instrumental reconstruction:
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    text(x = at-0.1, y = ylim[2]-0.2, labels = "CRUTEM5-based GMST", srt = 90, cex = cex)
    text(x = at+0.2, y = ylim[2]-0.2, labels = "HadSST4-based GMST", srt = 90, cex = cex)
    
    # add Obs. GMST datasets for comparison:
    HadCRUT5.trend = get.trend(x = HadCRUT5.global.annual$Anomaly, trend.years = trend.years, years = HadCRUT5.global.annual$Year)
    BEST.trend = get.trend(x = BEST.global.annual$Anomaly, trend.years = trend.years, years = BEST.global.annual$Year)
    # NOAA.period.diff = get.trend(x = NOAA.global.annual$Anomaly, trend.years = trend.years, years = NOAA.global.annual$Year)
    # JMA.global.annual
    # GISS.global.annual
    # ... only start in 1880 or later.
    
    points(x = at + 1 - 0.2, y =  HadCRUT5.trend[trend.ix], cex = 1, pch = 1)
    points(x = at + 1 , y =  BEST.trend[trend.ix], cex = 1, pch = 2)
    # points(x = at + 1 + 0.2, y =  NOAA.period.diff, cex = 1, pch = 3)
    
    text(x = at + 1 - 0.2, y = 0.3, labels = "HadCRUT5", srt = 90, cex = cex)
    text(x = at + 1, y = 0.3, labels = "Berkeley Earth", srt = 90, cex = cex)
    # text(x = at + 1 + 0.2, y = 0.3, labels = "NOAA Glob. Temp. v.5", srt = 90, cex = cex)
    
    # vioplot(x = CW14.period.diff, at = at + 1, wex = wex, col = make.transparent.color("grey", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = CW14.period.diff, col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    
    # Neukom Climate Field reconstructions:
    # vioplot(x = Neukom.periodmean$GMST_tas_land3 - Neukom.periodmean$GMST_tas_land1, at = at + 1, wex = wex, col = make.transparent.color("chartreuse4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tas_land3 - Neukom.periodmean$GMST_tas_land1, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tos3 - Neukom.periodmean$GMST_tos1, at = at + 1.1, wex = wex, col = make.transparent.color("blueviolet", alpha = 75),
    #        side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tos3 - Neukom.periodmean$GMST_tos1, col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 1.1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    # text(x = at+1-0.1, y = 0.3, labels = "Land CFR-based GMST", srt = 90, cex = cex)
    # text(x = at+1+0.2, y = 0.3, labels = "Ocean CFR-based GMST", srt = 90, cex = cex)
    # text(x = at+1+0.4, y = 0.3, labels = "(Neukom et al., 2019)", srt = 90, cex = cex)
    
    # Neukom GMST estimates:
    vioplot(x = c(sapply(neukom2019_trend_all, FUN = function(x) x[trend.ix,])), at = at + 2, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = c(sapply(neukom2019_trend_all, FUN = function(x) x[trend.ix,])), col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    at = 1
    neukom_BHM = neukom2019_trend_all[[1]][trend.ix,]
    vioplot(x = neukom_BHM, at = at + 2, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_BHM, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_CPS_new = neukom2019_trend_all[[2]][trend.ix,]
    vioplot(x = neukom_CPS_new, at = at + 2.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_CPS_new, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_DA = neukom2019_trend_all[[3]][trend.ix,]
    vioplot(x = neukom_DA, at = at + 3, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_DA, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 3, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_M08 = neukom2019_trend_all[[4]][trend.ix,]
    vioplot(x = neukom_M08, at = at + 3.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_M08, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 3.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_OIE = neukom2019_trend_all[[5]][trend.ix,]
    vioplot(x = neukom_OIE, at = at + 4, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_OIE, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 4, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_PAI = neukom2019_trend_all[[6]][trend.ix,]
    vioplot(x = neukom_PAI, at = at + 4.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_PAI, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 4.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_PCR = neukom2019_trend_all[[7]][trend.ix,]
    vioplot(x = neukom_PCR, at = at + 5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_PCR, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    at = 0.5
    text(x = at+2-0.1, y = 0.3, labels = "Full Ensemble", srt = 90, cex = cex)
    at = 1
    text(x = at+2-0.1, y = 0.3, labels = "BHM", srt = 90, cex = cex)
    text(x = at+2.5-0.1, y = 0.3, labels = "CPS", srt = 90, cex = cex)
    text(x = at+3-0.1, y = 0.3, labels = "DA", srt = 90, cex = cex)
    text(x = at+3.5-0.1, y = 0.3, labels = "M08", srt = 90, cex = cex)
    text(x = at+4-0.1, y = 0.3, labels = "OIE", srt = 90, cex = cex)
    text(x = at+4.5-0.1, y = 0.3, labels = "PAI", srt = 90, cex = cex)
    text(x = at+5-0.1, y = 0.3, labels = "PCR", srt = 90, cex = cex)
    text(x = at+3-0.1, y = 0.4, labels = substitute(paste(bold("PAGES2k GMST reconstructions \n (PAGES 2k Consortium, 2019)"))), cex = cex)
    
    
    # Percentile comparison:
    # length(which(HadCRUT5.period.diff > neukom_BHM[,3] - neukom_BHM[,1])) / 1000 * 100  # 0.4 percentile
    # length(which(HadCRUT5.period.diff > neukom_CPS_new[,3] - neukom_CPS_new[,1])) / 1000 * 100 # 0.6 percentile
    # length(which(HadCRUT5.period.diff > neukom_DA[,3] - neukom_DA[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > neukom_M08[,3] - neukom_M08[,1])) / 1000 * 100 # 13.3 percentile
    # length(which(HadCRUT5.period.diff > neukom_OIE[,3] - neukom_OIE[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > neukom_PAI[,3] - neukom_PAI[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > (neukom_PCR[,3] - neukom_PCR[,1]))) / 1000 * 100 # 2.4th percentile  
    
    # Percentile comparison of best estimates:
    # length(which(OBS.trends$GMST_tas_land > neukom2019_trend_all)) / 7000 * 100  # 48.9th percentile
    # length(which(OBS.trends$GMST_tos > neukom2019_trend_all)) / 7000 * 100  # 0th percentile
    
    # length(which(HadCRUT5.period.diff > neukom2019_trend_all)) / 7000 * 100  # 2.4th percentile
    # length(which(BEST.period.diff > neukom2019_trend_all)) / 7000 * 100  # 6.5th percentile
    # length(which(NOAA.period.diff > neukom2019_trend_all)) / 7000 * 100  # 0.3th percentile
    
    at = 0.5; cex = 0.5
    # text(x = at-0.1, y = -0.45, labels = round(length(which(mean(OBS.periodmean$GMST_tas_land3 - OBS.periodmean$GMST_tas_land1) > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at+0.2, y = -0.45, labels = round(length(which(mean(OBS.periodmean$GMST_tos3 - OBS.periodmean$GMST_tos1) > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    
    # text(x = at + 1 -0.25, y = -0.45, labels = round(length(which(HadCRUT5.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at + 1, y = -0.45, labels = round(length(which(BEST.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at + 1 + 0.25, y = -0.45, labels = round(length(which(NOAA.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    
    # text(x = at + 3, y = -0.45, labels = "Percentile of Best-Estimate Instrumental Dataset \n in Pages2k GMST reconstruction ensemble", srt = 0, cex = cex)
    
    
  }
  dev.off()







# ------------------------------------------------------------------------------------
# TREND FIGURES
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/04_paleo_data/")
library(vioplot)

trend.years[[2]]
trend.year.name = "1901-1940"
trend.ix = 3


plot.ocean2k <- function(file.name, trend.year.name, trend.ix, ylim = c(-0.9, 0.5), xlim = c(0,8)) {
  
  pdf(file = file.name, width = 5, height = 3)
  {
    par(mfrow=c(1, 1), mar=c(0.2,5,0.2,0.2))
    # ylim = c(-0.9, 0.5); xlim = c(0,8)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Regional sea surface temperature trend \n [° per n years], " , trend.year.name, sep = ""), xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    
    {
      # Tropics, full reconstruction:
      at = 0.7; wex = 0.6
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = OBS.Tropics_tos.trends[trend.ix,], at = at + 0.1, wex = 0.5, col = "lightblue", side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.Tropics_tos.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
      axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.Tropics.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.Tropics.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      points(x = at + 1, y = ocean2k.trends[trend.ix], pch = 25, bg = make.transparent.color("darkorchid4", 100), cex = 1.5)
      lines(x = c(2,2), y = c(-1, 1))
      
      # WPacific:
      at = 2.7
      vioplot(x = OBS.WPacific_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WPacific_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = OBS.WPacific_tos.trends[trend.ix,], at = at + 0.1, wex = 0.5, col = "lightblue", side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WPacific_tos.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WPacific[which(get.RE_score(cur.region = Tierney_regions$wpacific) > 0),trend.ix]
      paleo.trend.best = get.trend(x = Tierney_best$wpacific[[2]], trend.years = trend.years, years = Tierney_best$wpacific[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = wex, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(4,4), y = c(-1, 1))
      
      # Indian Ocean:  
      at = 4.7
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_IOcean[which(get.RE_score(cur.region = Tierney_regions$indian) > 0),trend.ix]
      paleo.trend.best = get.trend(x = Tierney_best$indian[[2]], trend.years = trend.years, years = Tierney_best$indian[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(6,6), y = c(-1, 1))
      
      # WAtlantic:
      at = 6.7
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WAtlantic[which(get.RE_score(cur.region = Tierney_regions$atlantic) > 0),trend.ix]
      paleo.trend.best = get.trend(x = Tierney_best$atlantic[[2]], trend.years = trend.years, years = Tierney_best$atlantic[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(8,8), y = c(-1, 1))
    }
    
    ### Legend:
    text(x = 1.4, y = ylim[2] - 0.1, labels = "Tropics", pos=2, cex = 0.7)
    text(x = 3.4, y = ylim[2] - 0.1, labels = "West \n Pacific", pos=2, cex = 0.7)
    text(x = 5.4, y = ylim[2] - 0.1, labels = "Indian \n Ocean", pos=2, cex = 0.7)
    text(x = 7.4, y = ylim[2] - 0.1, labels = "West \n Atlantic", pos=2, cex = 0.7)
    # text(x = 10, y = 0.4, labels = "East \n Pacific", pos=2)
    
    legend("bottomleft", c("Land-based reconstruction \n of ocean region", "HadSST4, region average", "Tierney et al., 2015, \n 'Best reconstruction'",
                           "Tierney et al., 2015, Range"), col = c(make.transparent.color("orange", 75), 
                                                                   make.transparent.color("blue", 75),
                                                                   rep(make.transparent.color("darkorchid4", 100), 2)), lwd = c(6,6, NA, 6), pch = c(NA, NA, 25, NA), pt.bg = c(NA, NA, "blueviolet", NA), cex = 0.6,
           inset = 0.01, ncol = 2, bg = "white")
  }
  dev.off()
  
}


## 'Turn-of-century cooling period'
plot.ocean2k(file.name = "SI_04a_ocean2k_trend_comparison_1871-1910.pdf", trend.year.name = "1871-1910", trend.ix = 2, 
             ylim = c(-1, 0.5), xlim = c(0,8))

plot.ocean2k(file.name = "SI_04a_ocean2k_trend_comparison_1871-1920.pdf", trend.year.name = "1871-1920", trend.ix = 6, 
             ylim = c(-1, 0.5), xlim = c(0,8))


## 'Early Twentieth Century Warming period' (ETCW)
plot.ocean2k(file.name = "SI_04b_ocean2k_trend_comparison_1901-1940.pdf", trend.year.name = "1901-1940", trend.ix = 3, 
             ylim = c(-0.5, 1), xlim = c(0,8))

plot.ocean2k(file.name = "SI_04b_ocean2k_trend_comparison_1901-1950.pdf", trend.year.name = "1901-1950", trend.ix = 7, 
             ylim = c(-0.5, 1), xlim = c(0,8))










# 2. Neukom global-scale climate field reconstructions + GMST estimates:
# ------------------------------------------------------------------------------------
file.name = "03a_pages2k_trend.pdf"
trend.year.name = "1871-1910"
trend.ix = 2


plot.pages2k <- function(file.name, trend.year.name, trend.ix, ylim = c(-0.5, 0.5), xlim = c(0,6.3)) {
  
  pdf(file = file.name, width = 5, height = 3)
  {
    par(mfrow=c(1, 1), mar=c(0.2,5,0.2,0.2))
  
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Global mean surface temperature trend \n [° per n years], " , trend.year.name, sep = ""), xlab = "", 
         main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    
    at = 0.5; wex = 0.6; cex = 0.6
    
    # Instrumental reconstruction:
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    text(x = at-0.1, y = ylim[2]-0.2, labels = "CRUTEM5-based GMST", srt = 90, cex = cex)
    text(x = at+0.2, y = ylim[2]-0.2, labels = "HadSST4-based GMST", srt = 90, cex = cex)
    
    # add Obs. GMST datasets for comparison:
    HadCRUT5.trend = get.trend(x = HadCRUT5.global.annual$Anomaly, trend.years = trend.years, years = HadCRUT5.global.annual$Year)
    BEST.trend = get.trend(x = BEST.global.annual$Anomaly, trend.years = trend.years, years = BEST.global.annual$Year)
    # NOAA.period.diff = get.trend(x = NOAA.global.annual$Anomaly, trend.years = trend.years, years = NOAA.global.annual$Year)
    # JMA.global.annual
    # GISS.global.annual
    # ... only start in 1880 or later.
    
    points(x = at + 1 - 0.2, y =  HadCRUT5.trend[trend.ix], cex = 1, pch = 1)
    points(x = at + 1 , y =  BEST.trend[trend.ix], cex = 1, pch = 2)
    # points(x = at + 1 + 0.2, y =  NOAA.period.diff, cex = 1, pch = 3)
    
    text(x = at + 1 - 0.2, y = 0.3, labels = "HadCRUT5", srt = 90, cex = cex)
    text(x = at + 1, y = 0.3, labels = "Berkeley Earth", srt = 90, cex = cex)
    # text(x = at + 1 + 0.2, y = 0.3, labels = "NOAA Glob. Temp. v.5", srt = 90, cex = cex)
    
    # vioplot(x = CW14.period.diff, at = at + 1, wex = wex, col = make.transparent.color("grey", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = CW14.period.diff, col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    
    # Neukom Climate Field reconstructions:
    # vioplot(x = Neukom.periodmean$GMST_tas_land3 - Neukom.periodmean$GMST_tas_land1, at = at + 1, wex = wex, col = make.transparent.color("chartreuse4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tas_land3 - Neukom.periodmean$GMST_tas_land1, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tos3 - Neukom.periodmean$GMST_tos1, at = at + 1.1, wex = wex, col = make.transparent.color("blueviolet", alpha = 75),
    #        side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    # vioplot(x = Neukom.periodmean$GMST_tos3 - Neukom.periodmean$GMST_tos1, col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 1.1, wex = 1.5, 
    #        axes = F, add = T, pchMed = 21, frame.plot = F)
    # text(x = at+1-0.1, y = 0.3, labels = "Land CFR-based GMST", srt = 90, cex = cex)
    # text(x = at+1+0.2, y = 0.3, labels = "Ocean CFR-based GMST", srt = 90, cex = cex)
    # text(x = at+1+0.4, y = 0.3, labels = "(Neukom et al., 2019)", srt = 90, cex = cex)
    
    # Neukom GMST estimates:
    vioplot(x = c(sapply(neukom2019_trend_all, FUN = function(x) x[trend.ix,])), at = at + 2, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = c(sapply(neukom2019_trend_all, FUN = function(x) x[trend.ix,])), col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    at = 1
    neukom_BHM = neukom2019_trend_all[[1]][trend.ix,]
    vioplot(x = neukom_BHM, at = at + 2, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_BHM, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_CPS_new = neukom2019_trend_all[[2]][trend.ix,]
    vioplot(x = neukom_CPS_new, at = at + 2.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_CPS_new, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 2.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_DA = neukom2019_trend_all[[3]][trend.ix,]
    vioplot(x = neukom_DA, at = at + 3, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_DA, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 3, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
  
    neukom_M08 = neukom2019_trend_all[[4]][trend.ix,]
    vioplot(x = neukom_M08, at = at + 3.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_M08, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 3.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_OIE = neukom2019_trend_all[[5]][trend.ix,]
    vioplot(x = neukom_OIE, at = at + 4, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_OIE, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 4, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_PAI = neukom2019_trend_all[[6]][trend.ix,]
    vioplot(x = neukom_PAI, at = at + 4.5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_PAI, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 4.5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    neukom_PCR = neukom2019_trend_all[[7]][trend.ix,]
    vioplot(x = neukom_PCR, at = at + 5, wex = wex, col = make.transparent.color("burlywood4", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = neukom_PCR, col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at + 5, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    
    at = 0.5
    text(x = at+2-0.1, y = 0.3, labels = "Full Ensemble", srt = 90, cex = cex)
    at = 1
    text(x = at+2-0.1, y = 0.3, labels = "BHM", srt = 90, cex = cex)
    text(x = at+2.5-0.1, y = 0.3, labels = "CPS", srt = 90, cex = cex)
    text(x = at+3-0.1, y = 0.3, labels = "DA", srt = 90, cex = cex)
    text(x = at+3.5-0.1, y = 0.3, labels = "M08", srt = 90, cex = cex)
    text(x = at+4-0.1, y = 0.3, labels = "OIE", srt = 90, cex = cex)
    text(x = at+4.5-0.1, y = 0.3, labels = "PAI", srt = 90, cex = cex)
    text(x = at+5-0.1, y = 0.3, labels = "PCR", srt = 90, cex = cex)
    text(x = at+3-0.1, y = 0.4, labels = substitute(paste(bold("PAGES2k GMST reconstructions \n (PAGES 2k Consortium, 2019)"))), cex = cex)
    
    # Percentile comparison:
    # length(which(HadCRUT5.period.diff > neukom_BHM[,3] - neukom_BHM[,1])) / 1000 * 100  # 0.4 percentile
    # length(which(HadCRUT5.period.diff > neukom_CPS_new[,3] - neukom_CPS_new[,1])) / 1000 * 100 # 0.6 percentile
    # length(which(HadCRUT5.period.diff > neukom_DA[,3] - neukom_DA[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > neukom_M08[,3] - neukom_M08[,1])) / 1000 * 100 # 13.3 percentile
    # length(which(HadCRUT5.period.diff > neukom_OIE[,3] - neukom_OIE[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > neukom_PAI[,3] - neukom_PAI[,1])) / 1000 * 100 # 0 percentile
    # length(which(HadCRUT5.period.diff > (neukom_PCR[,3] - neukom_PCR[,1]))) / 1000 * 100 # 2.4th percentile  
  
    # Percentile comparison of best estimates:
    # length(which(OBS.trends$GMST_tas_land > neukom2019_trend_all)) / 7000 * 100  # 48.9th percentile
    # length(which(OBS.trends$GMST_tos > neukom2019_trend_all)) / 7000 * 100  # 0th percentile
    
    # length(which(HadCRUT5.period.diff > neukom2019_trend_all)) / 7000 * 100  # 2.4th percentile
    # length(which(BEST.period.diff > neukom2019_trend_all)) / 7000 * 100  # 6.5th percentile
    # length(which(NOAA.period.diff > neukom2019_trend_all)) / 7000 * 100  # 0.3th percentile
    
    at = 0.5; cex = 0.5
    # text(x = at-0.1, y = -0.45, labels = round(length(which(mean(OBS.periodmean$GMST_tas_land3 - OBS.periodmean$GMST_tas_land1) > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at+0.2, y = -0.45, labels = round(length(which(mean(OBS.periodmean$GMST_tos3 - OBS.periodmean$GMST_tos1) > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    
    # text(x = at + 1 -0.25, y = -0.45, labels = round(length(which(HadCRUT5.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at + 1, y = -0.45, labels = round(length(which(BEST.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    # text(x = at + 1 + 0.25, y = -0.45, labels = round(length(which(NOAA.period.diff > neukom2019_diff_all)) / 7000 * 100, 1), srt = 0, cex = cex)
    
    # text(x = at + 3, y = -0.45, labels = "Percentile of Best-Estimate Instrumental Dataset \n in Pages2k GMST reconstruction ensemble", srt = 0, cex = cex)
  
  
  }
  dev.off()
}


plot.pages2k(file.name = "03a_pages2k_trend_comparison_1871-1910.pdf", trend.year.name = "1871-1910", trend.ix = 2, 
             ylim = c(-0.6, 0.6), xlim = c(0,8))
plot.pages2k(file.name = "03a_pages2k_trend_comparison_1871-1920.pdf", trend.year.name = "1871-1920", trend.ix = 6, 
             ylim = c(-0.6, 0.6), xlim = c(0,8))

plot.pages2k(file.name = "03b_pages2k_trend_comparison_1901-1940.pdf", trend.year.name = "1901-1940", trend.ix = 3, 
             ylim = c(0, 1.0), xlim = c(0,8))
plot.pages2k(file.name = "03b_pages2k_trend_comparison_1901-1950.pdf", trend.year.name = "1901-1950", trend.ix = 7, 
             ylim = c(0, 1.0), xlim = c(0,8))

# Hmm. Pretty strange.  Continue here...



##### END FOR NOW #####


##### THE END ########
  

# Cold ocean anomaly shows up clearly. Reasons why this cold ocean anomaly may be unrealistic:
# (1) literature: biases in ocean are prevalent; comparison to coastal stations appears to show that ocean anomaly is unrealistic.
# (2) "decoupling" of low-frequency variability (poor match between land- and ocean reconstruction) from high-frequency variability (very good match).
# (3) none of the cmip6 models shows a temporal "decoupling pattern" between low- and high-frequency variability.
# (4) none of the cmip6 models supports an ocean-land temperature difference close to that observed.





