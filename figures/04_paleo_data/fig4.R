
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021


## load all data for reconstructions:
source("scripts/04a_master_load_reconstructions.R")
source("scripts/04b_master_read_paleo_reconstructions.R")
source("scripts/04c_compute_trends_4paleo-comparison.R")



# ------------------------------------------------------------------------------------
# Period Difference Figures
# ------------------------------------------------------------------------------------
library(vioplot)


# 1. Ocean 2k:
# ------------------------------------------------------------------------------------

file.name = "figures/04_paleo_data/04b_ocean2k_period_diff.pdf"
ylim = c(-0.6, 0.4)
xlim = c(0,8)
trend.ix = 9


pdf(file = file.name, width = 3.5, height = 1.75, pointsize = 7)
{
    par(mfrow=c(1, 1), mar=c(1,5,1,1))
    # ylim = c(-0.9, 0.5); xlim = c(0,8)
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Regional sea surface temperature \n  change [°C], [1901-20] vs. [1871-90]", sep = ""), xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    lines(x = c(0, 1000), y = c(0,0), col = "lightgrey")
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    # mtext("a", side = 3, line = 0, at = -0.9, adj = 0, cex = 1.1, font = 2)
    mtext("b", side = 3, line = 0, at = -1.5, adj = 0, cex = 0.9, font = 2)
    
    {
      # Tropics, full reconstruction:
      at = 0.7; wex = 0.6
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.Tropics_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      #points(x = at + 0.2, y = BEST_Tropics[9], pch = 25, bg = make.transparent.color("darkblue", 250), cex = 1.5)
      
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
      #points(x = at + 0.2, y = BEST_WPacific.trend[9], pch = 25, bg = make.transparent.color("darkblue", 250), cex = 1.5)
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WPacific.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WPacific[which(get.RE_score(cur.region = Tierney_regions$wpacific) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$wpacific[[2]], trend.years = trend.years, years = Tierney_best$wpacific[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = wex, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(4,4), y = c(-1, 1))
      
      # Indian Ocean:  
      at = 4.7
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.IndianOcean_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      #points(x = at + 0.2, y = BEST_Indian_Ocean.trend[9], pch = 25, bg = make.transparent.color("darkblue", 250), cex = 1.5)
      
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.IndianOcean.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_IOcean[which(get.RE_score(cur.region = Tierney_regions$indian) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$indian[[2]], trend.years = trend.years, years = Tierney_best$indian[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(6,6), y = c(-1, 1))
      
      # WAtlantic:
      at = 6.7
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n")
      vioplot(x = OBS.WAtlantic_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      #points(x = at + 0.2, y = BEST_WAtlantic.trend[9], pch = 25, bg = make.transparent.color("darkblue", 250), cex = 1.5)
      
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n")
      vioplot(x = HadSST4.WAtlantic.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
              axes = F, add = T, pchMed = 21)
      
      paleo.trend = Tierney_WAtlantic[which(get.RE_score(cur.region = Tierney_regions$atlantic) > 0),trend.ix]
      paleo.trend.best = get.trend_perioddiff(x = Tierney_best$atlantic[[2]], trend.years = trend.years, years = Tierney_best$atlantic[[1]])[trend.ix]
      
      vioplot(x = paleo.trend, at = at + 0.8, wex = 0.5, col = make.transparent.color("darkorchid4", alpha = 100), side = "left", add = T, pchMed = 21, yaxt="n")
      points(x = at + 1, y = paleo.trend.best, pch = 25, bg = make.transparent.color("darkorchid4", alpha = 100), cex = 1)
      lines(x = c(8,8), y = c(-1, 1))
    }
    
    ### Legend:
    text(x = 1.4, y = ylim[2] - 0.1, labels = "Tropics", pos=2, cex = 0.71)
    text(x = 3.4, y = ylim[2] - 0.1, labels = "West \n Pacific", pos=2, cex = 0.71)
    text(x = 5.4, y = ylim[2] - 0.1, labels = "Indian \n Ocean", pos=2, cex = 0.71)
    text(x = 7.4, y = ylim[2] - 0.1, labels = "West \n Atlantic", pos=2, cex = 0.71)
    # text(x = 10, y = 0.4, labels = "East \n Pacific", pos=2)
    
    
    legend("bottomleft", c("Land-based reconstruction \n of ocean region", "HadSST4, region average", "Best reconstruction \n (Tierney et al., 2015)",
                           "Range (Tierney et al., 2015)"), col = c(make.transparent.color("orange", 75), 
                                                                   make.transparent.color("blue", 75),
                                                                   "black", make.transparent.color("darkorchid4", 100)), lwd = c(6,6, NA, 6), pch = c(NA, NA, 25, NA), pt.bg = c(NA, NA, make.transparent.color("darkorchid4", 100), NA), cex = 0.71,
           inset = 0.01, ncol = 2, bg = "white", pt.cex = 1)
  }
dev.off()







# 2. Neukom global-scale climate field reconstructions + GMST estimates:
# ------------------------------------------------------------------------------------
file.name = "figures/04_paleo_data/04a_Pages2k_period_diff.pdf"
ylim = c(-0.5, 0.5)
xlim = c(0,6.3)
trend.ix = 9


  pdf(file = file.name, width = 3.5, height = 1.75, pointsize = 7)
  {
    par(mfrow=c(1, 1), mar=c(0.2,5,0.2,0.2))
    
    plot(x = 1850:2020, y = 1850:2020, type="n", 
         ylab = paste("Global mean surface temperature \n change [°C], [1901-20] vs. [1871-90]", sep = ""), xlab = "", 
         main = "", ylim = ylim, xlim = xlim, las=1, yaxs="i", xaxt="n", xaxs="i", cex.axis = 0.8, cex.lab = 0.8)
    lines(x = c(0, 1000), y = c(0,0), col = "lightgrey")
    axis(side = 2, at = seq(-0.5, 0.5, 0.1), tcl=0.2, labels=F)
    mtext("a", side = 3, line = -0.5, at = -1.15, adj = 0, cex = 0.9, font = 2)
    
    at = 0.5; wex = 0.6; cex = 0.71
    
    # Instrumental reconstruction:
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], at = at, wex = wex, col = make.transparent.color("darkorange", alpha = 75), side = "left", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tas_land.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="left", at = at, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], at = at + 0.1, wex = wex, col = make.transparent.color("blue", alpha = 75), side = "right", add = T, pchMed = 21, yaxt="n", frame.plot = F)
    vioplot(x = OBS.GMST_tos.trends[trend.ix,], col = make.transparent.color("white", alpha = 0), border = NA, side="right", at = at + 0.1, wex = 1.5, 
            axes = F, add = T, pchMed = 21, frame.plot = F)
    text(x = at-0.1, y = ylim[2]-0.2, labels = c(expression(hat(T)[CRUTEM5]^{GMST})), srt = 90, cex = cex)
    text(x = at+0.2, y = ylim[2]-0.2, labels = c(expression(hat(T)[HadSST4]^{GMST})), srt = 90, cex = cex)
    
    # add Obs. GMST datasets for comparison:
    HadCRUT5.trend = get.trend_perioddiff(x = HadCRUT5.global.annual$Anomaly, trend.years = trend.years, years = HadCRUT5.global.annual$Year)
    BEST.trend = get.trend_perioddiff(x = BEST.global.annual$Anomaly, trend.years = trend.years, years = BEST.global.annual$Year)
    CW14.trend = get.trend_perioddiff(x = CW2014.global.annual_had4sst4_krig$Anomaly, trend.years = trend.years, years = CW2014.global.annual_had4sst4_krig$Year)
    CW14_COBE.trend = get.trend_perioddiff(x = CW2014.global.annual_cobe2cru_krig$Anomaly, trend.years = trend.years, years = CW2014.global.annual_cobe2cru_krig$Year)
    # NOAA.period.diff = get.trend_perioddiff(x = NOAA.global.annual$Anomaly, trend.years = trend.years, years = NOAA.global.annual$Year)
    # JMA.global.annual
    # GISS.global.annual
    # ... only start in 1880 or later.
    
    points(x = at + 1 - 0.2, y =  HadCRUT5.trend[trend.ix], cex = 1, pch = 8)
    points(x = at + 1 , y =  BEST.trend[trend.ix], cex = 1, pch = 2)
    points(x = at + 1 + 0.2, y =  CW14.trend[trend.ix], cex = 1, pch = 3)
    points(x = at + 1 + 0.4, y =  CW14_COBE.trend[trend.ix], cex = 1, pch = 4)
    
    text(x = at + 1 - 0.25, y = 0.48, labels = "HadCRUT5 (non-infilled)", srt = 90, cex = cex, adj = 1)
    text(x = at + 1, y = 0.48, labels = "Berkeley Earth", srt = 90, cex = cex, adj = 1)
    text(x = at + 1 + 0.25, y = 0.48, labels = "CW14 (HadSST4 based)", srt = 90, cex = cex, adj = 1)
    text(x = at + 1 + 0.53, y = 0.48, labels = "CW14 (COBE-SST2 based)", srt = 90, cex = cex, adj = 1)
    
    
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
    neukom2019_perioddiff_all = (sapply(X = neukom2019_trend_all, FUN=function(x) x[trend.ix,]))
    # length(which( median(OBS.GMST_tas_land.trends[trend.ix,]) > c(neukom2019_perioddiff_all))) / 7000 * 100  # 82.5th percentile
    # length(which( median(OBS.GMST_tos.trends[trend.ix,]) > c(neukom2019_perioddiff_all))) / 7000 * 100  # 0th percentile
    
    # length(which(HadCRUT5.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100  # 2.4th percentile
    # length(which(BEST.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100  # 6.5th percentile
    # length(which(CW14.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100  # 5.9th percentile
    # length(which(CW14_COBE.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100  # 24th percentile

    at = 0.5; cex = 0.71
    text(x = at-0.1, y = -0.45, labels = round(length(which( median(OBS.GMST_tas_land.trends[trend.ix,]) > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    text(x = at+0.2, y = -0.45, labels = round(length(which( median(OBS.GMST_tos.trends[trend.ix,]) > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    
    text(x = at + 1 -0.25, y = -0.45, labels = round(length(which(HadCRUT5.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    text(x = at + 1, y = -0.45, labels = round(length(which(BEST.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    text(x = at + 1 + 0.25, y = -0.45, labels = round(length(which(CW14.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    text(x = at + 1 + 0.55, y = -0.45, labels = round(length(which(CW14_COBE.trend[trend.ix] > c(neukom2019_perioddiff_all))) / 7000 * 100, 1), srt = 0, cex = cex)
    
    text(x = at + 3.5, y = -0.45, labels = "Percentile of instrumental dataset (median) \n in Pages2k GMST reconstruction ensemble", srt = 0, cex = cex)
    
  }
  dev.off()







##### END FOR NOW #####


##### THE END ########
  

# Cold ocean anomaly shows up clearly. Reasons why this cold ocean anomaly may be unrealistic:
# (1) literature: biases in ocean are prevalent; comparison to coastal stations appears to show that ocean anomaly is unrealistic.
# (2) "decoupling" of low-frequency variability (poor match between land- and ocean reconstruction) from high-frequency variability (very good match).
# (3) none of the cmip6 models shows a temporal "decoupling pattern" between low- and high-frequency variability.
# (4) none of the cmip6 models supports an ocean-land temperature difference close to that observed.





