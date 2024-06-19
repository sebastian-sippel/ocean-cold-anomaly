
# ------------------------------------------------------------------------------------
# Evaluate and plot reconstructions.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 03.01.2024


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")



## 01b. Get correction vector(s):
## ----------------------------------------
{
# get correction vector from Cowtan:
cor_vec_hybrid36 = GMST.hybrid36$mod_p1_min_50[1:167] - GMST.tos$mod_p1_min_50[1:167]
cor_vec_hadsst4 = GMST.tos$mod_p1_min_50 - GMST.tos_ua$mod_p1_min_50

mean(c(cor_vec_hybrid36+cor_vec_hadsst4[1:167])[51:90])
mean(c(cor_vec_hadsst4[1:167])[51:90])
}


### TEST IF CORRECTION VECTOR IS OK ?


plot(1850:2020, HadSST4.global.annual$Anomaly[1:171] - HadSST_ua.annual$Anomaly, type='l')
lines(x = c(1906, 1906), y = c(-1,1))
lines(x = c(1910, 1910), y = c(-1,1))


plot(1850:2020, y = colMedians(OBS.tos_$GMST_FM$ann$mod_p1_min[1:100,]) - GMST.tos_ua$mod_p1_min_50, type='l')
lines(1850:2020, y =  colMedians(OBS.tos_$GMST_FM$ann$mod_p1_min[101:200,]) - GMST.tos_ua$mod_p1_min_50, type='l', col = col.HadSST4)
lines(x = c(1906, 1906), y = c(-1,1))
lines(x = c(1910, 1910), y = c(-1,1))



plot(1850:2020, y = colMedians(OBS.tos_$GMST_FM$ann$mod_p1_min[1:100,]) - colMedians(OBS.tos_$GMST_FM$ann$mod_p1_min[101:200,]), type='l')

# plot(x= 1850:2020, y = HadSST_ua.annual$Anomaly, col = "violet", type='l')
# lines(x= 1850:2020, y = HadSST4.global.annual$Anomaly[1:171], col = col.HadSST4, type='l')

# SST Adjustment:
{
  c = -1+0.5*3; sc = 1
  
  plot(x = 1850:2020, y = cor_vec_hadsst4 * sc + c, type='l', ylim = c(-1, 1), col = "blue")
  #lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6) * sc + c, col = "blue", lty = 2)
  lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36) * sc + c, col = "blue", lty = 3)
}




# ------------------------------------------------------------------------------------
# Plot reconstruction(s):
# ------------------------------------------------------------------------------------


# 01_main_reconstruction: AGMT anomalies | the !main! reconstructions (no uncertainties/sensitivities!):
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures/01_reconstruction/")



library(RColorBrewer)
col = brewer.pal(n = 8, name = "Dark2")

library("dplR") # package for band-pass filtering




## ----------------------------------------------------------------------------------------
## Get Numbers for Abstract:
## ----------------------------------------------------------------------------------------


get.range.trends <- function(obs.data, cur.name = "GSAT", mod_type = "mod_p1", late.ix = 162:171, early.ix = 1:51) {
  
  trend.est = sapply(X = 1:200, FUN=function(ix) mean(obs.data[[cur.name]]$ann[[mod_type]][ix,late.ix]) - mean(obs.data[[cur.name]]$ann[[mod_type]][ix,early.ix]))
  trend.range = c(mean(trend.est), quantile(trend.est, probs = c(0.025, 0.5, 0.975)))
  return(trend.range)
}

get.OLS.trends <- function(obs.data, cur.name = "GSAT", mod_type = "mod_p1", years = 1880:2012) {
  
  year.ix = which(c(1850:2020) %in% years)
  trend.est = sapply(X = 1:200, FUN=function(ix) lm(c(obs.data[[cur.name]]$ann[[mod_type]][ix,year.ix]) ~ years)$coefficients[2] * length(years))
  trend.range = c(mean(trend.est), quantile(trend.est, probs = c(0.025, 0.5, 0.975)))
  return(trend.range)
}


## Warming trends:
## ----------------------------------------------------------------------------------------
# GMST reconstructions, 2011-2020 vs. 1850-1900:
get.range.trends(obs.data = OBS.tas_land_, cur.name = "GMST_FM", mod_type = "mod_p1_min", late.ix = 162:171, early.ix = 1:51)  # 1.06°C
get.range.trends(obs.data = OBS.tos_, cur.name = "GMST_FM", mod_type = "mod_p1_min", late.ix = 162:171, early.ix = 1:51)

get.OLS.trends(OBS.tos_, cur.name = "GMST_FM", mod_type = "mod_p1", years = 1901:1950) / 5
get.OLS.trends(OBS.tos_, cur.name = "GMST_FM", mod_type = "mod_p1", years = 1971:2020) / 5
get.OLS.trends(OBS.tas_land_, cur.name = "GMST_FM", mod_type = "mod_p1", years = 1971:2020) / 5



## Difference between land and ocean-based reconstruction:
## ----------------------------------------------------------------------------------------

# mod_p1_min: 1900-1930
mean((colMeans(OBS.tos_$GMST_FM$ann$mod_p1_min) - colMeans(OBS.tas_land_$GMST_FM$ann$mod_p1_min))[52:81])  # Ocean 0.26°C cooler than the land on average

# mod_p0: 1900-1930
mean((colMeans(OBS.tos_$GMST_FM$ann$mod_p0) - colMeans(OBS.tas_land_$GMST_FM$ann$mod_p0))[52:81])  # Ocean 0.20°C cooler than the land on average





## Correlations:
## ----------------------------------------------------------------------------------------
#### GMST:
### ANNUAL:
# 1850-1900
cor(CRUTEM5.global.annual$Anomaly[1:51], HadSST4.global.annual$Anomaly[1:51], use = "complete.obs")  # 0.47
cor(GMST.tas_land$mod_p1_min_50[1:51], GMST.tos$mod_p1_min_50[1:51]) # 0.70 ... hooray!

cor(CRUTEM5.global.annual$Anomaly[101:171], HadSST4.global.annual$Anomaly[101:171], use = "complete.obs") 

cor(GSAT.tas_land$hp_50[1:51], GSAT.tos$hp_50[1:51])
OBS_mod_p1_min = get.diff.window.cor_ens(x = OBS.tos_$GSAT$ann$mod_p1_min, y = OBS.tas_land_$GSAT$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 94:200)

OBS_mod_p1_min$cor_x.high_y.high_50


raw.data.scale$cor_x.high_y.high

cor(GMST.tas_land$mod_p1_min_50[101:171], GMST.tos$mod_p1_min_50[101:171])


mean(GMST.tas_land$mod_p1_min_50[51:(51+35)]) - mean(GMST.tos$mod_p1_min_50[51:(51+35)])



### MONTHLY:
cor(CRUTEM5.global.monthly$Anomaly[1:(51*12)], HadSST4.global.monthly$Anomaly[1:(51*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[(0*12):(51*12)], OBS.tos$GMST_FM$mon$mod_p1[(0*12):(51*12)])  # 0.37 ... hooray!
cor(OBS.tas_land$GMST_FM$mon$mod_p0[1:(51*12)], OBS.tos$GMST_FM$mon$mod_p0[1:(51*12)])  # 0.29 ... hooray!
cor(OBS.tas_land$GMST_FM$mon$mod_gta[1:(51*12)], OBS.tos$GMST_FM$mon$mod_gta[1:(51*12)])  # 0.30 ... hooray!

# 1850-2020
cor(CRUTEM5.global.monthly$Anomaly[1:(170*12)], HadSST4.global.monthly$Anomaly[1:(170*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[1:(170*12)], OBS.tos$GMST_FM$mon$mod_p1_min[1:(170*12)])





# Define colours for plotting: 
col.CRUTEM5 = "darkorange"
col.HadSST4 = "blue"
col.COBE2 = "deepskyblue1"
col.ERSST5 = "deepskyblue4"

col.BEST = "darkgoldenrod4"
col.hybrid36 = "brown4"
col.CLASSNMAT = "blueviolet"
col.HadSST4_ua = "magenta4"




## ----------------------------------------------------------------------------------------
## 01a. Plot GSAT filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GSAT_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)

  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")

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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }

  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$mod_p1_min_2.5, 
                                               rev(GSAT.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$mod_p1_min_2.5, 
                                               rev(GSAT.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GSAT.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$lp_2.5, 
                                               rev(GSAT.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$lp_2.5, 
                                               rev(GSAT.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GSAT.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
      c = 1+0.5*2; sc = 1
      polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$hp_2.5, 
                                                 rev(GSAT.tas_land$hp_97.5)) * sc + c, 
              col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
      
      polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$hp_2.5, 
                                                 rev(GSAT.tos$hp_97.5)) * sc + c, 
              col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
    }

  # Forced response:
  {
      c = -1 + 0.5*3; sc = 1

      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
      lines(x = 1850:2020, y = GSAT.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
  }
    
  # Unforced, LP:
  {
      c = -3 + 0.5*4; sc = 1
    
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tas_land$res_lp_2.5[1:165], 
                                                 rev(GSAT.tas_land$res_lp_97.5[1:165])) * sc + c, 
              col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
      
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tos$res_lp_2.5[1:165], 
                                                 rev(GSAT.tos$res_lp_97.5[1:165])) * sc + c, 
              col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
    }
    
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tas_land$res_hp_2.5[1:165], 
                                                 rev(GSAT.tas_land$res_hp_97.5[1:165])) * sc + c, 
              col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
      
      polygon(x = c(1850:2014, 2014:1850), y = c(GSAT.tos$res_hp_2.5[1:165], 
                                                 rev(GSAT.tos$res_hp_97.5[1:165])) * sc + c, 
              col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
      
      # Lines for main reconstructions:
      lines(x = 1850:2020, y = GSAT.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
      lines(x = 1850:2020, y = GSAT.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
      
      # Lines for other observational datasets:
      # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
      lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
      lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
      
      lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
      lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
      lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
      
      lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
    }
    
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
           lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
           cex = 0.75, ncol = 3, bg = "white")
    
    ### 
    # cor(scale(pass.filt(y = OBS.tas_land$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
    #    scale(pass.filt(y = OBS.tos$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
    
    # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
    #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
    
  

}
dev.off()



## ----------------------------------------------------------------------------------------
## 01b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GSAT_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$mod_p1_min_2.5, 
                                               rev(GSAT.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$mod_p1_min_2.5, 
                                               rev(GSAT.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1, lwd = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1, lwd = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1, lwd = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$lp_2.5, 
                                               rev(GSAT.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$lp_2.5, 
                                               rev(GSAT.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tas_land$hp_2.5, 
                                               rev(GSAT.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GSAT.tos$hp_2.5, 
                                               rev(GSAT.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GSAT.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GSAT.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GSAT.HadISST$Year, y = GSAT.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GSAT.COBE_SST2$Year, y = GSAT.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GSAT.ERSSTv5$Year, y = GSAT.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GSAT.BEST_Land$Year, y = GSAT.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GSAT.hybrid36$Year, y = GSAT.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GSAT.tas_sea$Year, y = GSAT.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GSAT.tos_ua$Year, y = GSAT.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 2, 4, 2, 2, 4), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.ERSST5, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GSAT$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()





## ----------------------------------------------------------------------------------------
## 02a. Plot GMST filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "d. Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "e. Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "f. Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
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
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMST.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 02b. Plot GMST filtered time series with hybrid correction vector:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1_exp.pdf", width = 8, height=11)
{
  par(mfrow=c(1, 1), mar=c(2,5,1,4))
  ylim = c(-4.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.5, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.5 + 0.5*3, labels = "d. Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.5 + 0.5*4, labels = "e. Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "f. Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-7.5, -6, 0.5) + 0.5*6, labels=seq(-0.5, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-7.5, -6, 0.1) + 0.5*6, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-7, -7) + 0.5*6, col = "darkgray", lty = 3)
    text(x = 1850, y = -6.1 + 0.5*6, labels = "g. Implied SST adjustments", col = "grey40", pos = 4)
    
    mtext(text = "Global implied \n SST Adjustment [°C]", side = 2, line = 3, at = -3.8, col = "red", cex = 0.9)
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
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMST.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }

  # SST Adjustment:
  {
    c = -7+0.5*6; sc = 1
    
    lines(x = 1850:2020, y = cor_vec_hadsst4 * sc + c, type='l', ylim = c(-1, 1), col = "red")
    #lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6) * sc + c, col = "blue", lty = 2)
    lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36) * sc + c, col = "red", lty = 3)
  }
  
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  legend("bottomright", c("HadSST4", "Cowtan-HybridSST"),
         lty = c(1, 3), col = c("red", "red"),  
         cex = 0.7, ncol = 2, bg = "white", inset = 0.01)

  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 02b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
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
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1, lwd = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1, lwd = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1, lwd = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()


## ----------------------------------------------------------------------------------------
## 02c. large version with corrections:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-4.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "d. Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "e. Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "f. Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-8, -6, 0.5) + 0.5*6, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-8, -6, 0.1) + 0.5*6, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-7, -7) + 0.5*6, col = "darkgray", lty = 3)
    text(x = 1850, y = -6 + 0.5*6, labels = "g. Implied SST adjustment", col = "pink3", pos = 4)
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
    lines(x = 1850:2020, y = GMST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMST.HadISST$Year, y = GMST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMST.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
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
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }

  # Implied SST Adjustment
  
  # SST Adjustment:
  {
    c = -7 + 0.5*6; sc = 1
    
    lines(x = 1850:2020, y = cor_vec_hadsst4 * sc + c, type='l', ylim = c(-1, 1), col = "pink2")
    #lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6) * sc + c, col = col.HadSST4, lty = 1)
    lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36) * sc + c, col = "pink3", lty = 3)
    
    #lines(x = 1850:2020, y = (HadSST4.global.annual$Anomaly[1:171] - HadSST_ua.annual$Anomaly) * sc + c, col = "black", lty = 3)
  }
  
  
  
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tas_land$res_hp_2.5[1:165], 
                                               rev(GMST.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMST.tos$res_hp_2.5[1:165], 
                                               rev(GMST.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMST.HadISST$Year, y = GMST.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMST.COBE_SST2$Year, y = GMST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMST.ERSSTv5$Year, y = GMST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMST.BEST_Land$Year, y = GMST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMST.hybrid36$Year, y = GMST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMST.tas_sea$Year, y = GMST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMST.tos_ua$Year, y = GMST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
    
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 03a. Plot GMSST filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMSST_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean sea surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$mod_p1_min_2.5, 
                                               rev(GMSST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$mod_p1_min_2.5, 
                                               rev(GMSST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$lp_2.5, 
                                               rev(GMSST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$lp_2.5, 
                                               rev(GMSST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$hp_2.5, 
                                               rev(GMSST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$hp_2.5, 
                                               rev(GMSST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMSST.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tas_land$res_lp_2.5[1:165], 
                                               rev(GMSST.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tos$res_lp_2.5[1:165], 
                                               rev(GMSST.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tas_land$res_hp_2.5[1:165], 
                                               rev(GMSST.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMSST.tos$res_hp_2.5[1:165], 
                                               rev(GMSST.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMSST.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 03b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMSST_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean sea surface temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$mod_p1_min_2.5, 
                                               rev(GMSST.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$mod_p1_min_2.5, 
                                               rev(GMSST.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1, lwd = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1, lwd = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1, lwd = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$lp_2.5, 
                                               rev(GMSST.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$lp_2.5, 
                                               rev(GMSST.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tas_land$hp_2.5, 
                                               rev(GMSST.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMSST.tos$hp_2.5, 
                                               rev(GMSST.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMSST.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMSST.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMSST.HadISST$Year, y = GMSST.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMSST.COBE_SST2$Year, y = GMSST.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMSST.ERSSTv5$Year, y = GMSST.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMSST.BEST_Land$Year, y = GMSST.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMSST.hybrid36$Year, y = GMSST.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMSST.tas_sea$Year, y = GMSST.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMSST.tos_ua$Year, y = GMSST.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMSST$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()




  



## ----------------------------------------------------------------------------------------
## 03a. Plot GMLSAT_NI filtered time series:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMLSAT_NI_land_vs_ocean_reconstruction_mod_p1.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(-3.6, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean land surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
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
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    text(x = 1850, y = -0.7 + 0.5*3, labels = "Forced response", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(-4, -2, 0.5) + 0.5*4, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(-4, -2, 0.1) + 0.5*4, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-3, -3) + 0.5*4, col = "darkgray", lty = 3)
    text(x = 1850, y = -2.7 + 0.5*4, labels = "Unforced, low-pass filtered (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    text(x = 1850, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$lp_2.5, 
                                               rev(GMLSAT_NI.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$lp_2.5, 
                                               rev(GMLSAT_NI.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    #lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$hp_2.5, 
                                               rev(GMLSAT_NI.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$hp_2.5, 
                                               rev(GMLSAT_NI.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Forced response:
  {
    c = -1 + 0.5*3; sc = 1
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$forced * sc + c, col = col.CRUTEM5, lwd = 3)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$forced * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$forced * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$forced * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$forced * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$forced * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$forced * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$forced * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$forced * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Unforced, LP:
  {
    c = -3 + 0.5*4; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tas_land$res_lp_2.5[1:165], 
                                               rev(GMLSAT_NI.tas_land$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tos$res_lp_2.5[1:165], 
                                               rev(GMLSAT_NI.tos$res_lp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$res_lp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$res_lp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$res_lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$res_lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$res_lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$res_lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$res_lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$res_lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # Unforced, HP:
  {
    c = -5 + 0.5*5; sc = 1
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tas_land$res_hp_2.5[1:165], 
                                               rev(GMLSAT_NI.tas_land$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2014, 2014:1850), y = c(GMLSAT_NI.tos$res_hp_2.5[1:165], 
                                               rev(GMLSAT_NI.tos$res_hp_97.5[1:165])) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$res_hp_50 * sc + c, col = col.CRUTEM5, lwd = 1)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$res_hp_50 * sc + c, col = col.HadSST4, lwd = 1)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$res_hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$res_hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$res_hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$res_hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$res_hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$res_hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()



## ----------------------------------------------------------------------------------------
## 03b. small version:
## ----------------------------------------------------------------------------------------

# mod_p1:
pdf(file = "01_GMLSAT_NI_land_vs_ocean_reconstruction_mod_p1_small.pdf", width = 8, height=5)
{
  par(mfrow=c(1, 1), mar=c(2,4,1,4))
  ylim = c(1, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "Global mean land surface air temperature anomaly [°C]", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  
  {
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 5.3, labels = "Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 3.3+0.5, labels = "Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.5 + 0.5*2, labels = "High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tas_land$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$mod_p1_min_2.5, 
                                               rev(GMLSAT_NI.tos$mod_p1_min_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$mod_p1_min_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$mod_p1_min_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$mod_p1_min_50 * sc + c, col = col.COBE2, lty = 1, lwd = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$mod_p1_min_50 * sc + c, col = col.ERSST5, lty = 1, lwd = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$mod_p1_min_50 * sc + c, col = col.BEST, lty = 1, lwd = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$mod_p1_min_50 * sc + c, col = col.hybrid36, lty = 1, lwd = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$mod_p1_min_50 * sc + c, col = col.CLASSNMAT, lty = 1, lwd = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$mod_p1_min_50 * sc + c, col = col.HadSST4_ua, lty = 1, lwd = 1)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$lp_2.5, 
                                               rev(GMLSAT_NI.tas_land$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$lp_2.5, 
                                               rev(GMLSAT_NI.tos$lp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$lp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$lp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$lp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$lp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$lp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$lp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$lp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$lp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 1
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tas_land$hp_2.5, 
                                               rev(GMLSAT_NI.tas_land$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.CRUTEM5, alpha = 30), border = make.transparent.color(col.CRUTEM5, 100))
    
    polygon(x = c(1850:2020, 2020:1850), y = c(GMLSAT_NI.tos$hp_2.5, 
                                               rev(GMLSAT_NI.tos$hp_97.5)) * sc + c, 
            col = make.transparent.color(col.HadSST4, alpha = 30), border = make.transparent.color(col.HadSST4, 75))
    
    # Lines for main reconstructions:
    lines(x = 1850:2020, y = GMLSAT_NI.tas_land$hp_50 * sc + c, col = col.CRUTEM5, lwd = 2)
    lines(x = 1850:2020, y = GMLSAT_NI.tos$hp_50 * sc + c, col = col.HadSST4, lwd = 2)
    
    # Lines for other observational datasets:
    # lines(x = GMLSAT_NI.HadISST$Year, y = GMLSAT_NI.HadISST$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    lines(x = GMLSAT_NI.COBE_SST2$Year, y = GMLSAT_NI.COBE_SST2$hp_50 * sc + c, col = col.COBE2, lty = 1)
    lines(x = GMLSAT_NI.ERSSTv5$Year, y = GMLSAT_NI.ERSSTv5$hp_50 * sc + c, col = col.ERSST5, lty = 1)
    
    lines(x = GMLSAT_NI.BEST_Land$Year, y = GMLSAT_NI.BEST_Land$hp_50 * sc + c, col = col.BEST, lty = 1)
    lines(x = GMLSAT_NI.hybrid36$Year, y = GMLSAT_NI.hybrid36$hp_50 * sc + c, col = col.hybrid36, lty = 1)
    lines(x = GMLSAT_NI.tas_sea$Year, y = GMLSAT_NI.tas_sea$hp_50 * sc + c, col = col.CLASSNMAT, lty = 1)
    
    lines(x = GMLSAT_NI.tos_ua$Year, y = GMLSAT_NI.tos_ua$hp_50 * sc + c, col = col.HadSST4_ua, lty = 1)
  }
  
  legend("top", c("CRUTEM5", "HadSST4", "HadSST4 unadj.", "COBE2-SST", "ERSSTv5", "ClassNMAT", "BEST-Land", "Cowtan-HybridSST"),
         lty = c(1, 1, 1, 1, 1, 1, 1, 1), lwd = c(2, 2, 1, 1, 1, 1, 1, 1), col = c(col.CRUTEM5, col.HadSST4, col.HadSST4_ua, col.COBE2, col.ERSST5, col.CLASSNMAT, col.BEST, col.hybrid36),  
         cex = 0.75, ncol = 3, bg = "white")
  
  ### 
  # cor(scale(pass.filt(y = OBS.tas_land$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = OBS.tos$GMLSAT_NI$ann$mod_p1_min, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  # cor(scale(pass.filt(y = tas_land.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165],
  #    scale(pass.filt(y = tos.mod$residuals, W = 20, type = "high", method = "Butterworth"), T, F)[1:165])
  
  
  
}
dev.off()




