
# ------------------------------------------------------------------------------------
# Calculate trends for comparisons and constraints based on proxy data:
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.02.2023
library(matrixStats)
library(dplR)
# requires to run 04a_master_load_reconstructions.R and 04b_master_read_paleo_reconstructions.R beforehand!


# 00. Define parameters for trend calculation
# ------------------------------------------------------------------------------------
probs = c(0.025, 0.5, 0.975) # Confidence intervals 
trend.years = list(1851:1890, 1871:1910, 1901:1940, 1975:2014, 1851:1900, 1871:1920, 1901:1950, 1965:2014)


# 01. Read HadSST4 raw regional estimates for comparison with trends from reconstructions (Figure 4):
# ------------------------------------------------------------------------------------

source("code/_convenience/convert.to.eurocentric.R")
raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw = c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]))
land_fraction = c(matrix(values(convert.to.pacificcentric(raster("/net/h2o/climphys1/sippels/_DATA/grid/5d00_static/cmip5_masks/sftlf_g025.nc"))), 72, 36)[,36:1])
lon = c(matrix(coordinates(raster.template)[,1], 72, 36)[,36:1])
lat = c(matrix(coordinates(raster.template)[,2], 72, 36)[,36:1])


# 01a. HadSST4 raw estimates for ocean regions:
# ------------------------------------------------------------------------------------
{
  IndianOcean_ <- WAtlantic_ <- WPacific_ <- EPacific_ <- list()
  
  for (mon in 1:12) {
    
    print(mon)
    
    load(paste("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/01_processed4train_CRU/HadSST4_mon", mon, ".RData", sep=""))
    
    # Indian Ocean (20◦N–15◦S, 40–100◦E):
    grid.ix = which(lat < 20 & lat > -15 & lon > 40 & lon < 100)
    w = areaw[grid.ix] * (1-land_fraction[grid.ix]) / sum(areaw[grid.ix] * (1-land_fraction[grid.ix]))  
    IndianOcean_[[mon]] = sapply(X = 1:dim(HadSST4_)[1], FUN=function(i) sapply(1:200, FUN=function(ens) weighted.mean(x = (HadSST4_[i,grid.ix] + HadSST4_ENS_anom_[[ens]][i,grid.ix]), w = w, na.rm = T)))
    
    # Western Pacific (25◦N–25◦S, 110–155◦E)
    grid.ix = which(lat < 25 & lat > -25 & lon > 110 & lon < 155)
    w = areaw[grid.ix] * (1-land_fraction[grid.ix]) / sum(areaw[grid.ix] * (1-land_fraction[grid.ix]))
    WPacific_[[mon]] = sapply(X = 1:dim(HadSST4_)[1], FUN=function(i) sapply(1:200, FUN=function(ens) weighted.mean(x = (HadSST4_[i,grid.ix] + HadSST4_ENS_anom_[[ens]][i,grid.ix]), w = w, na.rm = T)))
    
    # Eastern Pacific (10◦N–10◦S, 175◦E–85◦W)
    grid.ix = which(lat < 10 & lat > -10 & lon > 175 & lon < 275)
    # land_fraction1 = land_fraction; land_fraction1[grid.ix] = 5; image.plot(matrix(land_fraction1, 72, 36)); 
    w = areaw[grid.ix] * (1-land_fraction[grid.ix]) / sum(areaw[grid.ix] * (1-land_fraction[grid.ix]))
    EPacific_[[mon]] = sapply(X = 1:dim(HadSST4_)[1], FUN=function(i) sapply(1:200, FUN=function(ens) weighted.mean(x = (HadSST4_[i,grid.ix] + HadSST4_ENS_anom_[[ens]][i,grid.ix]), w = w, na.rm = T)))
    
    # Western Atlantic (15–30◦N, 60–90◦W)
    grid.ix = which(lat < 30 & lat > 15 & lon > 270 & lon < 300)
    # land_fraction1 = land_fraction; land_fraction1[grid.ix] = 5; image.plot(matrix(land_fraction1, 72, 36)); 
    w = areaw[grid.ix] * (1-land_fraction[grid.ix]) / sum(areaw[grid.ix] * (1-land_fraction[grid.ix]))
    WAtlantic_[[mon]] = sapply(X = 1:dim(HadSST4_)[1], FUN=function(i) sapply(1:200, FUN=function(ens) weighted.mean(x = (HadSST4_[i,grid.ix] + HadSST4_ENS_anom_[[ens]][i,grid.ix]), w = w, na.rm = T)))
  }
  
  ## Process into ANNUAL x ENSEMBLE format:
  IndianOcean_ann_ens = get.ens.avg_Apr_Mar(x = IndianOcean_)
  WAtlantic_ann_ens = get.ens.avg_Apr_Mar(x = WAtlantic_)
  WPacific_ann_ens = get.ens.avg_Apr_Mar(x = WPacific_)
  EPacific_ann_ens = get.ens.avg_Apr_Mar(x = EPacific_)
}
# get Tropics full ensemble with the same weights as Abram et al:
# str(IndianOcean_ann_ens)

Tropics_ann_ens = t(sapply(X = 1:200, FUN = function(ix) {
  return(sapply(X = 1:170, FUN=function(i) {
    weighted.mean(x = c(IndianOcean_ann_ens[ix,i], WAtlantic_ann_ens[ix,i], WPacific_ann_ens[ix,i]), w = c(w_IOcean, w_WAtlantic, w_WPacific), na.rm = T)
  }))
}))



# 03. Calculate trends in instrumental observations (Figure 4 and Figure 5): 
# ------------------------------------------------------------------------------------

## Correction of HadSST4 dataset with hybrid36-Cowtan et al. in the period 1891-1950:
hybrid.cor = rep(0, 171)
hybrid.cor[41:100] = -GMSST.tos$res_lp_50[41:100] + GMSST.hybrid36$res_lp_50[41:100]
OBS.tos_GMSST_cor = OBS.tos_$GMSST$ann$mod_p1_min + rep.row(x = hybrid.cor, n = 200)
hybrid.cor[1:167] = -GMSST.tos$res_lp_50[1:167] + GMSST.hybrid36$res_lp_50[1:167]
OBS.tos_GMSST_cor = OBS.tos_$GMSST$ann$mod_p1_min + rep.row(x = hybrid.cor, n = 200)


# Calculate trends in OBS from different reconstructions:
OBS.GMST_tos.trends = apply(X = OBS.tos_$GMST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMST_tas_land.trends = apply(X = OBS.tas_land_$GMST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)

OBS.GMSST_tos.trends = apply(X = OBS.tos_$GMSST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMSST_tos.trends_q = apply(X = OBS.GMSST_tos.trends, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

OBS.GMSST_tas_land.trends = apply(X = OBS.tas_land_$GMSST$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMSST_tas_land.trends_q = apply(X = OBS.GMSST_tas_land.trends, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

OBS.GMLSAT_tos.trends = apply(X = OBS.tos_$GMLSAT_NI$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMLSAT_tas_land.trends = apply(X = OBS.tas_land_$GMLSAT_NI$ann$mod_p1_min, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMLSAT_tas_land.trends_q = apply(X = OBS.GMLSAT_tas_land.trends, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

#OBS.GMSST_hybrid36.trends = get.trend(x = hybrid36.annual$Anomaly, trend.years = trend.years, years = 1850:2020)
OBS.GMSST_hybrid36cor.trends = apply(X = OBS.tos_GMSST_cor, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.GMSST_hybrid36cor.trends_q = apply(X = OBS.GMSST_hybrid36cor.trends, MARGIN = 1, FUN=function(x) quantile(x, probs=probs))

# Reconstructions processed to Apr.-Mar year for paleo-comparison:
OBS.Tropics_tas_land.trends = apply(X = Tropics_tas_land_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.Tropics_tos.trends = apply(X = Tropics_tos_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)

OBS.WAtlantic_tas_land.trends = apply(X = WAtlantic_tas_land_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.WAtlantic_tos.trends = apply(X = WAtlantic_tos_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)

OBS.WPacific_tas_land.trends = apply(X = WPacific_tas_land_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.WPacific_tos.trends = apply(X = WPacific_tos_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)

OBS.IndianOcean_tas_land.trends = apply(X = IndianOcean_tas_land_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
OBS.IndianOcean_tos.trends = apply(X = IndianOcean_tos_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)

OBS.EPacific_tas_land.trends = apply(X = EPacific_tas_land_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
#OBS.EPacific_tos.trends = apply(X = EPacific_tos_ens, MARGIN = 1, FUN=get.trend, trend.years = trend.years, years = 1850:2020)

# HadSST4 regional averages:
HadSST4.Tropics.trends = apply(X = Tropics_ann_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
HadSST4.WAtlantic.trends = apply(X = WAtlantic_ann_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
HadSST4.WPacific.trends = apply(X = WPacific_ann_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
HadSST4.IndianOcean.trends = apply(X = IndianOcean_ann_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)
HadSST4.EPacific.trends = apply(X = EPacific_ann_ens, MARGIN = 1, FUN=get.trend_perioddiff, trend.years = trend.years, years = 1850:2020)


## Trends in SST datasets:
# Calculate missing trends:
HadSST4.trends = get.trend_perioddiff(x = HadSST4.global.annual$Anomaly, trend.years = trend.years, years = 1850:2020)
COBE_SST2.trends = get.trend_perioddiff(x = COBE_SST2.annual$Anomaly, trend.years = trend.years, years = COBE_SST2.annual$Year)
ERSSTv5.trends = get.trend_perioddiff(x = ERSSTv5.annual$Anomaly, trend.years = trend.years, years = ERSSTv5.annual$Year)


# ------------------------------------------------------------------------------------
# 04. Trends in Paleo-reconstructions:
# ------------------------------------------------------------------------------------

# Ocean2k
ocean2k.trends = get.trend_perioddiff(x = ocean2k_$mod_p1_min_50, trend.years = trend.years, years = 1850:2001)

## Tierney et al. 2015:
Tierney_WAtlantic = t(sapply(X = 1:70, FUN=function(i) get.trend_perioddiff(x = Tierney_regions$atlantic[[7]][,i], trend.years =  trend.years, years = c(Tierney_regions$atlantic[[5]]))))
Tierney_IOcean = t(sapply(X = 1:98, FUN=function(i) get.trend_perioddiff(x = Tierney_regions$indian[[7]][,i], trend.years =  trend.years, years = c(Tierney_regions$indian[[5]]))))
Tierney_WPacific = t(sapply(X = 1:133, FUN=function(i) get.trend_perioddiff(x = Tierney_regions$wpacific[[7]][,i], trend.years =  trend.years, years = c(Tierney_regions$wpacific[[5]]))))
Tierney_EPacific = t(sapply(X = 1:63, FUN=function(i) get.trend_perioddiff(x = Tierney_regions$epacific[[7]][,i], trend.years =  trend.years, years = c(Tierney_regions$epacific[[5]]))))


### Neukom reconstructions:
neukom2019_trend_all = list(sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_BHM[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_CPS_new[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_DA[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_M08[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_OIE[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_PAI[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)),
                         sapply(2:1001, FUN=function(i) get.trend_perioddiff(x = neukom2019_PCR[,i], trend.years = (trend.years), years = neukom2019_BHM$Year_CE)))

# Crowley trends:
crowley2014.trend = get.trend_perioddiff(x = proxy.Crowley2014$`Global proxy (1782)`, trend.years = trend.years, years = proxy.Crowley2014$Year)






# ------------------------------------------------------------------------------------
# 05. Trends in CMIP6:
# ------------------------------------------------------------------------------------
w_WAtlantic = 5.5; w_WPacific = 26.9; w_IOcean = 25.5;  # -> same weights as Ocean2k.

CMIP6.Tropics = cbind(CMIP6.tos_all.df$ann$Y$IndianOcean, CMIP6.tos_all.df$ann$Y$WPacific, CMIP6.tos_all.df$ann$Y$WAtlantic)
CMIP6.Tropics_ = apply(X = CMIP6.Tropics, MARGIN = 1, FUN=function(x) weighted.mean(x = x, w = c(w_IOcean, w_WPacific, w_WAtlantic)))

{
  
  # select CMIP6 historical members:
  CMIP6.tas_land.all = data.frame(cbind(all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  CMIP6.tos.all = data.frame(cbind(all = paste(CMIP6.tos_all.df$M$mod, "_", CMIP6.tos_all.df$M$scen, "_", CMIP6.tos_all.df$M$ens.mem, sep="")))
  ens.mem = data.frame(cbind(mod=CMIP6.tas_land_all.df$M$mod, scen = CMIP6.tas_land_all.df$M$scen, ens.mem = CMIP6.tas_land_all.df$M$ens.mem, 
                             all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  ens.mem.un = unique(ens.mem)
  ens.mem.un = ens.mem.un[which(ens.mem.un$scen == "historical"),]
  # remove all ensemble members that contain NA's:
  na.mems = unique(CMIP6.tos.all$all[which(is.na(CMIP6.tos_all.df$ann$Yhat$GMST_FM))])
  omit.ix=na.omit(match(x = na.mems, table = ens.mem.un$all))
  if (length(omit.ix) > 0)  ens.mem.un = ens.mem.un[-omit.ix,]
  
  # Get trends for each ensemble member:
  CMIP6.TMLSAT_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.TMMSAT_tos.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMLSAT_NI_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMSST_tos.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMSST_tas_land.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMSST_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMLSAT_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.Tropics_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  CMIP6.GMST_true.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 9))
  
  
  for (en in 1:dim(ens.mem.un)[1]) {
    print(en)
    
    ix_tas_land = which(CMIP6.tas_land.all$all == ens.mem.un$all[en] & CMIP6.tas_land_all.df$M$year %in% 1850:2014)
    ix_tos = which(CMIP6.tos.all$all == ens.mem.un$all[en] & CMIP6.tos_all.df$M$year %in% 1850:2014)
    
    #CMIP6.tas_land_hist_mod_p1[en,] = CMIP6.tas_land_all.df$ann$Yhat$GSAT[ix_tas_land]
    #CMIP6.tos_hist_mod_p1[en,] = CMIP6.tos_all.df$ann$Yhat$GSAT[ix_tos]
    
    #CMIP6.tas_land_hist_mod_p1_pt1[en,] = CMIP6.tas_land_all.df$ann$Yhat_pt1$GSAT[ix_tas_land]
    #CMIP6.tos_hist_mod_p1_pt1[en,] = CMIP6.tos_all.df$ann$Yhat_pt1$GSAT[ix_tos]
    
    # Get trends: 
    # trend.years = list(1851:1890, 1871:1910, 1901:1940, 1975:2014)
    # CMIP6.TMLSAT_tas_land.trends[en,] = get.trend_(x = CMIP6.tas_land_all.df$ann$Yhat$TMLSAT_40S_40N[ix_tas_land], trend.length = 50, years = 1850:2014)
    CMIP6.TMLSAT_tas_land.trends[en,] = get.trend_perioddiff(x = CMIP6.tas_land_all.df$ann$Yhat$TMLSAT_40S_40N[ix_tas_land], trend.years = trend.years, years = 1850:2014)
    CMIP6.TMMSAT_tos.trends[en,] = get.trend_perioddiff(x = CMIP6.tos_all.df$ann$Yhat$TMMSAT_40S_40N[ix_tos], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMLSAT_NI_tas_land.trends[en,] = get.trend_perioddiff(x = CMIP6.tas_land_all.df$ann$Yhat$GMLSAT_NI[ix_tas_land], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMSST_tos.trends[en,] = get.trend_perioddiff(x = CMIP6.tos_all.df$ann$Yhat$GMSST[ix_tos], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMSST_tas_land.trends[en,] = get.trend_perioddiff(x = CMIP6.tas_land_all.df$ann$Yhat$GMSST[ix_tas_land], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMSST_true.trends[en,] = get.trend_perioddiff(x = CMIP6.tos_all.df$ann$Y$GMSST[ix_tos], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMLSAT_true.trends[en,] = get.trend_perioddiff(x = CMIP6.tos_all.df$ann$Y$GMLSAT_NI[ix_tas_land], trend.years = trend.years, years = 1850:2014)
    CMIP6.Tropics_true.trends[en,] = get.trend_perioddiff(x = CMIP6.Tropics_[ix_tas_land], trend.years = trend.years, years = 1850:2014)
    CMIP6.GMST_true.trends[en,] = get.trend_perioddiff(x = CMIP6.tos_all.df$ann$Y$GMST_FM[ix_tos], trend.years = trend.years, years = 1850:2014)
  }
}



