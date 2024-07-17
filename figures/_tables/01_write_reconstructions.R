

# ------------------------------------------------------------------------------------
# Write final reconstructions into output .txt files
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.06.2024
library(matrixStats)


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")



## ----------------------------------------------------------------------------------------
## Write reconstructions into one table:
## ----------------------------------------------------------------------------------------

df.out = data.frame(matrix(data = NA, nrow = 171, ncol = 13))
colnames(df.out) <- c("Year", "GMST_CRUTEM5", "GMST_CRUTEM5_p2.5", "GMST_CRUTEM5_p97.5",
                      "GMST_HadSST4", "GMST_HadSST4_p2.5", "GMST_HadSST4_p97.5",
                      "GMST_HadSST4_unadj", "GMST_ClassNMAT", "GMST_CoastalHybridSST",
                      "GMST_ERSST5", "GMST_COBE_SST2", "GMST_BESTLand")
df.out$Year = GMST.tas_land$Year
df.out[,1:8] = cbind(GMST.tas_land$Year, 
                   GMST.tas_land$mod_p1_min_50,
                   GMST.tas_land$mod_p1_min_2.5,
                   GMST.tas_land$mod_p1_min_97.5,
                   GMST.tos$mod_p1_min_50,
                   GMST.tos$mod_p1_min_2.5,
                   GMST.tos$mod_p1_min_97.5,    
                   GMST.tos_ua$mod_p1_min_50)
df.out$GMST_ClassNMAT[match(GMST.tas_sea$Year, df.out$Year)] = GMST.tas_sea$mod_p1_min_50
df.out$GMST_CoastalHybridSST[match(GMST.hybrid36$Year, df.out$Year)] = GMST.hybrid36$mod_p1_min_50
df.out$GMST_ERSST5[match(GMST.ERSSTv5$Year, df.out$Year)] = GMST.ERSSTv5$mod_p1_min_50
df.out$GMST_COBE_SST2[match(GMST.COBE_SST2$Year, df.out$Year)] = GMST.COBE_SST2$mod_p1_min_50
df.out$GMST_BESTLand[match(GMST.BEST_Land$Year, df.out$Year)] = GMST.BEST_Land$mod_p1_min_50

write.table(x = round(df.out, 3), file = "data/04_final/GMST_reconstructions.txt", quote = F, row.names = F, col.names = T, sep=",")




# 08. Attribution and Filtering of time series for mod_p0:
# ------------------------------------------------------------------------------------
{
  library("dplR") # package for band-pass filtering
  
  # GMST Prepare data for plotting:
  GMST.CRUTEM5 = get.df(Y = CRUTEM5.global.annual$Anomaly[8:171], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = c(1850:2020)[8:171], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tas_land = get.df(Y = OBS.tas_land_$GMST_FM$ann$mod_p0, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
  GMST.tos = get.df(Y = OBS.tos_$GMST$ann$mod_p0, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2014)
  GMST.COBE_SST2 = get.df(Y = COBE_SST2.annual$GMSST_mod_p0_1se, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = COBE_SST2.annual$Year, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.ERSSTv5 = get.df(Y = ERSSTv5.annual$GMST_FM_mod_p0_1se[1:167], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = ERSSTv5.annual$Year[1:167], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.BEST_Land = get.df(Y = BEST_Land.annual$GMST_FM_mod_p0_1se[101:271], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = BEST_Land.annual$Year[101:271], center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.hybrid36 = get.df(Y = colMedians(OBS_hybrid36.tos_$GMST_FM$ann$mod_p0), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2016, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tas_sea = get.df(Y = colMedians(OBS_CLASSNMAT.tas_$GMST_FM$ann$mod_p0), f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1880:2019, center = T, ens.ix = NULL, years.DA = 1850:2014)
  GMST.tos_ua = get.df(Y = HadSST_ua.annual$GMST_FM_mod_p0_1se, f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2014)
}



df.out = data.frame(matrix(data = NA, nrow = 171, ncol = 13))
colnames(df.out) <- c("Year", "GMST_CRUTEM5", "GMST_CRUTEM5_p2.5", "GMST_CRUTEM5_p97.5",
                      "GMST_HadSST4", "GMST_HadSST4_p2.5", "GMST_HadSST4_p97.5",
                      "GMST_HadSST4_unadj", "GMST_ClassNMAT", "GMST_CoastalHybridSST",
                      "GMST_ERSST5", "GMST_COBE_SST2", "GMST_BESTLand")
df.out$Year = GMST.tas_land$Year
df.out[,1:8] = cbind(GMST.tas_land$Year, 
                     GMST.tas_land$mod_p1_min_50,
                     GMST.tas_land$mod_p1_min_2.5,
                     GMST.tas_land$mod_p1_min_97.5,
                     GMST.tos$mod_p1_min_50,
                     GMST.tos$mod_p1_min_2.5,
                     GMST.tos$mod_p1_min_97.5,    
                     GMST.tos_ua$mod_p1_min_50)
df.out$GMST_ClassNMAT[match(GMST.tas_sea$Year, df.out$Year)] = GMST.tas_sea$mod_p1_min_50
df.out$GMST_CoastalHybridSST[match(GMST.hybrid36$Year, df.out$Year)] = GMST.hybrid36$mod_p1_min_50
df.out$GMST_ERSST5[match(GMST.ERSSTv5$Year, df.out$Year)] = GMST.ERSSTv5$mod_p1_min_50
df.out$GMST_COBE_SST2[match(GMST.COBE_SST2$Year, df.out$Year)] = GMST.COBE_SST2$mod_p1_min_50
df.out$GMST_BESTLand[match(GMST.BEST_Land$Year, df.out$Year)] = GMST.BEST_Land$mod_p1_min_50

write.table(x = round(df.out, 3), file = "data/04_final/GMST_reconstructions_no-training-bias-unc.txt", quote = F, row.names = F, col.names = T, sep=",")




  
  