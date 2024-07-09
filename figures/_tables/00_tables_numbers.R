
# ------------------------------------------------------------------------------------
# Tables and Numbers for Manuscript
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.06.2024
library(matrixStats)


## load all data for reconstructions:
source("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/scripts/04a_master_load_reconstructions.R")





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
cor(GMST.tas_land$mod_p1_min_50[1:51], GMST.tos$mod_p1_min_50[1:51]) # 0.71

### MONTHLY:
cor(CRUTEM5.global.monthly$Anomaly[1:(51*12)], HadSST4.global.monthly$Anomaly[1:(51*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[(0*12):(51*12)], OBS.tos$GMST_FM$mon$mod_p1[(0*12):(51*12)])  # 0.37
cor(OBS.tas_land$GMST_FM$mon$mod_p0[1:(51*12)], OBS.tos$GMST_FM$mon$mod_p0[1:(51*12)])  # 0.29
cor(OBS.tas_land$GMST_FM$mon$mod_gta[1:(51*12)], OBS.tos$GMST_FM$mon$mod_gta[1:(51*12)])  # 0.30

# 1850-2020
cor(CRUTEM5.global.monthly$Anomaly[1:(170*12)], HadSST4.global.monthly$Anomaly[1:(170*12)], use = "complete.obs")
cor(OBS.tas_land$GMST_FM$mon$mod_p1_min[1:(170*12)], OBS.tos$GMST_FM$mon$mod_p1_min[1:(170*12)])








########################################################################
## Differences in OBS periods:
########################################################################


# calculate period differences:
period.mean.years = list(1871:1890, 1901:1920, 1901:1930)

OBS.GMST_tos_mod_p1.pm = apply(X = OBS.tos_$GMST$ann$mod_p1_min, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)
OBS.GMST_tos_mod_p0.pm = apply(X = OBS.tos_$GMST$ann$mod_p0, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)

OBS.GMST_tas_land_mod_p1.pm = apply(X = OBS.tas_land_$GMST$ann$mod_p1_min, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)
OBS.GMST_tas_land_mod_p0.pm = apply(X = OBS.tas_land_$GMST$ann$mod_p0, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)

mean(OBS.GMST_tos_mod_p1.pm[2,] - OBS.GMST_tas_land_mod_p1.pm[2,])
mean(OBS.GMST_tos_mod_p0.pm[2,] - OBS.GMST_tas_land_mod_p0.pm[2,])

mean(OBS.GMST_tos_mod_p1.pm[3,] - OBS.GMST_tas_land_mod_p1.pm[3,])
mean(OBS.GMST_tos_mod_p0.pm[3,] - OBS.GMST_tas_land_mod_p0.pm[3,])


CMIP6.tos_hist_mod_p1.pm = apply(X = CMIP6_matrix$tos_mod_p1, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)
CMIP6.tas_land_hist_mod_p1.pm = apply(X = CMIP6_matrix$tas_land_mod_p1, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)
CMIP6.tos_hist_mod_p0.pm = apply(X = CMIP6_matrix$tos_mod_p0, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)
CMIP6.tas_land_hist_mod_p0.pm = apply(X = CMIP6_matrix$tas_land_mod_p0, MARGIN = 1, FUN=get.period.mean, period.years = period.mean.years, years = 1850:2020)



ix = 3
# mean(OBS.GMST_tos_mod_p1.pm[2,] - OBS.GMST_tas_land_mod_p1.pm[2,])
# min(OBS.GMST_tos_mod_p1.pm[2,] - OBS.GMST_tas_land_mod_p1.pm[2,])
quantile(OBS.GMST_tos_mod_p1.pm[ix,] - OBS.GMST_tas_land_mod_p1.pm[ix,], probs = c(0.95))
# mean(CMIP6.tos_hist_mod_p1.pm[2,] - CMIP6.tas_land_hist_mod_p1.pm[2,])
# quantile(CMIP6.tos_hist_mod_p1.pm[2,] - CMIP6.tas_land_hist_mod_p1.pm[2,], probs = c(0.975))
quantile(CMIP6.tos_hist_mod_p1.pm[ix,] - CMIP6.tas_land_hist_mod_p1.pm[ix,], probs = c(0.05))




########################################################################
## Differences in OBS periods:
########################################################################

## make table for ens.mem.un:
train.table = read.table(file = "/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/CMIP6-table-train.txt", header = F)
train.table[,3] = round(train.table[,3] / 165)
train.test = train.table
train.test[,4] <- rep(NA, 64)

for (m in 1:length(train.table[,1])) {
  print(train.table[m,1])
  train.test[m,4] = length(which(ens.mem.un$mod == train.table[m,1]))
}

write.table(x = train.test, file = "/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/data/CMIP6-table-test.txt", 
            quote = F, row.names = F, col.names = F, sep = "&")

