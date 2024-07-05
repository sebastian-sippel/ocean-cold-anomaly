
#### Read Palaeo temperature record:
# Sebastian Sippel
# 01.11.2022

source("code/_convenience/frenchcolormap.R")
library(R.matlab)

# setwd("/Users/sippels/Desktop/AbrametalPAGES2k2016_Supp2_inputdata/")
# setwd("/net/h2o/climphys1/sippels/_DATASET/Abram-et-al_PAGES2k/")
test = readMat("data/00_DATASET/paleoclimate/Abram-et-al_PAGES2k/archive_p2k_regional_reconstructions.mat")

# Western Atlantic:
WAtlantic = data.frame(cbind(test$archive.p2k.regional.reconstructions[,,5]$data[,,1]$year, test$archive.p2k.regional.reconstructions[,,5]$data[,,1]$rec, test$archive.p2k.regional.reconstructions[,,5]$data[,,1]$uncertainty))
names(WAtlantic) = c("year", "rec", "uncertainty")

# Western Pacific (best)
WPacific = data.frame(cbind(test$archive.p2k.regional.reconstructions[,,6]$data[,,1]$year, test$archive.p2k.regional.reconstructions[,,6]$data[,,1]$rec, test$archive.p2k.regional.reconstructions[,,6]$data[,,1]$uncertainty))
names(WPacific) = c("year", "rec", "uncertainty")

# Indian Ocean
IOcean = data.frame(cbind(test$archive.p2k.regional.reconstructions[,,7]$data[,,1]$year, test$archive.p2k.regional.reconstructions[,,7]$data[,,1]$rec, test$archive.p2k.regional.reconstructions[,,7]$data[,,1]$uncertainty))
names(IOcean) = c("year", "rec", "uncertainty")

# Indian Ocean, 25.5 × 106 km2; western Pacific, 26.9 × 106 km2; western Atlantic, 5.1 × 106 km2.
w_WAtlantic = 5.5; w_WPacific = 26.9; w_IOcean = 25.5;
w = c(w_WAtlantic, w_WPacific, w_IOcean)


# range of overlap years: 
ocean2k = data.frame(matrix(data = NA, nrow = length(1551:2009), ncol = 9))
colnames(ocean2k) = c("year", "WAtlantic", "WAtlantic_unc", "WPacific", "WPacific_unc", "IOcean", "IOcean_unc", "Tropics", "Tropics_unc")
ocean2k$year = 1551:2009
ocean2k$WAtlantic[match(x = WAtlantic$year, table = ocean2k$year)] = WAtlantic$rec
ocean2k$WAtlantic_unc[match(x = WAtlantic$year, table = ocean2k$year)] = WAtlantic$uncertainty
ocean2k$WPacific[match(x = WPacific$year, table = ocean2k$year)] = WPacific$rec
ocean2k$WPacific_unc[match(x = WPacific$year, table = ocean2k$year)] = WPacific$uncertainty
ocean2k$IOcean[match(x = IOcean$year, table = ocean2k$year)] = IOcean$rec
ocean2k$IOcean_unc[match(x = IOcean$year, table = ocean2k$year)] = IOcean$uncertainty

# range of overlap:
for (cur.year in 1621:2001) {
  year.ix = which(cur.year == ocean2k$year)
  ocean2k[year.ix,]$Tropics = weighted.mean(x = c(ocean2k[year.ix,]$WAtlantic, ocean2k[year.ix,]$WPacific, ocean2k[year.ix,]$IOcean), w = w)
}
# plot(ocean2k$year, ocean2k$Tropics)
# ocean2k$year[451]  # End: Year2001
ocean2k_ = get.df(Y = ocean2k$Tropics[match(x = 1850:2001, table = ocean2k$year)], f = rowMeans(CMIP6.GMST.f, na.rm=T), years = 1850:2001, center = T, ens.ix = NULL, years.DA = 1850:2001)




## ----------------------------------------------------------------------------------------------------------------
## Averages from April-March:
## ----------------------------------------------------------------------------------------------------------------
# OBS.tas_land$WAtlantic$mon$mod_p1_min
# x = OBS.tas_land_$WAtlantic$mon$mod_p1_min
get.ens.avg_Apr_Mar <- function(x, nrow = 200, ncol = 171) {
  ret.mat = matrix(data = NA, nrow = nrow, ncol = ncol)
  
  for ( i in 1:(ncol-1)) {
    ret.mat[,i] = rowMeans(sapply(X = 1:12, FUN=function(mon) {
      if(mon %in% 1:3) {
        x[[mon]][,i + 1]
      } else if (mon %in% 4:12) {
        x[[mon]][,i]
      }
    }), na.rm = T)
  }
  return(ret.mat)
}


## Get Tropics weighted average ensemble:
WAtlantic_tas_land_ens = get.ens.avg_Apr_Mar(x =  OBS.tas_land_$WAtlantic$mon$mod_p1_min)[,1:170]
WPacific_tas_land_ens = get.ens.avg_Apr_Mar(x =  OBS.tas_land_$WPacific$mon$mod_p1_min)[,1:170]
IndianOcean_tas_land_ens = get.ens.avg_Apr_Mar(x =  OBS.tas_land_$IndianOcean$mon$mod_p1_min)[,1:170]
EPacific_tas_land_ens = get.ens.avg_Apr_Mar(x =  OBS.tas_land_$EPacific$mon$mod_p1_min)[,1:170]
Tropics_tas_land_ens = matrix(NA, nrow = 200, ncol = 171)
WAtlantic_tos_ens = get.ens.avg_Apr_Mar(x =  OBS.tos_$WAtlantic$mon$mod_p1_min)[,1:170]
WPacific_tos_ens = get.ens.avg_Apr_Mar(x =  OBS.tos_$WPacific$mon$mod_p1_min)[,1:170]
# EPacific_tos_ens = get.ens.avg_Apr_Mar(x =  OBS.tos_$EPacific$mon$mod_p1_min)[,1:170]
IndianOcean_tos_ens = get.ens.avg_Apr_Mar(x =  OBS.tos_$IndianOcean$mon$mod_p1_min)[,1:170]
Tropics_tos_ens = matrix(NA, nrow = 200, ncol = 171)

for (i in 1:170) {
  print(i)
  Tropics_tas_land_ens[,i] = apply(X = rbind(WAtlantic_tas_land_ens[,i], WPacific_tas_land_ens[,i], IndianOcean_tas_land_ens[,i]), MARGIN = 2, FUN=function(x) weighted.mean(x = x, w = w))
  Tropics_tos_ens[,i] = apply(X = rbind(WAtlantic_tos_ens[,i], WPacific_tos_ens[,i], IndianOcean_tos_ens[,i]), MARGIN = 2, FUN=function(x) weighted.mean(x = x, w = w))
}
# plot(Tropics_ens[1,]); lines(Tropics_ens[2,])

Tropics.tas_land = get.df(Y = Tropics_tas_land_ens[,1:170], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2019, center = T, ens.ix = 1:200, years.DA = 1850:2014)
Tropics.tos = get.df(Y = Tropics_tos_ens[,1:170], f = rowMeans(CMIP6.GMSST.f, na.rm=T), years = 1850:2019, center = T, ens.ix = 1:200, years.DA = 1850:2014)


Tierney_best = readMat("data/00_DATASET/paleoclimate/Tierney-etal-ocean2k-Paleoceanography/data/rec_best.mat")
Tierney_regions = readMat("data/00_DATASET/paleoclimate/Tierney-etal-ocean2k-Paleoceanography/data/rec_alliters.mat")
# str(Tierney_regions$readme)


# str(Tierney_regions$atlantic[[7]])
# str(Tierney_best$atlantic[[5]])
# plot(x = WAtlantic$year, y = WAtlantic$rec, type='l')
# lines(x = Tierney_regions$atlantic[[5]], y = Tierney_best$atlantic[[2]], col = "red")




