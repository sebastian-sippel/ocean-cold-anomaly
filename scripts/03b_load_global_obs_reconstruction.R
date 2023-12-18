
# ------------------------------------------------------------------------------------
# Load tos and tas reconstruction(s) based on CMIP6
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 23.10.2022

# 00.(c) Load *new* reconstructions:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/OBS.tas_land_v4.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/OBS.tos_v4.RData")


# 00.(d) Derive reconstructions+percentiles:
# ---------------------------------------------------------------
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N")
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N", 
              "TMSST_40S_40N", "TMSST_25S_25N_", "IndianOcean", "WPacific", "EPacific", "WAtlantic")


YEARS.mon = seq(1850.042, 2020.958, length.out = 171*12)

# Observational predictions:
ens.ix = 94:200
OBS.tas_land = list(); OBS.tas_land[1:14] = sapply(1:14, FUN=function(i) OBS.tas_land[[i]] <- list() ); names(OBS.tas_land) = names.vec
OBS.tos = list(); OBS.tos[1:14] = sapply(1:14, FUN=function(i) OBS.tos[[i]] <- list() ); names(OBS.tos) = names.vec

for (i in 1:length(names.vec)) {
  cur.name = names.vec[i]
  print(cur.name)
  
  OBS.tas_land[[cur.name]]$mon = list(); OBS.tas_land[[cur.name]]$ann = list();
  
  OBS.tas_land[[cur.name]]$mon$mod_p0 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p0, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tas_land[[cur.name]]$mon$mod_p0_2.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p0, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tas_land[[cur.name]]$mon$mod_p0_97.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p0, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tas_land[[cur.name]]$mon$mod_p1 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tas_land[[cur.name]]$mon$mod_p1_2.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tas_land[[cur.name]]$mon$mod_p1_97.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tas_land[[cur.name]]$mon$mod_p1_min = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tas_land[[cur.name]]$mon$mod_p1_min_2.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tas_land[[cur.name]]$mon$mod_p1_min_97.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tas_land[[cur.name]]$mon$mod_gta = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_gta, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tas_land[[cur.name]]$mon$mod_gta_2.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_gta, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tas_land[[cur.name]]$mon$mod_gta_97.5 = c(t(sapply(X = OBS.tas_land_[[cur.name]]$mon$mod_gta, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tas_land[[cur.name]]$mon$YEARS = YEARS.mon
    
  OBS.tas_land[[cur.name]]$ann$mod_p0 = colMedians(OBS.tas_land_[[cur.name]]$ann$mod_p0[ens.ix,])
  OBS.tas_land[[cur.name]]$ann$mod_p0_2.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p0[ens.ix,], probs = 0.025)
  OBS.tas_land[[cur.name]]$ann$mod_p0_97.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p0[ens.ix,], probs = 0.975)
  OBS.tas_land[[cur.name]]$ann$mod_p1 = colMedians(OBS.tas_land_[[cur.name]]$ann$mod_p1[ens.ix,])
  OBS.tas_land[[cur.name]]$ann$mod_p1_2.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p1[ens.ix,], probs = 0.025)
  OBS.tas_land[[cur.name]]$ann$mod_p1_97.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p1[ens.ix,], probs = 0.975)
  OBS.tas_land[[cur.name]]$ann$mod_p1_min = colMedians(OBS.tas_land_[[cur.name]]$ann$mod_p1_min[ens.ix,])
  OBS.tas_land[[cur.name]]$ann$mod_p1_min_2.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p1_min[ens.ix,], probs = 0.025)
  OBS.tas_land[[cur.name]]$ann$mod_p1_min_97.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_p1_min[ens.ix,], probs = 0.975)
  OBS.tas_land[[cur.name]]$ann$mod_gta = colMedians(OBS.tas_land_[[cur.name]]$ann$mod_gta[ens.ix,])
  OBS.tas_land[[cur.name]]$ann$mod_gta_2.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_gta[ens.ix,], probs = 0.025)
  OBS.tas_land[[cur.name]]$ann$mod_gta_97.5 = colQuantiles(OBS.tas_land_[[cur.name]]$ann$mod_gta[ens.ix,], probs = 0.975)
  OBS.tas_land[[cur.name]]$ann$YEARS = 1850:2020
  
  
  
  OBS.tos[[cur.name]]$mon = list(); OBS.tos[[cur.name]]$ann = list();
  
  OBS.tos[[cur.name]]$mon$mod_p0 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p0, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tos[[cur.name]]$mon$mod_p0_2.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p0, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tos[[cur.name]]$mon$mod_p0_97.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p0, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tos[[cur.name]]$mon$mod_p1 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tos[[cur.name]]$mon$mod_p1_2.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tos[[cur.name]]$mon$mod_p1_97.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tos[[cur.name]]$mon$mod_p1_min = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tos[[cur.name]]$mon$mod_p1_min_2.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tos[[cur.name]]$mon$mod_p1_min_97.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_p1_min, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tos[[cur.name]]$mon$mod_gta = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_gta, FUN=function(x) colMedians(x[ens.ix,]))))
  OBS.tos[[cur.name]]$mon$mod_gta_2.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_gta, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.025))))
  OBS.tos[[cur.name]]$mon$mod_gta_97.5 = c(t(sapply(X = OBS.tos_[[cur.name]]$mon$mod_gta, FUN=function(x) colQuantiles(x[ens.ix,], probs = 0.975))))
  OBS.tos[[cur.name]]$mon$YEARS = YEARS.mon
  
  OBS.tos[[cur.name]]$ann$mod_p0 = colMedians(OBS.tos_[[cur.name]]$ann$mod_p0[ens.ix,])
  OBS.tos[[cur.name]]$ann$mod_p0_2.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p0[ens.ix,], probs = 0.025)
  OBS.tos[[cur.name]]$ann$mod_p0_97.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p0[ens.ix,], probs = 0.975)
  OBS.tos[[cur.name]]$ann$mod_p1 = colMedians(OBS.tos_[[cur.name]]$ann$mod_p1[ens.ix,])
  OBS.tos[[cur.name]]$ann$mod_p1_2.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p1[ens.ix,], probs = 0.025)
  OBS.tos[[cur.name]]$ann$mod_p1_97.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p1[ens.ix,], probs = 0.975)
  OBS.tos[[cur.name]]$ann$mod_p1_min = colMedians(OBS.tos_[[cur.name]]$ann$mod_p1_min[ens.ix,])
  OBS.tos[[cur.name]]$ann$mod_p1_min_2.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p1_min[ens.ix,], probs = 0.025)
  OBS.tos[[cur.name]]$ann$mod_p1_min_97.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_p1_min[ens.ix,], probs = 0.975)
  OBS.tos[[cur.name]]$ann$mod_gta = colMedians(OBS.tos_[[cur.name]]$ann$mod_gta[ens.ix,])
  OBS.tos[[cur.name]]$ann$mod_gta_2.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_gta[ens.ix,], probs = 0.025)
  OBS.tos[[cur.name]]$ann$mod_gta_97.5 = colQuantiles(OBS.tos_[[cur.name]]$ann$mod_gta[ens.ix,], probs = 0.975)
  OBS.tos[[cur.name]]$ann$YEARS = 1850:2020
}





# 00.(e) Derive blended reconstructions:
# ---------------------------------------------------------------
library(matrixStats)

# ocean-land blended "best estimate":
land_fraction = 0.29 # https://www.nationsonline.org/oneworld/earth.htm#:~:text=Surface%3A,total%20surface%20of%20the%20Earth.
ice_fraction =  0.04     # https://nsidc.org/learn/parts-cryosphere/sea-ice
land_ice_fraction = 0.33
  # Arctic average= (15.5+6.5) / 2 = 11 million square kilometers
  # Antarctic coverage = (18.5+2.5)/2 = 10.5 million square kilometers
  # -> 4% average sea ice. 21.5 * 10^6 / 509600000 # ( (2*6378)^2 * pi )
sea_fraction = 1 - land_ice_fraction

GSAT_blended_mon_mod_p1 = sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tos_$GMMSAT$mon$mod_p1_min[[mon]])) * sea_fraction + sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tas_land_$GMLSAT_MI$mon$mod_p1_min[[mon]])) * land_ice_fraction
GSAT_blended_mon_mod_p0 = sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tos_$GMMSAT$mon$mod_p0_min[[mon]])) * sea_fraction + sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tas_land_$GMLSAT_MI$mon$mod_p0_min[[mon]])) * land_ice_fraction
GSAT_blended_ann_mod_p1 = colMedians(OBS.tos_$GMMSAT$ann$mod_p1_min) * sea_fraction + colMedians(OBS.tas_land_$GMLSAT_MI$ann$mod_p1_min) * land_ice_fraction
GSAT_blended_ann_mod_p0 = colMedians(OBS.tos_$GMMSAT$ann$mod_p0_min) * sea_fraction + colMedians(OBS.tas_land_$GMLSAT_MI$ann$mod_p0_min) * land_ice_fraction

GMST_blended_mon_mod_p1 = sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tos_$GMSST$mon$mod_p1_min[[mon]])) * sea_fraction + sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tas_land_$GMLSAT_MI$mon$mod_p1_min[[mon]])) * land_ice_fraction
GMST_blended_mon_mod_p0 = sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tos_$GMSST$mon$mod_p1_min[[mon]])) * sea_fraction + sapply(X = 1:12, FUN = function(mon) colMedians(OBS.tas_land_$GMLSAT_MI$mon$mod_p1_min[[mon]])) * land_ice_fraction
GMST_blended_ann_mod_p1 = colMedians(OBS.tos_$GMSST$ann$mod_p1_min) * sea_fraction + colMedians(OBS.tas_land_$GMLSAT_MI$ann$mod_p1_min) * land_ice_fraction
GMST_blended_ann_mod_p0 = colMedians(OBS.tos_$GMSST$ann$mod_p1_min) * sea_fraction + colMedians(OBS.tas_land_$GMLSAT_MI$ann$mod_p1_min) * land_ice_fraction



# GMST.tas_land.mon.blended = sapply(X = 1:12, FUN = function(mon) colMedians(GMSST.tas_land_$mon$mod_p1[[mon]])) * 67/100 + sapply(X = 1:12, FUN = function(mon) colMedians(GMLSAT.tas_land_$mon$mod_p1[[mon]])) * 33/100
# GMST.tas_land.ann.blended_mod_p1 = colMeans(GMSST.tas_land_$ann$mod_p1) * 67/100 + colMeans(GMLSAT.tas_land_$ann$mod_p1) * 33/100
# GMST.tas_land.ann.blended_mod_p0 = colMeans(GMSST.tas_land_$ann$mod_p0) * 67/100 + colMeans(GMLSAT.tas_land_$ann$mod_p0) * 33/100
# GMST.tos.mon.blended = sapply(X = 1:12, FUN = function(mon) colMedians(GMSST.tos_$mon$mod_p1[[mon]])) * 67/100 + sapply(X = 1:12, FUN = function(mon) colMedians(GMLSAT.tos_$mon$mod_p1[[mon]])) * 33/100
# GMST.tos.ann.blended_mod_p1 = colMeans(GMSST.tos_$ann$mod_p1) * 67/100 + colMeans(GMLSAT.tos_$ann$mod_p1) * 33/100
# GMST.tos.ann.blended_mod_p0 = colMeans(GMSST.tos_$ann$mod_p1) * 67/100 + colMeans(GMLSAT.tos_$ann$mod_p0) * 33/100



# (1) Load other reconstructed datasets:
# ---------------------------------------------------------------

# (A) ClassNMAT reconstruction:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/OBS_CLASSNMAT.tas.RData")

# (B) Cowtan Coastal hybrid reconstruction:
# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/03_processed_small_data/AGMT.tos_hybrid36_.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/OBS_hybrid36.tos.RData")


# (2) Load other variability reconstruction datasets:
# ---------------------------------------------------------------

# (A) PSL reconstruction:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/_old/GSAT.HadSLP2_psl_v2.RData")
GSAT.psl_ = colMedians(GSAT.psl_$ann$mod_p0_1_05)

# (B) Mean Removed reconstruction:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/03_processed_small_data/AGMT.tos_MR_.RData")




