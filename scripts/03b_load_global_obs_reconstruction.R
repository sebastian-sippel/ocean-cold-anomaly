
# ------------------------------------------------------------------------------------
# Load tos and tas reconstruction(s) based on CMIP6
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 23.10.2022

# 00.(c) Load *new* reconstructions based on HadSST4 and CRUTEM5:
load("data/03_processedOBS_reconstr/OBS.tas_land_v5.RData")
load("data/03_processedOBS_reconstr/OBS.tos.RData")


# 00.(d) Derive reconstructions+percentiles:
# ---------------------------------------------------------------
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N")
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N", 
              "TMSST_40S_40N", "TMSST_25S_25N_", "IndianOcean", "WPacific", "EPacific", "WAtlantic")


YEARS.mon = seq(1850.042, 2020.958, length.out = 171*12)

# Observational predictions:
ens.ix = 1:200
OBS.tas_land = list(); OBS.tas_land[1:14] = sapply(1:14, FUN=function(i) OBS.tas_land[[i]] <- list() ); names(OBS.tas_land) = names.vec

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
}
  
# tos
names.vec = c("GSAT", "GMSST", "GMLSAT_NI", "GMST_FM", "GMMSAT", "IndianOcean", "WPacific", "WAtlantic")

OBS.tos = list(); OBS.tos[1:8] = sapply(1:8, FUN=function(i) OBS.tos[[i]] <- list() ); names(OBS.tos) = names.vec
  
for (i in 1:length(names.vec)) {
  cur.name = names.vec[i]
  print(cur.name)
  
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






# (1) Load other reconstructed datasets:
# ---------------------------------------------------------------

# (A) ClassNMAT reconstruction:
load("data/03_processedOBS_reconstr/OBS_CLASSNMAT.tas.RData")

# (B) Cowtan Coastal hybrid reconstruction:
load("data/03_processedOBS_reconstr/OBS_hybrid36.tos.RData")


# (2) Load other variability reconstruction datasets:
# ---------------------------------------------------------------

# (A) Mean Removed reconstruction:
load("data/03_processedOBS_reconstr/OBS.tos_MR.RData")




