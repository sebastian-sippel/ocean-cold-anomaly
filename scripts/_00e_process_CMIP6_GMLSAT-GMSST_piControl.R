
# ------------------------------------------------------------------------------------
# PROCESS GMLSAT-GMSST-blended temperatures for CMIP6
# ------------------------------------------------------------------------------------
## Needs to be run all on the same machine (n2o or CO2, not R-Studio Server, because list.files() produces different order).


# Sebastian Sippel
# 16.08.2022
# run on n2o-server with R version R-4.1.2

## cmip_split session on n2o, 06.04.2021:
# screen -S _process_GMrecon
# module load R/4.0.3-openblas 
# R (-> not R-3.6.1)

# screen -S _process_GMrecon_blended
# module load R/4.0.3-openblas 

library(ncdf4)
library(raster)
library(fields)

library(foreach)
library(doParallel)

# Set system locale to make sure file ordering is always the same:
# Sys.setlocale(category = "LC_COLLATE", locale = "en_US.UTF-8")

# 0. Load functions needed:
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP5_extr.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")



# 0a. Load/read files:
# --------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
CMIP6.tas_ann_piControl_ct$X = NA
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_ann_piControl_ct.RData")
CMIP6.tas_ann_piControl_ct$X = NA
str(CMIP6.tas_ann_piControl_ct$Y); str(CMIP6.tos_ann_piControl_ct$Y)

CMIP6.files_tas_mon = get.CMIP6.file.list(vari="tas", temp.res = "mon", scen = c("piControl"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/mon/g025/")
CMIP6.files_tos_mon = get.CMIP6.file.list(vari="tos", temp.res = "mon", scen = c("piControl"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tos/mon/g025/")
CMIP6.files_tas_mon = CMIP6.files_tas_mon[which(CMIP6.files_tas_mon$period.length >= 100*12),]
CMIP6.files_tos_mon = CMIP6.files_tos_mon[which(CMIP6.files_tos_mon$period.length >= 100*12),]

CMIP6.files_siconc_mon = get.CMIP6.file.list(vari="siconc", temp.res = "mon", scen = c("historical"), 
                                             CMIP6.dir = "/net/cfc/cmip6/Next_Generation/siconc/mon/native/")
CMIP6.files_sftlf = get.CMIP6.file.list(vari="sftlf", temp.res = "mon", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585", "1pctCO2", "piControl"), 
                                        CMIP6.dir = "/net/cfc/cmip6/Next_Generation/sftlf/fx/native/")


imodall = intersect(x = CMIP6.files_tas_mon$modall, y = CMIP6.files_tos_mon$modall)
CMIP6.files_tas_mon = CMIP6.files_tas_mon[which(CMIP6.files_tas_mon$modall %in% imodall),]
CMIP6.files_tos_mon = CMIP6.files_tos_mon[which(CMIP6.files_tos_mon$modall %in% imodall),]

CMIP6.files_tas_mon$modall == CMIP6.files_tos_mon$modall
ix_ann = apply(cbind(CMIP6.files_tas_mon$period.length/12, CMIP6.files_tos_mon$period.length/12), 1, min) # length index


# 0b. Check/merge files:
# --------------------------------------------------------
tos_file = rep(NA, 77); sftlf_file = rep(NA, 77); siconc_fxd_file = rep(NA, 77); siconc_var_file = rep(NA, 77); 

for (i in 1:length(CMIP6.files_tas_mon$file.name)) {
  print(i)
  
  cur.mod = CMIP6.files_tas_mon$mod[i]
  cur.modcl = CMIP6.files_tas_mon$modcl[i]
  cur.modall = CMIP6.files_tas_mon$modall[i]
  
  # get tos file:
  # tos is fine!
  # cur.ix = which(cur.mod == CMIP6.files_tos_mon$mod)[1]
  # if (!is.na(cur.ix)) tos_file[i] = paste("/net/cfc/cmip6/Next_Generation/tos/mon/g025/", CMIP6.files_tos_mon$file.name[cur.ix], sep="")
  
  
  # get sftlf_file:
  cur.ix = which(cur.mod == CMIP6.files_sftlf$mod)[1]
  if (!is.na(cur.ix)) sftlf_file[i] = paste("/net/cfc/cmip6/Next_Generation/sftlf/fx/native/", CMIP6.files_sftlf$file.name[cur.ix], sep="")
  if (is.na(cur.ix)) {
    print(paste(i, cur.mod, "modcl test", sep=" "))
    cur.ix = which(cur.modcl == CMIP6.files_sftlf$modcl)[1]
    if (!is.na(cur.ix) & cur.mod != "EC-Earth3-CC") {
      dim1 = dim(raster(paste("/net/cfc/cmip6/Next_Generation/sftlf/fx/native/", CMIP6.files_sftlf$file.name[cur.ix], sep="")))
      dim2 = dim(raster(paste("/net/cfc/cmip6/Next_Generation/tas/ann/native/", "tas_ann_", CMIP6.files_tas_mon$modall[i], "_native.nc", sep="")))
      if (all(dim1[1:2] == dim2[1:2])) {
        print("success")
        sftlf_file[i] = paste("/net/cfc/cmip6/Next_Generation/sftlf/fx/native/", CMIP6.files_sftlf$file.name[cur.ix], sep="")
      } 
    }
  }
  
  # assign siconc-fixed mask:
  if (cur.mod %in% c("GISS-E2-1-G" , "GISS-E2-2-G", "GISS-E2-1-G-CC")) cur.mod <- "GISS-E2-1-H"
  
  cur.ix = which(cur.mod == CMIP6.files_siconc_mon$mod & CMIP6.files_siconc_mon$scen == "historical")[1]
  if (!is.na(cur.ix)) {
    siconc_fxd_file[i] = paste("/net/cfc/cmip6/Next_Generation/siconc/mon/native/", CMIP6.files_siconc_mon$file.name[cur.ix], sep="")
    # siconc_fxd_file_fxd_mon = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/siconc_fxd_mon"
    #system(paste("cdo -O -splitmon -selvar,siconc", siconc_fxd_file, siconc_fxd_mon_file, sep=" "))
  }
}


length(which(is.na(sftlf_file)))
length(which(is.na(siconc_fxd_file)))
length(which(is.na(siconc_var_file)))

cbind(CMIP6.files_tas_mon$file.name, CMIP6.files_tos_mon$file.name, sftlf_file, siconc_fxd_file)


# 1. Calculate GMSST and GMLSAT and blended GMST.  + GMSSAT
# --------------------------------------------------------
# loosely follows: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2015GL064888
# blending only based on *absolute* values. Needs later preprocessing
registerDoParallel(cores=40)


  # co2: 1-500: CMIP6.blended.list = foreach(i=1:500) %dopar% {
  # xenon: 501 - length(CMIP6.files_tas_mon$file.name)  CMIP6.blended.list = foreach(i=501:length(CMIP6.files_tas_mon$file.name)) %dopar% {
    
## repeat calculation for E3S and FGO:
# mod.ix.rep = which(CMIP6.files_tas_mon$modcl %in% c("E3S", "FGO"))
# CMIP6.blended.list = foreach(i=mod.ix.rep) %dopar% {

CMIP6.blended.list = foreach(i=1:length(CMIP6.files_tas_mon$file.name)) %dopar% {
  
  print(i)
  # delete all previous files in /temp folder:
  # system("rm /net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/*.nc")
  if (is.na(sftlf_file[i]) | is.na(siconc_fxd_file)[i]) return(NULL);
  
  # Prepare tas and tos anomaly files:
  {
    tas_mon_file = paste("tas_mon_", CMIP6.files_tas_mon$modall[i], "_g025.nc", sep = "")
    tos_mon_file = paste("tos_mon_", CMIP6.files_tos_mon$modall[i], "_g025.nc", sep = "")
    
    tas_mon_file_ = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/tas_mon_file_", i, "_", sep="")
    tos_mon_file_ = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/tos_mon_file_", i, "_", sep="")

    ## regrid to 1° resolution, land and ocean: 
    system(paste("cdo -O -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tas/mon/g025/", tas_mon_file, sep=""),
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/", tas_mon_file, sep="")))

    system(paste("cdo -O -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tos/mon/g025/", tos_mon_file, sep=""),
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/", tos_mon_file, sep="")))

    # absolute temperature files:
    system(paste("cdo -O -splitmon", tas_mon_file, tas_mon_file_, sep=" "))
    system(paste("cdo -O -splitmon", tos_mon_file, tos_mon_file_, sep=" "))
  }
  
  # Prepare land fraction and sea ice files:
  {
    cur.sftlf = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/sftlf",i, ".nc", sep="")
    system(paste("cdo -O -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt -selvar,sftlf", sftlf_file[i], 
                 cur.sftlf, sep=" "))
    
    siconc_fxd_mon_file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/temp/siconc_fxd_mon", i, "_", sep="")
    system(paste("cdo -O -splitmon -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt -selvar,siconc", siconc_fxd_file[i], siconc_fxd_mon_file, sep=" "))
  }
  
  n = length(which(CMIP6.files_tas_mon$mod[i] == CMIP6.tas_ann_piControl_ct$M$mod & CMIP6.files_tas_mon$scen[i] == CMIP6.tas_ann_piControl_ct$M$scen & CMIP6.files_tas_mon$ens.mem[i] == CMIP6.tas_ann_piControl_ct$M$ens.mem))
  tas_ = apply(get_CMIP5_array(file = paste(tas_mon_file_, formatC(1, width = 2, format = "d", flag = "0"), ".nc", sep=""), var = "tas"), 3, c)[,1:ix_ann[i]] - 273.15
  print(n); print(dim(tas_))
  
  GSAT = matrix(data = NA, nrow = n, ncol = 12)  # NI = NO ice.
  GMLSAT_NI = matrix(data = NA, nrow = n, ncol = 12)  # NI = NO ice.
  GMLSAT_II = matrix(data = NA, nrow = n, ncol = 12)  # II = include ice. ice considered as land (!!), with 5% threshold (Cowtan blending paper).
  GMLSAT_MI = matrix(data = NA, nrow = n, ncol = 12)  # MI = Median ice mask of baseline period
  GMSST = matrix(data = NA, nrow = n, ncol = 12)   
  GMMSAT = matrix(data = NA, nrow = n, ncol = 12)  # global marine surface air temperature. -> same mask as GMSST
  GMST = matrix(data = NA, nrow = n, ncol = 12)   # properly blended surface temperature with varying sea ice fraction and absolute temperatures.
  GMST_FM = matrix(data = NA, nrow = n, ncol = 12)   # blending with fixed sea ice mask.
  TMMSAT_40S_40N = matrix(data = NA, nrow = n, ncol = 12)   # Tropical mean...
  TMSST_40S_40N = matrix(data = NA, nrow = n, ncol = 12)   
  TMLSAT_40S_40N = matrix(data = NA, nrow = n, ncol = 12)   
  wair_var = matrix(data = NA, nrow = n, ncol = 12)
  wair_fxd = matrix(data = NA, nrow = n, ncol = 12)
  wair_fxd_median = matrix(data = NA, nrow = n, ncol = 12)
  
  ## RUN MONTHLY CALCULATION OF GMLSAT / GMSST
  for (m in 1:12) {
    print(paste("Month", m))
    
    # READ FILES:
    sftlf_ = raster(cur.sftlf) + 0  
    sftlf = get_CMIP5_array(file = cur.sftlf, var = "sftlf", dims = 2) # image.plot(sftlf)
    land_dims = dim(sftlf)
    
    # convert raster to other format:
    # image.plot(matrix(values(sftlf_), land_dims[1], land_dims[2])[,land_dims[2]:1])
    aw_ = raster::area(sftlf_)
    aw = c(matrix(values(aw_), land_dims[1], land_dims[2])[,land_dims[2]:1]) # image.plot(matrix(values(aw_), land_dims[1], land_dims[2])[,land_dims[2]:1])
    lats = c(matrix(coordinates(sftlf_)[,2], land_dims[1], land_dims[2])[,land_dims[2]:1]) # image.plot(matrix(coordinates(sftlf_)[,2], land_dims[1], land_dims[2])[,land_dims[2]:1])
    
    # siconc_ = brick(paste(siconc_fxd_mon_file, formatC(m, width = 2, format = "d", flag = "0"), ".nc", sep="")) + 0
    siconc = get_CMIP5_array(file = paste(siconc_fxd_mon_file, formatC(m, width = 2, format = "d", flag = "0"), ".nc", sep=""), var = "siconc") 
    # generate sea ice mask:
    siconc_mask = apply(X = siconc, MARGIN = c(1,2), FUN = function(x) all(x>5)) # image.plot(siconc_mask)
    siconc_mask[which(is.na(siconc_mask))] = 0
    # median siconc in historical period:
    siconc_median = apply(X = siconc, MARGIN = c(1,2), FUN = function(x) median(x[1:50])) # image.plot(matrix(siconc_median, 360, 180))
    siconc_median[which(is.na(siconc_median))] = 0
    
    tas_ = apply(get_CMIP5_array(file = paste(tas_mon_file_, formatC(m, width = 2, format = "d", flag = "0"), ".nc", sep=""), var = "tas"), 3, c)[,1:ix_ann[i]] - 273.15 # image.plot(matrix(tas_anom[,20], land_dims[1], land_dims[2]))
    tos_ = apply(get_CMIP5_array(file = paste(tos_mon_file_, formatC(m, width = 2, format = "d", flag = "0"), ".nc", sep=""), var = "tos"), 3, c)[,1:ix_ann[i]]
    
    ## Calculate GMLSAT and GMSST:
    GSAT[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = aw, na.rm = T))
    land.weights = aw * c(sftlf)/100  # image.plot(matrix(land.weights, 360, 180))
    GMLSAT_NI[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = land.weights, na.rm = T))
    # include ice based on minimum mask:
    land.weights_II = aw * c(sftlf)/100 + c(siconc_mask) * aw * c(1 - sftlf/100)  # image.plot(matrix(land.weights_II, 360, 180))
    # image.plot(matrix(land.weights_II, 360, 180))
    GMLSAT_II[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = land.weights_II, na.rm = T))
    # include ice based on median mask:
    land.weights_MI = aw * c(sftlf)/100 + c(siconc_median)/100 * aw * c(1 - sftlf/100)  # image.plot(matrix(land.weights_MI, 360, 180))
    GMLSAT_MI[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = land.weights_MI, na.rm = T))
    
    ocean.weights = c(1 - siconc_mask) * aw * c(1 - sftlf/100) # image.plot(matrix(ocean.weights, 360, 180))
    GMSST[1:(dim(tos_)[2]),m] =  sapply(X = 1:dim(tos_)[2], FUN=function(ix) weighted.mean(x = tos_[,ix], w = ocean.weights, na.rm = T))
    GMMSAT[,m] =  sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = ocean.weights, na.rm = T))
    
    # get values for 40°S to 40°N:
    ocean.weights_40S_40N = ocean.weights; ocean.weights_40S_40N[which(lats > 40 | lats < -40)] = 0 # image.plot(matrix(ocean.weights_40S_40N, 360, 180))
    land.weights_40S_40N = land.weights; land.weights_40S_40N[which(lats > 40 | lats < -40)] = 0 # image.plot(matrix(land.weights_40S_40N, 360, 180))
    TMSST_40S_40N[1:dim(tos_)[2],m] =  sapply(X = 1:dim(tos_)[2], FUN=function(ix) weighted.mean(x = tos_[,ix], w = ocean.weights_40S_40N, na.rm = T))
    TMMSAT_40S_40N[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = ocean.weights_40S_40N, na.rm = T))
    TMLSAT_40S_40N[,m] = sapply(X = 1:dim(tas_)[2], FUN=function(ix) weighted.mean(x = tas_[,ix], w = land.weights_40S_40N, na.rm = T))
    
    # Produce blended dataset with fixed minimum mask: 
    wair = rep.col(c(sftlf) / 100 + ( 1 - c(sftlf) / 100 ) * c(siconc_mask), dim(tos_)[2]) # image.plot(matrix(wair[,80], 360, 180))
    wair[which(is.na(tos_) & wair < 1)] = 1
    Tblend = tas_[,1:dim(tos_)[2]] * wair  # image.plot(matrix(Tblend[,80], 360, 180))
    Tblend[which(wair < 1)] = Tblend[which(wair < 1)] + (tos_ * (1 - wair))[which(wair < 1)]  # image.plot(matrix(Tblend[,80], 360, 180))   # any(is.na(Tblend))
    GMST_FM[1:dim(tos_)[2],m] = sapply(X = 1:dim(Tblend)[2], FUN=function(ix) weighted.mean(x = Tblend[,ix], w = aw))
    wair_fxd[1:dim(tos_)[2],m] = sapply(X = 1:dim(wair)[2], FUN=function(ix) sum(wair[,ix] * aw) / sum(aw) )
    
    # median mask:
    wair_median = rep.col(c(sftlf) / 100 + ( 1 - c(sftlf) / 100 ) * c(siconc_median)/100, dim(tos_)[2]) # image.plot(matrix(wair_median[,80], 360, 180))
    wair_fxd_median[1:dim(tos_)[2],m] = sapply(X = 1:dim(wair_median)[2], FUN=function(ix) sum(wair_median[,ix] * aw) / sum(aw) )
    
  }
  
  # delete files:
  system(paste("rm", tas_mon_file, tos_mon_file, 
               paste(tas_mon_file_, "??.nc", sep=""), paste(tos_mon_file_, "??.nc", sep=""), 
               cur.sftlf, paste(siconc_fxd_mon_file, "??.nc", sep=""), sep=" "))
  
  ret.list = list( GSAT = GSAT, GMLSAT_NI = GMLSAT_NI, GMLSAT_II = GMLSAT_II, GMLSAT_MI = GMLSAT_MI, GMSST = GMSST, GMMSAT = GMMSAT,
                   TMSST_40S_40N = TMSST_40S_40N, TMMSAT_40S_40N = TMMSAT_40S_40N, TMLSAT_40S_40N = TMLSAT_40S_40N,
                   GMST = GMST, GMST_FM = GMST_FM, wair_var = wair_var, wair_fxd = wair_fxd, wair_fxd_median = wair_fxd_median,
                   M = CMIP6.files_tos_mon[i,] )
  save(list = c("ret.list"), 
       file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/out_piControl/out.list", i, ".RData", sep=""))
  return(ret.list)
}


# save CMIP6 blended global mean data files:
save(list = c("CMIP6.blended.list", "CMIP6.files_tas_mon"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/CMIP6.blended1_500.RData")
save(list = c("CMIP6.blended.list", "CMIP6.files_tas_mon"), 
     file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/CMIP6.blended501_959.RData")


# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/CMIP6.blended1_500.RData")
# str(CMIP6.blended.list[[211]])



# 2. Run Processing to merge with "CMIP6.tas_monX_ct" and "CMIP6.tas_monX_ct"
# --------------------------------------------------------
# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/")

# minimum ice mask is useful because otherwise SST may be confined to -1.8°C...
GSAT_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # NI = NO ice.
GMLSAT_NI_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # NI = NO ice.
GMLSAT_II_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # II = include ice. ice considered as land (!!)
GMLSAT_MI_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # MI = Median ice mask of baseline period
GMSST_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   
GMMSAT_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   
TMSST_40S_40N_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
TMMSAT_40S_40N_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
TMLSAT_40S_40N_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
GMST_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   # properly blended surface temperature with varying sea ice fraction and absolute temperatures.
GMST_FM_abs = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   # blending with fixed sea ice mask and anomalies.
wair_var = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
wair_fxd_max = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
wair_fxd_median = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)

GSAT = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # NI = NO ice.
GMLSAT_NI = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # NI = NO ice.
GMLSAT_II = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # II = include ice. ice considered as land (!!)
GMLSAT_MI = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)  # MI = Median ice mask of baseline period
GMSST = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   
GMMSAT = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   
TMSST_40S_40N = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
TMMSAT_40S_40N = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
TMLSAT_40S_40N = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)
GMST = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   # properly blended surface temperature with varying sea ice fraction and absolute temperatures.
GMST_FM = matrix(data = NA, nrow = length(CMIP6.tas_ann_piControl_ct$Y[,1]), ncol = 12)   # blending with fixed sea ice mask and anomalies.


# Run through absolute value blended files, and sort into CMIP6.tas_ann_ct structure
{
    control.vec = rep(NA, length(CMIP6.files_tas_mon$file.name))
    
    for(i in 1:length(CMIP6.files_tas_mon$file.name)) {
      print(i)
      if (is.na(sftlf_file[i]) | is.na(siconc_fxd_file)[i]) {
        print("no land file / siconc file")
        control.vec[i] = "no land file / siconc file"
        next;
      }
      
      if (!file.exists(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/out_piControl/out.list", i, ".RData", sep=""))) {
        print("no ret.list")
        control.vec[i] = "no ret.list"
        next;
      }
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_blended/out_piControl//out.list", i, ".RData", sep=""))
      
      ix = which(CMIP6.files_tas_mon$mod[i] == CMIP6.tas_ann_piControl_ct$M$mod & CMIP6.files_tas_mon$scen[i] == CMIP6.tas_ann_piControl_ct$M$scen & CMIP6.files_tas_mon$ens.mem[i] == CMIP6.tas_ann_piControl_ct$M$ens.mem)
      
        GSAT_abs[ix,] = ret.list$GSAT
        GMLSAT_NI_abs[ix,] = ret.list$GMLSAT_NI
        GMLSAT_II_abs[ix,] = ret.list$GMLSAT_II
        GMLSAT_MI_abs[ix,] = ret.list$GMLSAT_MI
        GMSST_abs[ix,] = ret.list$GMSST
        GMMSAT_abs[ix,] = ret.list$GMMSAT
        TMSST_40S_40N_abs[ix,] = ret.list$TMSST_40S_40N
        TMMSAT_40S_40N_abs[ix,] = ret.list$TMMSAT_40S_40N
        TMLSAT_40S_40N_abs[ix,] = ret.list$TMLSAT_40S_40N
        GMST_abs[ix,] = ret.list$GMST
        GMST_FM_abs[ix,] = ret.list$GMST_FM
        wair_var[ix,] = ret.list$wair_var
        wair_fxd_max[ix,] = ret.list$wair_fxd
        wair_fxd_median[ix,] = ret.list$wair_fxd_median
    }
  ## quick sanity check:
  {
    na.ix = which(GMST_abs[,1] > 20)
    GSAT_abs[na.ix,] = NA
    GMLSAT_NI_abs[na.ix,] = NA
    GMLSAT_II_abs[na.ix,] = NA
    GMLSAT_MI_abs[na.ix,] = NA
    GMSST_abs[na.ix,] = NA
    GMMSAT_abs[na.ix,] = NA
    TMSST_40S_40N_abs[na.ix,] = NA
    TMMSAT_40S_40N_abs[na.ix,] = NA
    TMLSAT_40S_40N_abs[na.ix,] = NA
    GMST_abs[na.ix,] = NA
    GMST_FM_abs[na.ix,] = NA
    wair_var[na.ix,] = NA
    wair_fxd_max[na.ix,] = NA
    wair_fxd_median[na.ix,] = NA
  }
}


 



# subtract reference period average across ensemble:
{
  ref.period.years = 1:2000
  ref.scen = "piControl"
  mod.un = unique(CMIP6.tas_ann_piControl_ct$M$file.name)
  
    for (m in 1:length(mod.un)) {
      print(paste(mod.un[m]))
      mod.ix = which(CMIP6.tas_ann_piControl_ct$M$file.name == mod.un[m])
      ref.mod.ix = which(CMIP6.tas_ann_piControl_ct$M$file.name == mod.un[m] & CMIP6.tas_ann_piControl_ct$M$year %in% ref.period.years & CMIP6.tas_ann_piControl_ct$M$scen == ref.scen)
      
      for (mon in 1:12) {
        GSAT[mod.ix,mon] = GSAT_abs[mod.ix,mon] - mean(GSAT_abs[ref.mod.ix,mon], na.rm = T)
        GMLSAT_NI[mod.ix,mon] = GMLSAT_NI_abs[mod.ix,mon] - mean(GMLSAT_NI_abs[ref.mod.ix,mon], na.rm = T)
        GMLSAT_II[mod.ix,mon] = GMLSAT_II_abs[mod.ix,mon] - mean(GMLSAT_II_abs[ref.mod.ix,mon], na.rm = T)
        GMLSAT_MI[mod.ix,mon] = GMLSAT_MI_abs[mod.ix,mon] - mean(GMLSAT_MI_abs[ref.mod.ix,mon], na.rm = T)
        GMSST[mod.ix,mon] = GMSST_abs[mod.ix,mon] - mean(GMSST_abs[ref.mod.ix,mon], na.rm = T)
        GMMSAT[mod.ix,mon] = GMMSAT_abs[mod.ix,mon] - mean(GMMSAT_abs[ref.mod.ix,mon], na.rm = T)
        TMSST_40S_40N[mod.ix,mon] = TMSST_40S_40N_abs[mod.ix,mon] - mean(TMSST_40S_40N_abs[ref.mod.ix,mon], na.rm = T)
        TMMSAT_40S_40N[mod.ix,mon] = TMMSAT_40S_40N_abs[mod.ix,mon] - mean(TMMSAT_40S_40N_abs[ref.mod.ix,mon], na.rm = T)
        TMLSAT_40S_40N[mod.ix,mon] = TMLSAT_40S_40N_abs[mod.ix,mon] - mean(TMLSAT_40S_40N_abs[ref.mod.ix,mon], na.rm = T)
        GMST[mod.ix,mon] = GMST_abs[mod.ix,mon] - mean(GMST_abs[ref.mod.ix,mon], na.rm = T)
        GMST_FM[mod.ix,mon] = GMST_FM_abs[mod.ix,mon] - mean(GMST_FM_abs[ref.mod.ix,mon], na.rm = T)
      }
    }
}
  





# Build into annual files and save (+ modify annual files to reflect GSAT naming convention):
{
  load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
  load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_ann_piControl_ct.RData")
  # CMIP6.tas_ann_piControl_ct$Y = CMIP6.tas_ann_piControl_ct$Y[,1]
  names(CMIP6.tas_ann_piControl_ct$Y) = c("GSAT")
  # CMIP6.tos_ann_ct$Y <- CMIP6.tas_ann_ct$Y
  names(CMIP6.tos_ann_piControl_ct$Y) = c("GSAT")
  
  df.ann = data.frame(cbind(GSAT_ = rowMeans(GSAT),
                            GMLSAT_NI = rowMeans(GMLSAT_NI), GMLSAT_II = rowMeans(GMLSAT_II), GMLSAT_MI = rowMeans(GMLSAT_MI),
                            GMSST = rowMeans(GMSST), GMMSAT = rowMeans(GMMSAT),
                            TMSST_40S_40N = rowMeans(TMSST_40S_40N), TMMSAT_40S_40N = rowMeans(TMMSAT_40S_40N), TMLSAT_40S_40N = rowMeans(TMLSAT_40S_40N),
                            GMST = rowMeans(GMST), GMST_FM = rowMeans(GMST_FM),
                            GSAT_abs = rowMeans(GSAT_abs),
                            GMLSAT_NI_abs = rowMeans(GMLSAT_NI_abs), GMLSAT_II_abs = rowMeans(GMLSAT_II_abs), GMLSAT_MI_abs = rowMeans(GMLSAT_MI_abs),
                            GMSST_abs = rowMeans(GMSST_abs), GMMSAT_abs = rowMeans(GMMSAT_abs),
                            TMSST_40S_40N_abs = rowMeans(TMSST_40S_40N_abs), TMMSAT_40S_40N_abs = rowMeans(TMMSAT_40S_40N_abs), TMLSAT_40S_40N_abs = rowMeans(TMLSAT_40S_40N_abs),
                            GMST_abs = rowMeans(GMST_abs), GMST_FM_abs = rowMeans(GMST_FM_abs), 
                            wair_var = rowMeans(wair_var), wair_fxd_max = rowMeans(wair_fxd_max), wair_fxd_median = rowMeans(wair_fxd_median)))
  CMIP6.tas_ann_piControl_ct$Y = data.frame(cbind(CMIP6.tas_ann_piControl_ct$Y$GSAT, df.ann))
  CMIP6.tos_ann_piControl_ct$Y = data.frame(cbind(CMIP6.tos_ann_piControl_ct$Y$GSAT, df.ann))
  save(list = "CMIP6.tas_ann_piControl_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData", sep=""))
  save(list = "CMIP6.tos_ann_piControl_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_ann_piControl_ct.RData", sep=""))
}

  
  # Build into monthly files and save:
  for (mon in 1:12) {
    print(mon)
    
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
    
    cur.df.mon = data.frame(cbind(GSAT_ = GSAT[,mon], GMLSAT_NI = GMLSAT_NI[,mon], GMLSAT_II = GMLSAT_II[,mon], GMLSAT_MI = GMLSAT_MI[,mon],
                                  GMSST = GMSST[,mon], GMMSAT = GMMSAT[,mon],
                                  TMSST_40S_40N = TMSST_40S_40N[,mon], TMMSAT_40S_40N = TMMSAT_40S_40N[,mon], TMLSAT_40S_40N = TMLSAT_40S_40N[,mon],
                                  GMST = GMST[,mon], GMST_FM = GMST_FM[,mon], 
                                  GSAT_abs = GSAT_abs[,mon], GMLSAT_NI_abs = GMLSAT_NI_abs[,mon], GMLSAT_II_abs = GMLSAT_II_abs[,mon], GMLSAT_MI_abs = GMLSAT_MI_abs[,mon],
                                  GMSST_abs = GMSST_abs[,mon], GMMSAT_abs = GMMSAT_abs[,mon],
                                  TMSST_40S_40N_abs = TMSST_40S_40N_abs[,mon], TMMSAT_40S_40N_abs = TMMSAT_40S_40N_abs[,mon], TMLSAT_40S_40N_abs = TMLSAT_40S_40N_abs[,mon],
                                  GMST_abs = GMST_abs[,mon], GMST_FM_abs = GMST_FM_abs[,mon],
                                  wair_var = wair_var[,mon], wair_fxd_max = wair_fxd_max[,mon], wair_fxd_median = wair_fxd_median[,mon]))
    
    # CMIP6.tas_monX_ct$Y = CMIP6.tas_monX_ct$Y[,1:7]
    # colnames(CMIP6.tas_monX_ct$Y) = c("GSAT", "GSAT_f", "GSAT_fl", "GSAT_f_noinclmem", "GSAT_ANN", "GSAT_ANN_f", "GSAT_ANN_fl")
    names(cmip6_piControl_tas_mon_5d00_XAX_ct$Y) = "GSAT"
    names(cmip6_piControl_tos_mon_5d00_XAX_ct$Y) = "AGMT"
    # CMIP6.tos_monX_ct$Y <- CMIP6.tas_monX_ct$Y
    
    ## make anomalies out of Y$GMST dataset?
    cmip6_piControl_tas_mon_5d00_XAX_ct$Y = data.frame(cbind(cmip6_piControl_tas_mon_5d00_XAX_ct$Y[,1], cur.df.mon))
    cmip6_piControl_tos_mon_5d00_XAX_ct$Y = data.frame(cbind(cmip6_piControl_tos_mon_5d00_XAX_ct$Y[,1], cur.df.mon))
    
    # save files:
    save(list = "cmip6_piControl_tas_mon_5d00_XAX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
    save(list = "cmip6_piControl_tos_mon_5d00_XAX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
  }
  





# 3. check if things look reasonable:
# --------------------------------------------------------

load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon8_piControl_ct.RData", sep=""))

plot(cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GSAT_[23000:24000])
unique(CMIP6.tas_ann_piControl_ct$M$mod[23500:23900])
unique(CMIP6.tas_ann_piControl_ct$M$ens.mem[23500:23900])


mod.un = unique(cmip6_piControl_tas_mon_5d00_XAX_ct$M$mod)
par(mar=c(5,5,1,1))

for (m in 1:length(mod.un)) {
  
  ix = which(cmip6_piControl_tas_mon_5d00_XAX_ct$M$mod == mod.un[m] & !(is.na(cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GMSST)) & !(is.na(cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GMLSAT_MI)))
  if(is.na(ix[1])) next;
  
  #plot(33/100 * CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] + 67 / 100 * CMIP6.tas_ann_ct$Y$GMSST[ix], CMIP6.tas_ann_ct$Y$GMST_FM[ix], main = mod.un[m])
  cur.cor = cor(33/100 * cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GMLSAT_MI[ix] + 67 / 100 * cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GMSST[ix], cmip6_piControl_tas_mon_5d00_XAX_ct$Y$GMST_FM[ix])
  #legend("topleft", paste("R = ", round(cur.cor, 5)))
  print(mod.un[m]); print(m); print(cur.cor)
  
  #plot(33/100 * CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] + 67 / 100 * CMIP6.tas_ann_ct$Y$GMSST[ix], CMIP6.tas_ann_ct$Y$GSAT[ix], main = mod.un[m])
  #abline(0,1 , col="red")
  
  
  # get land-ocean warming ratio:
  # plot(x = CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix], y = CMIP6.tas_ann_ct$Y$GMLSAT_NI[ix], main = mod.un[m])
  # plot(x = CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix], y = CMIP6.tas_ann_ct$Y$GMSST[ix], main = mod.un[m], xlim = c(-1, 7), ylim = c(-1, 4))
  # abline(0,1 , col="red")
  # print(round(lm(CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] ~ CMIP6.tas_ann_ct$Y$GMSST[ix])$coefficients[2], 3))

  # plot(x = CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix], y = CMIP6.tas_ann_ct$Y$GMSST[ix], main = mod.un[m], xlim = c(-1, 7), ylim = c(-1, 4))
  # abline(0,1 , col="red")
  # print(round(lm(CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] ~ CMIP6.tas_ann_ct$Y$GMSST[ix])$coefficients[2], 3))

  plot(x = cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMLSAT_40S_40N[ix], y = cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMMSAT_40S_40N[ix], main = mod.un[m], xlim = c(-1, 1), ylim = c(-1, 1))
  cor(cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMLSAT_40S_40N[ix], y = cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMMSAT_40S_40N[ix])
  abline(0,1 , col="red")
  print(round(lm(cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMLSAT_40S_40N[ix] ~ cmip6_piControl_tas_mon_5d00_XAX_ct$Y$TMMSAT_40S_40N[ix])$coefficients[2], 3))
  
  Sys.sleep(1)
  
}



## problems in fgoals?
## problems in E3SM?


plot(33/100 * CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] + 67 / 100 * CMIP6.tas_ann_ct$Y$GMSST[ix])
plot(CMIP6.tas_ann_ct$Y$GSAT[ix])
plot(33/100 * CMIP6.tas_ann_ct$Y$GMLSAT_MI[ix] + 67 / 100 * CMIP6.tas_ann_ct$Y$GMSST[ix], CMIP6.tas_ann_ct$Y$GSAT[ix], main = mod.un[m])

plot(CMIP6.tas_ann_ct$Y$GMST[ix], CMIP6.tas_ann_ct$Y$GSAT[ix], main = mod.un[m])
plot(CMIP6.tas_ann_ct$Y$GMST_FM[ix], CMIP6.tas_ann_ct$Y$GSAT[ix], main = mod.un[m])





