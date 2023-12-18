
## ------------------------------------------------------------------------
## Read DAMIP files 
## ------------------------------------------------------------------------

# Sebastian Sippel
# 26.04.2022

library(ncdf4)

# Read functions: 
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_functions_CMIP5_extr.R")


# Define directories:
CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/tas/ann/g025"
grid_5d00 = "/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt"
target_dir_5d00_monthly = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tas/mon/"
target_dir_5d00_annual = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tas/ann/"

target_dir_5d00_monthly_tos = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tos/mon/"
target_dir_5d00_annual_tos = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tos/ann/"



# ------------------------------------------------------------------
## Reprocess CMIP6 tas-ann historical+ssp245 files (for reference period mean-subtraction):
# ------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann.RData")

# CMIP6-tas-ann:
CMIP6.files_tas_ann = get.CMIP6.file.list(vari="tas", temp.res = "ann", scen = c("historical", "ssp245"), 
                                          CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann.RData")

# CMIP6.tas_ann = adjust.years.cmip6(XAX = CMIP6.tas_ann)
# CMIP6.tas_ann_ct = center.XAX_ann(XAX = CMIP6.tas_ann, fact = 1, areaw = areaw, ref.period.years = 1961:1990)



# 0.b) Process CMIP6 piControl files and regrid to 5x5°:
# ------------------------------------------------------------------
CMIP6.files_tas_ann = get.CMIP6.file.list(vari= "tas", temp.res = "ann", scen = c("piControl"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/ann/g025/")
CMIP6.files_tos_ann = get.CMIP6.file.list(vari="tos", temp.res = "ann", scen = c("piControl"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tos/ann/g025/")
CMIP6.files_tas_mon = get.CMIP6.file.list(vari= "tas", temp.res = "mon", scen = c("piControl"), 
                                                CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/mon/g025/")
CMIP6.files_tos_mon = get.CMIP6.file.list(vari="tos", temp.res = "mon", scen = c("piControl"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tos/mon/g025/")

CMIP6.files_tas_mon = CMIP6.files_tas_mon[which(CMIP6.files_tas_mon$period.length >= 100*12),]
CMIP6.files_tos_mon = CMIP6.files_tos_mon[which(CMIP6.files_tos_mon$period.length >= 100*12),]

imodall = intersect(x = CMIP6.files_tas_mon$modall, y = CMIP6.files_tos_mon$modall)
CMIP6.files_tas_ann = CMIP6.files_tas_ann[which(CMIP6.files_tas_ann$modall %in% imodall),]
CMIP6.files_tas_mon = CMIP6.files_tas_mon[which(CMIP6.files_tas_mon$modall %in% imodall),]
CMIP6.files_tos_ann = CMIP6.files_tos_ann[which(CMIP6.files_tos_ann$modall %in% imodall),]
CMIP6.files_tos_mon = CMIP6.files_tos_mon[which(CMIP6.files_tos_mon$modall %in% imodall),]

ix_ann = apply(cbind(CMIP6.files_tas_ann$period.length, CMIP6.files_tos_ann$period.length), 1, min) # length index
# apply(cbind(CMIP6.files_tas_mon$period.length, CMIP6.files_tos_mon$period.length), 1, min)
# which(CMIP6.files_tas_mon$period.length != CMIP6.files_tos_mon$period.length)
# which(CMIP6.files_tas_ann$period.length != CMIP6.files_tos_ann$period.length)
# CMIP6.files_tas_mon$file.name[2]
# CMIP6.files_tos_mon$file.name[2]

# 1. process for tas+tos:
CMIP6.orig_dir_ann = "/net/ch4/data/cmip6-Next_Generation/tas/ann/g025"
CMIP6.orig_dir = "/net/ch4/data/cmip6-Next_Generation/tas/mon/g025"
CMIP6.orig_dir_ann_tos = "/net/ch4/data/cmip6-Next_Generation/tos/ann/g025"
CMIP6.orig_dir_tos = "/net/ch4/data/cmip6-Next_Generation/tos/mon/g025"

for (i in 1:length(CMIP6.files_tas_ann$file.name)) { 
  print(i)

  ## tas files ANNUAL:
  tas.ifile = paste(CMIP6.orig_dir_ann, "/tas_ann_", CMIP6.files_tas_ann$modall[i], "_g025.nc", sep="")
  tas.ofile = paste(target_dir_5d00_annual, CMIP6.files_tas_ann$mod[i], "_", CMIP6.files_tas_ann$scen[i], "_", CMIP6.files_tas_ann$ens.mem[i], "_5d00.nc", sep="")
  ## get 5d00 annual file:
  system(paste("cdo -remapbil,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt ", 
               tas.ifile, " ", tas.ofile, sep=""))
  
  ## tas files MONTHLY:
  tas.ifile = paste(CMIP6.orig_dir, "/tas_mon_", CMIP6.files_tas_mon$modall[i], "_g025.nc", sep="")
  tas.ofile = paste(target_dir_5d00_monthly, CMIP6.files_tas_mon$mod[i], "_", CMIP6.files_tas_mon$scen[i], "_", CMIP6.files_tas_mon$ens.mem[i], "_5d00.nc", sep="")
  ## get 5d00 annual file:
  system(paste("cdo -remapbil,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt ", 
               tas.ifile, " ", tas.ofile, sep=""))
  
  
  
  ## tos files ANNUAL:
  tos.ifile = paste(CMIP6.orig_dir_ann_tos, "/tos_ann_", CMIP6.files_tos_ann$modall[i], "_g025.nc", sep="")
  tos.ofile = paste(target_dir_5d00_annual_tos, CMIP6.files_tos_ann$mod[i], "_", CMIP6.files_tos_ann$scen[i], "_", CMIP6.files_tos_ann$ens.mem[i], "_5d00.nc", sep="")
  ## get 5d00 annual file:
  system(paste("cdo -remapbil,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt ", 
               tos.ifile, " ", tos.ofile, sep=""))
  
  ## tos files MONTHLY:
  tos.ifile = paste(CMIP6.orig_dir_tos, "/tos_mon_", CMIP6.files_tos_mon$modall[i], "_g025.nc", sep="")
  tos.ofile = paste(target_dir_5d00_monthly_tos, CMIP6.files_tos_mon$mod[i], "_", CMIP6.files_tos_mon$scen[i], "_", CMIP6.files_tos_mon$ens.mem[i], "_5d00.nc", sep="")
  ## get 5d00 annual file:
  system(paste("cdo -remapbil,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt ", 
               tos.ifile, " ", tos.ofile, sep=""))
}


# --------------------------------
# 2. process for tas:
# --------------------------------

# READ CMIP6 piControl files:
cmip6_piControl_tas_ann_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files_tas_ann$mod, "_", CMIP6.files_tas_ann$scen, "_", CMIP6.files_tas_ann$ens.mem, "_5d00.nc", sep=""), 
                                                scen = CMIP6.files_tas_ann$scen, var="tas", res = "ann", 
                                                CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tas/ann/", 
                                                subtract.drift.piControl = F, time.count = ix_ann, rm.years = 0)
cmip6_piControl_tas_ann_5d00_XAX = get.XAX_ann(X = cmip6_piControl_tas_ann_5d00, M = CMIP6.files_tas_ann, start.year = 1, cmip = "cmip6", areaw = areaw)
CMIP6.tas_ann_piControl_ct = center.XAX_ann(XAX = cmip6_piControl_tas_ann_5d00_XAX, fact = 1, areaw = areaw, ref.period.years = 1:2000, ref.scen = "piControl")
save(list = c("CMIP6.tas_ann_piControl_ct"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")


# MON. save centered piControl dataset:
cmip6_piControl_tas_mon_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files_tas_mon$mod, "_", CMIP6.files_tas_mon$scen, "_", CMIP6.files_tas_mon$ens.mem, "_5d00.nc", sep=""), 
                                            scen = CMIP6.files_tas_mon$scen, var="tas", res = "mon", 
                                            CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tas/mon/", 
                                            subtract.drift.piControl = F, time.count = ix_ann*12, rm.years = 0)
for (mon in 1:12) {
  print(mon)
  
  cmip6_piControl_tas_mon_5d00_XAX = get.XAX_mon(X = cmip6_piControl_tas_mon_5d00, mon = mon, ncol = 72 * 36, M = CMIP6.files_tas_mon, start.year = 1, 
                                                 cmip = "CMIP6", areaw = areaw, rm.HIST.from.rcp26 = F)
  cmip6_piControl_tas_mon_5d00_XAX_ct = center.XAX_ann(XAX = cmip6_piControl_tas_mon_5d00_XAX, fact = 1, areaw = areaw, ref.period.years = 1:2000, ref.scen = "piControl")
  
  save(list = "cmip6_piControl_tas_mon_5d00_XAX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
}




# --------------------------------
# 3. process for tos:
# --------------------------------

# READ CMIP6 piControl files:
cmip6_piControl_tos_ann_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files_tos_ann$mod, "_", CMIP6.files_tos_ann$scen, "_", CMIP6.files_tos_ann$ens.mem, "_5d00.nc", sep=""), 
                                                scen = CMIP6.files_tos_ann$scen, var="tos", res = "ann", 
                                                CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tos/ann/", 
                                                subtract.drift.piControl = F, time.count = ix_ann, rm.years = 0)
cmip6_piControl_tos_ann_5d00_XAX = get.XAX_ann(X = cmip6_piControl_tos_ann_5d00, M = CMIP6.files_tos_ann, start.year = 1, cmip = "cmip6", areaw = areaw)
CMIP6.tos_ann_piControl_ct = center.XAX_ann(XAX = cmip6_piControl_tos_ann_5d00_XAX, fact = 1, areaw = areaw, ref.period.years = 1:2000, ref.scen = "piControl")
save(list = c("CMIP6.tos_ann_piControl_ct"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_ann_piControl_ct.RData")


# MON. save centered piControl datoset:
cmip6_piControl_tos_mon_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files_tos_mon$mod, "_", CMIP6.files_tos_mon$scen, "_", CMIP6.files_tos_mon$ens.mem, "_5d00.nc", sep=""), 
                                                scen = CMIP6.files_tos_mon$scen, var="tos", res = "mon", 
                                                CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tos/mon/", 
                                                subtract.drift.piControl = F, time.count = ix_ann * 12, rm.years = 0)
for (mon in 1:12) {
  print(mon)
  
  cmip6_piControl_tos_mon_5d00_XAX = get.XAX_mon(X = cmip6_piControl_tos_mon_5d00, mon = mon, ncol = 72 * 36, M = CMIP6.files_tos_mon, start.year = 1, 
                                                 cmip = "CMIP6", areaw = areaw, rm.HIST.from.rcp26 = F)
  cmip6_piControl_tos_mon_5d00_XAX_ct = center.XAX_ann(XAX = cmip6_piControl_tos_mon_5d00_XAX, fact = 1, areaw = areaw, ref.period.years = 1:2000, ref.scen = "piControl")
  
  save(list = "cmip6_piControl_tos_mon_5d00_XAX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
}




# --------------------------------
#### REPROCESSING FOR MONTHLY tos BASED ON NATIVE FILES (in order to get full coverage!!!):
# --------------------------------
str(CMIP6.files_tas_ann)
CMIP6.files_tas_ann_CURRENT = get.CMIP6.file.list(vari= "tas", temp.res = "ann", scen = c("piControl"), 
                                          CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tas/ann/")
cbind(CMIP6.files_tas_ann$file.name[1:77], CMIP6.files_tas_ann_CURRENT$file.name)


# 1. process for tas+tos:
CMIP6.orig_dir_tos = "/net/ch4/data/cmip6-Next_Generation/tos/mon/native"

for (i in 1:length(CMIP6.files_tas_ann$file.name)) { 
  print(i)
  
  ## tos files MONTHLY:
  tos.ifile = paste(CMIP6.orig_dir_tos, "/tos_mon_", CMIP6.files_tos_mon$modall[i], "_native.nc", sep="")
  tos.ofile = paste(target_dir_5d00_monthly_tos, CMIP6.files_tos_mon$mod[i], "_", CMIP6.files_tos_mon$scen[i], "_", CMIP6.files_tos_mon$ens.mem[i], "_5d00.nc", sep="")
  
  if (i == 2) {  ## tos_mon_ACCESS-ESM1-5_piControl_r1i1p1f1_native.nc
    tos.ifile = paste("/net/ch4/data/cmip6-Next_Generation/tos/mon/g025/", "/tos_mon_", CMIP6.files_tos_mon$modall[i], "_g025.nc", sep="")
  }
  ## get 5d00 annual file:
  system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt ", 
               tos.ifile, " ", tos.ofile, sep=""))
}


# MON. save centered piControl datoset:
CMIP6.files_tos_mon = CMIP6.files_tos_mon[-51,]  # exclude ICON which is new

cmip6_piControl_tos_mon_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files_tos_mon$mod, "_", CMIP6.files_tos_mon$scen, "_", CMIP6.files_tos_mon$ens.mem, "_5d00.nc", sep=""), 
                                                scen = CMIP6.files_tos_mon$scen, var="tos", res = "mon", 
                                                CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_piControl/tos/mon/", 
                                                subtract.drift.piControl = F, time.count = ix_ann[-51] * 12, rm.years = 0)


for (mon in 1:12) {
  print(mon)
  
  cmip6_piControl_tos_mon_5d00_XAX = get.XAX_mon(X = cmip6_piControl_tos_mon_5d00, mon = mon, ncol = 72 * 36, M = CMIP6.files_tos_mon, start.year = 1, 
                                                 cmip = "CMIP6", areaw = areaw, rm.HIST.from.rcp26 = F)
  cmip6_piControl_tos_mon_5d00_XAX_ct = center.XAX_ann(XAX = cmip6_piControl_tos_mon_5d00_XAX, fact = 1, areaw = areaw, ref.period.years = 1:2000, ref.scen = "piControl")
  
  save(list = "cmip6_piControl_tos_mon_5d00_XAX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
}


# write current list of files:
write.table(x = CMIP6.files_tos_mon, 
            file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/cmip6_piControl_files.txt", quote = F, sep=";", row.names = F)


load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
str(cmip6_piControl_tos_mon_5d00_XAX_ct)

## THE END ##



