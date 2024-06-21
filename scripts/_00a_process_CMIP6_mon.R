
# ------------------------------------------------------------------------------------
# READ DATA INTO FILES for cmip6 training for extremes
# ------------------------------------------------------------------------------------

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
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")

setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/")


# --------------------------------------------------------------------------
# 1. Derive file list with files to pick:
# --------------------------------------------------------------------------

# (a) get file list:
CMIP6.files_tas_ann = get.CMIP6.file.list(vari="tas", temp.res = "ann", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585"), 
                                         CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/ann/g025/")
CMIP6.files_tas_mon = get.CMIP6.file.list(vari="tas", temp.res = "mon", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/mon/g025/")

CMIP6.files_tos_ann = get.CMIP6.file.list(vari="tos", temp.res = "ann", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tos/ann/native/")
CMIP6.files_tos_mon = get.CMIP6.file.list(vari="tos", temp.res = "mon", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tos/mon/native/")

# intersect all files for GSAT reconstruction:
imodall = intersect(x = CMIP6.files_tas_ann$modall, y = CMIP6.files_tos_ann$modall)
CMIP6.files_tas_ann = CMIP6.files_tas_ann[which(CMIP6.files_tas_ann$modall %in% imodall),]
CMIP6.files_tas_mon = CMIP6.files_tas_mon[which(CMIP6.files_tas_mon$modall %in% imodall),]
CMIP6.files_tos_ann = CMIP6.files_tos_ann[which(CMIP6.files_tos_ann$modall %in% imodall),]
CMIP6.files_tos_mon = CMIP6.files_tos_mon[which(CMIP6.files_tos_mon$modall %in% imodall),]

# Save list of cmip5_files used for analysis:
write.table(x = CMIP6.files_tas_ann, file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/cmip6_historical_ssp_files.txt", quote = F, sep=";", row.names = F)
# read.table(file = "/net/h2o/climphys1/sippels/_projects/global_mean_prediction_v2/data/_metadata/cmip5_files.txt", sep = ";")
# ?read.table

# --------------------------------------------------------------------------
# 2. Process cmip6 data:
# --------------------------------------------------------------------------
registerDoParallel(cores=30)


# i. Process temperature files to 5x5° grid:
{
  foreach(i=1:length(imodall)) %dopar% {
    print(i)
    
    # (1) Regrid to 5d00 grid (-> CRUTEM5):
    tas_ann_file = paste("tas_ann_", CMIP6.files_tas_ann$modall[i], "_g025.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tas/ann/g025/", tas_ann_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/", tas_ann_file, sep="")))
    
    tas_mon_file = paste("tas_mon_", CMIP6.files_tas_mon$modall[i], "_g025.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tas/mon/g025/", tas_mon_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/mon/", tas_mon_file, sep="")))
    
    tos_ann_file = paste("tos_ann_", CMIP6.files_tos_ann$modall[i], "_native.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tos/ann/native/", tos_ann_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/ann/", tos_ann_file, sep="")))
    
    tos_mon_file = paste("tos_mon_", CMIP6.files_tos_mon$modall[i], "_native.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt -remapdis,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_1d00.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/tos/mon/native/", tos_mon_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/mon/", tos_mon_file, sep="")))
    }
}







# --------------------------------------------------------------------------
# 2. Read cmip6 data:
# --------------------------------------------------------------------------
raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw = c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]))
areaw_cos = sqrt(cos(c(matrix(raster::coordinates(raster.template)[,2], 72, 36)[,36:1])*pi/180))


# 2.(a) READ & Process CMIP6-tas files:
# --------------------------------------------------------------------------
{
   setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
  
  # CMIP6-tas-ann
   CMIP6.files_tas_ann = get.CMIP6.file.list(vari="tas", temp.res = "ann", scen = c("historical", "ssp245"), 
                                               CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
   CMIP6_tas_ann = read.CMIP6_novar(file.name = CMIP6.files_tas_ann$file.name, scen = CMIP6.files_tas_ann$scen, var="tas", res = "ann", 
                                              CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
   CMIP6.tas_ann = get.XAX_ann(X = CMIP6_tas_ann, ncol = 72 * 36, M = CMIP6.files_tas_ann, start.year = 1850, cmip = "CMIP6", areaw = areaw)
   CMIP6.tas_ann = adjust.years.cmip6(XAX = CMIP6.tas_ann)
   CMIP6.tas_ann_ct = center.XAX_ann(XAX = CMIP6.tas_ann, fact = 1, areaw = areaw, ref.period.years = 1961:1990)
   ## generate the mod.phys.df to know the number of members:  # -> in a large number of cases, only about 3 members. 
   # mod.df = extract.mod.df(XAX_scen = CMIP6.tas_ann) # hist(mod.df$nmem[-which(mod.df$nmem == 0)])
   CMIP6.tas_ann_ct = extract.fraw2(XAX_scen = CMIP6.tas_ann_ct, nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
   # ix = which(CMIP6.tas_ann_ct$M$mod == "KACE-1-0-G"); plot(CMIP6.tas_ann_ct$Y$AGMT_f[ix]); lines(CMIP6.tas_ann_ct$Y$AGMT_fl[ix], col = "blue")
   # plot(CMIP6.tas_ann_ct$Y$AGMT_fl[ix], CMIP6.tas_ann_ct$Y$AGMT_f[ix]); abline(0, 1, col = "red")
   # give number of members as variable in M?
   save(list = "CMIP6.tas_ann", file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann.RData")
   save(list = "CMIP6.tas_ann_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_ct.RData", sep=""))
   
   
  # CMIP6-tas-mon
   setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/mon/")
   
   CMIP6.files_tas_mon = get.CMIP6.file.list(vari="tas", temp.res = "mon", scen = c("historical", "ssp245"), 
                                             CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/mon/")
   CMIP6_tas_mon = read.CMIP6_novar(file.name = CMIP6.files_tas_mon$file.name, scen = CMIP6.files_tas_mon$scen, var="tas", res = "mon", 
                                    CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/mon/")
   CMIP6.tas_mon = list()
   CMIP6.tas_mon_ct = list()
   for (mon in 1:12) {
      CMIP6.tas_mon[[mon]] = get.XAX_mon(X = CMIP6_tas_mon, mon = mon, ncol = 72 * 36, M = CMIP6.files_tas_mon, start.year = 1850, cmip = "CMIP6", areaw = areaw, rm.HIST.from.rcp26 = F)
      CMIP6.tas_mon[[mon]] = adjust.years.cmip6(XAX = CMIP6.tas_mon[[mon]])
      CMIP6.tas_mon_ct[[mon]] = center.XAX_ann(XAX = CMIP6.tas_mon[[mon]], fact = 1, areaw = areaw, ref.period.years = 1961:1990)
      CMIP6.tas_mon_ct[[mon]] = extract.fraw2(XAX_scen = CMIP6.tas_mon_ct[[mon]], nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
      CMIP6.tas_mon_ct[[mon]]$Y$AGMT_ANN = CMIP6.tas_ann_ct$Y$AGMT
      CMIP6.tas_mon_ct[[mon]]$Y$AGMT_ANN_f = CMIP6.tas_ann_ct$Y$AGMT_f
      CMIP6.tas_mon_ct[[mon]]$Y$AGMT_ANN_fl = CMIP6.tas_ann_ct$Y$AGMT_fl
         
      # plot(CMIP6.tas_mon_ct[[mon]]$Y$AGMT_fl[1:1000], type='l')
      CMIP6.tas_monX_ct = CMIP6.tas_mon_ct[[mon]]
      save(list = "CMIP6.tas_monX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_ct.RData", sep=""))
   }
}




# 2.(b) READ & Process CMIP6-tos files:
# --------------------------------------------------------------------------
{
   setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/ann/")
   load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_ct.RData")
   
   # CMIP6-tos-ann
   CMIP6.files_tos_ann = get.CMIP6.file.list(vari="tos", temp.res = "ann", scen = c("historical", "ssp245"), 
                                             CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/ann/")
   CMIP6_tos_ann = read.CMIP6_novar(file.name = CMIP6.files_tos_ann$file.name, scen = CMIP6.files_tos_ann$scen, var="tos", res = "ann", 
                                    CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/ann/")
   CMIP6.tos_ann = get.XAX_ann(X = CMIP6_tos_ann, ncol = 72 * 36, M = CMIP6.files_tos_ann, start.year = 1850, cmip = "CMIP6", areaw = areaw)
   CMIP6.tos_ann = adjust.years.cmip6(XAX = CMIP6.tos_ann)
   CMIP6.tos_ann_ct = center.XAX_ann(XAX = CMIP6.tos_ann, fact = 1, areaw = areaw, ref.period.years = 1961:1990)
   ## generate the mod.phys.df to know the number of members:  # -> in a large number of cases, only about 3 members. 
   # mod.df = extract.mod.df(XAX_scen = CMIP6.tos_ann) # hist(mod.df$nmem[-which(mod.df$nmem == 0)])
   CMIP6.tos_ann_ct = extract.fraw2(XAX_scen = CMIP6.tos_ann_ct, nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
   save(list = "CMIP6.tos_ann_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_ann_ct.RData", sep=""))
   
   
   # CMIP6-tos-mon:
   setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/mon/")
   
   CMIP6.files_tos_mon = get.CMIP6.file.list(vari="tos", temp.res = "mon", scen = c("historical", "ssp245"), 
                                             CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/mon/")
   CMIP6_tos_mon = read.CMIP6_novar(file.name = CMIP6.files_tos_mon$file.name, scen = CMIP6.files_tos_mon$scen, var="tos", res = "mon", 
                                    CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tos/mon/")
   CMIP6.tos_mon = list()
   CMIP6.tos_mon_ct = list()
   for (mon in 1:12) {
      CMIP6.tos_mon[[mon]] = get.XAX_mon(X = CMIP6_tos_mon, mon = mon, ncol = 72 * 36, M = CMIP6.files_tos_mon, start.year = 1850, cmip = "CMIP6", areaw = areaw, rm.HIST.from.rcp26 = F)
      CMIP6.tos_mon[[mon]] = adjust.years.cmip6(XAX = CMIP6.tos_mon[[mon]])
      CMIP6.tos_mon_ct[[mon]] = center.XAX_ann(XAX = CMIP6.tos_mon[[mon]], fact = 1, areaw = areaw, ref.period.years = 1961:1990)
      CMIP6.tos_mon_ct[[mon]] = extract.fraw2(XAX_scen = CMIP6.tos_mon_ct[[mon]], nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
      CMIP6.tos_mon_ct[[mon]]$Y$AGMT_ANN = CMIP6.tas_ann_ct$Y$AGMT
      CMIP6.tos_mon_ct[[mon]]$Y$AGMT_ANN_f = CMIP6.tas_ann_ct$Y$AGMT_f
      CMIP6.tos_mon_ct[[mon]]$Y$AGMT_ANN_fl = CMIP6.tas_ann_ct$Y$AGMT_fl
      
      # plot(CMIP6.tos_mon_ct[[mon]]$Y$AGMT_fl[1:1000], type='l')
      CMIP6.tos_monX_ct = CMIP6.tos_mon_ct[[mon]]
      save(list = "CMIP6.tos_monX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_ct.RData", sep=""))
   }
}







# 2.(c) READ & Process CMIP6-psl files:
# --------------------------------------------------------------------------

CMIP6.files_tas_ann = read.table("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/cmip6_historical_ssp_files.txt", sep=";", header = T)
imodall = CMIP6.files_tas_ann$modall

CMIP6.files_tas_ann = get.CMIP6.file.list(vari="tas", temp.res = "ann", scen = c("historical", "ssp245"), 
                                          CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
CMIP6.files_tas_mon = get.CMIP6.file.list(vari="tas", temp.res = "mon", scen = c("historical", "ssp245"), 
                                          CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/mon/")

CMIP6.files_psl_mon = get.CMIP6.file.list(vari="psl", temp.res = "mon", scen = c("historical", "ssp245"), 
                                          CMIP6.dir = "/net/cfc/cmip6/Next_Generation/psl/mon/g025/")


## check if all CMIP6 files are available -> Yes they are!
any(!(CMIP6.files_tas_mon$modall %in% CMIP6.files_psl_mon$modall)) ## all CMIP6 psl files needed are available !!


# i. Process psl files to 5x5° grid:
registerDoParallel(cores=30)
{
  foreach(i=1:length(CMIP6.files_tas_mon$modall)) %dopar% {
    print(i)
    
    # (1) Regrid to 5d00-72x37 grid (-> HadSLP2 -> uninterpolated):
    psl_ann_file = paste("psl_ann_", CMIP6.files_tas_mon$modall[i], "_g025.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00_72x37.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/psl/ann/g025/", psl_ann_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/ann/g72x37/", psl_ann_file, sep="")))
    
    psl_mon_file = paste("psl_mon_", CMIP6.files_tas_mon$modall[i], "_g025.nc", sep="")
    system(paste("cdo -O -remapcon,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00_72x37.txt", 
                 paste("/net/cfc/cmip6/Next_Generation/psl/mon/g025/", psl_mon_file, sep=""), 
                 paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/", psl_mon_file, sep="")))
  }
}



## ii. Read psl files and attach the "correct" target variables:
{
  setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/ann/g72x37/")
  load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_ct.RData")
  
  # CMIP6-psl-ann
  CMIP6.files_psl_ann = get.CMIP6.file.list(vari="psl", temp.res = "ann", scen = c("historical", "ssp245"), 
                                            CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/ann/g72x37/")
  CMIP6_psl_ann = read.CMIP6_novar(file.name = CMIP6.files_psl_ann$file.name, scen = CMIP6.files_psl_ann$scen, var="psl", res = "ann", 
                                   CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/ann/g72x37/")
  CMIP6.psl_ann = get.XAX_ann(X = CMIP6_psl_ann, ncol = 72 * 37, M = CMIP6.files_psl_ann, start.year = 1850, cmip = "CMIP6", areaw = rep(1, 72*37))
  CMIP6.psl_ann = adjust.years.cmip6(XAX = CMIP6.psl_ann)
  CMIP6.psl_ann_ct = center.XAX_ann(XAX = CMIP6.psl_ann, fact = 1, areaw = rep(1, 72*37), ref.period.years = 1961:1990)
  CMIP6.psl_ann_ct = extract.fraw2(XAX_scen = CMIP6.psl_ann_ct, nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
  ## CHECK WHETHER SAME FILE ORDER:
  all(CMIP6.psl_ann_ct$M$mod == CMIP6.tas_ann_ct$M$mod & CMIP6.psl_ann_ct$M$ens.mem == CMIP6.tas_ann_ct$M$ens.mem)
  CMIP6.psl_ann_ct$Y = CMIP6.tas_ann_ct$Y
  save(list = "CMIP6.psl_ann_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.psl_ann_ct.RData", sep=""))
  
  
  # CMIP6-tos-mon:
  setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/")
  CMIP6.files_psl_mon = get.CMIP6.file.list(vari="psl", temp.res = "mon", scen = c("historical", "ssp245"), 
                                            CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/")
  CMIP6_psl_mon = read.CMIP6_novar(file.name = CMIP6.files_psl_mon$file.name[1:10], scen = CMIP6.files_psl_mon$scen, var="psl", res = "mon", 
                                   CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/")
  CMIP6.psl_mon = list()
  CMIP6.psl_mon_ct = list()
  for (mon in 1:12) {
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_ct.RData", sep=""))
    
    CMIP6.psl_mon[[mon]] = get.XAX_mon(X = CMIP6_psl_mon, mon = mon, ncol = 72 * 37, M = CMIP6.files_psl_mon, start.year = 1850, cmip = "CMIP6", areaw = rep(1, 72*37), rm.HIST.from.rcp26 = F)
    CMIP6.psl_mon[[mon]] = adjust.years.cmip6(XAX = CMIP6.psl_mon[[mon]])
    CMIP6.psl_mon_ct[[mon]] = center.XAX_ann(XAX = CMIP6.psl_mon[[mon]], fact = 1, areaw = rep(1, 72*37), ref.period.years = 1961:1990)
    CMIP6.psl_mon_ct[[mon]] = extract.fraw2(XAX_scen = CMIP6.psl_mon_ct[[mon]], nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)
    
    if (!(all(CMIP6.psl_mon_ct$M$mod == CMIP6.tas_monX_ct$M$mod & CMIP6.psl_mon_ct$M$ens.mem == CMIP6.tas_monX_ct$M$ens.mem))) break;
    CMIP6.psl_mon_ct[[mon]]$Y = CMIP6.tas_monX_ct$Y
    CMIP6.psl_monX_ct = CMIP6.psl_mon_ct[[mon]]
    save(list = "CMIP6.psl_monX_ct", file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.psl_mon", mon, "_ct.RData", sep=""))
  }
}




### check whether files are reasonable:
# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.psl_ann_ct.RData")
# hist(CMIP6.psl_ann_ct$X[,1000])
# areaw= rep(1, 72*37)
# test = CMIP6.psl_ann_ct$X %*% areaw
# plot(test)
# hist(test)

  


