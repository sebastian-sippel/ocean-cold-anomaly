
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
target_dir_5d00_monthly = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/tas/ann/"



# ------------------------------------------------------------------
## Reprocess CMIP6 tas-ann historical+ssp245 files (for reference period mean-subtraction):
# ------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")

# CMIP6-tas-ann:
CMIP6.files_tas_ann = get.CMIP6.file.list(vari="tas", temp.res = "ann", scen = c("historical", "ssp245"), 
                                          CMIP6.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
# CMIP6_tas_ann = read.CMIP6_novar(file.name = CMIP6.files_tas_ann$file.name, scen = CMIP6.files_tas_ann$scen, var="tas", res = "ann", 
#                                 CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/tas/ann/")
# CMIP6.tas_ann = get.XAX_ann(X = CMIP6_tas_ann, ncol = 72 * 36, M = CMIP6.files_tas_ann, start.year = 1850, cmip = "CMIP6", areaw = areaw)
# CMIP6.tas_ann = adjust.years.cmip6(XAX = CMIP6.tas_ann)
# save(list = "CMIP6.tas_ann", file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann.RData")
# CMIP6.tas_ann = adjust.years.cmip6(XAX = CMIP6.tas_ann)
# CMIP6.tas_ann_ct = center.XAX_ann(XAX = CMIP6.tas_ann, fact = 1, areaw = areaw, ref.period.years = 1961:1990)



# 0.b) Process CMIP6 DAMIP files and regrid to 5x5°:
# ------------------------------------------------------------------

CMIP6.files_tas       = get.CMIP6.file.list(vari= "tas", temp.res = "ann", scen = c("hist-GHG", "hist-aer", "hist-nat"), CMIP6.dir = CMIP6.orig_dir)
# CMIP6.files_psl       = get.CMIP6.file.list(vari= "psl", temp.res = "ann", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp534", "ssp585", "hist-GHG", "hist-aer", "hist-nat"), CMIP6.dir = "/net/ch4/data/cmip6-Next_Generation/psl/ann/g025")
CMIP6.files = CMIP6.files_tas # check.same.2files(CMIP6.files_tas, CMIP6.files_psl)

for (i in 1:length(CMIP6.files$file.name)) { 
  print(i)
  
  ## tas files:
  tas.ifile = paste(CMIP6.orig_dir, "/tas_ann_", CMIP6.files$modall[i], "_g025.nc", sep="")
  tas.ofile = paste(target_dir_5d00_monthly, CMIP6.files$mod[i], "_", CMIP6.files$scen[i], "_", CMIP6.files$ens.mem[i], "_5d00.nc", sep="")
  ## get 5d00 annual file:
  system(paste("cdo -remapbil,/net/h2o/climphys1/sippels/_code/spatial_data/cdo_grids/Global/global_5d00.txt ", 
               tas.ifile, " ", tas.ofile, sep=""))
}


# READ CMIP6 DAMIP files:
cmip6_DAMIP_tas_ann_5d00 = read.CMIP6_novar(file.name = paste(CMIP6.files$mod, "_", CMIP6.files$scen, "_", CMIP6.files$ens.mem, "_5d00.nc", sep=""), 
                                      scen = CMIP6.files$scen, var="tas", res = "ann", 
                                      CMIP5.dir = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/tas/ann/", time.count = rep(-1, 1579))
cmip6_DAMIP_tas_ann_5d00_XAX = get.XAX_ann(X = cmip6_DAMIP_tas_ann_5d00, M = CMIP6.files, start.year = 1, cmip = "cmip6", areaw = areaw)



## MERGE historical+ssp245 + DAMIP simulations:
CMIP6.tas_ann_ALL = list()
CMIP6.tas_ann_ALL$X = rbind(CMIP6.tas_ann$X, cmip6_DAMIP_tas_ann_5d00_XAX$X)
CMIP6.tas_ann_ALL$Y = rbind(CMIP6.tas_ann$Y, cmip6_DAMIP_tas_ann_5d00_XAX$Y)
CMIP6.tas_ann_ALL$M = rbind(CMIP6.tas_ann$M, cmip6_DAMIP_tas_ann_5d00_XAX$M)

CMIP6.tas_ann_ALL = adjust.years.cmip6(XAX = CMIP6.tas_ann_ALL)
CMIP6.tas_ann_ALL_ct = center.XAX_ann(XAX = CMIP6.tas_ann_ALL, fact = 1, areaw = areaw, ref.period.years = 1961:1990)
# plot(CMIP6.tas_ann_ALL_ct$Y$AGMT)

# extract forced response:
CMIP6.tas_ann_ALL_ct = extract.fraw2(XAX_scen = CMIP6.tas_ann_ALL_ct, nmem = 3, f.var = "AGMT", incl.scen.for.hist = T)

## produce MMM-forced responses:
CMIP6.tas_ann_ALL_ct$X = NA

# save forced response dataset:
save(list = c("CMIP6.tas_ann_ALL_ct"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_ALL_ct.RData")



# --------------------------------------------------------
# produce MMM forced response:
# --------------------------------------------------------

CMIP6.tas_ann_ALL_ct_f = list()
CMIP6.tas_ann_HIST_ct_f = list()
CMIP6.tas_ann_HISTssp245_ct_f = list()
#CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist <- CMIP6.tas_ann_ALL_ct_f$AGMT_fl_hist <- CMIP6.tas_ann_ALL_ct_f$AGMT_f_histGHG <- CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histGHG <-
#  CMIP6.tas_ann_ALL_ct_f$AGMT_f_histAER <- CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histAER <- CMIP6.tas_ann_ALL_ct_f$AGMT_f_histNAT <- CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histNAT <- numeric()
#CMIP6.tas_ann_ALL_ct_f$M <- numeric()


## which models have forced response for hist-ghg, hist-nat, hist-aer, and historical?
mod.un = unique(CMIP6.tas_ann_ALL_ct$M$mod)

for (m in 1:length(mod.un)) {
  print(m)
  
  # plot(CMIP6.tas_ann_ALL_ct$Y$AGMT_f)
  mod.ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & !is.na(CMIP6.tas_ann_ALL_ct$Y$AGMT_f))
  # print(c(mod.un[m], unique(CMIP6.tas_ann_ALL_ct$M$scen[ix])))
  
  # ALL DAMIP: 
  if (all(c("historical", "hist-GHG", "hist-aer", "hist-nat") %in% unique(CMIP6.tas_ann_ALL_ct$M$scen[mod.ix]))) {
    
    # produce vectors of respective forced responses:
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "historical" )
    CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist, CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165]))
    CMIP6.tas_ann_ALL_ct_f$AGMT_fl_hist = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_hist, CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165]))
    
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "hist-GHG" )
    CMIP6.tas_ann_ALL_ct_f$AGMT_f_histGHG = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histGHG, CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165]))
    CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histGHG = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histGHG, CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165]))
    
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "hist-aer" )
    CMIP6.tas_ann_ALL_ct_f$AGMT_f_histAER = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histAER, CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165]))
    CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histAER = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histAER, CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165]))
    
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "hist-nat" )
    CMIP6.tas_ann_ALL_ct_f$AGMT_f_histNAT = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histNAT, CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165]))
    CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histNAT = data.frame(cbind(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histNAT, CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165]))
    
    # get model names:
    CMIP6.tas_ann_ALL_ct_f$M = rbind(CMIP6.tas_ann_ALL_ct_f$M, CMIP6.tas_ann_ALL_ct$M[ix[1],])
  }
  
  
  # only historical:
  if (all(c("historical") %in% unique(CMIP6.tas_ann_ALL_ct$M$scen[mod.ix]))) {
    
    # produce vectors of respective forced responses:
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "historical" )
    CMIP6.tas_ann_HIST_ct_f$AGMT_f_hist = data.frame(cbind(CMIP6.tas_ann_HIST_ct_f$AGMT_f_hist, CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165]))
    CMIP6.tas_ann_HIST_ct_f$AGMT_fl_hist = data.frame(cbind(CMIP6.tas_ann_HIST_ct_f$AGMT_fl_hist, CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165]))
  
    # get model names:
    CMIP6.tas_ann_HIST_ct_f$M = rbind(CMIP6.tas_ann_HIST_ct_f$M, CMIP6.tas_ann_ALL_ct$M[ix[1],])
  }
  
  
  # historical+ssp245:
  if (all(c("historical", "ssp245") %in% unique(CMIP6.tas_ann_ALL_ct$M$scen[mod.ix]))) {
    
    # produce vectors of respective forced responses:
    ix = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "historical" )
    ix1 = which(mod.un[m] == CMIP6.tas_ann_ALL_ct$M$mod & CMIP6.tas_ann_ALL_ct$M$scen == "ssp245" )
    
    CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist = data.frame(cbind(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist, 
                                                                 c(CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix][1:165], CMIP6.tas_ann_ALL_ct$Y$AGMT_f[ix1][1:86])))
    
    CMIP6.tas_ann_HISTssp245_ct_f$AGMT_fl_hist = data.frame(cbind(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_fl_hist, 
                                                                  c(CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix][1:165], CMIP6.tas_ann_ALL_ct$Y$AGMT_fl[ix1][1:86])))
    
    # get model names:
    CMIP6.tas_ann_HISTssp245_ct_f$M = rbind(CMIP6.tas_ann_HISTssp245_ct_f$M, CMIP6.tas_ann_ALL_ct$M[ix[1],])
  }
}


names(CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_hist) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histGHG) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histGHG) <-
  names(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histAER) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histAER) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_f_histNAT) <- names(CMIP6.tas_ann_ALL_ct_f$AGMT_fl_histNAT) <-
  CMIP6.tas_ann_ALL_ct_f$M$mod

names(CMIP6.tas_ann_HIST_ct_f$AGMT_f_hist) <- names(CMIP6.tas_ann_HIST_ct_f$AGMT_fl_hist) <- CMIP6.tas_ann_HIST_ct_f$M$mod

names(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist) <- names(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_fl_hist) <- CMIP6.tas_ann_HISTssp245_ct_f$M$mod

# save forced response estimates:
save(list = c("CMIP6.tas_ann_ALL_ct_f"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_ALL_ct_f.RData")
save(list = c("CMIP6.tas_ann_HIST_ct_f"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HIST_ct_f.RData")
save(list = c("CMIP6.tas_ann_HISTssp245_ct_f"), file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HISTssp245_ct_f.RData")


