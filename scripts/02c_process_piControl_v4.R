
# ------------------------------------------------------------------------------------
# Evaluate tas_land and tos reconstruction based on land/ocean grid cells
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 08.2022
#library(hydroGOF)
library(mvnfast)
library(ncdf4)

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_ridge_fun_v2.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/scripts/02a_load_global_observations.R")



# 00.(b) Define names.vec for for-loop processing:
names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N", 
              "TMSST_40S_40N", "TMSST_25S_25N_", "IndianOcean", "WPacific", "EPacific", "WAtlantic")

file.names = c("tas_land_predGSAT_v3", "tas_land_predGMSST_v3", "tas_land_predGMLSAT_MI_v3", "tas_land_predGMLSAT_NI_v3",
               "tas_land_predGMST_v3", "tas_land_predGMMSAT_v3", "tas_land_predTMLSAT_v3", "tas_land_predTMMSAT_v3",
               "tas_land_predTMSST_v4", "tas_land_predTMSST25_v4",
               "tas_land_predIndianOcean_v4", "tas_land_predWPacific_v4", "tas_land_predEPacific_v4", "tas_land_predWAtlantic_v4")

# Define functions:
# For filling NAs for the tos predictions:
fill.na_from_outlier <- function(Yhat, outlier.ix, transform.Yhat = F) {
  
  if (transform.Yhat == T) {
    Yhat = data.frame(cbind(Yhat$pt0[,1], Yhat$pt0[,2], Yhat$pt1[,1], Yhat$pt1[,2], Yhat$pt2[,1], Yhat$pt2[,2]))
    names(Yhat) = c("pt0.min", "pt0.1se", "pt1.min", "pt1.1se", "pt2.min", "pt2.1se")
  }
  
  new_Yhat = data.frame(matrix(data = NA, nrow = dim(Yhat)[1] + length(outlier.ix), ncol = dim(Yhat)[2]))
  names(new_Yhat) = names(Yhat)
  
  if (length(CMIP6.df$outlier.ix) == 0) {
    new_Yhat = Yhat
  } else {
  new_Yhat[-outlier.ix,] = Yhat
  }
  return(new_Yhat)
}

# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
# XAX = CMIP6.tas_ann_piControl_ct
# str(CMIP6.tas_ann_piControl_ct)


extract.nyear.trend_from_reconstruction <- function(XAX = CMIP6.tas_ann_piControl_ct, nyears = 40, trend.sep = 10, 
                     var.names = c("TMLSAT_40S_40N"), 
                     cur.file.name = "tas_land_predTMLSAT_v3", start.year = NULL, start.year.mask = 1900) {

  # run through each model:
  mod.un = unique(paste(XAX$M$mod, "_", XAX$M$scen, "_", XAX$M$ens.mem, sep=""))
  XAX.out = list()
  XAX.out$X = matrix(nrow=1,ncol=dim(XAX$X)[2])
  XAX.out$M = XAX$M[1,]
  names(XAX.out$M) = names(XAX$M)
  XAX.out$Y = data.frame(matrix(nrow=1, ncol=length(var.names)*2))
  
  s=1
  start.ix = numeric(length=1)
  
  for (m in 1:length(mod.un)) {
    print(m)
    cur.mod = strsplit(x = mod.un[m], split = "_")[[1]][1]
    cur.scen = strsplit(x = mod.un[m], split = "_")[[1]][2]
    cur.ens = strsplit(x = mod.un[m], split = "_")[[1]][3]
    
    ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens)
    if (is.null(start.year)) {
      start.year_ = seq(head(as.numeric(XAX$M$year[ix]), 1), tail(as.numeric(XAX$M$year[ix])-nyears, 1), by = trend.sep)
    } else {
      start.year_ = start.year
    }
    cur.ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens & XAX$M$year %in% start.year_)
    if (length(cur.ix)==0) next;
    start.ix[s:(s+length(cur.ix)-1)] = cur.ix
    s = length(start.ix) + 1
  }
  
  ## Project reconstruction for each year in start.ix and for nyears, with start.year.mask:
  our.arr = array(data = NA, dim = c(12, length(start.ix), nyears))
  
  for (mon in 1:12) {
    print(mon)
    
    if (substr(cur.file.name, 1, 3) == "tas") {
    ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/CRUTEM5_mon", mon, ".RData", sep=""))
    ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
    } else if (substr(cur.file.name, 1, 3) == "tos") {
      ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/HadSST4_mon", mon, ".RData", sep=""))
      ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
    }

    # str(cmip6_piControl_tas_mon_5d00_XAX_ct$X)
    date.iix = 0
    for (date.ix in (start.year.mask-1850+1):(start.year.mask-1850+nyears)) {
      date.iix = date.iix+1
      print(date.ix)
      
        # Load trained model & Prepare OBS dataset:
      if (substr(cur.file.name, 1, 3) == "tas") {
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tas_mon_5d00_XAX_ct$X
      } else if (substr(cur.file.name, 1, 3) == "tos") {
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(HadSST4_calendar_[date.ix], "%Y"), "-", format(HadSST4_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tos_mon_5d00_XAX_ct$X
      }
      # project 
      our.arr[mon,,date.iix] = X[start.ix+date.iix,CMIP6.df$grid.ix] %*% mod_p1$beta[,1] + mod_p1$a0[1]
    }
  }
  
  # apply: 
  our.arr_ = apply(X = our.arr, MARGIN=c(2,3), FUN=mean)
  ## Now run and calculate trends based on start.ix:
  for (i in 1:(length(start.ix))) {
    print(i)
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    
    # NA check:
    if(any(is.na(XAX$Y[[var.names]][ix]))) {
      XAX.out$Y[i,] = rep(NA, length(var.names))
      XAX.out$M[i,] = XAX$M[ix[1],]
      next;
    }
    
    # calculate trends:
    XAX.out$Y[i,1] = sapply(X = var.names, FUN=function(cur.name) lm(XAX$Y[[cur.name]][ix] ~ c(1:nyears))$coefficients[2]) # c(lm(c(XAX$Y$AGMT[ix]) ~ c(1:nyears))$coefficients[2], NA, NA)
    XAX.out$Y[i,2] = lm(our.arr_[i,] ~ c(1:nyears))$coefficients[2]
    XAX.out$M[i,] = XAX$M[ix[1],]
  }
  
  names(XAX.out$Y) <- c(var.names, paste(var.names, "_pred", sep=""))
  return(XAX.out)
}



## Run from 1880 to 1950.
## Is there decadal variability in the difference in the mean-included and mean-removed simulation? -> this decadal variability would be the same pattern as the observed one...

extract.nyear.sequence_from_reconstruction <- function(XAX = CMIP6.tas_ann_piControl_ct, nyears = 40, trend.sep = 20, 
                                                    var.names = c("TMLSAT_40S_40N"), 
                                                    cur.file.name = "tas_land_predTMLSAT_v3", start.year = NULL, start.year.mask = 1900) {
  
  # run through each model:
  mod.un = unique(paste(XAX$M$mod, "_", XAX$M$scen, "_", XAX$M$ens.mem, sep=""))
  XAX.out = list()
  XAX.out$X = matrix(nrow=1,ncol=dim(XAX$X)[2])
  XAX.out$M = XAX$M[1,]
  names(XAX.out$M) = names(XAX$M)
  XAX.out$Y = data.frame(matrix(nrow=1, ncol=length(var.names)*2))
  XAX.out$Yhat = data.frame(matrix(nrow=1, ncol=length(var.names)*2))
  
  s=1
  start.ix = numeric(length=1)
  
  for (m in 1:length(mod.un)) {
    print(m)
    cur.mod = strsplit(x = mod.un[m], split = "_")[[1]][1]
    cur.scen = strsplit(x = mod.un[m], split = "_")[[1]][2]
    cur.ens = strsplit(x = mod.un[m], split = "_")[[1]][3]
    
    ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens)
    if (is.null(start.year)) {
      start.year_ = seq(head(as.numeric(XAX$M$year[ix]), 1), tail(as.numeric(XAX$M$year[ix])-nyears, 1), by = trend.sep)
    } else {
      start.year_ = start.year
    }
    cur.ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens & XAX$M$year %in% start.year_)
    if (length(cur.ix)==0) next;
    start.ix[s:(s+length(cur.ix)-1)] = cur.ix
    s = length(start.ix) + 1
  }
  
  ## Project reconstruction for each year in start.ix and for nyears, with start.year.mask:
  our.arr = array(data = NA, dim = c(12, length(start.ix), nyears))
  y = array(data = NA, dim = c(12, length(start.ix), nyears))
  
  for (mon in 1:12) {
    print(mon)
    
    if (substr(cur.file.name, 1, 3) == "tas") {
      ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/CRUTEM5_mon", mon, ".RData", sep=""))
      ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
    } else if (substr(cur.file.name, 1, 3) == "tos") {
      ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/HadSST4_mon", mon, ".RData", sep=""))
      ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
    }
    
    # str(cmip6_piControl_tas_mon_5d00_XAX_ct$X)
    date.iix = 0
    for (date.ix in (start.year.mask-1850+1):(start.year.mask-1850+nyears)) {
      date.iix = date.iix+1
      print(date.ix)
      
      # Load trained model & Prepare OBS dataset:
      if (substr(cur.file.name, 1, 3) == "tas") {
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tas_mon_5d00_XAX_ct$X
      } else if (substr(cur.file.name, 1, 3) == "tos") {
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(HadSST4_calendar_[date.ix], "%Y"), "-", format(HadSST4_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tos_mon_5d00_XAX_ct$X
      }
      # project 
      our.arr[mon,,date.iix] = X[start.ix+date.iix,CMIP6.df$grid.ix] %*% mod_p1$beta[,1] + mod_p1$a0[1]
      y[mon,,date.iix] = XAX$Y[[var.names]][start.ix+date.iix]
    }
  }
  
  # apply: 
  our.arr_ = apply(X = our.arr, MARGIN=c(2,3), FUN=mean)
  y_ = apply(X = y, MARGIN=c(2,3), FUN=mean)
  
  ## Now run and calculate trends based on start.ix:
  for (i in 1:(length(start.ix))) {
    print(i)
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    
    # calculate trends:
    XAX.out$M[i,] = XAX$M[ix[1],]
  }
  
  XAX.out$Yhat = our.arr_
  XAX.out$Y = y_
  
  return(XAX.out)
}




extract.nyear.sequence_from_reconstruction_4MR <- function(XAX = CMIP6.tas_ann_piControl_ct, nyears = 70, trend.sep = 20, 
                                                       var.names = c("AGMT"), 
                                                       cur.file.name = "tas_land_predTMLSAT_v3", start.year = NULL, start.year.mask = 1900) {
  
  # run through each model:
  mod.un = unique(paste(XAX$M$mod, "_", XAX$M$scen, "_", XAX$M$ens.mem, sep=""))
  XAX.out = list()
  XAX.out$X = matrix(nrow=1,ncol=dim(XAX$X)[2])
  XAX.out$M = XAX$M[1,]
  names(XAX.out$M) = names(XAX$M)
  XAX.out$Y = data.frame(matrix(nrow=1, ncol=length(var.names)*2))
  XAX.out$Yhat = data.frame(matrix(nrow=1, ncol=length(var.names)*2))
  
  s=1
  start.ix = numeric(length=1)
  
  for (m in 1:length(mod.un)) {
    print(m)
    cur.mod = strsplit(x = mod.un[m], split = "_")[[1]][1]
    cur.scen = strsplit(x = mod.un[m], split = "_")[[1]][2]
    cur.ens = strsplit(x = mod.un[m], split = "_")[[1]][3]
    
    ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens)
    if (is.null(start.year)) {
      start.year_ = seq(head(as.numeric(XAX$M$year[ix]), 1), tail(as.numeric(XAX$M$year[ix])-nyears, 1), by = trend.sep)
    } else {
      start.year_ = start.year
    }
    cur.ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens & XAX$M$year %in% start.year_)
    if (length(cur.ix)==0) next;
    start.ix[s:(s+length(cur.ix)-1)] = cur.ix
    s = length(start.ix) + 1
  }
  
  ## Project reconstruction for each year in start.ix and for nyears, with start.year.mask:
  our.arr = array(data = NA, dim = c(12, length(start.ix), nyears))
  y = array(data = NA, dim = c(12, length(start.ix), nyears))
  
  for (mon in 1:12) {
    print(mon)
    
    if (substr(cur.file.name, 1, 3) == "tas") {
      ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/CRUTEM5_mon", mon, ".RData", sep=""))
      ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_piControl_ct.RData", sep=""))
    } else if (substr(cur.file.name, 1, 3) == "tos") {
      ## 0.1 load CRU data:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/HadSST4_mon", mon, ".RData", sep=""))
      ## 0.1 Load climate model monthly data & select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_piControl_ct.RData", sep=""))
    }
    
    # str(cmip6_piControl_tas_mon_5d00_XAX_ct$X)
    date.iix = 0
    for (date.ix in (start.year.mask-1850+1):(start.year.mask-1850+nyears)) {
      date.iix = date.iix+1
      print(date.ix)
      
      # Load trained model & Prepare OBS dataset:
      if (substr(cur.file.name, 1, 3) == "tas") {
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tas_mon_5d00_XAX_ct$X
      } else if (substr(cur.file.name, 1, 3) == "tos") {
        # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/", cur.file.name,"/", format(HadSST4_calendar_[date.ix], "%Y"), "-", format(HadSST4_calendar_[date.ix], "%m"), ".RData", sep=""))
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/02_trained_models/tos_sea.AGMT_MR/", format(HadSST4_calendar_[date.ix], "%Y"), "-", format(HadSST4_calendar_[date.ix], "%m"), ".RData", sep=""))
        X = cmip6_piControl_tos_mon_5d00_XAX_ct$X
        X_MR = X[,CMIP6.df$grid.ix] - rep.col(rowMeans(X[,CMIP6.df$grid.ix]), n = length(CMIP6.df$grid.ix))
      }
      # project 
      our.arr[mon,,date.iix] = X_MR[start.ix+date.iix,] %*% mod_p1$beta[,mod_p1$lambda.min] + mod_p1$a0[mod_p1$lambda.min]
      y[mon,,date.iix] = XAX$Y[[var.names]][start.ix+date.iix]
    }
  }
  
  # apply: 
  our.arr_ = apply(X = our.arr, MARGIN=c(2,3), FUN=mean)
  y_ = apply(X = y, MARGIN=c(2,3), FUN=mean)
  
  ## Now run and calculate trends based on start.ix:
  for (i in 1:(length(start.ix))) {
    print(i)
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    
    # calculate trends:
    XAX.out$M[i,] = XAX$M[ix[1],]
  }
  
  XAX.out$Yhat = our.arr_
  XAX.out$Y = y_
  
  return(XAX.out)
}






# remove mean:



