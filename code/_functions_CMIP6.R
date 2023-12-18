




rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}







get.CMIP6.file.list <- function(vari= "tas", temp.res = "ann", scen = c("historical", "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp534", "ssp585"), 
                                CMIP6.dir = "/net/cfc/cmip6/Next_Generation/tas/ann/g025") {
  
  # 1. get  list for each scenario:
  scen.list = list()
  for (scen.idx in 1:length(scen)) {
    scen.list[[scen.idx]] = list.files(path = paste(CMIP6.dir, sep=""), pattern = scen[scen.idx])
  }
  file.name=unlist(scen.list)
  
  # get data.frame with overview of variables:
  vari=rep(vari, length(file.name))
  res=rep(temp.res, length(file.name))
  scen=sapply(strsplit(file.name,"_"),function(x) paste(x[4],collapse="_"))
  mod=sapply(strsplit(file.name,"_"),function(x) paste(x[3],collapse="_"))
  ens.mem=sapply(strsplit(file.name,"_"),function(x) paste(x[5],collapse="_"))
  modall=sapply(strsplit(file.name,"_"),function(x) paste(x[3:5],collapse="_"))
  modcl=sapply(as.character(mod),function(x) substr(x,1,3)[[1]][1])
  
  if ( scen[1] == "piControl" ) {
    period.length=rep(NA, length(file.name))
    
    for(i in 1:length(file.name)) {
      print(i)
      test=nc_open(paste(CMIP6.dir, "/", file.name, sep="")[i])
      # brick(paste(CMIP6.dir, file.name[i], sep=""), varname="tas")
      period.length[i]=dim(ncvar_get(test, varid = vari[1], start=c(1,1,1), count=c(1,1,-1))) # [3]
    }
    return(data.frame(file.name, vari, res, mod, modcl, scen, ens.mem, modall, period.length, stringsAsFactors=F))
  }

  return(data.frame(file.name, vari, res, mod, modcl, scen, ens.mem, modall, stringsAsFactors=F))
}



# get_CMIP5_array <- function(file, var="tas", time.count=-1) {
#  ncfile=nc_open(file)
#  ncdata=ncvar_get(nc = ncfile, varid=var, start=c(1,1,1), count=c(-1,-1,time.count))
#  nc_close(ncfile)
#  return(ncdata)
#}

get_CMIP5_array <- function(file, var="tas", dims = 3, time.count=-1) {
  ncfile=nc_open(file)
  ncdata=ncvar_get(nc = ncfile, varid=var, start=rep(1, dims) , count=c(rep(-1, dims-1), time.count))
  nc_close(ncfile)
  return(ncdata)
}




check.same.2files <- function(f.var1, f.var2) {
  
  ix.var1 = which(f.var1$modall %in% f.var2$modall)
  ix.var2 = which(f.var2$modall %in% f.var1$modall)
  print(all(f.var1$modall[ix.var1] == f.var2$modall[ix.var2]))
  
  # ret.list = list()
  # ret.list[[1]] = f.var1[ix.var1,]
  # ret.list[[2]] = f.var2[ix.var2,]
  return(f.var1[ix.var1,])
}






read.CMIP6_novar <- function(file.name, var, res="ann", scen, CMIP5.dir, time.count = NA, subtract.drift.piControl = F, rm.years = 100) {
  X = list()
  for (i in 1:length(file.name)) {
    print(i)
    if (scen[i] == "piControl") {  # piControl runs are shortened to 200 years...:
      if(res=="mon") { res.fact<-12 } else res.fact<-1;
      
        if (subtract.drift.piControl == F) {
          temp = apply(get_CMIP5_array(file = paste(CMIP5.dir, "/", file.name[i], sep=""), time.count = time.count[i], var=var), 3, c)
          X[[i]] = temp[,(rm.years+1) :(dim(temp)[2])]
        } else if (subtract.drift.piControl == T) {
          temp = apply(get_CMIP5_array(file = paste(CMIP5.dir, "/", file.name[i], sep=""), time.count = time.count[i], var=var), 3, c)
          X[[i]] = t(apply(X = temp[,(rm.years+1):(dim(temp)[2])], MARGIN=1, FUN=function(x) x - lm(x ~ c(1:(dim(temp)[2]-rm.years)))$fitted)) + rep.col(rowMeans(temp[,(rm.years+1):(dim(temp)[2])]), dim(temp)[2]-rm.years)
        }
    } else {
      X[[i]] = apply(get_CMIP5_array(file = paste(CMIP5.dir, "/", file.name[i], sep=""), var=var), 3, c)
      # get_CMIP5_vec(file = paste(CMIP5.dir, "/", file.name[i], sep=""), var="time")
    }
  }
  return(X)
}




# X = cmip6_tas_ann_scen_5d00
# M = CMIP5.files_tas_scen

library(raster)
raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw=c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72,36)[,36:1]))
areaw_cos = sqrt(cos(c(matrix(raster::coordinates(raster.template)[,2], 72, 36)[,36:1])*pi/180))



get.XAX_ann <- function(X, ncol = 72*36, M, start.year = 1850, cmip, areaw) {
  
  # OVERALL DIMENSION
  s = sum(sapply(X = X, FUN=function(x) dim(x)[2]))
  
  MAX <- data.frame(matrix(nrow=s,ncol=9))
  names(MAX) = c("vari", "res", "file.name", "cmip", "mod", "modcl", "scen", "ens.mem", "year")
  YAX <- data.frame(matrix(nrow=s, ncol=1))
  names(YAX) <- "AGMT"
  
  XAX <- matrix(nrow=s,ncol=ncol)
  
  cc <- 0
  for (k in 1:length(X)){
    cat("\r ",k)
    Xsc <- t(X[[k]])
    # Ysc <- Y[[k]][seq(mon, dim(X[[k]])[3], 12),]
    
    for (pc in 1:(dim(X[[k]])[2])){
      cc <- cc+1
      MAX[cc,] <- as.character(c(M$vari[k], M$res[k], M$file.name[k], cmip, M$mod[k], M$modcl[k], M$scen[k], M$ens.mem[k], start.year-1+as.numeric(pc)))
      YAX[cc,] <- Xsc[pc,] %*% areaw
      XAX[cc,]  <- Xsc[pc,] 
    }
    print(MAX[cc,])
  }
  
  # save to .RData file:
  return(list(X = XAX, Y=YAX, M=MAX))
}






## FUNCTION TO DEFINE ARRAYS:
# ------------------------------
# X = cmip5_tas_mon_5d00_monthly
# Y = cmip5_tas_mon_GLmean_anom
# M = CMIP5.files
# mon = 1


get.XAX_mon <- function(X, mon, ncol = 72*36, M, rm.HIST.from.rcp26 = F, start.year = 1850, cmip = "cmip6", areaw) {
  
  # Y = data.frame(AGMT = t(X) %*% areaw)
  
  # OVERALL DIMENSION
  s = sum(sapply(X = X, FUN=function(x) dim(x)[2]))/12
  
  MAX <- data.frame(matrix(nrow=s,ncol=10))
  names(MAX) = c("vari", "res", "file.name", "cmip", "mod", "modcl", "scen", "ens.mem", "year", "mon")
  # YAX <- data.frame(matrix(nrow=s, ncol=6))
  # names(YAX) <- names(Y)
  
  XAX <- matrix(nrow=s,ncol=ncol)
  YAX <- matrix(nrow=s,ncol=1)
  
  cc <- 0
  for (k in 1:length(X)){
    cat("\r ",k)
    Xsc <- t(X[[k]][,seq(mon, dim(X[[k]])[2], 12)])
    # Xsc <- t(apply(X[[k]][,seq(mon, dim(X[[k]])[2], 12)],3,as.vector))
    Ysc <- Xsc %*% areaw
    
    for (pc in 1:(dim(X[[k]])[2]/12)){
      cc <- cc+1
      MAX[cc,] <- as.character(c(M$vari[k], M$res[k], M$file.name[k], cmip, M$mod[k], M$modcl[k], M$scen[k], M$ens.mem[k], start.year-1+as.numeric(pc), mon))
      YAX[cc,] <- Ysc[pc,]
      XAX[cc,]  <- Xsc[pc,]
    }
    print(MAX[cc,])
  }
  
  if (rm.HIST.from.rcp26 == T) {
    rm.idx = which(MAX$year < 2006 & MAX$scen == "rcp26")
    MAX = MAX[-rm.idx,]
    YAX = YAX[-rm.idx,]
    XAX = XAX[-rm.idx,]
  }
  
  # XAX = as.big.matrix(x = XAX, type = "double", separated = F, shared = T, backingfile = paste("XAX_", mon,".bin", sep=""), descriptorfile = paste("XAX_", mon,".desc", sep=""))
  # X.bm = as.big.matrix(x = t(values(X.RB)), type = "double", separated = FALSE, shared = TRUE,
  #                     backingfile = "X.bin", descriptorfile = "X.desc")
  
  # save to .RData file:
  return(list(X = XAX, Y=data.frame(AGMT=YAX), M=MAX))
}




adjust.years.cmip6 <- function(XAX) {
  file.un = unique(XAX$M$file.name)
  for (i in 1:length(file.un)) {
    # print(i)
    ix = which(XAX$M$file.name == file.un[i])
    if (XAX$M$scen[ix[1]] == "historical") {
      XAX$M$year[ix] = 1850:2014
    } else if (XAX$M$scen[ix[1]] %in% c("hist-GHG", "hist-aer", "hist-nat")) {
      XAX$M$year[ix] = 1850:(1850-1+length(ix))
    } else if (XAX$M$scen[ix[1]] %in% c("ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp534", "ssp585")) {
      XAX$M$year[ix] = 2015:2100
      if (length(XAX$M$year[ix]) != 86) print(i)
    }
  }
  return(XAX)
}


adjust.years.LENS <- function(XAX) {
  file.un = unique(XAX$M$file.name)
  for (i in 1:length(file.un)) {
    # print(i)
    ix = which(XAX$M$file.name == file.un[i])
    if (XAX$M$mod[ix[1]] == "MPI-ESM") {
      XAX$M$year[ix] = 1850:2099
    } else if (XAX$M$mod[ix[1]] == c("CESM1-CAM5")) {
      XAX$M$year[ix] = 1920:2100
    } else if (XAX$M$mod[ix[1]] == c("CSIRO-Mk3-6-0")) {
      XAX$M$year[ix] = 1850:2100
    } else if (XAX$M$mod[ix[1]] == c("CanESM2")) {
      XAX$M$year[ix] = 1950:2100
    } else if (XAX$M$mod[ix[1]] == c("EC-EARTH")) {
      XAX$M$year[ix] = 1860:2100
    } else if (XAX$M$mod[ix[1]] == c("GFDL-CM3")) {
      XAX$M$year[ix] = 1920:2100
    } else if (XAX$M$mod[ix[1]] == c("GFDL-ESM2M")) {
      XAX$M$year[ix] = 1950:2100
    }
  }
  return(XAX)
}


## Center CMIP5 scenarios with control run:
center.cmip.scen_with_piC <- function(XAX_piC, XAX_hist, rm.yr = 100, rm.na = F) {
  
  hist.mod.un = unique(XAX_hist$M$mod)
  
  for (i in 1:length(hist.mod.un)) {
    # print(i)
    
    hist.cur.mod = hist.mod.un[i]
    piC.ix = which(XAX_piC$M$mod == hist.cur.mod)[-c(1:rm.yr)]
    hist.ix = which(XAX_hist$M$mod == hist.cur.mod)
    
    
    print(paste(hist.cur.mod, " piC-length: ", length(piC.ix), " Number of hist-scen: ", length(unique(XAX_hist$M$ens.mem[hist.ix])), sep=""))
    
    AGMT.avg = mean(XAX_piC$Y$AGMT[piC.ix])
    pattern.avg = colMeans(XAX_piC$X[piC.ix,])
    # image.plot(matrix(pattern.avg, 72, 36))
    
    ## Process pattern in historical files:
    XAX_hist$X[hist.ix,] = XAX_hist$X[hist.ix,] - rep.row(x = pattern.avg, n = dim(XAX_hist$X[hist.ix,])[1])
    XAX_hist$Y[hist.ix,] = XAX_hist$Y[hist.ix,] - rep.row(x = AGMT.avg, n = length(XAX_hist$Y$AGMT[hist.ix]))
    
  }
  
  if (rm.na == T) {
      na.ix=which(is.na(XAX_hist$Y$AGMT))
      XAX_hist$Y <- data.frame(XAX_hist$Y[-na.ix,1])
      colnames(XAX_hist$Y) <- "AGMT"
      XAX_hist$X = XAX_hist$X[-na.ix,]
      XAX_hist$M = XAX_hist$M[-na.ix,]
  }
  
  return(XAX_hist)
}





## Center each ensemble member with historical period but over all realizations:
# XAX_scen = cmip6_ann_scen_5d00_XAX
# cmip5_ann_scen_5d00_XAX_ct$M$mod.phys.bool[which(cmip5_ann_scen_5d00_XAX_ct$Y$AGMT > 200)]
center.ensemble.member.cmip <- function(XAX_scen, ref.scen = "historical") {
  
  mod.phys.un = unique(XAX_scen$M$mod.phys)
  XAX_scen.out = XAX_scen
  mod.phys.bool = rep(NA, length(XAX_scen$M$mod.phys))
  
  # run through each mod.phys.un:
  for (m in 1:length(mod.phys.un)) {
    scen.un = unique(XAX_scen$M$scen[which(XAX_scen$M$mod.phys == mod.phys.un[m])])
    ens.mem.un = unique(XAX_scen$M$ens.mem[which(XAX_scen$M$mod.phys == mod.phys.un[m])])
    print(paste(m, mod.phys.un[m], sep=" "))
    print(scen.un)
    print(ens.mem.un)
    
    # for (em in 1:length(ens.mem.un)) {
      all.ix = which(XAX_scen$M$mod.phys == mod.phys.un[m])
      ref.ix = which(XAX_scen$M$mod.phys == mod.phys.un[m] & XAX_scen$M$scen %in% ref.scen & XAX_scen$M$year %in% 1870:1920)
      
      if (length(ref.ix) == 0) {
        mod.phys.bool[all.ix] = F
        ref.ix = which(XAX_scen$M$mod == substring(mod.phys.un[m], 1, nchar(mod.phys.un[m])-3)  & XAX_scen$M$scen %in% ref.scen & XAX_scen$M$year %in% 1870:1920)
        # if (length(ref.ix) == 0) ref.ix = which(XAX_scen$M$mod == substring(mod.phys.un[m], 1, nchar(mod.phys.un[m])-3) & XAX_scen$M$scen %in% ref.scen & XAX_scen$M$year %in% 1870:1920)
      } else {
        mod.phys.bool[all.ix] = T  # physics preserved in mean-subtraction
      }
      # Center based on ref. period:
      print(length(ref.ix))
      XAX_scen.out$Y$AGMT[all.ix] = XAX_scen$Y$AGMT[all.ix] - mean(XAX_scen$Y$AGMT[ref.ix])
      XAX_scen.out$tas[all.ix,] = XAX_scen$tas[all.ix,] - rep.row(x = colMeans(XAX_scen$tas[ref.ix,]), n = length(all.ix))
      # XAX_scen.out$psl[all.ix,] = XAX_scen$psl[all.ix,] - rep.row(x = colMeans(XAX_scen$psl[ref.ix,]), n = length(all.ix))
  }
  print(any(is.na(mod.phys.bool)))
  
  XAX_scen.out$M$mod.phys.bool = mod.phys.bool
  return(XAX_scen.out)
  
  # plot(cmip5_ann_piControl_5d00_XAX_ct$Y$AGMT[which(cmip5_ann_piControl_5d00_XAX_ct$M$mod == "GISS-E2-H" & cmip5_ann_piControl_5d00_XAX_ct$M$scen == "piControl_DT")])
  # unique(XAX_scen$M$mod)
  # unique(XAX_scen$M$mod[which(XAX_scen$M$scen == "historicalANT")])
  # unique(XAX_scen$M$mod[which(XAX_scen$M$scen == "historicalNat")])
  # unique(XAX_scen$M$mod[which(XAX_scen$M$scen == "historicalGHG")])
}



## Center each model with historical period
# 06.04.2021
center.XAX_ann <- function(XAX, fact = 24 * 3600, areaw, ref.period.years = 1950:2000, ref.scen = "historical") {
  XAX_norm = XAX
  XAX_norm$X[,] = NA
  XAX_norm$Y$AGMT = rep(NA, length(XAX_norm$Y$AGMT))
  
  # subtract reference period average:
  mod.un = unique(XAX$M$mod)
  
  for (m in 1:length(mod.un)) {
    print(paste(mod.un[m]))
    mod.ix = which(XAX$M$mod == mod.un[m])
    ref.mod.ix = which(XAX$M$mod == mod.un[m] & XAX$M$year %in% ref.period.years & XAX$M$scen == ref.scen)
    test = apply(X = XAX$X[ref.mod.ix,], MARGIN = 2, FUN = mean)
    XAX_norm$X[mod.ix,] = t(t(XAX$X[mod.ix,]) - test) * fact
    XAX_norm$Y$AGMT[mod.ix] = c(XAX_norm$X[mod.ix,] %*% areaw)
  }
  return(XAX_norm)
}





## get Santer et al. fingerprint.
# 06.04.2021
get.EOF.FP <- function(XAX, FP.scen = "ssp245", areaw_cos) {
  
  ssp245.modcl = unique(XAX$M$modcl[which(XAX$M$scen == FP.scen)])
  year.un = as.numeric(unique(XAX$M$year))
  
  ssp245.fingerprint = list()
  for (m in 1:length(ssp245.modcl)) {
    print(ssp245.modcl[m])
    
    hist.years = as.numeric(unique(XAX$M$year[which(ssp245.modcl[m] == XAX$M$modcl & XAX$M$scen %in% c("historical"))]))
    ssp.years = as.numeric(unique(XAX$M$year[which(ssp245.modcl[m] == XAX$M$modcl & XAX$M$scen %in% c(FP.scen))]))
    
    ## Fingerprint historical years:
    if ( length(which(ssp245.modcl[m] == XAX$M$modcl & 1900 == XAX$M$year & XAX$M$scen %in% c("historical"))) > 1 ) {
      hist.fingerprint = sapply(X = 1:length(hist.years), FUN=function(y) colMeans(XAX$X[which(ssp245.modcl[m] == XAX$M$modcl & hist.years[y] == XAX$M$year & XAX$M$scen %in% c("historical")),]))
    } else {
      hist.fingerprint = sapply(X = 1:length(hist.years), FUN=function(y) XAX$X[which(ssp245.modcl[m] == XAX$M$modcl & hist.years[y] == XAX$M$year & XAX$M$scen %in% c("historical")),])
    }
    ## Fingerprint ssp years:
    if ( length(which(ssp245.modcl[m] == XAX$M$modcl & 2020 == XAX$M$year & XAX$M$scen %in% c(FP.scen))) > 1 ) {
      ssp.fingerprint = sapply(X = 1:length(ssp.years), FUN=function(y) colMeans(XAX$X[which(ssp245.modcl[m] == XAX$M$modcl & ssp.years[y] == XAX$M$year & XAX$M$scen %in% c(FP.scen)),]))
    } else {
      ssp.fingerprint = sapply(X = 1:length(ssp.years), FUN=function(y) XAX$X[which(ssp245.modcl[m] == XAX$M$modcl & ssp.years[y] == XAX$M$year & XAX$M$scen %in% c(FP.scen)),])
    }
    ssp245.fingerprint[[m]] = cbind(hist.fingerprint, ssp.fingerprint)
    names(ssp245.modcl[m]) = ssp245.modcl[m]
  }
  
  ssp245 = sapply(X = 1:length(c(hist.years, ssp.years)), FUN=function(y) rowMeans(sapply(X = ssp245.fingerprint, FUN=function(x) x[,y]), na.rm=T))
  test.aw = t(scale(t(ssp245), T, F)) * rep.col(areaw_cos, n = 251)
  ssp245.svd = svd(test.aw)   # rcp85.svd_ = svd(test.aw_orig)
  plot(ssp245.svd$v[,1])
  image.plot(matrix(ssp245.svd$u[,1], 144,73))
  
  
  ## OLD AND WRONG:  # test.aw = (ssp585)  * rep.col(areaw, 213)  # ssp585.svd = svd(test.aw)
  # image.plot(matrix(ssp585.svd$u[,1] / areaw, 72, 36))
  y_Fssp245 = c( as.matrix(XAX$X) %*% ( ssp245.svd$u[,1] * areaw_cos * (-1)) )
  XAX$Y$EOF245 = y_Fssp245
  XAX$EOF_fingerprint_svd = ssp245.svd
  XAX$forced_fingerprint_by_mod = ssp245.fingerprint
  XAX$forced_fingerprint_all_mod = ssp245
  XAX$areaw_cos = areaw_cos
  return(XAX)
}



calculate.trend.slope.XAX <- function(XAX, var = "AGMT_f", years = 1950:2014, scen = c("historical", "ssp245"), modcl = F) {
  mod.un = unique(XAX$M$mod[which(!is.na(XAX$Y[[var]]) & XAX$M$scen %in% scen)])
  
  trends.df = data.frame(t(rep(NA, 8)))
  names(trends.df) = c("modcl", "mod", "scen", "ens.mem", "start.year", "end.year", "slope", "intercept")
  
  i = 0
  for (m in 1:length(mod.un)) {
    print(mod.un[m])
    ens.mem.un = unique(XAX$M$ens.mem[which(!is.na(XAX$Y[[var]]) & XAX$M$scen %in% scen & XAX$M$mod == mod.un[m])])
    
    if (modcl == T) ens.mem.un = ens.mem.un[1]
    
    for (e in 1:length(ens.mem.un)) {
      i = i+1
      ix = which(!is.na(XAX$Y[[var]]) & XAX$M$scen %in% scen & XAX$M$mod == mod.un[m] & 
                   XAX$M$ens.mem == ens.mem.un[e] & XAX$M$year %in% years)
      lm.mod = lm(XAX$Y[[var]][ix] ~ years)
      trends.df[i,1:4] = c(substr(mod.un[m], 1, 3), mod.un[m], unique(XAX$M$scen[ix])[1], ens.mem.un[e]) 
      trends.df[i,5:8] = c(years[1], tail(years,1), lm.mod$coefficients[2], lm.mod$coefficients[1])
    }
  }
  return(trends.df)
}



## Extract simplified forced response
# 28.08.2020
# XAX_scen = Rx1day_ann_ct
extract.fraw2 <- function(XAX_scen, nmem = 3, f.var = "EOF245", incl.scen.for.hist = F) {
  
  scen.un = unique(XAX_scen$M$scen)
  l = dim(XAX_scen$M)[1]
  
  # Extract "raw" forced response:
  XAX_scen$Y[[paste(f.var, "_f", sep="")]] = rep(NA, l)
  XAX_scen$Y[[paste(f.var, "_fl", sep="")]] = rep(NA, l)
  XAX_scen$Y[[paste(f.var, "_f_noinclmem", sep="")]] = rep(NA, l)
  # XAX_scen$Y$nmem = rep(NA, l)
  
  
  # XAX_scen$Y$ftot_l = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.10y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.20y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.30y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.50y = rep(NA, l)
  
  
  for (s in 1:length(scen.un)) {
    print(paste(s, scen.un[s], sep=" "))
    cur.scen = scen.un[s]
    
    # Determine number of members per model:
    no.mem.mod = sapply(X = unique(XAX_scen$M$mod), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod & XAX_scen$M$scen == cur.scen)])))
    # no.mem = sapply(X = unique(XAX_scen$M$mod.phys), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod.phys & XAX_scen$M$scen == cur.scen)])))
    # length(which(no.mem.mod >= 3))
    # length(which(no.mem >= 3))
    
    cur.no.mem = no.mem.mod[which(no.mem.mod >= nmem)]
    cur.mod.phys.un = names(which(no.mem.mod >= nmem))

    for (m in 1:length(cur.mod.phys.un)) {
      print(paste(m, cur.mod.phys.un[m], sep=" "))
      
      # all scenarios:
      ix_scen = which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen)
      ens.mat = sapply(X = unique(XAX_scen$M$ens.mem[ix_scen]), FUN=function(cur.ens) XAX_scen$Y[[f.var]][which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen & XAX_scen$M$ens.mem == cur.ens)])
      if (any(is.na(ens.mat))) next 
      
      if (cur.no.mem[m] > 2) {
        ftot_IV = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat[,-k], na.rm=T))
      } else {
        ftot_IV = matrix(NA, nrow = dim(ens.mat)[1], ncol = dim(ens.mat)[2])
      }
      
      ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
      # ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
      # print(dim(ftot))
      # Fit LOWESS: fr.loess_0.75 = loess(fr.raw ~ c(1:231), span = 0.75, degree = 2)$fittedn=dim(ens.mat)[1]
      # equivalent to 10-year running mean: 15 / 251 # 15 / 165
      n = dim(ens.mat)[1]
      ftot_l = rep.col(loess(rowMeans(ftot) ~ c(1:n), span = 0.0909, degree = 1, surface = "direct")$fitted, dim(ens.mat)[2])
      if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
      
      if (cur.scen %in% c("ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585")) {
        ix_hist = c(which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical"))
        ens.mat.hist = sapply(X = unique(XAX_scen$M$ens.mem[ix_hist]), FUN=function(cur.ens) XAX_scen$Y[[f.var]][which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical" & XAX_scen$M$ens.mem == cur.ens)])
        ftot_hist = c(sapply(X = 1:dim(ens.mat.hist)[2], FUN=function(k) rowMeans(ens.mat.hist, na.rm=T))[,1])
        ftot_long = c(ftot_hist, rowMeans(ftot))
        ftot_l <- ftot_l_all <- loess(ftot_long ~ c(1:(165+n)), span = 0.0597, degree = 1, surface = "direct")$fitted
        if (is.na(ftot_long[1])) ftot_l = c(NA, ftot_l)
        # if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
        ftot_l = rep.col(ftot_l[-c(1:165)], dim(ens.mat)[2])
        
        
        # change historical with respect to ssp scenario:
        if (incl.scen.for.hist == T) {
        ftot_l_hist = rep.col(ftot_l_all[c(1:165)], dim(ens.mat.hist)[2])
        XAX_scen$Y[[paste(f.var, "_fl", sep="")]][ix_hist] = c(ftot_l_hist)
        }
      }
      
      ## Fill in forced responses and smoothed versions:
      XAX_scen$Y[[paste(f.var, "_f", sep="")]][ix_scen] = c(ftot)
      XAX_scen$Y[[paste(f.var, "_f_noinclmem", sep="")]][ix_scen] = c(ftot_IV)
      XAX_scen$Y[[paste(f.var, "_fl", sep="")]][ix_scen] = c(ftot_l)
      
      # IV.fraw = ens.mat - ftot_IV
      # XAX_scen$Y$IV.fraw.GMT.10y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 10, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.20y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 20, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.30y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 30, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.50y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 50, fill = NA))
    }
  }
  
  # which(is.na(XAX_scen$Y$EOF245_f))
  # plot(XAX_scen$Y$EOF245_fl[50000:100000], type='l')
  # XAX_scen$M$scen[20001:21000]
  
  return(XAX_scen)
}



extract.mod.df <- function(XAX_scen) {
  
  # data.frame of number of ensemble members per model and scenario combination:
  mod.un = unique(XAX_scen$M$mod)
  scen.un = unique(XAX_scen$M$scen)
  mod.df = data.frame(mod = NA, modcl = NA, scen = NA, nmem = NA)
  
  
  for (s in 1:length(scen.un)) {
    print(paste(s, scen.un[s], sep=" "))
    cur.scen = scen.un[s]
    
    # Determine number of members per model:
    no.mem.mod = sapply(X = unique(XAX_scen$M$mod), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod & XAX_scen$M$scen == cur.scen)])))
    # no.mem = sapply(X = unique(XAX_scen$M$mod.phys), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod.phys & XAX_scen$M$scen == cur.scen)])))
    # length(which(no.mem.mod >= 3))
    # length(which(no.mem >= 3))
    
    cur.no.mem = no.mem.mod
    cur.mod.phys.un = names((no.mem.mod))
  
    mod.df = rbind(mod.df, data.frame(mod = cur.mod.phys.un, modcl = substr(x = cur.mod.phys.un, 1, 3), scen = cur.scen, nmem = cur.no.mem))
  }
  return(mod.df[-1,])
}
  
  
  
  # data.frame of individual model+ensemble member combination:
  # train.mod.ens = unique(paste(XAX_scen$M$mod, XAX_scen$M$ens.mem, sep="_"))
  # mod = XAX_scen$M$mod
  # ens.mem = XAX_scen$M$ens.mem

  # for (me in 1:length(train.mod.ens)) {
#    print(me)
#    cur.mod = strsplit(x = train.mod.ens[me], split = "_")[[1]][1]
#    cur.ens.mem = strsplit(x = train.mod.ens[me], split = "_")[[1]][2]
#    cur.modcl = substr(cur.mod, 1, 3)
#    n = length(which(cur.mod == mod & cur.ens.mem == ens.mem))
#  }
  
  # strsplit(x = train.mod.ens, split = "_")[[1]][1]
  




## Extract forced response from CMIP6 file:
# XAX_scen = cmip6_ann_scen_5d00_XAX_ct
# ERF.df = ERF_FAIR

extract.fraw <- function(XAX_scen, ERF.df, nmem = 3) {
  
  scen.un = unique(XAX_scen$M$scen)
  l = dim(XAX_scen$M)[1]
  
  # Extract "raw" forced response:
  XAX_scen$Y$ftot = rep(NA, l)
  XAX_scen$Y$ftot_l = rep(NA, l)
  XAX_scen$Y$IV.fraw.GMT.10y = rep(NA, l)
  XAX_scen$Y$IV.fraw.GMT.20y = rep(NA, l)
  XAX_scen$Y$IV.fraw.GMT.30y = rep(NA, l)
  XAX_scen$Y$IV.fraw.GMT.50y = rep(NA, l)
  
  # Add contribution to forced response from D&A runs:
  XAX_scen$Y$fant = rep(NA, l)
  XAX_scen$Y$fghg = rep(NA, l)
  XAX_scen$Y$faer = rep(NA, l)
  XAX_scen$Y$fnat = rep(NA, l)
  
  # Add Radiative Forcing from FAIR-model.
  # XAX_scen$Y$ERF_nat  = rep(NA, l)

  XAX_scen$Y$ERF_tot_r3  = rep(NA, l)
  XAX_scen$Y$ERF_ant_r3  = rep(NA, l)
  XAX_scen$Y$ERF_nat_r3  = rep(NA, l)
  XAX_scen$Y$ERF_aer_r3  = rep(NA, l)
  XAX_scen$Y$ERF_ghg_r3  = rep(NA, l)
  
  # XAX_scen$Y$fant_natERF  = rep(NA, l)
  XAX_scen$Y$ftot_ERF  = rep(NA, l)
  XAX_scen$Y$fant_ERF  = rep(NA, l)
  XAX_scen$Y$fghg_ERF  = rep(NA, l)
  XAX_scen$Y$faer_ERF  = rep(NA, l)
  XAX_scen$Y$fnat_ERF  = rep(NA, l)

    
  for (s in 1:length(scen.un)) {
    print(paste(s, scen.un[s], sep=" "))
    cur.scen = scen.un[s]
    
    # Determine number of members per model:
    no.mem.mod = sapply(X = unique(XAX_scen$M$mod), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod & XAX_scen$M$scen == cur.scen)])))
    no.mem = sapply(X = unique(XAX_scen$M$mod.phys), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod.phys & XAX_scen$M$scen == cur.scen)])))
    # length(which(no.mem.mod >= 3))
    # length(which(no.mem >= 3))
    
    cur.mod.phys.un = names(which(no.mem >= nmem))
    cur.no.mem = no.mem[which(no.mem >= nmem)]
    if (cur.scen %in% c("hist-aer", "hist-GHG", "hist-nat")) {
      cur.mod.phys.un = names(which(no.mem >= 1))
      cur.no.mem = no.mem[which(no.mem >= 1)]
    } 
    
    for (m in 1:length(cur.mod.phys.un)) {
      print(paste(m, cur.mod.phys.un[m], sep=" "))
      
      # all scenarios:
        ix_scen = which(XAX_scen$M$mod.phys == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen)
        ens.mat = sapply(X = unique(XAX_scen$M$ens.mem[ix_scen]), FUN=function(cur.ens) XAX_scen$Y$AGMT[which(XAX_scen$M$mod.phys == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen & XAX_scen$M$ens.mem == cur.ens)])
          if (cur.no.mem[m] > 2) {
            ftot_IV = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat[,-k], na.rm=T))
          } else {
            ftot_IV = matrix(NA, nrow = dim(ens.mat)[1], ncol = dim(ens.mat)[2])
          }

        ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
        # ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
        # print(dim(ftot))
        # Fit LOWESS: fr.loess_0.75 = loess(fr.raw ~ c(1:231), span = 0.75, degree = 2)$fittedn=dim(ens.mat)[1]
        # 231 * 0.25 / n
        n = dim(ens.mat)[1]
        ftot_l = rep.col(loess(rowMeans(ftot) ~ c(1:n), span = 231 * 0.25 / n, degree = 2)$fitted, dim(ens.mat)[2])
        if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
        
      if (cur.scen %in% c("ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp534-over", "ssp585")) {
        ix_hist = c(which(XAX_scen$M$mod.phys == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical"))
        ens.mat.hist = sapply(X = unique(XAX_scen$M$ens.mem[ix_hist]), FUN=function(cur.ens) XAX_scen$Y$AGMT[which(XAX_scen$M$mod.phys == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical" & XAX_scen$M$ens.mem == cur.ens)])
        ftot_hist = c(sapply(X = 1:dim(ens.mat.hist)[2], FUN=function(k) rowMeans(ens.mat.hist, na.rm=T))[,1])
        ftot_long = c(ftot_hist, rowMeans(ftot))
        ftot_l = loess(ftot_long ~ c(1:(165+n)), span = 231 * 0.25 / n, degree = 2)$fitted
        if (is.na(ftot_long[1])) ftot_l = c(NA, ftot_l)
        # if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
        ftot_l = rep.col(ftot_l[-c(1:165)], dim(ens.mat)[2])
      }
        
        ## Fill in forced responses and smoothed versions:
        XAX_scen$Y$ftot[ix_scen] = c(ftot)
        XAX_scen$Y$ftot_l[ix_scen] = c(ftot_l)
        
        IV.fraw = ens.mat - ftot_IV
        XAX_scen$Y$IV.fraw.GMT.10y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 10, fill = NA))
        XAX_scen$Y$IV.fraw.GMT.20y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 20, fill = NA))
        XAX_scen$Y$IV.fraw.GMT.30y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 30, fill = NA))
        XAX_scen$Y$IV.fraw.GMT.50y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 50, fill = NA))
        
        ## Fill in ERF from FAIR dataset:
        cur.scen.erf = cur.scen
        
        if (cur.scen == "hist-GHG") {
          cur.scen.erf <- "historical"
          ERF.ix = match(x = as.numeric(XAX_scen$M$year[ix_scen]), table = ERF.df[[cur.scen.erf]]$year)
          XAX_scen$Y$ERF_nat[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_nat_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_aer_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_ghg_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$co2_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$ch4_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$n2o_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$other_wmghg_r3[ERF.ix] + 
            ERF.df[[cur.scen.erf]]$o3_tropospheric_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$o3_stratospheric_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$h2o_stratospheric_r3[ERF.ix]
          XAX_scen$Y$ERF_tot_r3[ix_scen]  = XAX_scen$Y$ERF_ghg_r3[ix_scen]
          XAX_scen$Y$ERF_ant_r3[ix_scen]  = XAX_scen$Y$ERF_ghg_r3[ix_scen]
        } else if (cur.scen == "hist-aer") {
          cur.scen.erf <- "historical"
          ERF.ix = match(x = as.numeric(XAX_scen$M$year[ix_scen]), table = ERF.df[[cur.scen.erf]]$year)
          XAX_scen$Y$ERF_nat[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_nat_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_aer_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$aerosol.radiation_interactions_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$aerosol.cloud_interactions_r3[ERF.ix]
          XAX_scen$Y$ERF_ghg_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_tot_r3[ix_scen]  = XAX_scen$Y$ERF_aer_r3[ix_scen]
          XAX_scen$Y$ERF_ant_r3[ix_scen]  = XAX_scen$Y$ERF_aer_r3[ix_scen]
        } else if (cur.scen == "hist-nat") {
          cur.scen.erf <- "historical"
          ERF.ix = match(x = as.numeric(XAX_scen$M$year[ix_scen]), table = ERF.df[[cur.scen.erf]]$year)
          XAX_scen$Y$ERF_nat[ix_scen]  = ERF.df[[cur.scen.erf]]$total_natural[ERF.ix]
          XAX_scen$Y$ERF_nat_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$total_natural_r3[ERF.ix]
          XAX_scen$Y$ERF_aer_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_ghg_r3[ix_scen]  = rep(0, length(ix_scen))
          XAX_scen$Y$ERF_tot_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$total_natural_r3[ERF.ix]
          XAX_scen$Y$ERF_ant_r3[ix_scen]  = rep(0, length(ix_scen))
        }  else {
        ERF.ix = match(x = as.numeric(XAX_scen$M$year[ix_scen]), table = ERF.df[[cur.scen.erf]]$year)
        XAX_scen$Y$ERF_nat[ix_scen]  = ERF.df[[cur.scen.erf]]$total_natural[ERF.ix]
        XAX_scen$Y$ERF_tot_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$total_r3[ERF.ix]
        XAX_scen$Y$ERF_ant_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$total_anthropogenic_r3[ERF.ix]
        XAX_scen$Y$ERF_nat_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$total_natural_r3[ERF.ix]
        XAX_scen$Y$ERF_aer_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$aerosol.radiation_interactions_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$aerosol.cloud_interactions_r3[ERF.ix]
        XAX_scen$Y$ERF_ghg_r3[ix_scen]  = ERF.df[[cur.scen.erf]]$co2_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$ch4_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$n2o_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$other_wmghg_r3[ERF.ix] + 
          ERF.df[[cur.scen.erf]]$o3_tropospheric_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$o3_stratospheric_r3[ERF.ix] + ERF.df[[cur.scen.erf]]$h2o_stratospheric_r3[ERF.ix]
        }
      }
    }    
    

  # ---------------------------------------------------------------------------------
  # ADD D&A runs: T_ant, T_nat, T_aer, T_GHG based on D&A runs... 
  # Forcing Regression to predict ftot_ERF, fant_ERF, fghg_ERF, faer_ERF: 
  DA.mod.phys.un = unique(XAX_scen$M$mod.phys[which(XAX_scen$M$scen %in% c("hist-aer"))])
  # unique(XAX_scen$M$mod.phys)
  
  for  (m in 1:length(DA.mod.phys.un)) {
    
    ix_aer = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & XAX_scen$M$scen == "hist-aer")
    ix_nat = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & XAX_scen$M$scen == "hist-nat")
    ix_ghg = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & XAX_scen$M$scen == "hist-GHG")
    ix_historical = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & XAX_scen$M$scen == "historical")
    
    # plot(XAX_scen$Y$ftot[ix_aer])

    ## Forced Aerosol Response:
    XAX_scen$Y$faer[ix_aer] = XAX_scen$Y$ftot[ix_aer]
    XAX_scen$Y$faer[ix_nat] = rep(0, length(ix_nat))
    XAX_scen$Y$faer[ix_ghg] = rep(0, length(ix_ghg))
    XAX_scen$Y$faer[ix_historical] = c(rep.col(x = XAX_scen$Y$ftot[ix_aer][1:165], length(ix_historical)/165))
    
    ## Forced natural Response:
    XAX_scen$Y$fnat[ix_aer] = rep(0, length(ix_aer))
    XAX_scen$Y$fnat[ix_nat] = XAX_scen$Y$ftot[ix_nat]
    XAX_scen$Y$fnat[ix_ghg] = rep(0, length(ix_ghg))
    XAX_scen$Y$fnat[ix_historical] = c(rep.col(x = XAX_scen$Y$ftot[ix_nat][1:165], length(ix_historical)/165))
    
    ## Forced GHG Response:
    XAX_scen$Y$fghg[ix_aer] = rep(0, length(ix_aer))
    XAX_scen$Y$fghg[ix_nat] = rep(0, length(ix_nat)) # XAX_scen$Y$ftot[ix_nat]
    XAX_scen$Y$fghg[ix_ghg] = XAX_scen$Y$ftot[ix_ghg]
    XAX_scen$Y$fghg[ix_historical] = c(rep.col(x = XAX_scen$Y$ftot[ix_ghg][1:165], length(ix_historical)/165))
    
    ## Forced ANTHROPOGENIC Response:
    XAX_scen$Y$fant[ix_aer] = XAX_scen$Y$faer[ix_aer]
    XAX_scen$Y$fant[ix_nat] = rep(0, length(ix_nat))
    XAX_scen$Y$fant[ix_ghg] = XAX_scen$Y$ftot[ix_ghg]
    XAX_scen$Y$fant[ix_historical] = c(rep.col(x = XAX_scen$Y$ftot[ix_historical][1:165] - XAX_scen$Y$ftot[ix_nat][1:165], length(ix_historical)/165))
    
    ## Regression based on radiative forcings:
    # ----------------------------------------
    cur.ix = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & !is.na(XAX_scen$Y$ftot) & !is.na(XAX_scen$Y$ERF_tot_r3))
    # cur.ix = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & !is.na(XAX_scen$Y$ftot) & !is.na(XAX_scen$Y$ERF_tot_r3) & 
    #                 XAX_scen$M$scen %in% c("hist-aer", "hist-GHG", "hist-nat"))
    # cur.ix = which(XAX_scen$M$mod.phys == DA.mod.phys.un[m] & !is.na(XAX_scen$Y$ftot) & !is.na(XAX_scen$Y$ERF_tot_r3) & 
    #                 !(XAX_scen$M$scen %in% c("hist-aer", "hist-GHG", "hist-nat")))
    print(c(paste(m, DA.mod.phys.un[m], sep =" "), unique(XAX_scen$M$scen[cur.ix])))
    
    y = XAX_scen$Y$ftot[cur.ix]; # plot(y)
    ERF_ant = XAX_scen$Y$ERF_ant_r3[cur.ix]
    ERF_nat = XAX_scen$Y$ERF_nat_r3[cur.ix]
    ERF_aer = XAX_scen$Y$ERF_aer_r3[cur.ix]
    ERF_ghg = XAX_scen$Y$ERF_ghg_r3[cur.ix]
    
    ## Forcing regression: 
    # lm1 = lm(y ~ ERF_ant + ERF_nat)
    # plot(lm1$fitted, y); abline(c(0, 1), col="red"); cor(lm1$fitted, y)^2
    # plot(y); lines(lm1$fitted, col = "darkblue", lwd = 2)
    lm2 = lm(y ~ ERF_aer + ERF_ghg + ERF_nat)
    print(coefficients(lm2))
    # plot(lm2$fitted, y); abline(c(0, 1), col="red"); legend("topleft", DA.mod.phys.un[m], cex=0.7); cor(lm2$fitted, y) ^ 2
    # plot(y); lines(lm2$fitted, col = "darkblue", lwd = 2); legend("topright", DA.mod.phys.un[m], cex=0.7)
    
    ## Predict for AER, GHG, and NAT:
    XAX_scen$Y$ftot_ERF[cur.ix] = predict(lm2, newdata=data.frame(ERF_aer = ERF_aer, ERF_ghg = ERF_ghg, ERF_nat = ERF_nat))
    XAX_scen$Y$fghg_ERF[cur.ix] = predict(lm2, newdata=data.frame(ERF_aer = rep(0, length(cur.ix)), ERF_ghg = ERF_ghg, ERF_nat = rep(0, length(cur.ix))))
    XAX_scen$Y$faer_ERF[cur.ix] = predict(lm2, newdata=data.frame(ERF_aer = ERF_aer, ERF_ghg = rep(0, length(cur.ix)), ERF_nat = rep(0, length(cur.ix))))
    XAX_scen$Y$fnat_ERF[cur.ix] = predict(lm2, newdata=data.frame(ERF_aer = rep(0, length(cur.ix)), ERF_ghg = rep(0, length(cur.ix)), ERF_nat = ERF_nat))
    XAX_scen$Y$fant_ERF[cur.ix] = XAX_scen$Y$ftot_ERF[cur.ix] - XAX_scen$Y$fnat_ERF[cur.ix]
  }
  
  ## Evaluation over all models:
  # all.ix = which(XAX_scen$M$mod.phys %in% DA.mod.phys.un & !is.na(XAX_scen$Y$ftot) & !is.na(XAX_scen$Y$ERF_tot_r3))
  # plot(XAX_scen$Y$ftot_ERF[all.ix], XAX_scen$Y$ftot[all.ix]); abline(c(0, 1), col="red", lwd = 2); legend("topleft", c("All Forc."), cex=0.7); cor(lm2$fitted, y) ^ 2
  # plot(XAX_scen$Y$faer_ERF[all.ix], XAX_scen$Y$faer[all.ix]); abline(c(0, 1), col="red", lwd = 2); cor(lm2$fitted, y) ^ 2
  # plot(XAX_scen$Y$fghg_ERF[all.ix], XAX_scen$Y$fghg[all.ix]); abline(c(0, 1), col="red", lwd = 2); cor(lm2$fitted, y) ^ 2
  # plot(XAX_scen$Y$fnat_ERF[all.ix], XAX_scen$Y$fnat[all.ix]); abline(c(0, 1), col="red", lwd = 2); cor(lm2$fitted, y) ^ 2
  
  return(XAX_scen)
}


get.XAX_ann_piC <- function(X, ncol = 72*36, M, start.year = 1, center = T, scale = F, cmip) {
  
  # OVERALL DIMENSION
  s = sum(sapply(X = X, FUN=function(x) dim(x)[2]))
  
  MAX <- data.frame(matrix(nrow=s,ncol=9))
  names(MAX) = c("vari", "res", "file.name", "cmip", "mod", "modcl", "scen", "ens.mem", "year")
  YAX <- data.frame(matrix(nrow=s, ncol=5))
  names(YAX) <- c("AGMT", "GMT.10y", "GMT.20y", "GMT.30y", "GMT.50y")
  
  XAX <- matrix(nrow=s,ncol=ncol)
  
  cc <- 0
  for (k in 1:length(X)){
    cat("\r ",k)
    Xsc <- t(X[[k]])
    
    if (center == T & scale == F) {
      Xsc = scale(Xsc, T, F)
    } else if (scale == T & scale == T) {
      Xsc = scale(Xsc, T, T)
    }
    # Ysc <- Y[[k]][seq(mon, dim(X[[k]])[3], 12),]
    
    # do slow averages:
    AGMT = c(Xsc %*% areaw)
    AGMT.10y = rollmean(x = AGMT, k = 10, fill = NA)
    AGMT.20y = rollmean(x = AGMT, k = 20, fill = NA)
    AGMT.30y = rollmean(x = AGMT, k = 30, fill = NA)
    AGMT.50y = rollmean(x = AGMT, k = 50, fill = NA)
    
    for (pc in 1:(dim(X[[k]])[2])){
      cc <- cc+1
      MAX[cc,] <- as.character(c(M$vari[k], M$res[k], M$file.name[k], cmip, M$mod[k], M$modcl[k], M$scen[k], M$ens.mem[k], start.year-1+as.numeric(pc)))
      YAX[cc,1] <- AGMT[pc]
      YAX[cc,2] <- AGMT.10y[pc]
      YAX[cc,3] <- AGMT.20y[pc]
      YAX[cc,4] <- AGMT.30y[pc]
      YAX[cc,5] <- AGMT.50y[pc]
      XAX[cc,]  <- Xsc[pc,] 
    }
    print(MAX[cc,])
  }
  
  # save to .RData file:
  return(list(X = XAX, Y=YAX, M=MAX))
}





library(foreach)
library(doParallel)
registerDoParallel(cores=46)



# get average over ensemble members:
get.ens.avg.trends <- function(cmip.trends) {
  
  mod.un=unique(cmip.trends$M$mod)
  CMIP6.files_tas_historical_trend_ensavg = list()
  CMIP6.files_tas_historical_trend_ensavg$X = matrix(data=NA, nrow = length(mod.un), ncol = 72*36)
  CMIP6.files_tas_historical_trend_ensavg$Y = cmip.trends$Y[1,]
  CMIP6.files_tas_historical_trend_ensavg$M = cmip.trends$M[1,]
  
  for (i in 1:length(mod.un)) {
    ix=which(mod.un[i] == cmip.trends$M$mod)
    print(cmip.trends$M$mod[ix])
    
    if (length(ix) == 1) {
      CMIP6.files_tas_historical_trend_ensavg$X[i,] = cmip.trends$X[ix,] 
    } else {
      CMIP6.files_tas_historical_trend_ensavg$X[i,] = apply(cmip.trends$X[ix,], MARGIN=2, FUN=mean)
    }
    CMIP6.files_tas_historical_trend_ensavg$M[i,] = cmip.trends$M[ix[1],] 
    CMIP6.files_tas_historical_trend_ensavg$Y[i,] = cmip.trends$Y[ix[1],] 
    CMIP6.files_tas_historical_trend_ensavg$Y[i,1] = c(CMIP6.files_tas_historical_trend_ensavg$X[i,] %*% areaw)
  }
  
  return(CMIP6.files_tas_historical_trend_ensavg)
}






extract.nyear.trend <- function(XAX, TCR, nyears = 34) {
  
  # run through each model:
  mod.un = unique(paste(XAX$M$mod, "_", XAX$M$scen, "_", XAX$M$ens.mem, sep=""))
  XAX.out = list()
  XAX.out$X = matrix(nrow=1,ncol=dim(XAX$X)[2])
  XAX.out$M = XAX$M[1,]
  names(XAX.out$M) = names(XAX$M)
  XAX.out$Y = data.frame(matrix(nrow=1, ncol=3))
  names(XAX.out$Y) <- c("HIST.TREND", "TCR", "ECS")
  s=1
  start.ix = numeric(length=1)
  
  for (m in 1:length(mod.un)) {
    print(m)
    cur.mod = strsplit(x = mod.un[m], split = "_")[[1]][1]
    cur.scen = strsplit(x = mod.un[m], split = "_")[[1]][2]
    cur.ens = strsplit(x = mod.un[m], split = "_")[[1]][3]
    
    ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens)
    start.year = seq(head(as.numeric(XAX$M$year[ix]), 1), tail(as.numeric(XAX$M$year[ix])-nyears, 1), by = nyears)
    cur.ix = which(XAX$M$mod == cur.mod & XAX$M$scen == cur.scen & XAX$M$ens.mem == cur.ens & XAX$M$year %in% start.year)
    start.ix[s:(s+length(cur.ix)-1)] = cur.ix
    s = length(start.ix) + 1
  }
  
  ## Now run and calculate trends based on start.ix:
  for (i in 1:(length(start.ix))) {
    print(i)
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    
    # calculate trends / get TCR+ECS:
    XAX.out$Y[i,] = c(lm(c(XAX$Y$AGMT[ix]) ~ c(1:nyears))$coefficients[2], NA, NA)
          tcr.ix = which(XAX$M$mod[ix[1]] == TCR$Model)
          if (length(TCR$TCR[tcr.ix]) == 1) XAX.out$Y[i,2] = TCR$TCR[tcr.ix]
          if (length(TCR$ECS[tcr.ix]) == 1) XAX.out$Y[i,3] = TCR$ECS[tcr.ix]
    XAX.out$M[i,] = XAX$M[ix[1],]
  }
  
  # Run trend calculation in parallel:
  test = foreach(i=1:(length(start.ix))) %dopar% {
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    test = apply(X = XAX$X[ix,], MARGIN=2, FUN=function(x) lm(x ~ c(1:nyears))$coefficients[2])
    return(test)
  }
  XAX.out$X = do.call(what = rbind, args = test)
  
  return(XAX.out)
}




# start.year = NULL -> start years taken with trend-sep difference from each other
extract.nyear.trend_ <- function(XAX, nyears = 34, start.year = NULL, trend.sep = 10, var.names = "GSAT") {
  
  # run through each model:
  mod.un = unique(paste(XAX$M$mod, "_", XAX$M$scen, "_", XAX$M$ens.mem, sep=""))
  XAX.out = list()
  XAX.out$X = matrix(nrow=1,ncol=dim(XAX$X)[2])
  XAX.out$M = XAX$M[1,]
  names(XAX.out$M) = names(XAX$M)
  XAX.out$Y = data.frame(matrix(nrow=1, ncol=length(var.names)))
  
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
  
  ## Now run and calculate trends based on start.ix:
  for (i in 1:(length(start.ix))) {
    print(i)
    ix = start.ix[i]:(start.ix[i]+nyears-1)
    
    # NA check:
    if(any(apply(XAX$Y[var.names][ix,], MARGIN=2, FUN = function( x) any(is.na(x))))) {
      XAX.out$Y[i,] = rep(NA, length(var.names))
      XAX.out$M[i,] = XAX$M[ix[1],]
      next;
    }
    
    # calculate trends / get TCR+ECS:
    XAX.out$Y[i,] = sapply(X = var.names, FUN=function(cur.name) lm(XAX$Y[[cur.name]][ix] ~ c(1:nyears))$coefficients[2]) # c(lm(c(XAX$Y$AGMT[ix]) ~ c(1:nyears))$coefficients[2], NA, NA)
    XAX.out$M[i,] = XAX$M[ix[1],]
  }
  
  # Run trend calculation in parallel:
  # test = foreach(i=1:(length(start.ix))) %dopar% {
  #  ix = start.ix[i]:(start.ix[i]+nyears-1)
  #  test = apply(X = XAX$X[ix,], MARGIN=2, FUN=function(x) lm(x ~ c(1:nyears))$coefficients[2])
  #  return(test)
  #}
  # XAX.out$X = do.call(what = rbind, args = test)
  names(XAX.out$Y) = var.names
  return(XAX.out)
}





extract.fraw_map2 <- function(XAX_scen, nmem = 3, f.var = "EOF245") {
  
  scen.un = unique(XAX_scen$M$scen)
  l = dim(XAX_scen$M)[1]
  
  # Extract "raw" forced response:
  XAX_scen$Y[[paste(f.var, "_f", sep="")]] = rep(NA, l)
  XAX_scen$Y[[paste(f.var, "_fl", sep="")]] = rep(NA, l)
  # XAX_scen$Y$ftot_l = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.10y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.20y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.30y = rep(NA, l)
  # XAX_scen$Y$IV.fraw.GMT.50y = rep(NA, l)
  
  XAX_scen$X_f = list() # matrix(NA, nrow = nrow(XAX_scen$X), ncol = ncol(XAX_scen$X))
  # XAX_scen$X_f_noincl = matrix(NA, nrow = nrow(XAX_scen$X), ncol = ncol(XAX_scen$X))
  XAX_scen$X_f_r31 = list() # matrix(NA, nrow = nrow(XAX_scen$X), ncol = ncol(XAX_scen$X))
  
  for (s in 1:length(scen.un)) {
    print(paste(s, scen.un[s], sep=" "))
    cur.scen = scen.un[s]
    
    # Determine number of members per model:
    no.mem.mod = sapply(X = unique(XAX_scen$M$mod), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod & XAX_scen$M$scen == cur.scen)])))
    # no.mem = sapply(X = unique(XAX_scen$M$mod.phys), FUN=function(x) length(unique(XAX_scen$M$ens.mem[which(x == XAX_scen$M$mod.phys & XAX_scen$M$scen == cur.scen)])))
    # length(which(no.mem.mod >= 3))
    # length(which(no.mem >= 3))
    
    cur.no.mem = no.mem.mod[which(no.mem.mod >= nmem)]
    cur.mod.phys.un = names(which(no.mem.mod >= nmem))
    
    for (m in 1:length(cur.mod.phys.un)) {
      print(paste(m, cur.mod.phys.un[m], sep=" "))
      
      # all scenarios:
      ix_scen = which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen)
      ens.mat = sapply(X = unique(XAX_scen$M$ens.mem[ix_scen]), FUN=function(cur.ens) XAX_scen$Y[[f.var]][which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen & XAX_scen$M$ens.mem == cur.ens)])
      year.un = unique(as.numeric(XAX_scen$M$year[ix_scen]))
      
      XAX_scen$X_f[[m]] = t(sapply(X = 1:length(year.un), FUN=function(y) {
        ix = which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == cur.scen & XAX_scen$M$year == year.un[y])
        return(colMeans(XAX_scen$X[ix,]))
      }))  
        
      if (any(is.na(ens.mat))) next 
      if (cur.no.mem[m] > 2) {
        ftot_IV = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat[,-k], na.rm=T))
      } else {
        ftot_IV = matrix(NA, nrow = dim(ens.mat)[1], ncol = dim(ens.mat)[2])
      }
      ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
      # ftot = sapply(X = 1:dim(ens.mat)[2], FUN=function(k) rowMeans(ens.mat, na.rm=T))
      # print(dim(ftot))
      # Fit LOWESS: fr.loess_0.75 = loess(fr.raw ~ c(1:231), span = 0.75, degree = 2)$fittedn=dim(ens.mat)[1]
      # 231 * 0.25 / n
      n = dim(ens.mat)[1]
      ftot_l = rep.col(loess(rowMeans(ftot) ~ c(1:n), span = 0.99, degree = 2)$fitted, dim(ens.mat)[2])
      if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
      
      # get two smoothed versions of the forced responses:
      XAX_scen$X_f_r31[[m]] = apply(X = XAX_scen$X_f[[m]], MARGIN = 2, FUN=function(x) rollmean(x = x, k = 31, fill = NA))
      
      if (cur.scen %in% c("ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp585")) {
        ix_hist = c(which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical"))
        ens.mat.hist = sapply(X = unique(XAX_scen$M$ens.mem[ix_hist]), FUN=function(cur.ens) XAX_scen$Y[[f.var]][which(XAX_scen$M$mod == cur.mod.phys.un[m] & XAX_scen$M$scen == "historical" & XAX_scen$M$ens.mem == cur.ens)])
        ftot_hist = c(sapply(X = 1:dim(ens.mat.hist)[2], FUN=function(k) rowMeans(ens.mat.hist, na.rm=T))[,1])
        ftot_long = c(ftot_hist, rowMeans(ftot))
        ftot_l = loess(ftot_long ~ c(1:(165+n)), span = 0.99, degree = 2)$fitted
        if (is.na(ftot_long[1])) ftot_l = c(NA, ftot_l)
        # if (all(is.na(ftot[1,]))) ftot_l = rbind(rep(NA, dim(ens.mat)[2]), ftot_l)
        ftot_l = rep.col(ftot_l[-c(1:165)], dim(ens.mat)[2])
      }
      
      ## Fill in forced responses and smoothed versions:
      XAX_scen$Y[[paste(f.var, "_f", sep="")]][ix_scen] = c(ftot)
      XAX_scen$Y[[paste(f.var, "_fl", sep="")]][ix_scen] = c(ftot_l)
      
      # IV.fraw = ens.mat - ftot_IV
      # XAX_scen$Y$IV.fraw.GMT.10y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 10, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.20y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 20, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.30y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 30, fill = NA))
      # XAX_scen$Y$IV.fraw.GMT.50y[ix_scen] = c(apply(X = IV.fraw, MARGIN=2, FUN=rollmean, k = 50, fill = NA))
    }
  }
  
  # which(is.na(XAX_scen$Y$EOF245_f))
  # plot(XAX_scen$Y$EOF245_fl[50000:100000], type='l')
  # XAX_scen$M$scen[20001:21000]
  
  return(XAX_scen)
}

