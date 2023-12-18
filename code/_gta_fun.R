

# TEST GTA1:

library(ncdf4);
library(ncdf4.helpers);


# prepare covariance/correlation matrix
prepare_cov = function( tmap, dist ) {
  xs = (((1:dim(tmap)[2])-0.5)*180/dim(tmap)[2]- 90.0)*pi/180;
  ys = (((1:dim(tmap)[1])-0.5)*360/dim(tmap)[1]-180.0)*pi/180;
  # switch each and times if necessary
  las = rep(xs,each=length(ys));
  lns = rep(ys,times=length(xs));
  dists = matrix(0,nrow=length(las),ncol=length(las));
  for ( i in 1:length(las) ) {
    dists[i,] = 6371.0*acos(pmin(pmax(sin(las[i])*sin(las) + cos(las[i])*cos(las)*cos(lns[i]-lns),-1.0),1.0));
  }
  cov = exp(-dists/dist);
  return(cov);
}

prepare_cov_Pacific_centered = function( tmap, dist ) {
  xs = (((1:dim(tmap)[2])-0.5)*180/dim(tmap)[2]- 90.0)*pi/180;
  ys = (((1:dim(tmap)[1])-0.5)*360/dim(tmap)[1]-0)*pi/180;
  # switch each and times if necessary
  las = rep(xs,each=length(ys));
  lns = rep(ys,times=length(xs));
  dists = matrix(0,nrow=length(las),ncol=length(las));
  for ( i in 1:length(las) ) {
    dists[i,] = 6371.0*acos(pmin(pmax(sin(las[i])*sin(las) + cos(las[i])*cos(las)*cos(lns[i]-lns),-1.0),1.0));
  }
  cov = exp(-dists/dist);
  return(cov);
}

# calculate GTA1 estimator for a given map using the correlation matrix
gta1 = function( t, cov ) {
  data = as.vector(t);
  # set up matrices
  unobsflag = is.na(data);
  obsflag = !unobsflag;
  y = data[obsflag];
  w = cov[obsflag,][,obsflag];
  # solve for gls mean
  wi = solve(w);
  swi = rowSums(wi);
  return(sum(swi*y)/sum(swi));
}




# MAIN PROGRAM
# read data: available from https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/gridded_fields/HadCRUT.4.6.0.0.median_netcdf.zip
# args = commandArgs(trailingOnly=TRUE);
# setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadCRUT4/5d00_monthly/")
# nc = nc_open("HadCRUT.4.6.0.0.median.nc");
# tobs = ncvar_get(nc,"temperature_anomaly");
# ts <- ncvar_get(nc,"time");

# do rough coversion of time to year/month
# ts = ts/365.25+1850.0;
# ts = floor(ts)+floor(12*ts%%1)/12.0+1/24.0;

# flag missings
# tobs[t<-90] = NA;
# tobs[t>490] = NA;

# prepare correlation matrix
# cov = prepare_cov(tobs[,,0],800.0);
# cov2 = prepare_cov_Pacific_centered(tobs[,,0],800.0); 

# image(cov)
# image(cov2)
# min(cov - cov2)  # cov and cov2 are virtually identical.


# calculate temperatures
# T_GTA1 = rep(NA, length(ts))

# for (m in 1:(dim(tobs)[3])) {
#  t = gta1(tobs[,,m],cov);
#  T_GTA1[m] = t
#  cat(ts[m]," ",round(t,4),"\n");
# }



validate.pert.dataset_gta <- function(X_pert, Y, X_obs, beta.gta, ens.mem.split, ens.ix = NULL) {
  
  # 1. Prepare regression/cross-validation:
  n.ens = length(unique(ens.mem.split))
  
  beta = beta.gta

  # generate list of Yhat:
  Yhat = list()
  lm.cal = list()
  Yhat_raw = list()
  MSE.by.mod = list()
  cor.by.mod = list()
  me.by.mod = list()
  MSE.by.ensreal = list()
  cor.by.ensreal = list()
  me.by.ensreal = list()
  
  for (i in 1:length(X_pert)) {
    Yhat_raw[[i]] = rep(NA, length(Y))
    Yhat[[i]] = rep(NA, length(Y))
    MSE.by.mod[[i]] = rep(NA, length(Y))
    cor.by.mod[[i]] = rep(NA, length(Y))
    me.by.mod[[i]] = rep(NA, length(Y))
    
    MSE.by.ensreal[[i]] = rep(NA, length(unique(ens.mem.split)))
    cor.by.ensreal[[i]] = rep(NA, length(unique(ens.mem.split)))
    me.by.ensreal[[i]] = rep(NA, length(unique(ens.mem.split)))
  }
  
  # 2. Predict Yhat and compute model-specific MSE and other metrics:
  # for (sim in 1:nsim) {
    # print(paste("\r ***", sim))
  #    ix = which(foldid == sim)
    
    for (i in 1:length(X_pert)) {
      # (i) predict out of sample observations:
      Yhat_raw[[i]] = unname(c(X_pert[[i]] %*% beta))
      Yhat[[i]] = unname(lm(Y ~ Yhat_raw[[i]])$fitted)
      lm.cal[[i]] = lm(Y ~ Yhat_raw[[i]])
      
      # (ii) Calculate error metrics per model:
      MSE.by.mod[[i]] = mse(sim = Yhat[[i]], obs = Y)
      cor.by.mod[[i]] = cor(Yhat[[i]], Y) # sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) cor(Yhat[[i]][ix,lambda.ix], Y[ix]))
      me.by.mod[[i]] = me(sim = Yhat[[i]], obs = Y)
    }
  # }
  
  # 3. Calculate error metrics for each bias realization ensemble member:
  for (en in 1:n.ens) {
    # print(paste("\r ***", en))
    ix = which(ens.mem.split == en)
    
    for (i in 1:length(X_pert)) {
      MSE.by.ensreal[[i]][en] = mse(sim = Yhat[[i]][ix], obs = Y[ix])
      cor.by.ensreal[[i]][en] = cor(Yhat[[i]][ix], Y[ix])
      me.by.ensreal[[i]][en] = me(sim = Yhat[[i]][ix], obs = Y[ix])
    }
  }
  
  # Generate averages:
  beta_ = beta
  # a0_ = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mean(sapply(X = a0, FUN=function(x) x[lambda.ix]), na.rm=T))
  MSE = MSE.by.mod
  
  # Evaluate \lambda.min and \lambda.min+0.05
  # pt.lambda.min = lapply(X = MSE, FUN=function(x) which.min(x))
  # pt.lambda.1.05 = lapply(X = MSE, FUN=function(x) which(x < (min(x) * 1.05))[1])
  
  # Generate prediction on observations:
  obs_pred = c(( X_obs$X_pert %*% beta_ ) * lm.cal[[2]]$coefficients[2] + lm.cal[[2]]$coefficients[1])
  obs_pred_pert2 = c(( X_obs$X_pert2 %*% beta_ ) * lm.cal[[3]]$coefficients[2] + lm.cal[[3]]$coefficients[1])
  
  # obs_pred_min = obs_pred[,lambda.min]
  # obs_pred_1_05 = obs_pred[,lambda.1_05]
  
  
  names(MSE.by.mod) <- names(cor.by.mod) <- names(me.by.mod) <- names(MSE.by.ensreal) <- names(cor.by.ensreal) <- names(me.by.ensreal) <- names(MSE) <- c("pt0", "pt1", "pt2")
  names(Yhat_raw) <- names(Yhat) <- c("pt0", "pt1", "pt2")
  
  ret.list = list(MSE.by.mod, cor.by.mod, me.by.mod, MSE.by.ensreal, cor.by.ensreal, me.by.ensreal, 
                  beta_, MSE, 
                  obs_pred, obs_pred_pert2, ens.ix,
                  Yhat, Yhat_raw)
  names(ret.list) = c("MSE.by.mod", "cor.by.mod", "me.by.mod", "MSE.by.ensreal", "cor.by.ensreal", "me.by.ensreal", 
                      "beta", "MSE", 
                      "obs_pred", "obs_pred_pert2", "ens.ix",
                      "Yhat", "Yhat_raw")
  return(ret.list)
}






