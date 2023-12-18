
# ------------------------------------------------------------------------------------
# Train tas-reconstruction based on land-only grid cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 17.06.2021
# run on xenon-server:
# screen -S train_psl
# module load R/4.0.3-openblas
# R (-> not R-3.6.1)

library(ncdf4)
library(raster)
library(fields)
library(foreach)
library(doParallel)
library(glmnet)
library(hydroGOF)
library(openblasctl) # install.packages("openblasctl", repos="https://hpcran.org")
openblas_set_num_threads(num_threads = 3)


# Change things: 
# 1. run with interpolated grid cells but only where error estimates are available AND exclude very high-latitude grid cells. OK
# 2. run with original error estimates (Had_SLP2_obs_error_) -> OK | as well as the filled error estimates (Had_SLP2_obs_error_f_) -> TODO evtl.
# 3. pick \lambda such that the correlation with tas_land/tos reconstruction is optimized. -> later.
# 



# --------------------------------------------------------------------------
# 0.a Load functions needed:
# --------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_ridge_fun.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_gta_fun.R")


# --------------------------------------------------------------------------
# 1. Train to predict forced response based on land-only record.
# --------------------------------------------------------------------------

setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/")

## check if spatial format is correct:
# image.plot(matrix(HadSLP2_uninterpol_[100,], 72, 37))
# image.plot(matrix(cur.CMIP6.psl_monX_ct$X[13020,], 72, 37)[,37:1])
# test = subset(brick("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/psl_mon_ACCESS-CM2_historical_r1i1p1f1_g025.nc"), seq(1, 360*2, 12)) + 0
# plot(mean(test))
# test0 = subset(brick("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/psl/mon/g72x37/psl_mon_ACCESS-CM2_historical_r1i1p1f1_g025.nc"), 1) - mean(test)
# image.plot(matrix(cur.CMIP6.psl_monX_ct$X[1,], 72, 37))
# image.plot(matrix(values(test0), 72, 37)[,37:1])
# plot(c(matrix(cur.CMIP6.psl_monX_ct$X[1,], 72, 37)), c(matrix(values(test0), 72, 37)[,37:1]))


for (mon in 1:12) {
  print(mon)
  
  ## 0.1 Load climate model monthly data & select training model indices:
  {
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.psl_mon", mon, "_ct.RData", sep=""))
    
    # select training model indices:
    CMIP6.psl_monX_ct$M$modcl[which(CMIP6.psl_monX_ct$M$modcl == "UKE")] <- "Had"
      modcl.un = unique(CMIP6.psl_monX_ct$M$modcl)
      train.ix = list(); train.ix1 = list();
    for (m in 1:length(modcl.un)) { 
      print(m)
      cur.modcl = modcl.un[m]
        ix = which(cur.modcl == CMIP6.psl_monX_ct$M$modcl & CMIP6.psl_monX_ct$M$year %in% c(1850:2020)  & 
                     !is.na(CMIP6.psl_monX_ct$Y$GMLSAT_MI) & !is.na(CMIP6.psl_monX_ct$Y$GMSST))
          ens.mem = unique(CMIP6.psl_monX_ct$M$ens.mem[ix])[1:5]
      train.ix[[m]] = which(cur.modcl == CMIP6.psl_monX_ct$M$modcl & CMIP6.psl_monX_ct$M$year %in% c(1850:2020) & CMIP6.psl_monX_ct$M$ens.mem %in% ens.mem & 
                              !is.na(CMIP6.psl_monX_ct$Y$GMLSAT_MI) & !is.na(CMIP6.psl_monX_ct$Y$GMSST))
      if (length(train.ix[[m]]) < 495) train.ix[[m]] = integer(0)
      # train.ix1[[m]] = which(cur.modcl == CMIP6.psl_monX_ct$M$modcl & !is.na(CMIP6.psl_monX_ct$Y$GSAT_f) & CMIP6.psl_monX_ct$M$ens.mem %in% paste("r", 1:5,"i1p1f1", sep="") & CMIP6.psl_monX_ct$M$year %in% 1850:2020 )
    }
      train.mod = cbind(1:30, modcl.un, sapply(X = train.ix, FUN = function(x) length(x)))
      train.ix = unlist(train.ix) # length(train.ix)
    # old version:
      # train.ix1 = which(!is.na(CMIP6.psl_monX_ct$Y$GSAT_f) & CMIP6.psl_monX_ct$M$ens.mem %in% paste("r", 1:5,"i1p1f1", sep="") & CMIP6.psl_monX_ct$M$year %in% 1850:2020 )
      # train.mod1 = cbind(1:30, modcl.un, sapply(X = train.ix1, FUN = function(x) length(x)))
      # train.ix1 = unlist(train.ix1)
      # cbind(train.mod, train.mod1)
    
    cur.CMIP6.psl_monX_ct = list()
    cur.CMIP6.psl_monX_ct$X = CMIP6.psl_monX_ct$X[train.ix,]
    cur.CMIP6.psl_monX_ct$Y = CMIP6.psl_monX_ct$Y[train.ix,]
    cur.CMIP6.psl_monX_ct$M = CMIP6.psl_monX_ct$M[train.ix,]
    rm(CMIP6.psl_monX_ct)
  }
  
  ## 0.2 load CRU data:
  {
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_HadSLP2/HadSLP2_mon", mon, ".RData", sep=""))
    # correct units to [Pa]:
    HadSLP2_interpol_ = HadSLP2_interpol_ * 100; HadSLP2_uninterpol_ = HadSLP2_uninterpol_ * 100; 
    HadSLP2_obs_error_ = HadSLP2_obs_error_ * 100; HadSLP2_obs_error_f_ = HadSLP2_obs_error_f_ * 100; HadSLP2_obs_error_f2x_ = HadSLP2_obs_error_f2x_ * 100; 
  }
  
  
  ## Run through all time steps from 1850 through 2020
  registerDoParallel(cores = 6)
  cv.nr.cores = 4
  
  foreach(date.ix=1:155) %dopar% {
    
    print(HadSLP2_calendar_[date.ix])
      
    # 1.1 Select grid cell indices:
    {
      raster.template = raster(xmn = 0, xmx=360, ymn = -90, ymx=90, nrows = 37, ncols = 72)
      lats = c(matrix(coordinates(raster.template)[,2], 72, 37)[,37:1])

      # grid.ix with original error estimates:
      grid.ix = which(!is.na(HadSLP2_interpol_[date.ix,]) & !is.na(HadSLP2_obs_error_[date.ix,]) & lats < 80 & lats > -80)
      unc = HadSLP2_obs_error_[date.ix,grid.ix]
      # image.plot(matrix(HadSLP2_obs_error_f_[date.ix,], 72, 37))
      # image.plot(matrix(HadSLP2_obs_error_[date.ix,], 72, 37))
      # image.plot(matrix(cur.CMIP6.tas_monX_ct$X[10300,], 72, 36))
      
      # grid.ix with filled error estimates:
      # grid.ix2 = which(!is.na(HadSLP2_uninterpol_[date.ix,]) & !is.na(HadSLP2_obs_error_f_[date.ix,]) & lats < 80 & lats > -80)
      # unc2 = HadSLP2_obs_error_f2x_[date.ix,grid.ix]
    }
      

    # 1.2 Combine CMIP6 training data into one data.frame and save removed outliers:
      {
        CMIP6.df <- list()
        CMIP6.df$Y <- cur.CMIP6.psl_monX_ct$Y
        CMIP6.df$M <- cur.CMIP6.psl_monX_ct$M
        # CMIP6.df$cur.land.weights = cur.land.weights
        CMIP6.df$grid.ix = grid.ix
        CMIP6.df$train.ix = train.ix
      }
      
      
    # 1.3 Generate perturbed training dataset: str(X)
      {
        X = list()
        X$X_orig = cur.CMIP6.psl_monX_ct$X[,grid.ix]
        X$X_pert = gen.pert.traintest.data_noENS(X = cur.CMIP6.psl_monX_ct$X[,grid.ix], M = cur.CMIP6.psl_monX_ct$M, 
                                           bias.ens = NULL, unc = unc, fact = 1)
        X$X_pert2 = gen.pert.traintest.data_noENS(X = cur.CMIP6.psl_monX_ct$X[,grid.ix], M = cur.CMIP6.psl_monX_ct$M, 
                                           bias.ens = NULL, unc = unc, fact = 2)
        #ens.mem.split = gen.pert.traintest.data(X = cur.CMIP6.psl_monX_ct$X[,grid.ix], M = cur.CMIP6.psl_monX_ct$M, 
        #                                        bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1, ret.ens.mem.split = T)
        
        train.ix1 =  which(!is.na(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl))
        CMIP6.df$train.ix1 = train.ix1
        X_IV = list()
        X_IV$X_orig = X$X_orig[train.ix1,]
        X_IV$X_pert = X$X_pert[train.ix1,]
        X_IV$X_pert2 = X$X_pert2[train.ix1,]
      }
      
    # 1.4 Generate observational input data with uncertainties:
    {
        cur.HadSLP2_ENS = list()
        cur.HadSLP2_ENS$X_pert0 = rep.row(HadSLP2_interpol_[date.ix,], n = 200)[,grid.ix]
        cur.HadSLP2_ENS$X_pert1 = rep.row(HadSLP2_interpol_[date.ix,], n = 200)[,grid.ix] + 
                  sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
        cur.HadSLP2_ENS$X_pert2 = rep.row(HadSLP2_interpol_[date.ix,], n = 200)[,grid.ix] + 
                  2 * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
    }
     
    # 1.5 Define model training parameters:
    {
      # Define training parameters:
      modcl.un = unique(cur.CMIP6.psl_monX_ct$M$modcl)
      crossclass = sapply(X = cur.CMIP6.psl_monX_ct$M$modcl, FUN=function(x) which(x == modcl.un))
      lambda.seq = rev(10 ^ seq(-5,5,length.out=100))
      penalty.factor1 = rep(1, length(grid.ix))
      train.weights = unlist(sapply(X = unique(crossclass), FUN=function(i) rep(1 / length(which(crossclass == unique(crossclass)[i])) / length(unique(crossclass)), length(which(crossclass == unique(crossclass)[i])))))
    }
      
    # 4. TRAIN MODELS WITH DIFFERENT TARGETS:  
    # --------------------------------------------------
      
    # 4.1 Predict GSAT:  
      {
        
      # Start the clock!
      # ptm <- proc.time()
    mod0.GSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.psl_monX_ct$Y$GSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T)
    mod_p0 = validate.pert.dataset_4psl(X_pert = X, Y = cur.CMIP6.psl_monX_ct$Y$GSAT, X_obs = cur.HadSLP2_ENS,
                               mod.glmnet = mod0.GSAT, foldid = crossclass)
    
    mod1.GSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.psl_monX_ct$Y$GSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T)
    mod_p1 = validate.pert.dataset_4psl(X_pert = X, Y = cur.CMIP6.psl_monX_ct$Y$GSAT, X_obs = cur.HadSLP2_ENS,
                               mod.glmnet = mod1.GSAT, foldid = crossclass)
    
    
    ## try prediction of variability:
    mod0.GSAT_IV = cv.glmnet2(x = X_IV$X_orig, y = c(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl)[train.ix1], alpha = 0,
                           lambda = lambda.seq, foldid = crossclass[train.ix1], 
                           penalty.factor = penalty.factor1[train.ix1], weights = train.weights[train.ix1],
                           nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T)
    mod_p0_IV = validate.pert.dataset_4psl(X_pert = X_IV, Y = c(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl)[train.ix1], X_obs = cur.HadSLP2_ENS,
                                        mod.glmnet = mod0.GSAT_IV, foldid = crossclass[train.ix1])
    # plot(mod_p0_IV$Yhat$pt0[,2], c(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl)[train.ix1])
    
    mod1.GSAT_IV = cv.glmnet2(x = X_IV$X_pert, y = c(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl)[train.ix1], alpha = 0,
                           lambda = lambda.seq, foldid = crossclass[train.ix1], 
                           penalty.factor = penalty.factor1[train.ix1], weights = train.weights[train.ix1],
                           nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T)
    mod_p1_IV = validate.pert.dataset_4psl(X_pert = X_IV, Y = c(cur.CMIP6.psl_monX_ct$Y$GSAT - cur.CMIP6.psl_monX_ct$Y$GSAT_fl)[train.ix1], X_obs = cur.HadSLP2_ENS,
                                        mod.glmnet = mod1.GSAT_IV, foldid = crossclass[train.ix1])
    
    # Stop the clock
    # proc.time() - ptm
    ## runs roughly 4 minutes. 
    ## 4 * 171 * 12 ~ 136.8 h 
    
    # save models for given time step:
    save(list = c("mod_p0", "mod_p1", "mod_p0_IV", "mod_p1_IV", "CMIP6.df"), 
         file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/HadSLP2_predGSAT_v2/", 
                      format(HadSLP2_calendar_[date.ix], "%Y"), "-", format(HadSLP2_calendar_[date.ix], "%m"), ".RData", sep=""))
      }
  }
}




# -THE END-

## Check different estimates:
# plot(mod_p0$Yhat$pt0[,1], cur.CMIP6.psl_monX_ct$Y$GSAT)
# cor(mod_p0$Yhat$pt0[,1], cur.CMIP6.psl_monX_ct$Y$GSAT)

# plot(mod_p1$Yhat$pt1[,1], cur.CMIP6.psl_monX_ct$Y$GSAT)
# cor(mod_p1$Yhat$pt1[,1], cur.CMIP6.psl_monX_ct$Y$GSAT)


## compare fingerprint estimates (very different!!):
# test = raster(xmn = 0, xmx=360, ymn = -90, ymx=90, nrows = 37, ncols = 72)
# beta = rep(NA, 72*37)
# beta[CMIP6.df$grid.ix] <- mod_p0$beta[,mod_p0$lambda.1_05]
# image.plot(matrix(beta, 72, 37))

# test = raster(xmn = 0, xmx=360, ymn = -90, ymx=90, nrows = 37, ncols = 72)
# beta = rep(NA, 72*37)
# beta[CMIP6.df$grid.ix] <- mod_p1$beta[,mod_p1$lambda.1_05]
# image.plot(matrix(beta, 72, 37))


