
# ------------------------------------------------------------------------------------
# Train tas-reconstruction based on land-only grid cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 17.06.2021
# run on co2-server:
# screen -S train_tas_land
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



# --------------------------------------------------------------------------
# 0.a Load functions needed:
# --------------------------------------------------------------------------
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_ridge_fun_v2.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_gta_fun.R")


raster.template = raster(res = 5, xmn = 0, xmx=360, ymn = -90, ymx=90)
areaw = c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]) / sum(c(matrix(values(raster::area(raster.template)), 72, 36)[,36:1]))
areaw_cos = sqrt(cos(c(matrix(raster::coordinates(raster.template)[,2], 72, 36)[,36:1])*pi/180))

# --------------------------------------------------------------------------
# 1. Train to predict forced response based on land-only record.
# --------------------------------------------------------------------------

setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/")

# New paleo reconstructions: tas_land_predTMSST_v4; tas_land_predTMSST25_v4; tas_land_predWAtlantic_v4; tas_land_predWPacific_v4; tas_land_predEPacific_v4; tas_land_predIndianOcean_v4; tas_land_predTropics_v4; 

for (mon in 1:12) {
  print(mon)
  
  ## 0.1 Load climate model monthly data & select training model indices:
  {
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_mon", mon, "_ct.RData", sep=""))
    
    # select training model indices:
    CMIP6.tas_monX_ct$M$modcl[which(CMIP6.tas_monX_ct$M$modcl == "UKE")] <- "Had"
      mod.un = unique(CMIP6.tas_monX_ct$M$mod)
      train.ix = list(); train.ix1 = list();
      
    for (m in 1:length(mod.un)) { 
      print(m)
      cur.mod = mod.un[m]
      
      if (substr(mod.un, 1, 3)[m] %in% c("CES", "EC-", "Had", "MIR", "MPI", "Nor")) {   # models with more than 3 versions... -> choose only two ensemble members per version
        ix = which(cur.mod == CMIP6.tas_monX_ct$M$mod & CMIP6.tas_monX_ct$M$year %in% c(1850:2020) & 
                     !is.na(CMIP6.tas_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tas_monX_ct$Y$GMSST) & !is.na(CMIP6.tas_monX_ct$Y$GMST_FM))
        ens.mem = unique(CMIP6.tas_monX_ct$M$ens.mem[ix])[1:2]
        train.ix[[m]] = which(cur.mod == CMIP6.tas_monX_ct$M$mod & CMIP6.tas_monX_ct$M$year %in% c(1850:2020) & CMIP6.tas_monX_ct$M$ens.mem %in% ens.mem & 
                                !is.na(CMIP6.tas_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tas_monX_ct$Y$GMSST) & !is.na(CMIP6.tas_monX_ct$Y$GMST_FM))
        if (length(train.ix[[m]]) < 330) train.ix[[m]] = integer(0)
      } else {
        ix = which(cur.mod == CMIP6.tas_monX_ct$M$mod & CMIP6.tas_monX_ct$M$year %in% c(1850:2020) & 
                     !is.na(CMIP6.tas_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tas_monX_ct$Y$GMSST) & !is.na(CMIP6.tas_monX_ct$Y$GMST_FM))
          ens.mem = unique(CMIP6.tas_monX_ct$M$ens.mem[ix])[1:3]
        train.ix[[m]] = which(cur.mod == CMIP6.tas_monX_ct$M$mod & CMIP6.tas_monX_ct$M$year %in% c(1850:2020) & CMIP6.tas_monX_ct$M$ens.mem %in% ens.mem & 
                              !is.na(CMIP6.tas_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tas_monX_ct$Y$GMSST) & !is.na(CMIP6.tas_monX_ct$Y$GMST_FM))
        if (length(train.ix[[m]]) < 495) train.ix[[m]] = integer(0)
      }
      # train.ix1[[m]] = which(cur.modcl == CMIP6.tas_monX_ct$M$modcl & !is.na(CMIP6.tas_monX_ct$Y$GSAT_f) & CMIP6.tas_monX_ct$M$ens.mem %in% paste("r", 1:5,"i1p1f1", sep="") & CMIP6.tas_monX_ct$M$year %in% 1850:2020 )
    }
      
      train.mod = cbind(1:64, mod.un, sapply(X = train.ix, FUN = function(x) length(x)))
      train.ix = unlist(train.ix) # length(train.ix)
    # old version:
      # train.ix1 = which(!is.na(CMIP6.tas_monX_ct$Y$GSAT_f) & CMIP6.tas_monX_ct$M$ens.mem %in% paste("r", 1:5,"i1p1f1", sep="") & CMIP6.tas_monX_ct$M$year %in% 1850:2020 )
      # train.mod1 = cbind(1:30, modcl.un, sapply(X = train.ix1, FUN = function(x) length(x)))
      # train.ix1 = unlist(train.ix1)
      # cbind(train.mod, train.mod1)
      
    cur.CMIP6.tas_monX_ct = list()
    cur.CMIP6.tas_monX_ct$X = CMIP6.tas_monX_ct$X[train.ix,]
    cur.CMIP6.tas_monX_ct$Y = CMIP6.tas_monX_ct$Y[train.ix,]
    cur.CMIP6.tas_monX_ct$M = CMIP6.tas_monX_ct$M[train.ix,]
    rm(CMIP6.tas_monX_ct)
  }
  
  ## 0.2 load CRU data:
  load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/CRUTEM5_mon", mon, ".RData", sep=""))
  
  
  ## Run through all time steps from 1850 through 2020
  registerDoParallel(cores = 6)
  cv.nr.cores = 4
  
  foreach(date.ix=1:171) %dopar% {
    
    print(CRUTEM5_calendar_[date.ix])
      
    # 1.1 Select grid cell indices:
    {
      cur.land.weights = land.weights[date.ix,]
      cur.land.weights[which(cur.land.weights == 0)] = NA
      grid.ix = which(cur.land.weights > 0)
      ## define bias and uncertainties according to CRU dataset:
      bias.ens = lapply(X = CRUTEM5_ENS_anom_, FUN=function(x) x[date.ix,grid.ix])
      unc = CRUTEM5_sampling_unc_[date.ix,grid.ix]
      # image.plot(matrix(cur.land.weights, 72, 36))
      # image.plot(matrix(cur.CMIP6.tas_monX_ct$X[10300,], 72, 36))
    }
      

    # 1.2 Combine CMIP6 training data into one data.frame and save removed outliers:
    {
        CMIP6.df <- list()
        #CMIP6.df$Y <- cur.CMIP6.tas_monX_ct$Y
        #CMIP6.df$M <- cur.CMIP6.tas_monX_ct$M
        CMIP6.df$cur.land.weights = cur.land.weights
        CMIP6.df$grid.ix = grid.ix
        CMIP6.df$train.ix = train.ix
        CMIP6.df$train.mod = train.mod
      }
      
      
    # 1.3 Generate perturbed training dataset. str(X)
    {
        X = list()
        X$X_orig = cur.CMIP6.tas_monX_ct$X[,grid.ix]
        X$X_pert = gen.pert.traintest.data(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                           bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1)
        X$X_pert2 = gen.pert.traintest.data(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                           bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 2)
        ens.mem.split = gen.pert.traintest.data(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                                bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1, ret.ens.mem.split = T)
      }
      
    # 1.4 Generate observational input data with uncertainties:
    {
        cur.CRUTEM5_ENS = list()
        cur.CRUTEM5_ENS$X_pert = rep.row(CRUTEM5_[date.ix,], n = 200)[,grid.ix] + 
                  t(sapply(X = bias.ens, FUN=function(x) x)) +
                  sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
        cur.CRUTEM5_ENS$X_pert2 = rep.row(CRUTEM5_[date.ix,], n = 200)[,grid.ix] + 
                  2 * t(sapply(X = bias.ens, FUN=function(x) x)) +
                  2 * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
    }
     
    # 1.5 Define model training parameters:
    {
      # Define training parameters:
      modcl.un = unique(cur.CMIP6.tas_monX_ct$M$modcl)
      crossclass = sapply(X = cur.CMIP6.tas_monX_ct$M$modcl, FUN=function(x) which(x == modcl.un))
      lambda.seq = rev(10 ^ seq(-5,5,length.out=100))
      penalty.factor1 = rep(1, length(grid.ix))
      penalty.factor2 = 1 / (cur.land.weights[grid.ix])
      train.weights = unlist(sapply(X = unique(crossclass), FUN=function(i) rep(1 / length(which(crossclass == unique(crossclass)[i])) / length(unique(crossclass)), length(which(crossclass == unique(crossclass)[i])))))
    }
      
    # 4. TRAIN MODELS WITH DIFFERENT TARGETS:  
    # --------------------------------------------------
      
    # 2.1 Get SVD computations for each setup:
    {
      # mod_p0:
      svd_mod_p0 = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", ret.svd = T)
      # mod_p1:
      svd_mod_p1 = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", ret.svd = T)
      
      # mod_gta:
      cov2 = prepare_cov_Pacific_centered(matrix(cur.land.weights, 72, 36),800.0); 
      w = cov2[grid.ix,][,grid.ix]; wi = solve(w); beta.gta = rowSums(wi) / sum(rowSums(wi));
    }
    
    
    
    # 4.1 Predict GSAT:  
    {
      
      mod0.GSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GSAT, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GSAT, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GSAT, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGSAT_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.2 Predict GMST:  
    {
      
      mod0.GMST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GMST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GMST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.3 Predict GMLSAT_MI:  
    {
      mod0.GMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMLSAT_MI_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.4 Predict GMLSAT_NI:  
    {
      
      mod0.GMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMLSAT_NI_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.5 Predict GMSST:  
    {
      mod0.GMSST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMSST, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMSST, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMSST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMSST, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMSST, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMSST, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMSST_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.6 Predict GMMSAT:  
    {
      mod0.GMMSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMMSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMMSAT_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.7 Predict TMLSAT:  
    {
      mod0.TMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$TMLSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.TMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.TMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$TMLSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.TMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predTMLSAT_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.8 Predict TMMSAT:  
    {
      mod0.TMMSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$TMMSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.TMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.TMMSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$TMMSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.TMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predTMMSAT_v3/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.9 Predict TMSST_40S_40N:  
    {
      mod0.TMSST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$TMSST_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.TMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.TMSST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$TMSST_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_40S_40N, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.TMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_40S_40N, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predTMSST_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.10 Predict TMSST_25S_25N:  
    {
      mod0.TMSST25 = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$TMSST_25S_25N_, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_25S_25N_, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0.TMSST25, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.TMSST25 = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$TMSST_25S_25N_, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_25S_25N_, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1.TMSST25, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$TMSST_25S_25N_, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predTMSST25_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    # 4.11 Predict IndianOcean:  
    {
      mod0. = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$IndianOcean, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$IndianOcean, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1. = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$IndianOcean, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$IndianOcean, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$IndianOcean, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predIndianOcean_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.12 Predict WPacific:  
    {
      mod0. = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$WPacific, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WPacific, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1. = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$WPacific, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WPacific, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WPacific, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predWPacific_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.13 Predict EPacific:  
    {
      mod0. = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$EPacific, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$EPacific, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1. = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$EPacific, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$EPacific, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$EPacific, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predEPacific_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.14 Predict WAtlantic:  
    {
      mod0. = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$WAtlantic, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WAtlantic, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod0., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1. = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$WAtlantic, alpha = 0,
                         lambda = lambda.seq, foldid = crossclass, 
                         penalty.factor = penalty.factor1, weights = train.weights,
                         nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WAtlantic, X_obs = cur.CRUTEM5_ENS,
                                     mod.glmnet = mod1., foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$WAtlantic, X_obs = cur.CRUTEM5_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predWAtlantic_v4/", 
                        format(CRUTEM5_calendar_[date.ix], "%Y"), "-", format(CRUTEM5_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
  }
}


# -THE END-


