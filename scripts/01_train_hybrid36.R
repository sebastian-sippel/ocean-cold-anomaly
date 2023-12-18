
# ------------------------------------------------------------------------------------
# Train tos-reconstruction based on tos grid cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 17.06.2021
# run on xenon-server:
# screen -S train_tos
# module load R/4.0.3-openblas
# R (-> not R-3.6.1)

library(ncdf4)
library(raster)
library(fields)

library(foreach)
library(doParallel)
library(glmnet)
library(MASS)
library(mvnfast)
library(hydroGOF)

# install.packages("openblasctl", repos="https://hpcran.org")
library(openblasctl)
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
# 1. Train to predict forced response on sea-only record.
# --------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/")

for (mon in 1:12) {
  print(mon)
  
  ## 0.1 Load climate model monthly data & select training model indices:
  {
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_ct.RData", sep=""))
    # load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/_processed/CMIP6.tas_mon", mon, "_ct.RData", sep=""))
      
    # select training model indices:
    CMIP6.tos_monX_ct$M$modcl[which(CMIP6.tos_monX_ct$M$modcl == "UKE")] <- "Had"
    mod.un = unique(CMIP6.tos_monX_ct$M$mod)
    train.ix = list()
    
    for (m in 1:length(mod.un)) {
      print(m)
      cur.mod = mod.un[m]
      
      if (substr(mod.un, 1, 3)[m] %in% c("CES", "EC-", "Had", "MIR", "MPI", "Nor")) {   # models with more than 3 versions... -> choose only two ensemble members per version
        ix = which(cur.mod == CMIP6.tos_monX_ct$M$mod & CMIP6.tos_monX_ct$M$year %in% c(1850:2020) & 
                     !is.na(CMIP6.tos_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tos_monX_ct$Y$GMSST) & !is.na(CMIP6.tos_monX_ct$Y$GMST_FM))
        ens.mem = unique(CMIP6.tos_monX_ct$M$ens.mem[ix])[1:2]
        train.ix[[m]] = which(cur.mod == CMIP6.tos_monX_ct$M$mod & CMIP6.tos_monX_ct$M$year %in% c(1850:2020) & CMIP6.tos_monX_ct$M$ens.mem %in% ens.mem & 
                                !is.na(CMIP6.tos_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tos_monX_ct$Y$GMSST) & !is.na(CMIP6.tos_monX_ct$Y$GMST_FM))
        if (length(train.ix[[m]]) < 330) train.ix[[m]] = integer(0)
      } else {
        ix = which(cur.mod == CMIP6.tos_monX_ct$M$mod & CMIP6.tos_monX_ct$M$year %in% c(1850:2020)  & 
                     !is.na(CMIP6.tos_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tos_monX_ct$Y$GMSST)  & !is.na(CMIP6.tos_monX_ct$Y$GMST_FM) )
        ens.mem = unique(CMIP6.tos_monX_ct$M$ens.mem[ix])[1:3]
        train.ix[[m]] = which(cur.mod == CMIP6.tos_monX_ct$M$mod & CMIP6.tos_monX_ct$M$year %in% c(1850:2020) & CMIP6.tos_monX_ct$M$ens.mem %in% ens.mem & 
                                !is.na(CMIP6.tos_monX_ct$Y$GMLSAT_NI) & !is.na(CMIP6.tos_monX_ct$Y$GMSST)  & !is.na(CMIP6.tos_monX_ct$Y$GMST_FM))
        if (length(train.ix[[m]]) < 495) train.ix[[m]] = integer(0)
      }
    }
    train.mod = cbind(1:64, mod.un, sapply(X = train.ix, FUN = function(x) length(x)))
    train.ix = unlist(train.ix) # length(train.ix)
    # old version:
    # train.ix = which(!is.na(CMIP6.tas_monX_ct$Y$AGMT) & !is.na(CMIP6.tos_monX_ct$Y$AGMT_ANN_fl) & CMIP6.tos_monX_ct$M$ens.mem %in% paste("r", 1:5,"i1p1f1", sep="") & CMIP6.tos_monX_ct$M$year %in% 1850:2020 )
  
      cur.CMIP6.tos_monX_ct = list()
      cur.CMIP6.tos_monX_ct$X = CMIP6.tos_monX_ct$X[train.ix,]
      cur.CMIP6.tos_monX_ct$Y = CMIP6.tos_monX_ct$Y[train.ix,]
      cur.CMIP6.tos_monX_ct$M = CMIP6.tos_monX_ct$M[train.ix,]
      rm(CMIP6.tos_monX_ct)
  }
  
  
  ## load hybrid36/HadSST3 data:
  load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/_processed_Cowtan_hybrid36/hybrid_mon", mon, ".RData", sep=""))
  
  
  ## Run through all time steps from 1850 through 2020
  registerDoParallel(cores = 7)  # used to be 9
  cv.nr.cores = 4 # used to be 5
  
  foreach(date.ix=1:167) %dopar% {   # date.ix = 122
    
    print(HadSST3_calendar_[date.ix])
    
    # 1.1 Select grid cell indices & Check for outliers in data:
    {
      cur.land.weights = land.weights[date.ix,]
      cur.land.weights[which(cur.land.weights == 1)] = NA
      
      # grid.ix_HadSST3 = which(cur.land.weights < 1)
      grid.ix_HadSST3 = which(!is.na(hybrid_36m_[date.ix,]))
      na.val = unique(which(is.na(cur.CMIP6.tos_monX_ct$X), arr.ind = T)[,2])
      grid.ix_CMIP6 = c(1:2592)[-na.val]
      grid.ix_unc = which(HadSST3_meas_sampling_unc_[date.ix,] < 1.5)   # SD of uncertainty < 1°C. results in about 1.25% grid cells removed due to uncertainty.
      grid.ix_ = intersect(intersect(grid.ix_CMIP6, grid.ix_HadSST3), grid.ix_unc)

      ### (2) Check for outliers in CMIP6 data. Both on (i) grid.ix and (ii) train.ix:
      
      ## (i) check where the inter-model standard deviation diverges extremely:
      modcl.un = unique(cur.CMIP6.tos_monX_ct$M$modcl)
      mod.sd = matrix(data = NA, nrow = length(modcl.un), ncol = 72*36)
      for (m in 1:length(modcl.un)) {
        ix = which(cur.CMIP6.tos_monX_ct$M$modcl == modcl.un[m] & cur.CMIP6.tos_monX_ct$M$year %in% 1900:2000)
        mod.sd[m,] = colSds(cur.CMIP6.tos_monX_ct$X[ix,])
      }
      # minimum SD:
      # image.plot(matrix(apply(X = mod.sd, MARGIN = 2, min), 72, 36))
      # get fifth-smallest SD value:
      test = apply(X = mod.sd, MARGIN = 2, FUN=function(x) head(sort(x, decreasing = F))[5])
      # test[c(1:(72*36))[-grid.ix]] = NA
      # image.plot(matrix(test, 72, 36) < 0.05)
      grid.ix = intersect(which(test > 0.05), grid.ix_)
      
      # uncertainty out overall:
      # (1 - length(grid.ix) / length(grid.ix_HadSST4)) * 100  # 4.9 % of data not used:
      # (1 - length(grid.ix_unc) / length(grid.ix_HadSST4)) * 100 # 1.9 % because of large uncertainty.
      # (1 - length(intersect(grid.ix_HadSST4, grid.ix_CMIP6)) / length(grid.ix_HadSST4)) * 100 # 3.17% because of CMIP6
      # (1 - length(grid.ix) / length(grid.ix_)) * 100 # 0.5 % of data because of large CMIP6 disagreement.
      
      # (ii) refine train.ix to remove outliers in training data:
      X_sc = scale(cur.CMIP6.tos_monX_ct$X[,grid.ix], center = T, scale = T)
      # check for extreme outliers:
      outlier.ix = sort(unique(which(X_sc > 8 | X_sc < -8, arr.ind = T)[,1]))
        # hist(as.numeric(cur.CMIP6.tos_monX_ct$M$year[outlier.ix]))
        # cur.CMIP6.tos_monX_ct$M$mod[outlier.ix]
      ## Remove any remaining outliers:
      if ( length(outlier.ix) > 0 ) {
        cur.CMIP6.tos_monX_ct$X = cur.CMIP6.tos_monX_ct$X[-outlier.ix,]
        cur.CMIP6.tos_monX_ct$Y = cur.CMIP6.tos_monX_ct$Y[-outlier.ix,]
        cur.CMIP6.tos_monX_ct$M = cur.CMIP6.tos_monX_ct$M[-outlier.ix,]
      }
    }
    
    # 1.2 Combine CMIP6 training data into one data.frame and save removed outliers:
    {
      CMIP6.df <- list()
      # CMIP6.df$Y <- cur.CMIP6.tos_monX_ct$Y
      # CMIP6.df$M <- cur.CMIP6.tos_monX_ct$M
      CMIP6.df$outlier.ix = outlier.ix
      CMIP6.df$cur.land.weights = cur.land.weights
      CMIP6.df$grid.ix_HadSST3 = grid.ix_HadSST3
      CMIP6.df$grid.ix_CMIP6 = grid.ix_CMIP6
      CMIP6.df$grid.ix_unc = grid.ix_unc
      CMIP6.df$mod.sd = mod.sd
      CMIP6.df$grid.ix = grid.ix
      CMIP6.df$train.ix = train.ix
      CMIP6.df$train.mod = train.mod
    }
      
    ## (3) Define biases and uncertainties according to CRU dataset:
    HadSST3_ENS_anom_[101:200] = HadSST3_ENS_anom_[1:100]
    bias.ens = lapply(X = HadSST3_ENS_anom_, FUN=function(x) x[date.ix,grid.ix])
    unc = HadSST3_meas_sampling_unc_[date.ix,grid.ix]
    
    if (any(is.na(bias.ens[[1]]))) {
      print("NA in bias ensemble!")
      bias.ens = lapply(X = bias.ens, FUN=function(x) {
        na.ix=which(is.na(x))
        x[na.ix] = mean(x, na.rm=T)
        return(x)
        })
    }
    
    
    
    # 1. generate perturbed training dataset. str(X)
    {
      X = list()
      X$X_orig = cur.CMIP6.tos_monX_ct$X[,grid.ix]
      X$X_pert = gen.pert.traintest.data(X = cur.CMIP6.tos_monX_ct$X[,grid.ix], M = cur.CMIP6.tos_monX_ct$M, 
                                         bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1)
      X$X_pert2 = gen.pert.traintest.data(X = cur.CMIP6.tos_monX_ct$X[,grid.ix], M = cur.CMIP6.tos_monX_ct$M, 
                                          bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 2)
      ens.mem.split = gen.pert.traintest.data(X = cur.CMIP6.tos_monX_ct$X[,grid.ix], M = cur.CMIP6.tos_monX_ct$M, 
                                              bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1, ret.ens.mem.split = T)
      ens.ix = gen.pert.traintest.data(X = cur.CMIP6.tos_monX_ct$X[,grid.ix], M = cur.CMIP6.tos_monX_ct$M, 
                                       bias.ens = bias.ens, unc = unc, randomize.ens = F, fact = 1, ret.ens.mem.split = F, ret.ens.ix = T)
    }
    
    # 2. Generate observational input data with uncertainties:
    {
      cur.HadSST3_ENS = list()
      cur.HadSST3_ENS$X_pert = rep.row(hybrid_36m_[date.ix,], n = 200)[,grid.ix] + 
        t(sapply(X = bias.ens, FUN=function(x) x)) +
        sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
      # rmvn(n = 200, mu = rep(0, dim(unc_cov)[1]), sigma = unc_cov)
      cur.HadSST3_ENS$X_pert2 = rep.row(hybrid_36m_[date.ix,], n = 200)[,grid.ix] + 
        2 * t(sapply(X = bias.ens, FUN=function(x) x)) +
        2 * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
      #3 * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
    }
    
    
    # 1.5 Define model training parameters:
    {
        modcl.un = unique(cur.CMIP6.tos_monX_ct$M$modcl)
        crossclass = sapply(X = cur.CMIP6.tos_monX_ct$M$modcl, FUN=function(x) which(x == modcl.un))
        lambda.seq = rev(10 ^ seq(-5,5,length.out=100))
        penalty.factor1 = rep(1, length(grid.ix))
        train.weights = unlist(sapply(X = unique(crossclass), FUN=function(i) rep(1 / length(which(crossclass == unique(crossclass)[i])) / length(unique(crossclass)), length(which(crossclass == unique(crossclass)[i])))))
    }
    
    # 2. TRAIN MODELS WITH DIFFERENT TARGETS:  
    # --------------------------------------------------
    
    # 2.1 Get SVD computations for each setup:
    {
    # mod_p0:
    svd_mod_p0 = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GSAT, alpha = 0,
                          lambda = lambda.seq, foldid = crossclass, 
                          penalty.factor = penalty.factor1, weights = train.weights,
                          nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", ret.svd = T)
    # mod_p1:
    svd_mod_p1 = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GSAT, alpha = 0,
                            lambda = lambda.seq, foldid = crossclass, 
                            penalty.factor = penalty.factor1, weights = train.weights,
                            nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", ret.svd = T)
    
    # mod_gta:
    cov2 = prepare_cov_Pacific_centered(matrix(cur.land.weights, 72, 36),800.0); 
    w = cov2[grid.ix,][,grid.ix]; wi = solve(w); beta.gta = rowSums(wi) / sum(rowSums(wi));
    }

    # Start the clock!
    # ptm <- proc.time()    
    # Stop the clock
    # proc.time() - ptm
    ## runs roughly 4 minutes. 
    ## 4 * 171 * 12 ~ 136.8 h 
    
    
    # 4.1 Predict GSAT:  
    {

        mod0.GSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GSAT, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
        mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GSAT, X_obs = cur.HadSST3_ENS,
                                   mod.glmnet = mod0.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
        
        mod1.GSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GSAT, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
        mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GSAT, X_obs = cur.HadSST3_ENS,
                                   mod.glmnet = mod1.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
        
        # get gta model prediction:
        mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GSAT, X_obs = cur.HadSST3_ENS, 
                                              beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)

        # save models for given time step:
        save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
             file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGSAT/", 
                          format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }

    # 4.2 Predict GMST:  
    {
      
      mod0.GMST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMST_FM, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.GMST, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.GMST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMST_FM, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.GMST, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMST_FM, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGMST/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.3 Predict GMLSAT_MI:  
    {
      mod0.GMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_MI, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_MI, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.GMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_MI, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_MI, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_MI, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGMLSAT_MI/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }

    # 4.4 Predict GMLSAT_NI:  
    {
      
      mod0.GMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_NI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_NI, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.GMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_NI, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_NI, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.GMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMLSAT_NI, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)

      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGMLSAT_NI/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }

    # 4.5 Predict GMSST:  
    {
      mod0.GMSST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GMSST, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMSST, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.GMSST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GMSST, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMSST, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMSST, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)

      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGMSST/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.6 Predict GMMSAT:  
    {
      mod0.GMMSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$GMMSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMMSAT, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.GMMSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$GMMSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMMSAT, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$GMMSAT, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)

      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGMMSAT/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    # 4.7 Predict TMLSAT:  
    {
      mod0.TMLSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$TMLSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.TMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.TMLSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$TMLSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.TMLSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMLSAT_40S_40N, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)

      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predTMLSAT/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }

    # 4.8 Predict TMMSAT:  
    {
      mod0.TMMSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$TMMSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0.TMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1.TMMSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$TMMSAT_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1.TMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMMSAT_40S_40N, X_obs = cur.HadSST3_ENS, 
                                            beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predTMMSAT/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.9 Predict TMSST_40S_40N:  
    {
      mod0. = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tos_monX_ct$Y$TMSST_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMSST_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod0., foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      mod1. = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tos_monX_ct$Y$TMSST_40S_40N, alpha = 0,
                               lambda = lambda.seq, foldid = crossclass, 
                               penalty.factor = penalty.factor1, weights = train.weights,
                               nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMSST_40S_40N, X_obs = cur.HadSST3_ENS,
                                     mod.glmnet = mod1., foldid = crossclass, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # get gta model prediction:
      mod_gta = validate.pert.dataset_gta(X_pert = X, Y = cur.CMIP6.tos_monX_ct$Y$TMSST_40S_40N, X_obs = cur.HadSST3_ENS, 
                                          beta.gta = beta.gta, ens.mem.split = ens.mem.split, ens.ix = ens.ix)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "mod_gta", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predTMSST/", 
                        format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
      rm(mod0., mod_p0, mod1., mod_p1, mod_gta)
    }
    
  }
}




