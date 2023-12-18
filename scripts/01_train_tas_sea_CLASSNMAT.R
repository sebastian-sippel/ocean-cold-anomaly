
# ------------------------------------------------------------------------------------
# Train tas-reconstruction based on tas-air maring grid cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 17.06.2021
# run on co2-server:
# screen -S train_CLASSnmat
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
openblas_set_num_threads(num_threads = 5)


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
# 1. Train to predict forced response and AGMT based on sea-only record.
# --------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/")
# setwd("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/")


for (mon in 1:12) {
  print(mon) # mon = 11
  
  ## load climate model monthly data:
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
  load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/_processed_CRU/CLASSnmat_mon", mon, ".RData", sep=""))
  
  
  ## Run through all time steps from 1850 through 2020
  registerDoParallel(cores = 6)  # used to be 9
  cv.nr.cores = 4 # used to be 5
  
  foreach(date.ix=60:78) %dopar% {   # date.ix = 122; CLASSnmat_calendar_[140]
    
    print(CLASSnmat_calendar_[date.ix])
    # date.ix = 56 # Jan 1920: 71
    
    
    ### (1) look at sea mask for Jan 1920:
    grid.ix_CLASSnmat = which(!is.na(CLASSnmat_[date.ix,]))
    # image.plot(matrix(CLASSnmat_total_unc_[date.ix,], 72, 36))
    grid.ix_unc = which(CLASSnmat_total_unc_[date.ix,] < 3)   # SD of uncertainty < 1°C. results in about 1.25% grid cells removed due to uncertainty.
    grid.ix_ = intersect(grid.ix_CLASSnmat, grid.ix_unc)
    
    
    ## (i) check where the inter-model standard deviation diverges extremely:
    modcl.un = unique(cur.CMIP6.tas_monX_ct$M$modcl)
    mod.sd = matrix(data = NA, nrow = length(modcl.un), ncol = 72*36)
    for (m in 1:length(modcl.un)) {
      ix = which(cur.CMIP6.tas_monX_ct$M$modcl == modcl.un[m] & cur.CMIP6.tas_monX_ct$M$year %in% 1900:2000)
      mod.sd[m,] = colSds(cur.CMIP6.tas_monX_ct$X[ix,])
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
    X_sc = scale(cur.CMIP6.tas_monX_ct$X[,grid.ix], center = T, scale = T)
    # check for extreme outliers:
    outlier.ix = sort(unique(which(X_sc > 8 | X_sc < -8, arr.ind = T)[,1]))
    # hist(as.numeric(cur.CMIP6.tas_monX_ct$M$year[outlier.ix]))
    # cur.CMIP6.tas_monX_ct$M$mod[outlier.ix]
    ## Remove any remaining outliers:
    if ( length(outlier.ix) > 0 ) {
      cur.CMIP6.tas_monX_ct$X = cur.CMIP6.tas_monX_ct$X[-outlier.ix,]
      cur.CMIP6.tas_monX_ct$Y = cur.CMIP6.tas_monX_ct$Y[-outlier.ix,]
      cur.CMIP6.tas_monX_ct$M = cur.CMIP6.tas_monX_ct$M[-outlier.ix,]
    }
    
    # 1.2 Combine CMIP6 training data into one data.frame and save removed outliers:
    {
      CMIP6.df <- list()
      CMIP6.df$Y <- cur.CMIP6.tas_monX_ct$Y
      CMIP6.df$M <- cur.CMIP6.tas_monX_ct$M
      CMIP6.df$train.ix = train.ix
      CMIP6.df$outlier.ix = outlier.ix
      CMIP6.df$grid.ix_unc = grid.ix_unc
      CMIP6.df$grid.ix = grid.ix
      CMIP6.df$mod.sd = mod.sd
    }
      
    ## (3) Define biases and uncertainties according to CRU dataset:
    # bias.ens = lapply(X = HadSST4_ENS_anom_, FUN=function(x) x[date.ix,grid.ix])
    
    # HadSST4 - Covariance Matrices: 
    # load and process HadSST.4.0.1.0_error_covariance matrices:
    {
      # f.name = paste("/net/h2o/climphys1/sippels/_DATASET/CRU/HadSST4/HadSST.4.0.1.0/HadSST.4.0.1.0_error_covariance/HadSST.4.0.1.0_error_covariance_", 
      #               substr(HadSST4_calendar_[date.ix], 1,4), "/HadSST.4.0.1.0_error_covariance_",
      #               substr(HadSST4_calendar_[date.ix], 1, 4), substr(HadSST4_calendar_[date.ix], 6, 7), ".nc", sep="")
      # tos_cov = ncvar_get(nc = nc_open(filename = f.name), varid = "tos_cov")[transfer.CRU.to.CMIP6.grid,transfer.CRU.to.CMIP6.grid][grid.ix,grid.ix]
      
      # test = ncvar_get(nc = nc_open(filename = f.name), varid = "tos_cov")
      # test the equivalence of matrix reshaping:
      # test[248, 246]
      # test[transfer.CRU.to.CMIP6.grid,transfer.CRU.to.CMIP6.grid][transfer.CRU.to.CMIP6.grid[248],transfer.CRU.to.CMIP6.grid[246]]
      # image.plot(matrix(diag(test[transfer.CRU.to.CMIP6.grid,transfer.CRU.to.CMIP6.grid]), 72, 36))
      # image.plot(matrix(HadSST4_sampling_unc_[date.ix,], 72, 36))
      # 
      # unc_cov = tos_cov  # + diag(HadSST4_sampling_unc_[date.ix,grid.ix]) ^ 2 + diag(HadSST4_measurement_unc_[date.ix,grid.ix]) ^ 2
      # image(unc_cov, zlim = c(0, 0.1))
      unc = CLASSnmat_total_unc_[date.ix, grid.ix]
    }
    
    ## CHECK IF ALL VARIABLES UNCERTAINTY LOOKS REASONABLE:
    # bias.ens_ = sapply(X = HadSST4_ENS_anom_, FUN=function(x) x[date.ix,grid.ix])
    # par(mfrow=c(1,1))
    # hist(bias.ens_)
    # hist(sqrt(diag(tos_cov)))
    # hist(HadSST4_measurement_unc_[date.ix,grid.ix])
    # hist(HadSST4_sampling_unc_[date.ix,grid.ix])
      
    # image.plot(matrix(cur.land.weights, 72, 36))
    # image.plot(matrix(cur.CMIP6.tos_monX_ct$X[10300,], 72, 36))
    
    # Generate training/testing dataset with perturbations:
    # X = cur.CMIP6.tas_monX_ct$X[,grid.ix]
    # Y = cur.CMIP6.tas_monX_ct$Y$AGMT_f
    # M = cur.CMIP6.tas_monX_ct$M
    # bias.ens = lapply(X = CRUTEM5_ENS_anom_, FUN=function(x) x[date.ix,grid.ix])
    # unc = CRUTEM5_sampling_unc_[date.ix,grid.ix]
    
    
    # 1.3 Generate perturbed training dataset. str(X)
    {
        X = list()
        X$X_orig = cur.CMIP6.tas_monX_ct$X[,grid.ix]
        X$X_pert = gen.pert.traintest.data_noENS(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                           unc = unc, fact = 1)
        X$X_pert2 = gen.pert.traintest.data_noENS(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                            unc = unc, fact = 2)
        ens.mem.split = gen.pert.traintest.data_noENS(X = cur.CMIP6.tas_monX_ct$X[,grid.ix], M = cur.CMIP6.tas_monX_ct$M, 
                                                unc = unc, fact = 1, ret.ens.mem.split = T)
    }

        
    # 1.4 Generate observational input data with uncertainties:
    {
        cur.CLASSnmat = list()
        # cur.CLASSnmat$X_orig = CLASSnmat_[date.ix,grid.ix]
        cur.CLASSnmat$X_pert = rep.row(CLASSnmat_[date.ix,], n = 200)[,grid.ix] + 
          # CLASSnmat_[date.ix,grid.ix] + 
          sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
        cur.CLASSnmat$X_pert2 = rep.row(CLASSnmat_[date.ix,], n = 200)[,grid.ix] + 
          # CLASSnmat_[date.ix,grid.ix] + 
          2 * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = 200, mean = 0, sd = unc[i]))
    }

        
    # 1.5 Define model training parameters:
    {
      # Define training parameters:
      modcl.un = unique(cur.CMIP6.tas_monX_ct$M$modcl)
      crossclass = sapply(X = cur.CMIP6.tas_monX_ct$M$modcl, FUN=function(x) which(x == modcl.un))
      lambda.seq = rev(10 ^ seq(-5,5,length.out=100))
      penalty.factor = rep(1, length(grid.ix))
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
      
    }
    
    # 4.1 Predict GSAT:  
    {
      
      mod0.GSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GSAT, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GSAT, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GSAT, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGSAT/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    # 4.2 Predict GMST:  
    {
      mod0.GMST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GMST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMST_FM, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GMST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGMST/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    # 4.3 Predict GMLSAT_MI:  
    {
      mod0.GMLSAT_MI = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GMLSAT_MI, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMLSAT_MI = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, alpha = 0,
                             lambda = lambda.seq, foldid = crossclass, 
                             penalty.factor = penalty.factor1, weights = train.weights,
                             nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_MI, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GMLSAT_MI, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGMLSAT_MI/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    
    # 4.4 Predict GMLSAT_NI:  
    {
      mod0.GMLSAT_NI = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GMLSAT_NI, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMLSAT_NI = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMLSAT_NI, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GMLSAT_NI, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGMLSAT_NI/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    
    # 4.5 Predict GMSST:  
    {
      mod0.GMSST = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMSST, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMSST, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMSST = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMSST, alpha = 0,
                                  lambda = lambda.seq, foldid = crossclass, 
                                  penalty.factor = penalty.factor1, weights = train.weights,
                                  nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMSST, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GMSST, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGMSST/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    # 4.6 Predict GMMSAT:  
    {
      mod0.GMMSAT = cv.glmnet2(x = X$X_orig, y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p0)
      mod_p0 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod0.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      mod1.GMMSAT = cv.glmnet2(x = X$X_pert, y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, alpha = 0,
                              lambda = lambda.seq, foldid = crossclass, 
                              penalty.factor = penalty.factor1, weights = train.weights,
                              nr.cores = cv.nr.cores, nr.subsample = NULL, cv = "leave.model.out", nsim = NULL, include.all.mod.fit = F, ret.beta.only = T, fit.fun = "weighted.ridge", svd.list = svd_mod_p1)
      mod_p1 = validate.pert.dataset(X_pert = X, Y = cur.CMIP6.tas_monX_ct$Y$GMMSAT, X_obs = cur.CLASSnmat,
                                     mod.glmnet = mod1.GMMSAT, foldid = crossclass, ens.mem.split = ens.mem.split)
      
      
      # save models for given time step:
      save(list = c("mod_p0", "mod_p1", "CMIP6.df"), 
           file = paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/CLASSNMAT/tas_predGMMSAT/", 
                        format(CLASSnmat_calendar_[date.ix], "%Y"), "-", format(CLASSnmat_calendar_[date.ix], "%m"), ".RData", sep=""))
    }
    
    
    
    

  }
  
}




