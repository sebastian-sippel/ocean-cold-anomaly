
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 08.2022
library(hydroGOF)
library(mvnfast)
library(ncdf4)

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_ridge_fun_v2.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/scripts/02a_load_global_observations.R")



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

# Generate perturbed dataset:
gen.pert.data_diff_grid <- function(X, M, bias.ens, unc, grid.ix, ens.ix = 94:200, fact = 1, random.seed = 4) {
  
  X_pert <- matrix(NA, nrow = dim(X)[1], ncol = dim(X)[2]) 
  train.mod.ens = unique(paste(M$mod, M$scen, M$ens.mem, sep="_"))
  mod = M$mod
  scen = M$scen
  ens.mem = M$ens.mem
  year = M$year
  
  train.mod.ens_hist_un = train.mod.ens[which(sapply(X = strsplit(x = train.mod.ens, split = "_"), FUN=function(x) x[2]) == "historical")]
  mod_hist_un = sapply(X = strsplit(x = train.mod.ens_hist_un, split = "_"), FUN=function(x) x[1])
  ens.mem_hist_un = sapply(X = strsplit(x = train.mod.ens_hist_un, split = "_"), FUN=function(x) x[3])
  
  # sample ensemble members:
  set.seed(seed = random.seed)
  sample.ix = sample(x = ens.ix, size = length(train.mod.ens_hist_un), replace = T)
  
  for (date.ix in 1:165) {
    # print(date.ix)
    n = length(grid.ix[[date.ix]])
    
    cur.year = c(1850:2014)[date.ix]
    # ix = sapply(X = 1:length(train.mod.ens_hist_un), FUN=function(me) {
    #  which(mod_hist_un[me] == mod & ens.mem_hist_un[me] == ens.mem & "historical" == scen & cur.year == year)
    #}) ## ask to make sure for each model, ...
    ix = which(cur.year == year)
    
    # get bias ensemble:
    cur.bias.ens = t(sapply(X = sample.ix, FUN=function(i) fact * bias.ens[[i]][date.ix,grid.ix[[date.ix]]]))
    # get uncertainty ensemble:
    if (is.numeric(unc)) {
      cur.unc = t(sapply(X = sample.ix, FUN=function(i) fact * rnorm(n = n, mean = rep(0, n), sd = unc[date.ix,grid.ix[[date.ix]]])))
    } else if (is.list(unc)) {
      cur.unc = fact * rmvn(n = length(sample.ix), mu = rep(0, dim(unc[[date.ix]])[1]), sigma = unc[[date.ix]], ncores = 1)
      # fact * rmvn(n = length(ix), mu = rep(0, n), sigma = unc[[date.ix]], ncores = 1)
    }
    
    X_pert[ix,grid.ix[[date.ix]]] = X[ix,grid.ix[[date.ix]]] + cur.bias.ens + cur.unc
  }
  return(X_pert)
}




# 00.(A-2) load tos data from statistical models:
# -------------------------------------------------------
{
  names.vec = c("GSAT", "GMSST", "GMLSAT_MI", "GMLSAT_NI", "GMST_FM", "GMMSAT", "TMLSAT_40S_40N", "TMMSAT_40S_40N", 
                "TMSST_40S_40N")
  
  file.names = c("tos_predGSAT", "tos_predGMSST", "tos_predGMLSAT_MI", "tos_predGMLSAT_NI",
                 "tos_predGMST", "tos_predGMMSAT", "tos_predTMLSAT", "tos_predTMMSAT",
                 "tos_predTMSST")
  
  
  # OBS_hybrid36_hybrid36ervational predictions:
  OBS_hybrid36.tos = list(); OBS_hybrid36.tos[1:9] = sapply(1:9, FUN=function(i) OBS_hybrid36.tos[[i]] <- list() ); names(OBS_hybrid36.tos) = names.vec
  OBS_hybrid36.tos_ = list(); OBS_hybrid36.tos_[1:9] = sapply(1:9, FUN=function(i) OBS_hybrid36.tos_[[i]] <- list() ); names(OBS_hybrid36.tos_) = names.vec
  
  # 4EVALUATION:
  # Data frame that contains target values (for evaluation) and is the same for all different targets:
  CMIP6_hybrid36.tos.df = list(); CMIP6_hybrid36.tos.df$mon = list(); CMIP6_hybrid36.tos.df$ann = list()   
  
  # 4FULL_LENGTH:
  CMIP6_hybrid36.tos_all.df = list(); CMIP6_hybrid36.tos_all.df$mon = list(); CMIP6_hybrid36.tos_all.df$ann = list()
  
  for (mon in 1:12) {
    print(mon)
    
    ## 0.1 load CRU data:
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/_processed_Cowtan_hybrid36/hybrid_mon", mon, ".RData", sep=""))
    
    ## 0.1 Load climate model monthly data & select training model indices:
    {
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tos_mon", mon, "_ct.RData", sep=""))
      
      # select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/tos_predGSAT/1850-01.RData", sep=""))
      modcl.un = unique(CMIP6.tos_monX_ct$M$modcl)
      train.ix = CMIP6.df$train.ix
    }
    
    ## 0.2 Set up lists for respective month:
    {
      # OBS_hybrid36_hybrid36ervations: Set up list for respective month:
      {
        for (i in 1:9) { OBS_hybrid36.tos[[i]][[mon]] = list() }
      }
      
      # 4evaluation (training data!):
      {
        CMIP6_hybrid36.tos.df$mon[[mon]] = list()
        CMIP6_hybrid36.tos.df$mon[[mon]]$Y = CMIP6.tos_monX_ct$Y[train.ix,]
        CMIP6_hybrid36.tos.df$mon[[mon]]$M = CMIP6.tos_monX_ct$M[train.ix,]
        CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat = list();
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta = list();
        for (l in 1:9) {
          CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[l]] = lapply(X = 1:2, FUN=function(i) { ret.df = data.frame(matrix(NA, nrow = length(train.ix), ncol = 6)); names(ret.df) = c("pt0.min", "pt0.1se", "pt1.min", "pt1.1se", "pt2.min", "pt2.1se"); return(ret.df) })
          CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[l]]$mod_gta = data.frame(matrix(NA, nrow = length(train.ix), ncol = 6)); names(CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[l]]$mod_gta) = c("pt0", "pt1", "pt2", "pt0.raw", "pt1.raw", "pt2.raw")
          names(CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[l]]) = c("mod_p0", "mod_p1", "mod_gta")
          CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[l]] = list();
        }
        names(CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat) <- names(CMIP6_hybrid36.tos.df$mon[[mon]]$beta) <- names.vec
        CMIP6_hybrid36.tos.df$mon[[mon]]$train.ix = train.ix
      }
      
      # 4full-length analysis:
      {
        CMIP6_hybrid36.tos_all.df$mon[[mon]] = list()
        CMIP6_hybrid36.tos_all.df$mon[[mon]]$Y = CMIP6.tos_monX_ct$Y[names.vec]
        # sapply(X = names.vec, FUN=function(x) CMIP6_hybrid36.tas_monX_ct$Y[x])
        CMIP6_hybrid36.tos_all.df$M = CMIP6.tos_monX_ct$M
        CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat = data.frame(matrix(NA, nrow = dim(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Y)[1], ncol = dim(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Y)[2]))
        names(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat) = names.vec
        CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat_pt1 = data.frame(matrix(NA, nrow = dim(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Y)[1], ncol = dim(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Y)[2]))
        names(CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat_pt1) = names.vec
      }
      # rm(CMIP6_hybrid36.tas_monX_ct); rm(CMIP6_hybrid36.df); rm(mod_gta); rm(mod_p0); rm(mod_p1); rm(mod_p2) 
    }
    
    ## 0.3 Generate dataset with bias/uncertainty perturbations for the whole of CMIP6_hybrid36:
    {
      grid.ix = list(); unc_cov = list()
      for (date.ix in 1:167) {
        print(date.ix)
        
        # Generate uncertainty estimates and bias ensemble:
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/", "tos_predGSAT","/", format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))

        ## (3) Define biases and uncertainties according to CRU dataset:
        HadSST3_ENS_anom_[101:200] = HadSST3_ENS_anom_[1:100]
        bias.ens = lapply(X = HadSST3_ENS_anom_, FUN=function(x) x[date.ix,CMIP6.df$grid.ix])
        unc = HadSST3_meas_sampling_unc_[date.ix,CMIP6.df$grid.ix]
        
        if (any(is.na(bias.ens[[1]]))) {
          print("NA in bias ensemble!")
          bias.ens = lapply(X = bias.ens, FUN=function(x) {
            na.ix=which(is.na(x))
            x[na.ix] = mean(x, na.rm=T)
            return(x)
          })
        }
        
        grid.ix[[date.ix]] = CMIP6.df$grid.ix
        unc_cov[[date.ix]] = diag(unc) ^ 2 # + diag(HadSST4_sampling_unc_[date.ix,grid.ix[[date.ix]]]) ^ 2 + diag(HadSST4_measurement_unc_[date.ix,grid.ix[[date.ix]]]) ^ 2
      }
      X_pert = gen.pert.data_diff_grid(X = CMIP6.tos_monX_ct$X, M = CMIP6.tos_monX_ct$M, bias.ens = HadSST3_ENS_anom_, unc = unc_cov, grid.ix = grid.ix, ens.ix = 94:200, fact = 1, random.seed = 4)
    }
    
    for (date.ix in 1:167) {
      print(date.ix)
      
      # For LOOP OVER ALL variables
      # ---------------------------------------------
      for (i in 1:length(names.vec)) {
        cur.name = names.vec[i]
        print(cur.name)
        
        # Load trained model & Prepare OBS_hybrid36 dataset:
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/hybrid36/", file.names[i],"/", format(HadSST3_calendar_[date.ix], "%Y"), "-", format(HadSST3_calendar_[date.ix], "%m"), ".RData", sep=""))
        OBS_hybrid36.tos[[cur.name]][[mon]][[date.ix]] = list()
        OBS_hybrid36.tos[[cur.name]][[mon]][[date.ix]]$mod_p0 = mod_p0; OBS_hybrid36.tos[[cur.name]][[mon]][[date.ix]]$mod_p1 = mod_p1; OBS_hybrid36.tos[[cur.name]][[mon]][[date.ix]]$mod_gta = mod_gta
        
        # select model coefficients:
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]] = list()
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$mod_p0 = mod_p0$beta
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$mod_p0_a0 = mod_p0$a0
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$mod_p1 = mod_p1$beta
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$mod_p1_a0 = mod_p1$a0
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$mod_gta = mod_gta$beta
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$grid.ix = CMIP6.df$grid.ix
        CMIP6_hybrid36.tos.df$mon[[mon]]$beta[[cur.name]][[date.ix]]$cur.land.weights = CMIP6.df$cur.land.weights
        
        # reconstruct Yhat for full data.frame:
        # only select those dates/times that are from the mon/date.ix:
        cur.ix = which(CMIP6_hybrid36.tos.df$mon[[mon]]$M$year == format(HadSST3_calendar_[date.ix], "%Y") & CMIP6_hybrid36.tos.df$mon[[mon]]$M$mon == as.numeric(format(HadSST3_calendar_[date.ix], "%m")))
        
        CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[cur.name]]$mod_p0[cur.ix,] = fill.na_from_outlier(Yhat = mod_p0$Yhat, outlier.ix = CMIP6.df$outlier.ix, transform.Yhat = T)[cur.ix,]
        CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[cur.name]]$mod_p1[cur.ix,] = fill.na_from_outlier(Yhat = mod_p1$Yhat, outlier.ix = CMIP6.df$outlier.ix, transform.Yhat = T)[cur.ix,]
        CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[cur.name]]$mod_gta[cur.ix,1:3] = fill.na_from_outlier(Yhat = sapply(mod_gta$Yhat, FUN=function(x) x), outlier.ix = CMIP6.df$outlier.ix, transform.Yhat = F)[cur.ix,]
        CMIP6_hybrid36.tos.df$mon[[mon]]$Yhat[[cur.name]]$mod_gta[cur.ix,4:6] = fill.na_from_outlier(Yhat = sapply(mod_gta$Yhat_raw, FUN=function(x) x), outlier.ix = CMIP6.df$outlier.ix, transform.Yhat = F)[cur.ix,]
        
        # Predict mod_p1 for the long time series:
        cur.ix_all = which(CMIP6_hybrid36.tos_all.df$M$year == format(HadSST3_calendar_[date.ix], "%Y") & CMIP6_hybrid36.tos_all.df$M$mon == as.numeric(format(HadSST3_calendar_[date.ix], "%m")))
        CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat[[cur.name]][cur.ix_all] = CMIP6.tos_monX_ct$X[cur.ix_all,CMIP6.df$grid.ix] %*% mod_p1$beta[,1] + mod_p1$a0[1]
        CMIP6_hybrid36.tos_all.df$mon[[mon]]$Yhat_pt1[[cur.name]][cur.ix_all] = X_pert[cur.ix_all,CMIP6.df$grid.ix] %*% mod_p1$beta[,1] + mod_p1$a0[1]
      }
    }
  }
  
  
  # 2. Process into time series of annual and monthly OBS_hybrid36ervations:
  # ---------------------------------------------------------------
  
  # 2.1 Process monthly / annual OBS_hybrid36ervational predictions:
  # ---------------------------------------------------------------
  
  {      
    # for loop over 8 different prediction metrics:
    for (i in 1:length(names.vec)) {
      print(i)
      cur.name = names.vec[i]
      
      OBS_hybrid36.tos_[[cur.name]]$mon = list(); OBS_hybrid36.tos_[[cur.name]]$ann = list();
      OBS_hybrid36.tos_[[cur.name]]$mon$mod_p0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_1_05))
      OBS_hybrid36.tos_[[cur.name]]$mon$mod_p1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_1_05))
      OBS_hybrid36.tos_[[cur.name]]$mon$mod_p0_min = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_min))
      OBS_hybrid36.tos_[[cur.name]]$mon$mod_p1_min = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_min))
      OBS_hybrid36.tos_[[cur.name]]$mon$mod_gta = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_gta$obs_pred))
      
      OBS_hybrid36.tos_[[cur.name]]$ann$mod_p0 = sapply(X = 1:167, FUN=function(i) rowMeans(sapply(X = OBS_hybrid36.tos_[[cur.name]]$mon$mod_p0, FUN=function(x) x[1:200,i])))
      OBS_hybrid36.tos_[[cur.name]]$ann$mod_p1 = sapply(X = 1:167, FUN=function(i) rowMeans(sapply(X = OBS_hybrid36.tos_[[cur.name]]$mon$mod_p1, FUN=function(x) x[1:200,i])))
      OBS_hybrid36.tos_[[cur.name]]$ann$mod_p0_min = sapply(X = 1:167, FUN=function(i) rowMeans(sapply(X = OBS_hybrid36.tos_[[cur.name]]$mon$mod_p0_min, FUN=function(x) x[1:200,i])))
      OBS_hybrid36.tos_[[cur.name]]$ann$mod_p1_min = sapply(X = 1:167, FUN=function(i) rowMeans(sapply(X = OBS_hybrid36.tos_[[cur.name]]$mon$mod_p1_min, FUN=function(x) x[1:200,i])))
      OBS_hybrid36.tos_[[cur.name]]$ann$mod_gta = sapply(X = 1:167, FUN=function(i) rowMeans(sapply(X = OBS_hybrid36.tos_[[cur.name]]$mon$mod_gta, FUN=function(x) x[1:200,i])))
      OBS_hybrid36.tos_[[cur.name]]$ann$HadSST3_calendar_ = HadSST3_calendar_
      
    }
    save(list = c("OBS_hybrid36.tos_"), 
         file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/OBS_hybrid36.tos.RData")
  }
  
  # 2.2 4EVALUATION: Process MSE and annual files:
  # ---------------------------------------------------------------
  {
    CMIP6_hybrid36.tos.df$mon_MSE = list(); CMIP6_hybrid36.tos.df$ann_MSE = list(); CMIP6_hybrid36.tos.df$ann_RMSE = list(); CMIP6_hybrid36.tos.df$mon_ME_sq = list();
    for (i in 1:length(names.vec)) {
      print(i)
      cur.name = names.vec[i]
      
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]] = list()
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_p0_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$MSE$pt0[y$mod_p0$lambda.1_05]))
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_p0_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$MSE$pt1[y$mod_p0$lambda.1_05]))
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_p1_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$MSE$pt0[y$mod_p1$lambda.1_05]))
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_p1_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$MSE$pt1[y$mod_p1$lambda.1_05]))
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_gta_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_gta$MSE$pt0))
      CMIP6_hybrid36.tos.df$mon_MSE[[cur.name]]$mod_gta_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) y$mod_gta$MSE$pt1))
      
      CMIP6_hybrid36.tos.df$mon_ME_sq[[cur.name]] = list()
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_p0_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p0$me.by.ensreal$pt0[,y$mod_p0$lambda.1_05]^2)))
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_p0_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p0$me.by.ensreal$pt1[,y$mod_p0$lambda.1_05]^2)))
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_p1_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p1$me.by.ensreal$pt0[,y$mod_p0$lambda.1_05]^2)))
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_p1_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p1$me.by.ensreal$pt1[,y$mod_p0$lambda.1_05]^2)))
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_gta_pt0 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_gta$me.by.ensreal$pt0^2)))
      CMIP6_hybrid36.tos.df$mon_ME_sq$mod_gta_pt1 = lapply(X = OBS_hybrid36.tos[[cur.name]], FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_gta$me.by.ensreal$pt1^2)))
      
      # Process CMIP6_hybrid36 predictions into annual files (training eval):
      CMIP6_hybrid36.tos.df$ann[[cur.name]] = list()
      CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_p0 = data.frame(matrix(data = rowMeans(sapply(X = CMIP6_hybrid36.tos.df$mon, FUN=function(x) c(as.matrix(x$Yhat[[cur.name]]$mod_p0)))), ncol = 6))
      CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_p1 = data.frame(matrix(data = rowMeans(sapply(X = CMIP6_hybrid36.tos.df$mon, FUN=function(x) c(as.matrix(x$Yhat[[cur.name]]$mod_p1)))), ncol = 6))
      CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_gta = data.frame(matrix(data = rowMeans(sapply(X = CMIP6_hybrid36.tos.df$mon, FUN=function(x) c(as.matrix(x$Yhat[[cur.name]]$mod_gta)))), ncol = 6))
      
      names(CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_p0) <- names(CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_p1) <- names(CMIP6_hybrid36.tos.df$mon[[1]]$Yhat$GSAT$mod_p0)
      names(CMIP6_hybrid36.tos.df$ann[[cur.name]]$mod_gta) <- names(CMIP6_hybrid36.tos.df$mon[[1]]$Yhat$GSAT$mod_gta)
    }
    save(list = c("CMIP6_hybrid36.tos.df"), 
         file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6_hybrid36.tos.df.RData")
  }
  
  # 2.3 LONG RECONSTRUCTION: 
  # ---------------------------------------------------------------
  {
    CMIP6_hybrid36.tos_all.df$ann$Yhat = list(); CMIP6_hybrid36.tos_all.df$ann$Y = list(); 
    for (i in 1:length(names.vec)) {
      print(i)
      cur.name = names.vec[i]
      CMIP6_hybrid36.tos_all.df$ann$Yhat[[cur.name]] = rowMeans(sapply(X = CMIP6_hybrid36.tos_all.df$mon, FUN=function(x) x$Yhat[[cur.name]]))
      CMIP6_hybrid36.tos_all.df$ann$Yhat_pt1[[cur.name]] = rowMeans(sapply(X = CMIP6_hybrid36.tos_all.df$mon, FUN=function(x) x$Yhat_pt1[[cur.name]]))
      CMIP6_hybrid36.tos_all.df$ann$Y[[cur.name]] = rowMeans(sapply(X = CMIP6_hybrid36.tos_all.df$mon, FUN=function(x) x$Y[[cur.name]]))
    }
    save(list = c("CMIP6_hybrid36.tos_all.df"), 
         file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6_hybrid36.tos_all.df.RData")
  }
}






