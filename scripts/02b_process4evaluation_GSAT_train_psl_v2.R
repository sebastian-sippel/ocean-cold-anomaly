
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 08.2022
library(hydroGOF)

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")
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




# 00.(A-1) load psl data from statistical models:
# -------------------------------------------------------
{
  GSAT.psl = list(); 
  # CMIP6.psl.df. Data frame contains target values (for evaluation) and is the same for all different targets:
  CMIP6.psl.df = list(); CMIP6.psl.df$mon = list(); CMIP6.psl.df$ann = list()   
  # CMIP6.psl.GSAT.Yhat.df Data frame contains predictions for GSAT:
  CMIP6.psl.GSAT.Yhat.df = list(); CMIP6.psl.GSAT.Yhat.df$mon = list(); CMIP6.psl.GSAT.Yhat.df$ann = list(); 
  CMIP6.psl.GSAT.beta = list(); 
  
  
  for (mon in 1:12) {
    print(mon)
    
    
    ## 0.1 load CRU data:
    load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_HadSLP2/HadSLP2_mon", mon, ".RData", sep=""))
    
    ## 0.1 Load climate model monthly data & select training model indices:
    {
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.psl_mon", mon, "_ct.RData", sep=""))
      
      # select training model indices:
      load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/HadSLP2_predGSAT_v2/1850-01.RData", sep=""))
      modcl.un = unique(CMIP6.psl_monX_ct$M$modcl)
      train.ix = CMIP6.df$train.ix
      
      CMIP6.psl.df$mon[[mon]] = list()
      CMIP6.psl.df$mon[[mon]]$Y = CMIP6.psl_monX_ct$Y[train.ix,]
      CMIP6.psl.df$mon[[mon]]$M = CMIP6.psl_monX_ct$M[train.ix,]
      rm(CMIP6.tas_monX_ct); rm(CMIP6.df); rm(mod_gta); rm(mod_p0); rm(mod_p1);  
    }
    
    # 0.2 initialize CMIP6 Yhat lists:
    {
      # GSAT:
      CMIP6.psl.GSAT.Yhat.df$mon[[mon]] = lapply(X = 1:4, FUN=function(i) { ret.df = data.frame(matrix(NA, nrow = length(train.ix), ncol = 4)); names(ret.df) = c("pt0.min", "pt0.1se", "pt1.min", "pt1.1se"); return(ret.df) })
      names(CMIP6.psl.GSAT.Yhat.df$mon[[mon]]) = c("mod_p0", "mod_p1", "mod_p0_IV", "mod_p1_IV")
    }
    
    # Set up list for respective month:
    GSAT.psl[[mon]] = list(); 
    CMIP6.psl.GSAT.beta[[mon]] = list(); 
    
    for (date.ix in 1:155) {
      print(date.ix)
      
      # 1.1 GSAT: 
      # ---------------------------------------------
      {
        # Load trained model:
        load(paste("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/HadSLP2_predGSAT_v2/", format(HadSLP2_calendar_[date.ix], "%Y"), "-", format(HadSLP2_calendar_[date.ix], "%m"), ".RData", sep=""))

        # select model coefficients:
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]] = list()
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p0 = list(); CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p1 = list(); CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p0_IV = list(); CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p1_IV = list(); 
        
        # get higher \lambda values:
        mod_p0$lambda.1_05 = which(mod_p0$MSE$pt0 < min(mod_p0$MSE$pt0) * 1.05)[1]
        mod_p0$lambda.1_1 = which(mod_p0$MSE$pt0 < min(mod_p0$MSE$pt0) * 1.1)[1]
        mod_p0$lambda.1_15 = which(mod_p0$MSE$pt0 < min(mod_p0$MSE$pt0) * 1.15)[1]
        mod_p0$lambda.1_2 = which(mod_p0$MSE$pt0 < min(mod_p0$MSE$pt0) * 1.2)[1]
        mod_p0$lambda.1_25 = which(mod_p0$MSE$pt0 < min(mod_p0$MSE$pt0) * 1.25)[1]
        mod_p1$lambda.1_05 = which(mod_p1$MSE$pt1 < min(mod_p1$MSE$pt1) * 1.05)[1]
        mod_p1$lambda.1_1 = which(mod_p1$MSE$pt1 < min(mod_p1$MSE$pt1) * 1.1)[1]
        mod_p1$lambda.1_15 = which(mod_p1$MSE$pt1 < min(mod_p1$MSE$pt1) * 1.15)[1]
        mod_p1$lambda.1_2 = which(mod_p1$MSE$pt1 < min(mod_p1$MSE$pt1) * 1.2)[1]
        mod_p1$lambda.1_25 = which(mod_p1$MSE$pt1 < min(mod_p1$MSE$pt1) * 1.25)[1]
        
        mod_p0_IV$lambda.1_05 = which(mod_p0_IV$MSE$pt0 < min(mod_p0_IV$MSE$pt0) * 1.05)[1]
        mod_p0_IV$lambda.1_1 = which(mod_p0_IV$MSE$pt0 < min(mod_p0_IV$MSE$pt0) * 1.1)[1]
        mod_p0_IV$lambda.1_15 = which(mod_p0_IV$MSE$pt0 < min(mod_p0_IV$MSE$pt0) * 1.15)[1]
        mod_p0_IV$lambda.1_2 = which(mod_p0_IV$MSE$pt0 < min(mod_p0_IV$MSE$pt0) * 1.2)[1]
        mod_p0_IV$lambda.1_25 = which(mod_p0_IV$MSE$pt0 < min(mod_p0_IV$MSE$pt0) * 1.25)[1]
        
        mod_p1_IV$lambda.1_05 = which(mod_p1_IV$MSE$pt1 < min(mod_p1_IV$MSE$pt1) * 1.05)[1]
        mod_p1_IV$lambda.1_1 = which(mod_p1_IV$MSE$pt1 < min(mod_p1_IV$MSE$pt1) * 1.1)[1]
        mod_p1_IV$lambda.1_15 = which(mod_p1_IV$MSE$pt1 < min(mod_p1_IV$MSE$pt1) * 1.15)[1]
        mod_p1_IV$lambda.1_2 = which(mod_p1_IV$MSE$pt1 < min(mod_p1_IV$MSE$pt1) * 1.2)[1]
        mod_p1_IV$lambda.1_25 = which(mod_p1_IV$MSE$pt1 < min(mod_p1_IV$MSE$pt1) * 1.25)[1]
        
        # Identify new predictions based on different \lambda choices:
        mod_p0$obs_pred_pert0_lamvar = mod_p0$obs_pred_pert0[,c(mod_p0$lambda.min, mod_p0$lambda.1_05, mod_p0$lambda.1_1, mod_p0$lambda.1_15, mod_p0$lambda.1_2, mod_p0$lambda.1_25)]
        mod_p0$obs_pred_pert1_lamvar = mod_p0$obs_pred_pert1[,c(mod_p0$lambda.min, mod_p0$lambda.1_05, mod_p0$lambda.1_1, mod_p0$lambda.1_15, mod_p0$lambda.1_2, mod_p0$lambda.1_25)]
        mod_p1$obs_pred_pert0_lamvar = mod_p1$obs_pred_pert0[,c(mod_p1$lambda.min, mod_p1$lambda.1_05, mod_p1$lambda.1_1, mod_p1$lambda.1_15, mod_p1$lambda.1_2, mod_p1$lambda.1_25)]
        mod_p1$obs_pred_pert1_lamvar = mod_p1$obs_pred_pert1[,c(mod_p1$lambda.min, mod_p1$lambda.1_05, mod_p1$lambda.1_1, mod_p1$lambda.1_15, mod_p1$lambda.1_2, mod_p1$lambda.1_25)]

        mod_p0_IV$obs_pred_pert0_lamvar = mod_p0_IV$obs_pred_pert0[,c(mod_p0_IV$lambda.min, mod_p0_IV$lambda.1_05, mod_p0_IV$lambda.1_1, mod_p0_IV$lambda.1_15, mod_p0_IV$lambda.1_2, mod_p0_IV$lambda.1_25)]
        mod_p0_IV$obs_pred_pert1_lamvar = mod_p0_IV$obs_pred_pert1[,c(mod_p0_IV$lambda.min, mod_p0_IV$lambda.1_05, mod_p0_IV$lambda.1_1, mod_p0_IV$lambda.1_15, mod_p0_IV$lambda.1_2, mod_p0_IV$lambda.1_25)]
        mod_p1_IV$obs_pred_pert0_lamvar = mod_p1_IV$obs_pred_pert0[,c(mod_p1_IV$lambda.min, mod_p1_IV$lambda.1_05, mod_p1_IV$lambda.1_1, mod_p1_IV$lambda.1_15, mod_p1_IV$lambda.1_2, mod_p1_IV$lambda.1_25)]
        mod_p1_IV$obs_pred_pert1_lamvar = mod_p1_IV$obs_pred_pert1[,c(mod_p1_IV$lambda.min, mod_p1_IV$lambda.1_05, mod_p1_IV$lambda.1_1, mod_p1_IV$lambda.1_15, mod_p1_IV$lambda.1_2, mod_p1_IV$lambda.1_25)]
        
        GSAT.psl[[mon]][[date.ix]] = list()
        GSAT.psl[[mon]][[date.ix]]$mod_p0 = mod_p0; GSAT.psl[[mon]][[date.ix]]$mod_p1 = mod_p1; GSAT.psl[[mon]][[date.ix]]$mod_p0_IV = mod_p0_IV; GSAT.psl[[mon]][[date.ix]]$mod_p1_IV = mod_p1_IV; 
        
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p0$beta = mod_p0$beta[,c(mod_p0$lambda.min, mod_p0$lambda.1_05, mod_p0$lambda.1_1, mod_p0$lambda.1_15, mod_p0$lambda.1_2, mod_p0$lambda.1_25)]
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p1$beta = mod_p1$beta[,c(mod_p1$lambda.min, mod_p1$lambda.1_05, mod_p1$lambda.1_1, mod_p1$lambda.1_15, mod_p1$lambda.1_2, mod_p1$lambda.1_25)]
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p0_IV$beta = mod_p0_IV$beta[,c(mod_p0_IV$lambda.min, mod_p0_IV$lambda.1_05, mod_p0_IV$lambda.1_1, mod_p0_IV$lambda.1_15, mod_p0_IV$lambda.1_2, mod_p0_IV$lambda.1_25)]
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$mod_p1_IV$beta = mod_p1_IV$beta[,c(mod_p1_IV$lambda.min, mod_p1_IV$lambda.1_05, mod_p1_IV$lambda.1_1, mod_p1_IV$lambda.1_15, mod_p1_IV$lambda.1_2, mod_p1_IV$lambda.1_25)]
        
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$grid.ix = CMIP6.df$grid.ix
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$train.ix = CMIP6.df$train.ix
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$train.ix1 = CMIP6.df$train.ix1
        CMIP6.psl.GSAT.beta[[mon]][[date.ix]]$cur.land.weights = CMIP6.df$cur.land.weights
        
        # reconstruct Yhat for full data.frame:
        # only select those dates/times that are from the mon/date.ix:
        cur.ix = which(CMIP6.psl.df$mon[[mon]]$M$year == format(HadSLP2_calendar_[date.ix], "%Y") & CMIP6.psl.df$mon[[mon]]$M$mon == as.numeric(format(HadSLP2_calendar_[date.ix], "%m")))
        CMIP6.psl.GSAT.Yhat.df$mon[[mon]]$mod_p0[cur.ix,] = cbind(mod_p0$Yhat$pt0, mod_p0$Yhat$pt1)[cur.ix,]
        CMIP6.psl.GSAT.Yhat.df$mon[[mon]]$mod_p1[cur.ix,] = cbind(mod_p1$Yhat$pt0, mod_p1$Yhat$pt1)[cur.ix,]
      }
    }
  }
  
  
  # 2. Process into time series of annual and monthly observations:
  # ---------------------------------------------------------------
  
    # 2.1 GSAT:
    # ---------------------------------------------
    {
    # Process monthly / annual observational predictions:
    GSAT.psl_ = list()
    GSAT.psl_$mon = list(); GSAT.psl_$ann = list();
    
    GSAT.psl_$mon$mod_p0_min = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,1]))
    GSAT.psl_$mon$mod_p0_1_05 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,2]))
    GSAT.psl_$mon$mod_p0_1_1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,3]))
    GSAT.psl_$mon$mod_p0_1_15 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,4]))
    GSAT.psl_$mon$mod_p0_1_2 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,5]))
    GSAT.psl_$mon$mod_p0_1_25 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$obs_pred_pert0_lamvar[,6]))
    GSAT.psl_$mon$mod_p1_min = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,1]))
    GSAT.psl_$mon$mod_p1_1_05 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,2]))
    GSAT.psl_$mon$mod_p1_1_1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,3]))
    GSAT.psl_$mon$mod_p1_1_15 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,4]))
    GSAT.psl_$mon$mod_p1_1_2 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,5]))
    GSAT.psl_$mon$mod_p1_1_25 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$obs_pred_pert1_lamvar[,6]))
    
    GSAT.psl_$mon$mod_p0_IV_min = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,1]))
    GSAT.psl_$mon$mod_p0_IV_1_05 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,2]))
    GSAT.psl_$mon$mod_p0_IV_1_1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,3]))
    GSAT.psl_$mon$mod_p0_IV_1_15 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,4]))
    GSAT.psl_$mon$mod_p0_IV_1_2 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,5]))
    GSAT.psl_$mon$mod_p0_IV_1_25 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0_IV$obs_pred_pert0_lamvar[,6]))
    GSAT.psl_$mon$mod_p1_IV_min = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,1]))
    GSAT.psl_$mon$mod_p1_IV_1_05 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,2]))
    GSAT.psl_$mon$mod_p1_IV_1_1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,3]))
    GSAT.psl_$mon$mod_p1_IV_1_15 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,4]))
    GSAT.psl_$mon$mod_p1_IV_1_2 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,5]))
    GSAT.psl_$mon$mod_p1_IV_1_25 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1_IV$obs_pred_pert1_lamvar[,6]))
    
    GSAT.psl_$ann$mod_p0_min = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_min, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_1_05 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_1_05, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_1_1 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_1_1, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_1_15 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_1_15, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_1_2 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_1_2, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_1_25 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_1_25, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_min = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_min, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_1_05 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_1_05, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_1_1 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_1_1, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_1_15 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_1_15, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_1_2 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_1_2, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_1_25 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_1_25, FUN=function(x) x[1:200,i])))
    
    GSAT.psl_$ann$mod_p0_IV_min = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_min, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_IV_1_05 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_1_05, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_IV_1_1 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_1_1, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_IV_1_15 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_1_15, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_IV_1_2 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_1_2, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p0_IV_1_25 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p0_IV_1_25, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_min = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_min, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_1_05 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_1_05, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_1_1 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_1_1, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_1_15 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_1_15, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_1_2 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_1_2, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$mod_p1_IV_1_25 = sapply(X = 1:155, FUN=function(i) rowMeans(sapply(X = GSAT.psl_$mon$mod_p1_IV_1_25, FUN=function(x) x[1:200,i])))
    GSAT.psl_$ann$HadSLP2_calendar_ = HadSLP2_calendar_
    
    # Calculate / Process MSE:
    CMIP6.psl.GSAT.Yhat.df$mon_MSE = list(); CMIP6.psl.GSAT.Yhat.df$ann_MSE = list(); CMIP6.psl.GSAT.Yhat.df$ann_RMSE = list(); CMIP6.psl.GSAT.Yhat.df$mon_ME_sq = list();
    CMIP6.psl.GSAT.Yhat.df$mon_MSE$mod_p0_pt0 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$MSE$pt0[y$mod_p0$lambda.1_05]))
    CMIP6.psl.GSAT.Yhat.df$mon_MSE$mod_p0_pt1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p0$MSE$pt1[y$mod_p0$lambda.1_05]))
    CMIP6.psl.GSAT.Yhat.df$mon_MSE$mod_p1_pt0 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$MSE$pt0[y$mod_p1$lambda.1_05]))
    CMIP6.psl.GSAT.Yhat.df$mon_MSE$mod_p1_pt1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) y$mod_p1$MSE$pt1[y$mod_p1$lambda.1_05]))

    CMIP6.psl.GSAT.Yhat.df$mon_ME_sq$mod_p0_pt0 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p0$me.by.ensreal$pt0[,y$mod_p0$lambda.1_05]^2)))
    CMIP6.psl.GSAT.Yhat.df$mon_ME_sq$mod_p0_pt1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p0$me.by.ensreal$pt1[,y$mod_p0$lambda.1_05]^2)))
    CMIP6.psl.GSAT.Yhat.df$mon_ME_sq$mod_p1_pt0 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p1$me.by.ensreal$pt0[,y$mod_p0$lambda.1_05]^2)))
    CMIP6.psl.GSAT.Yhat.df$mon_ME_sq$mod_p1_pt1 = lapply(X = GSAT.psl, FUN=function(x) sapply(X = x, FUN=function(y) mean(y$mod_p1$me.by.ensreal$pt1[,y$mod_p0$lambda.1_05]^2)))

    # Process CMIP6 predictions into annual files:
    CMIP6.psl.GSAT.Yhat.df$ann$mod_p0 = data.frame(matrix(data = rowMeans(sapply(X = CMIP6.psl.GSAT.Yhat.df$mon, FUN=function(x) c(as.matrix(x$mod_p0)))), ncol = 4))
    CMIP6.psl.GSAT.Yhat.df$ann$mod_p1 = data.frame(matrix(data = rowMeans(sapply(X = CMIP6.psl.GSAT.Yhat.df$mon, FUN=function(x) c(as.matrix(x$mod_p1)))), ncol = 4))
    names(CMIP6.psl.GSAT.Yhat.df$ann$mod_p0) <- names(CMIP6.psl.GSAT.Yhat.df$ann$mod_p1) <- names(CMIP6.psl.GSAT.Yhat.df$mon[[1]]$mod_p0)

  
    # for (i in 1:length(names(CMIP6.psl.GSAT.Yhat.df$mon_MSE))) {
    #  CMIP6.psl.GSAT.Yhat.df$ann_MSE[[i]] = rowMeans(sapply(X = CMIP6.psl.GSAT.Yhat.df$mon_MSE[[i]], FUN=function(x) x))
    #  CMIP6.psl.GSAT.Yhat.df$ann_RMSE[[i]] = rowMeans(sqrt(sapply(X = CMIP6.psl.GSAT.Yhat.df$mon_MSE[[i]], FUN=function(x) x)))
    # }
    # names(CMIP6.psl.GSAT.Yhat.df$ann_MSE) = names(CMIP6.psl.GSAT.Yhat.df$mon_MSE)
    # names(CMIP6.psl.GSAT.Yhat.df$ann_RMSE) = names(CMIP6.psl.GSAT.Yhat.df$mon_MSE)
  }
  
    
  
  # Save observational reconstruction:
  save(list = c("GSAT.psl_"), 
       file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/GSAT.HadSLP2_psl_v2.RData")

  save(list = c("CMIP6.psl.df", "CMIP6.psl.GSAT.Yhat.df", "CMIP6.psl.GSAT.beta"), 
       file = "/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.GSAT.HadSLP2_psl_4eval_v2.RData")
}




