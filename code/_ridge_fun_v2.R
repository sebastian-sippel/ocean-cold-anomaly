
# ---------------------------------------------------
# My own functions for ridge regression:
# ---------------------------------------------------

# Sebastian Sippel
# 08.12.2019
require(matrixStats)
require(weights)


## Function to repeat vector either by row or column:
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


# Error metrics / weighted root mean squared error:
# https://stats.stackexchange.com/questions/230517/weighted-root-mean-square-error
weighted.rmse <- function(actual, predicted, weight){
  sqrt(sum((predicted-actual)^2*weight)/sum(weight))
}





# Implementation of weighted ridge regression via SVD
# Sebastian Sippel, 30.11.2021
# 
# following p.9 'Lecture notes on ridge regression': https://arxiv.org/pdf/1509.09169.pdf.
# and uses MASS::lm.ridge code and https://gist.github.com/MGallow/6f2763e9a36bb200a3cb4b264d144dae
# Note: Weights are scaled internally to sum to n. This is identical to glmnet: 
# 'weights is for the observation weights. Default is 1 for each observation. 
# (Note: glmnet rescales the weights to sum to N, the sample size.)'
weighted_ridge = function(X, y, lam = 0, weights = NULL, svd.list = NULL, ret.svd = F) {
  
  # checks
  n = dim(X)[1]
  p = dim(X)[2]
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  if (length(weights) != n)
    stop("weights must be length ", n)
  # if (length(lam) > 1)
  #  stop("lam must be a scalar!")
  if (any(lam < 0))
    stop("lam must be nonnegative!")
  if (!is.null(svd.list)) {
    svd = svd.list
  }
  
  # initialization
  X. = as.matrix(X)
  y. = as.matrix(y)
  
  # rescale weights:
  weights = weights / sum(weights) * n
  
  # No intercept calculation is being performed (removed from original function)
  #
  
  root = sqrt(weights)
  X. = root * X.
  y. = root * y.
  
  
  # SVD
  if (is.null(svd.list)) {
    svd = svd(X.)
  }
  if(ret.svd == T) return(svd)
  rhs <- t(svd$u) %*% y.
  # d <- svd$d
  
  # adjust d vector for regularization and diagonalize for the entire lambda sequence 'lam':
  k <- length(lam)
  dx <- length(svd$d)
  
  div = svd$d^2 + rep(lam, rep(dx, k))
  a <- drop(svd$d * rhs)/div
  dim(a) <- c(dx, k)
  
  betas <- svd$v %*% a
  
  # calculate beta estimates
  # old implementation for length(lam)==0
  # d_adj = diag(svd$d/(svd$d^2 + lam))
  # betas = svd$v %*% d_adj %*% t(svd$u) %*% y.
  # rownames(betas) = NULL
  
  
  # add intercept, if needed
  # not implemented...
  
  returns = list(coef = betas)
  return(returns)
}




## Fit (fast) ridge regression for sequence of \lambda values:
# Sebastian Sippel, 30.11.2021
# This is a convenient function that performs (weighted) standardization before ridge regression.
# This is an adjustment to the previous anchor regression function, but without \gamma
# standardization is the default
# (Much) faster than glmnet with the openblas libraries on co2, etc. Some benchmarks for n>25000 and p~1400 in benchmark script.
fit.ridge <- function(y, X, lambda = 100,  standardize = T, weights = NULL, svd.list = NULL, ret.svd = F) {
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  if (is.null(weights)) {
    weights = rep(1, n)
  }
  # rescale weights:
  weights = weights / sum(weights) * n
  
  if (standardize == T) {
    # Standardize x and y before ridge regression:
    mu_x = colWeightedMeans(x = X, w = weights)
    sd_x = colWeightedSds(x = X, w = weights) * sqrt((n-1) / n)
    mu_y = weighted.mean(x = y, w = weights); 
    sd_y = weightedSd(x = y, w = weights) * sqrt((n-1) / n)
    X_sc = (X - rep.row(mu_x, n = n)) / rep.row(sd_x, n = n)
    Y_sc = (y - mu_y) / sd_y
    
    # Fit (fast) weighted ridge regression for standardized data:
    if (ret.svd == T) return(weighted_ridge(X = X_sc, y = Y_sc, lam = lambda, weights = weights, ret.svd = T))
    beta = weighted_ridge(X = X_sc, y = Y_sc, lam = lambda, weights = weights, svd.list = svd.list)$coef
    
    # convert to un-standardized coefficients:
    # https://stats.stackexchange.com/questions/155362/glmnet-unstandardizing-linear-regression-coefficients
    a0 = mu_y - colSums(beta * mu_x / sd_x) * sd_y
    beta0 = beta / sd_x * sd_y
    return(list(a0 = a0, beta = beta0, lambda = lambda, beta_std = beta))
  } else if (standardize == F) {
    # not implemented currently!!
    return(NULL)
  }
}





## Fit ridge regression with cross-validation.
# Function implements cross-validation and calls fit_ridge or glmnet  .
# Calculations in parallel over number of simulations (with cv = "model.subagging") or number of models (with cv = "leave.model.out")
# Function calculates several error metrics.
cv.glmnet2 <- function(x, y, alpha, lambda, foldid, 
                       penalty.factor, weights,
                       nr.cores = NULL, nr.subsample = 3000, cv = "leave.model.out", nsim = 20, 
                       keep = F, adj.mean.by.mod = F, include.all.mod.fit = F, ret.beta.only = T,
                       fit.fun = "weighted.ridge", svd.list = NULL, ret.svd = F) {
  
  # Purpose of function: 
  # Similar to cv.glmnet, but including options for different \gamma values
  ## Standardization of regression coefficients is based on: https://stats.stackexchange.com/questions/155362/glmnet-unstandardizing-linear-regression-coefficients
  
  # 1. Prepare regression/cross-validation:
  foldid.un = na.omit(unique(foldid))
  
  # Different cross-validation strategy to use:
  if (cv == "model.subagging") {
    nuse=0.5
    design=sapply(X = 1:nsim, FUN=function(ix) {
      set.seed(ix+6) 
      sample(x = 1:length(foldid.un), size = ceiling(length(foldid.un)*nuse), replace = F) })
  } else if (cv == "leave.model.out") {
    nr.subsample = NULL
    nsim = length(foldid.un)
    design = sapply(X = 1:nsim, FUN=function(ix) c(1:nsim)[-ix])
  } else if (cv == "model.by.model") {     # -> need to set nr.subsample == NULL
    # nsim = length(um)
    # design <- matrix(data = c(1:nsim), nrow = 1, ncol = nsim)
  }
  
  require(doParallel)
  registerDoParallel(cores = nr.cores)
  # ptm <- proc.time()
  glmnet2.list = foreach(sim=1:nsim) %dopar% {
    
    # print(paste("\r ***", sim))
    train <- numeric(0)
    for (ucc in 1:length(foldid.un)){ if( ucc %in% design[,sim]) train <- c(train, foldid.un[ucc])}
    test <- foldid.un[-which(foldid.un %in% train)]
    
    itrainX <- which(foldid %in% train)
    itestX <- which(!(foldid %in% train))
    
    ## Subsample training indices and standardize data:
    # Subsample itrainX to a train ix of equal weight:
    set.seed(sim)
    
    if (is.null(nr.subsample)) {
      train.ix = itrainX
    } else {
      train.ix = c(sapply(X = design[,sim], FUN=function(cc) sort(sample(x = which(cc == foldid), size = nr.subsample, replace=T))))
    }
    cv.glmnet2.out = list()
    # start_time <- Sys.time()
    # cv.glmnet2.out$anchor.fit = fit.anchor_t(x = x[train.ix,], y = y[train.ix], A = A[train.ix], lambda = lambda, gamma = gamma, standardize = T)
    # end_time <- Sys.time()
    # end_time - start_time
    if (fit.fun == "glmnet") {
      cv.glmnet2.out$fit = glmnet(x = x[train.ix,], y = y[train.ix], family = "gaussian", weights = weights[train.ix], 
                                  alpha = alpha, lambda = lambda, penalty.factor = penalty.factor, standardize = T)
    } else if (fit.fun == "weighted.ridge") {
      lambda_new = lambda * length(train.ix) / weightedSd(x = y[train.ix], w = weights[train.ix] / sum(weights[train.ix]) * length(train.ix))
      if(ret.svd == T) return(fit.ridge(X = x[train.ix,], y = y[train.ix], lambda = lambda_new, weights = weights[train.ix] / sum(weights[train.ix]) * length(train.ix), standardize = T, ret.svd = T))
      
      if (is.null(svd.list)) {
        cv.glmnet2.out$fit = fit.ridge(X = x[train.ix,], y = y[train.ix], lambda = lambda_new, weights = weights[train.ix] / sum(weights[train.ix]) * length(train.ix), standardize = T)
      } else {
        cv.glmnet2.out$fit = fit.ridge(X = x[train.ix,], y = y[train.ix], lambda = lambda_new, weights = weights[train.ix] / sum(weights[train.ix]) * length(train.ix), standardize = T, svd.list = svd.list[[sim]])
      }
    }
    # start_time <- Sys.time()
    # end_time <- Sys.time()
    # end_time - start_time
    # plot(cv.glmnet2.out$weighted_ridge$beta[,50], cv.glmnet2.out$glmnet.fit$beta[,50])
    # cor(cv.glmnet2.out$weighted_ridge$beta[,50], cv.glmnet2.out$glmnet.fit$beta[,50])
    cv.glmnet2.out$Yhat = matrix(data = NA, nrow = length(y), ncol = length(lambda))
    if (adj.mean.by.mod == F) {
      cv.glmnet2.out$Yhat[itestX,] = as.matrix(x[itestX,] %*% cv.glmnet2.out$fit$beta + rep.row(cv.glmnet2.out$fit$a0, n = length(itestX)))
    } else if (adj.mean.by.mod == T) {
      # pred = x[itestX,] %*% cv.glmnet2.out$glmnet.fit$beta # + rep.row(cv.anchor.out$anchor.fit$a0, n = length(itestX))
      # cv.anchor.out$Yhat[itestX,] = pred - rep.row(colMeans(pred), length(itestX)) + mean(y[itestX])
    }
    return(cv.glmnet2.out)
  }
  # proc.time() - ptm
  if(ret.svd==T) return(glmnet2.list)
  
  ## SUMMARY STATISTICS:
  # Weighted RMSE + Weighted correlation
  Yhat = sapply(X = 1:length(lambda), FUN=function(lambda.ix) rowMeans(sapply(X = glmnet2.list, FUN=function(x) x$Yhat[,lambda.ix]), na.rm=T))
  Yhat.sd = sapply(X = 1:length(lambda), FUN=function(lambda.ix) rowSds(sapply(X = glmnet2.list, FUN=function(x) x$Yhat[,lambda.ix]), na.rm=T))
  beta = sapply(X = 1:length(lambda), FUN=function(lambda.ix) rowMeans(sapply(X = glmnet2.list, FUN=function(x) x$fit$beta[,lambda.ix]), na.rm=T))
  # beta_std = sapply(X = 1:length(lambda), FUN=function(lambda.ix) rowMeans(sapply(X = glmnet2.list, FUN=function(x) x$fit$beta_std[,lambda.ix]), na.rm=T))
  a0 = sapply(X = 1:length(lambda), FUN=function(lambda.ix) mean(sapply(X = glmnet2.list, FUN=function(x) x$fit$a0[lambda.ix]), na.rm=T))
  
  # Save prediction from each simulation for later lambda selection:
  allsim = list()
  allsim$design = design
  allsim$beta = lapply(X = glmnet2.list, FUN=function(x) x$fit$beta)
  allsim$a0 = lapply(X = glmnet2.list, FUN=function(x) x$fit$a0)
  if (keep == T) allsim$Yhat = lapply(X = glmnet2.list, FUN=function(x) x$Yhat)
  
  w=rep(NA, length(y)); for (cc in 1:length(foldid.un)) w[which(cc == foldid)] = 1/length(which(cc == foldid)) / length(foldid.un)
  MSE = sapply(X = 1:length(lambda), FUN=function(lambda.ix) weighted.rmse(actual = y, predicted = Yhat[,lambda.ix], weight = w)^2)
  cur.resid = Yhat - rep.col(y, n = length(lambda))
  # res.cor = sapply(X = 1:length(lambda), FUN=function(lambda.ix) wtd.cor(x = cur.resid[,lambda.ix], y = A, weight = w)[1])
  
  ## MSE and residual correlation by model:
  MSE.by.mod = matrix(data = NA, nrow = length(foldid.un), ncol = length(lambda))
  MSE.by.mod[foldid.un,] = t(sapply(foldid.un, FUN=function(cc) {
    c.ix = which(cc == foldid)
    sapply(X = 1:length(lambda), FUN=function(lambda.ix) mse(sim = Yhat[c.ix,lambda.ix], obs = y[c.ix]))
  }))
  # MSE standard deviation across models:
  ## CONTINUE HERE: STANDARD ERROR OF THE MEAN (i.e. variation of the MEAN given folds...)
  SE.MSE =  colSds(x = MSE.by.mod) / sqrt(length(foldid.un))
  #plot(MSE)
  #lines(MSE + SE.MSE, col="red")  # -> 1SE criterion...
  cor.by.mod = matrix(data = NA, nrow = length(foldid.un), ncol = length(lambda))
  cor.by.mod[foldid.un,] = t(sapply(foldid.un, FUN=function(cc) {
    c.ix = which(cc == foldid)
    sapply(X = 1:length(lambda), FUN=function(lambda.ix) cor(x = Yhat[c.ix,lambda.ix], y = y[c.ix]))
  }))
  
  # get lambda.min and lambda.1.05
  lambda.min = which.min(colMeans(MSE.by.mod))
  lambda.1_05 = which(colMeans(MSE.by.mod) < (min(colMeans(MSE.by.mod)) * 1.05))[1]
  
  
  # Return list:
  ret.list = list()
  ret.list$beta = beta
  ret.list$a0 = a0
  
  ret.list$Yhat = Yhat
  ret.list$Yhat.sd = Yhat.sd
  ret.list$Y = y
  ret.list$foldid = foldid
  ret.list$MSE = MSE
  ret.list$SE.MSE = SE.MSE
  # ret.list$res.cor = res.cor
  ret.list$MSE.by.mod = MSE.by.mod
  ret.list$cor.by.mod = cor.by.mod
  ret.list$allsim = allsim
  ret.list$train.mod = unique(names(foldid))
  
  # include fit on all input data?
  if (include.all.mod.fit == T) {
    ## get glmnet on all data (FITTING ON ALL DATA JOINTLY IS EQUAL TO AVERAGING??):
    glmnet.fit.all = glmnet(x = x, y = y, family = "gaussian", weights = weights, 
                            alpha = alpha, lambda = lambda, penalty.factor = penalty.factor, standardize = T)
    # plot(glmnet.fit.all$beta[,80], beta[,80])
    # cor(glmnet.fit$beta[,80], beta[,80])
    # abline(0, 1, col = "red")
    ret.list$beta_all = as.matrix(glmnet.fit.all$beta)
    ret.list$a0_all = as.matrix(glmnet.fit.all$a0)
  }
  
  if (ret.beta.only == T) {
    ret.list1 = list()
    ret.list1$beta = allsim$beta
    ret.list1$a0 = allsim$a0
    ret.list1$beta_avg = beta
    ret.list1$a0_avg = a0
    ret.list1$lambda.min = lambda.min
    ret.list1$lambda.1_05 = lambda.1_05
    return(ret.list1)
  }
  
  
  return(ret.list)
}






## Generate perturbed train/test data.
gen.pert.traintest.data <- function(X, M, bias.ens, unc, 
                                    randomize.ens = F, fact = 1, ret.ens.mem.split = F, ret.ens.ix = F) {
  
  X_pert <- X
  train.mod.ens = unique(paste(M$mod, M$ens.mem, sep="_"))
  mod = M$mod
  ens.mem = M$ens.mem
  
  # keep first 2 members for testing: 
  if (randomize.ens == T) {
    set.seed(seed = 4)
    test.ens = c(sort(sample(x = 1:200, size = 50, replace = F)))
    set.seed(seed = 5)
    ens.ix = sample(x = c(1:200)[-test.ens], size = length(train.mod.ens), replace = T)
  } else {
    ens.ix = 1:200
  }
  if(ret.ens.ix == T) return(ens.ix)
  
  ens.mem.split = rep(NA, dim(X_pert)[1])
  
  # Start the clock!
  # ptm <- proc.time()
  for (me in 1:length(train.mod.ens)) {
    # print(me)
    cur.mod = strsplit(x = train.mod.ens[me], split = "_")[[1]][1]
    cur.ens.mem = strsplit(x = train.mod.ens[me], split = "_")[[1]][2]
    ix = which(cur.mod == mod & cur.ens.mem == ens.mem)
    
    if(ret.ens.mem.split == T) {  ens.mem.split[ix] = me; next }
    
    if (length(dim(unc)) == 1 | length(dim(unc)) == 0) {  # in case of one-dimensional uncertainties
      X_pert[ix,] = X_pert[ix,] + 
        fact * rep.row(x = bias.ens[[ens.ix[me]]], n = length(ix)) +    # biases from ensemble anomalies
        fact * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = length(ix), mean = 0, sd = unc[i]))   # realize uncorrelated uncertainties...
    } else if (length(dim(unc)) == 2) {
      X_pert[ix,] = X_pert[ix,] + 
        fact * rep.row(x = bias.ens[[ens.ix[me]]], n = length(ix)) +    # biases from ensemble anomalies
        # fact * mvrnorm(n = length(ix), mu = rep(0, dim(unc)[1]), Sigma = unc)
        fact * rmvn(n = length(ix), mu = rep(0, dim(unc)[1]), sigma = unc, ncores = 1)
    }
    
    # X_pert3[ix,] = X_pert3[ix,] + 
    #  3 * rep.row(x = CRUTEM5_ENS_anom[[me]][date.ix,grid.ix], n = length(ix)) +    # biases from ensemble anomalies
    #  3 * sapply(X = grid.ix, FUN=function(i) rnorm(n = length(ix), mean = 0, sd = CRUTEM5_sampling_unc[date.ix,i]))   # realize uncorrelated uncertainties...
  }
  # ptm <- proc.time()
  # proc.time() - ptm
  
  
  if(ret.ens.mem.split == T) {  return(ens.mem.split) }
  
  #ret.list = list()
  #ret.list$X_pert = X_pert
  #ret.list$ens.mem.split = ens.mem.split
  
  return(X_pert)
}




## Generate perturbed train/test data.
gen.pert.traintest.data_noENS <- function(X, M, bias.ens, unc, fact = 1, ret.ens.mem.split = F) {
  
  X_pert <- X
  train.mod.ens = unique(paste(M$mod, M$ens.mem, sep="_"))
  mod = M$mod
  ens.mem = M$ens.mem
  ens.mem.split = rep(NA, dim(X_pert)[1])
  
  # Start the clock!
  # ptm <- proc.time()
  for (me in 1:length(train.mod.ens)) {
    # print(me)
    cur.mod = strsplit(x = train.mod.ens[me], split = "_")[[1]][1]
    cur.ens.mem = strsplit(x = train.mod.ens[me], split = "_")[[1]][2]
    ix = which(cur.mod == mod & cur.ens.mem == ens.mem)
    
    if(ret.ens.mem.split == T) {  ens.mem.split[ix] = me; next }
    
    if (length(dim(unc)) == 1 | length(dim(unc)) == 0 ) {  # in case of one-dimensional uncertainties
      X_pert[ix,] = X_pert[ix,] + 
        # fact * rep.row(x = bias.ens[[ens.ix[me]]], n = length(ix)) +    # biases from ensemble anomalies
        fact * sapply(X = 1:length(unc), FUN=function(i) rnorm(n = length(ix), mean = 0, sd = unc[i]))   # realize uncorrelated uncertainties...
    } else if (length(dim(unc)) == 2) {
      X_pert[ix,] = X_pert[ix,] + 
        # fact * rep.row(x = bias.ens[[ens.ix[me]]], n = length(ix)) +    # biases from ensemble anomalies
        # fact * mvrnorm(n = length(ix), mu = rep(0, dim(unc)[1]), Sigma = unc)
        fact * rmvn(n = length(ix), mu = rep(0, dim(unc)[1]), sigma = unc, ncores = 1)
    }
    
    # X_pert3[ix,] = X_pert3[ix,] + 
    #  3 * rep.row(x = CRUTEM5_ENS_anom[[me]][date.ix,grid.ix], n = length(ix)) +    # biases from ensemble anomalies
    #  3 * sapply(X = grid.ix, FUN=function(i) rnorm(n = length(ix), mean = 0, sd = CRUTEM5_sampling_unc[date.ix,i]))   # realize uncorrelated uncertainties...
  }
  # ptm <- proc.time()
  # proc.time() - ptm
  
  if(ret.ens.mem.split == T) {  return(ens.mem.split) }
  
  #ret.list = list()
  #ret.list$X_pert = X_pert
  #ret.list$ens.mem.split = ens.mem.split
  
  return(X_pert)
}





# Cross-validate on perturbed dataset(s):
validate.pert.dataset <- function(X_pert, Y, X_obs,
                                  mod.glmnet, lambda.min, lambda.ix,
                                  foldid, ens.mem.split, ens.ix = NULL) {
  
  # 1. Prepare regression/cross-validation:
  foldid.un = na.omit(unique(foldid))
  nsim = length(foldid.un)
  n.ens = length(unique(ens.mem.split))
  
  beta = mod.glmnet$beta
  a0 = mod.glmnet$a0
  lambda.min = mod.glmnet$lambda.min
  lambda.1_05 = mod.glmnet$lambda.1_05
  
  # generate list of Yhat:
  Yhat = list()
  MSE.by.mod = list()
  cor.by.mod = list()
  me.by.mod = list()
  MSE.by.ensreal = list()
  cor.by.ensreal = list()
  me.by.ensreal = list()
  
  for (i in 1:length(X_pert)) {
    Yhat[[i]] = matrix(data = NA, nrow = length(Y), ncol = dim(beta[[1]])[2])
    MSE.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
    cor.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
    me.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
    
    MSE.by.ensreal[[i]] = matrix(data = NA, nrow = length(unique(ens.mem.split)), ncol = dim(beta[[1]])[2])
    cor.by.ensreal[[i]] = matrix(data = NA, nrow = length(unique(ens.mem.split)), ncol = dim(beta[[1]])[2])
    me.by.ensreal[[i]] = matrix(data = NA, nrow = length(unique(ens.mem.split)), ncol = dim(beta[[1]])[2])
  }
  
  # 2. Predict Yhat and compute model-specific MSE and other metrics:
  for (sim in 1:nsim) {
    # print(paste("\r ***", sim))
    ix = which(foldid == sim)
    
    for (i in 1:length(X_pert)) {
      # (i) predict out of sample observations:
      Yhat[[i]][ix,] = as.matrix(X_pert[[i]][ix,] %*% beta[[sim]] + rep.row(a0[[sim]], n = length(ix)))
      # (ii) Calculate error metrics per model:
      MSE.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mse(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
      cor.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) cor(Yhat[[i]][ix,lambda.ix], Y[ix]))
      me.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) me(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
    }
  }
  
  # 3. Calculate error metrics for each bias realization ensemble member:
  for (en in 1:n.ens) {
    # print(paste("\r ***", en))
    ix = which(ens.mem.split == en)
    
    for (i in 1:length(X_pert)) {
      MSE.by.ensreal[[i]][en,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mse(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
      cor.by.ensreal[[i]][en,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) cor(Yhat[[i]][ix,lambda.ix], Y[ix]))
      me.by.ensreal[[i]][en,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) me(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
    }
  }
  
  # Generate averages:
  beta_ = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) rowMeans(sapply(X = beta, FUN=function(x) x[,lambda.ix]), na.rm=T))
  a0_ = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mean(sapply(X = a0, FUN=function(x) x[lambda.ix]), na.rm=T))
  MSE = lapply(X = MSE.by.mod, FUN=function(x) colMeans(x))
  
  # Evaluate \lambda.min and \lambda.min+0.05
  pt.lambda.min = lapply(X = MSE, FUN=function(x) which.min(x))
  pt.lambda.1.05 = lapply(X = MSE, FUN=function(x) which(x < (min(x) * 1.05))[1])
  
  # Generate prediction on observations:
  obs_pred = X_obs$X_pert %*% beta_ + rep.row(a0_, n = 200)
  obs_pred_pert2 = X_obs$X_pert2 %*% beta_ + rep.row(a0_, n = 200)
  
  obs_pred_min = obs_pred[,lambda.min]
  obs_pred_1_05 = obs_pred[,lambda.1_05]
  
  names(MSE.by.mod) <- names(cor.by.mod) <- names(me.by.mod) <- names(MSE.by.ensreal) <- names(cor.by.ensreal) <- names(me.by.ensreal) <- names(MSE) <-
    names(pt.lambda.min) <- names(pt.lambda.1.05) <- c("pt0", "pt1", "pt2")
  
  Yhat.ret = lapply(Yhat, FUN = function(x) x[,c(lambda.min, lambda.1_05)])
  names(Yhat.ret) <- c("pt0", "pt1", "pt2")
  
  ret.list = list(MSE.by.mod, cor.by.mod, me.by.mod, MSE.by.ensreal, cor.by.ensreal, me.by.ensreal, 
                  as.matrix(beta_)[,c(lambda.min, lambda.1_05)], a0_[c(lambda.min, lambda.1_05)],
                  MSE, lambda.min, lambda.1_05, pt.lambda.min, pt.lambda.1.05, 
                  obs_pred, obs_pred_pert2, obs_pred_min, obs_pred_1_05, ens.ix,
                  Yhat.ret)
  names(ret.list) = c("MSE.by.mod", "cor.by.mod", "me.by.mod", "MSE.by.ensreal", "cor.by.ensreal", "me.by.ensreal", 
                      "beta", "a0",
                      "MSE", "lambda.min", "lambda.1_05", "pt.lambda.min", "pt.lambda.1_05",
                      "obs_pred", "obs_pred_pert2", "obs_pred_min", "obs_pred_1_05", "ens.ix",
                      "Yhat")
  return(ret.list)
}



# Cross-validate on perturbed dataset(s):
validate.pert.dataset_4psl <- function(X_pert, Y, X_obs,
                                       mod.glmnet, lambda.min, lambda.ix,
                                       foldid) {
  
  # 1. Prepare regression/cross-validation:
  foldid.un = na.omit(unique(foldid))
  nsim = length(foldid.un)
  
  beta = mod.glmnet$beta
  a0 = mod.glmnet$a0
  lambda.min = mod.glmnet$lambda.min
  lambda.1_05 = mod.glmnet$lambda.1_05
  
  # generate list of Yhat:
  Yhat = list()
  MSE.by.mod = list()
  cor.by.mod = list()
  me.by.mod = list()
  MSE.by.ensreal = list()
  cor.by.ensreal = list()
  me.by.ensreal = list()
  
  for (i in 1:length(X_pert)) {
    Yhat[[i]] = matrix(data = NA, nrow = length(Y), ncol = dim(beta[[1]])[2])
    MSE.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
    cor.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
    me.by.mod[[i]] = matrix(data = NA, nrow = length(foldid.un), ncol = dim(beta[[1]])[2])
  }
  
  # 2. Predict Yhat and compute model-specific MSE and other metrics:
  for (sim in 1:nsim) {
    # print(paste("\r ***", sim))
    ix = which(foldid == sim)
    
    for (i in 1:length(X_pert)) {
      # (i) predict out of sample observations:
      Yhat[[i]][ix,] = as.matrix(X_pert[[i]][ix,] %*% beta[[sim]] + rep.row(a0[[sim]], n = length(ix)))
      # (ii) Calculate error metrics per model:
      MSE.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mse(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
      cor.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) cor(Yhat[[i]][ix,lambda.ix], Y[ix]))
      me.by.mod[[i]][sim,] = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) me(sim = Yhat[[i]][ix,lambda.ix], obs = Y[ix]))
    }
  }
  
  # Generate averages:
  beta_ = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) rowMeans(sapply(X = beta, FUN=function(x) x[,lambda.ix]), na.rm=T))
  a0_ = sapply(X = 1:(dim(beta[[1]])[2]), FUN=function(lambda.ix) mean(sapply(X = a0, FUN=function(x) x[lambda.ix]), na.rm=T))
  MSE = lapply(X = MSE.by.mod, FUN=function(x) colMeans(x))
  
  # Evaluate \lambda.min and \lambda.min+0.05
  pt.lambda.min = lapply(X = MSE, FUN=function(x) which.min(x))
  pt.lambda.1.05 = lapply(X = MSE, FUN=function(x) which(x < (min(x) * 1.05))[1])
  
  # Generate prediction on observations:
  obs_pred_pert0 = X_obs$X_pert0 %*% beta_ + rep.row(a0_, n = 200)
  obs_pred_pert1 = X_obs$X_pert1 %*% beta_ + rep.row(a0_, n = 200)
  obs_pred_pert2 = X_obs$X_pert2 %*% beta_ + rep.row(a0_, n = 200)
  
  obs_pred_pert0_min = obs_pred_pert0[,lambda.min]
  obs_pred_pert0_1_05 = obs_pred_pert0[,lambda.1_05]
  
  obs_pred_pert1_min = obs_pred_pert1[,lambda.min]
  obs_pred_pert1_1_05 = obs_pred_pert1[,lambda.1_05]
  
  names(MSE.by.mod) <- names(cor.by.mod) <- names(me.by.mod) <- names(MSE) <-
    names(pt.lambda.min) <- names(pt.lambda.1.05) <- c("pt0", "pt1", "pt2")
  
  Yhat.ret = lapply(Yhat, FUN = function(x) x[,c(lambda.min, lambda.1_05)])
  names(Yhat.ret) <- c("pt0", "pt1", "pt2")
  
  ret.list = list(MSE.by.mod, cor.by.mod, me.by.mod, 
                  as.matrix(beta_), a0_,
                  MSE, lambda.min, lambda.1_05, pt.lambda.min, pt.lambda.1.05, 
                  obs_pred_pert0, obs_pred_pert1, obs_pred_pert2, 
                  obs_pred_pert0_min, obs_pred_pert0_1_05, obs_pred_pert1_min, obs_pred_pert1_1_05,
                  Yhat.ret)
  names(ret.list) = c("MSE.by.mod", "cor.by.mod", "me.by.mod", 
                      "beta", "a0",
                      "MSE", "lambda.min", "lambda.1_05", "pt.lambda.min", "pt.lambda.1_05",
                      "obs_pred_pert0", "obs_pred_pert1", "obs_pred_pert2", 
                      "obs_pred_pert0_min", "obs_pred_pert0_1_05", "obs_pred_pert1_min", "obs_pred_pert1_1_05",
                      "Yhat")
  return(ret.list)
}


