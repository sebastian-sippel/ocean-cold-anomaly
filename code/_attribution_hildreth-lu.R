
## ----------------------------------------------------------------------
## Attribution functions for 1D time series with different noise models:
## ----------------------------------------------------------------------

# Sebastian Sippel
# 26.04.2022
library(HoRM)


## Set up attribution function using AR(1) and FD covariance matrix with Hildreth-Lu method
# Following Jara Imbers et al., Journal of Climate
# https://journals.ametsoc.org/view/journals/clim/27/10/jcli-d-12-00622.1.xml
# https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1002/jgrd.50296


# X = cbind(rnorm(n = 1000), rnorm(n = 1000))
# beta = c(0.9, 0.5)
# epsilon = arima.sim(model = list(ar = 0.7), n = 1000, n.start = 10000)
# Y = X %*% beta + epsilon
# acf(lm(Y ~ X)$residuals)
# rho = lm(lm(Y ~ X)$residuals[2:100] ~ lm(Y ~ X)$residuals[1:99] + 0)$coefficients
# rho.seq = seq(0, 0.9, 0.02)
# HL_mod = lapply(X = rho.seq, FUN=function(rho) hildreth_lu_AR1(Y = Y, X = X, rho = rho, transform = "Young"))
# plot(rho.seq, sapply(X = HL_mod, FUN=function(x) x$mod.param$RSS))
# rho.seq[which.min(sapply(X = HL_mod, FUN=function(x) x$mod.param$RSS))]
# plot(seq(0, 0.9, 0.1), sapply(X = seq(0, 0.9, 0.1), FUN=function(rho) hildreth_lu_AR1(Y = Y, X = X, rho = rho, transform = "Young")$X_2))
  

### continue here with implementation (!!!)
hildreth_lu_AR1 <- function(Y, X, rho, transform = "Young", invert = T) {
  
  # https://www.taylorfrancis.com/books/mono/10.1201/9781315154701/handbook-regression-methods-derek-scott-young
  n = length(Y)
  
  # get X with constant predictor column:
  X = t(t(X))
  X__ = cbind(rep(1, n), X)
  
  na.ix = which(is.na(cbind(X__, Y)), arr.ind = T)
  if (length(na.ix) > 0) {
    nna.ix = c(1:165)[-na.ix[,1]]
    Y = Y[-na.ix[,1]]  
    X = t(t(X[-na.ix[,1],]))
    X__ = X__[-na.ix[,1],]   
    n = length(Y)
  }
  
  # construct covariance matrix for AR(1) model:
  delta <- rep(1, n)
  phi <- rho
  sigma <- 1
  t <- cumsum(delta)
  
  Sigma <- sigma^2/(1-phi^2)*phi^abs(outer(t,t,"-"))
  
  if (invert == T) {
    Sigma_inv = solve(Sigma)
    # image(Sigma)
    # image(Sigma_inv)
    A = chol(Sigma) 
  }
  
  # transformation following Young:
  Y_ = Y[2:n] - rho * (Y[1:(n-1)])
  X_ = apply(X = X, MARGIN = 2, FUN = function(x) x[2:n] - rho * (x[1:(n-1)]))
  
  if (transform == "Young") {
    
    mod = lm(Y_ ~ X_)
    # X_ = X[2:n] - rho * (X[1:(n-1)])
    # get coefficients:
    beta_hat = mod$coefficients[-1]
    # get intercept through transformation following Young 2017 Handbook of Regression Methods:
    beta0_hat = mean(Y) - colMeans(X) %*% beta_hat
    
  } else if (transform == "matrix.projection") {
    Y_ = solve(t(A)) %*% Y
    X_ = solve(t(A)) %*% X
    
    mod = lm(Y_[2:(n)] ~ X_[2:(n),])
    # sum(lm(Y_[1:(n)] ~ X_[1:(n),])$residuals^2)
    
    # get coefficients:
    beta_hat = mod$coefficients[-1]
    # get intercept through transformation following Young 2017 Handbook of Regression Methods:
    beta0_hat = mean(Y) - colMeans(X) %*% beta_hat
    
  } else if (transform == "gls.estimate") {
    
    beta_hat = c(solve((t(X__) %*% Sigma_inv %*% X__)) %*% t(X__) %*% Sigma_inv %*% Y)[-1]
    # get intercept:
    # beta0_hat = c(solve((t(X__) %*% Sigma_inv %*% X__)) %*% t(X__) %*% Sigma_inv %*% Y)[1]
    beta0_hat = mean(Y) - colMeans(X) %*% beta_hat
  } else if (transform == "HoRM") {
    cur.mod = hildreth.lu(y = Y, x = X, rho = rho)
    beta_hat = cur.mod$coefficients[2]
    beta0_hat = cur.mod$coefficients[1]
  }
  
  
  # Return fitted values and reisudals:
  fitted = (X %*% beta_hat + rep(beta0_hat, n))
  res = Y - fitted
  
    # calculate RSS for transformed model for model selection according to Hildreth-Lu:
    mod_TF = lm(Y_ ~ X_)
    RSS = sum(mod_TF$residuals^2)
    ## get variance of the innovations following AR(1) model:
    sigma_sq = var(res[2:n] - rho * res[1:(n-1)])
    # get variance around beta_hat following formula in Imbers et al., 2014:
    # beta_hat_var = diag(solve(t(X__) %*% solve(sigma_sq * Sigma) %*% X__))
    
    if (invert == T) {
      beta_hat_var = diag(solve(t(X__) %*% ( (1 / sigma_sq) * Sigma_inv) %*% X__))  # instead of beta_hat_var could also use std error of 'mod'.    
    } else {
      beta_hat_var = NA
    }
    
    
    # test = (((solve(sigma_sq * Sigma)) -     (1 / sigma_sq * Sigma_inv)))
    # min(test)

    # Return model parameters:
    ret.list = list()
    ret.list$mod.param = data.frame(t(c(beta_hat, beta0_hat, RSS, sigma_sq, beta_hat_var, rho)))
    names(ret.list$mod.param) = c(names(coefficients(mod_TF)[-1]), "beta0_hat", "RSS", "sigma_sq", paste(names(coefficients(mod_TF)), "_var", sep=""), "rho")
    
    ret.list$fitted = c(fitted)
    ret.list$residuals = c(res)
  
    if (length(na.ix) > 0) {
      n0 = length(nna.ix) + length(na.ix[,1])
      res = rep(NA, n0); res[nna.ix] = ret.list$residuals; ret.list$residuals = res
      fitted = rep(NA, n0); fitted[nna.ix] = ret.list$fitted; ret.list$fitted = fitted
    }
    return(ret.list)
}



# run Hildreth-Lu regression in sequence:
hildreth_lu_AR1_ <- function(Y, X, rho.seq = seq(0, 0.95, 0.05), transform = "Young", invert = T) {
  
  test.mod = lapply(X = rho.seq, FUN=function(rho) hildreth_lu_AR1(Y, X, rho, transform = transform, invert = invert))
  # plot(c(unlist(lapply(X = test.mod, FUN = function(x) x$mod.param$RSS))))
  ix=which.min(c(unlist(lapply(X = test.mod, FUN = function(x) x$mod.param$RSS))))
  
  return(test.mod[[ix]])
}






## Run band-pass filter:
run.pass.filt <- function(x, years = 1850:2020, W = 20, type = "low", center = T, ens.ix = NULL) {
  if (is.null(ens.ix)) {
    test  = pass.filt(y = x, W = W, type = type, method = "Butterworth")
    if (center == T) test = test - mean(test[match(x = 1961:1990, table = years)])
    return(test)
  } else {
    test = t(sapply(X = ens.ix, FUN = function(i) pass.filt(y = x[i,], W = W, type = type, method = "Butterworth")))
    if (center == T) test = test - mean(test[,match(x = 1961:1990, table = years)])
    test_q = data.frame(colQuantiles(x = test, probs = c(0.025, 0.5, 0.975)))
    
    names(test_q) = c(paste(substr(type, 1, 1), "p_2.5", sep=""), paste(substr(type, 1, 1), "p_50", sep=""), paste(substr(type, 1, 1), "p_97.5", sep=""))
    return(test_q)
  }
}



# Y = OBS.tas_land$GSAT$ann$mod_p1_min
# f = rowMeans(CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist)
# years = 1850:2020
# ens.ix = NULL

get.df <- function(Y, f, years, center = T, ens.ix = NULL, rho.seq = seq(0, 0.9, 0.05), years.DA = 1850:2014) {
  
  ret.df = data.frame(matrix(NA, nrow = length(years), ncol = 21))
  names(ret.df) = c("Year", "mod_p1_min_2.5", "mod_p1_min_50", "mod_p1_min_97.5", "lp_2.5", "lp_50", "lp_97.5", "hp_2.5", "hp_50", "hp_97.5",
                           "forced", "residual", 
                    "res_2.5", "res_50", "res_97.5", "res_lp_2.5", "res_lp_50", "res_lp_97.5", "res_hp_2.5", "res_hp_50", "res_hp_97.5")
  ret.df$Year = years
  if (is.null(ens.ix)) {
    ret.df$mod_p1_min_50 = Y
    ret.df$lp_50 = run.pass.filt(x = Y, years = years, W = 20, type = "low", center = T, ens.ix = ens.ix)
    ret.df$hp_50 = run.pass.filt(x = Y, years = years, W = 20, type = "high", center = T, ens.ix = ens.ix)
    Y_ATT = Y
  } else {
    ret.df[2:4] = colQuantiles(Y[ens.ix,], probs = c(0.025, 0.5, 0.975))
    ret.df[5:7] = cbind(run.pass.filt(x = Y, years = years, W = 20, type = "low", center = T, ens.ix = ens.ix))
    ret.df[8:10] = cbind(run.pass.filt(x = Y, years = years, W = 20, type = "high", center = T, ens.ix = ens.ix))
    Y_ATT = colMedians(Y[ens.ix,])
  }

  # Run Attribution:
    year.ix = na.omit(match(x =  years.DA, table = years))
    f.year.ix = na.omit(match(x =  years, table = years.DA))
    if (is.null(dim(f))) {
      test = hildreth_lu_AR1_(Y = Y_ATT[year.ix], X = f[f.year.ix], rho.seq = rho.seq, transform = "Young")
    } else {
      test = hildreth_lu_AR1_(Y = Y_ATT[year.ix], X = f[f.year.ix,], rho.seq = rho.seq, transform = "Young")
    }
    print(test$mod.param)
  if (center == T) {
    test$fitted = test$fitted - mean(test$fitted[match(x = 1961:1990, table = c(years.DA)[f.year.ix])])
    test$residuals = test$residuals - mean(test$residuals[match(x = 1961:1990, table = c(years.DA)[f.year.ix])])
  }
    ret.df$forced[year.ix] = test$fitted
    ret.df$residual[year.ix] = test$residuals
  # Subtract from each OBS ensemble member:
    if (is.null(ens.ix)) {
      Y_res = Y[year.ix] - test$fitted
      ret.df$res_50[year.ix] = Y_res
      ret.df$res_lp_50[year.ix] = run.pass.filt(x = Y_res, years = years, W = 20, type = "low", center = center, ens.ix = ens.ix)
      ret.df$res_hp_50[year.ix] = run.pass.filt(x = Y_res, years = years, W = 20, type = "high", center = center, ens.ix = ens.ix)
    } else {
      Y_res = Y[,year.ix] - rep.row(x = test$fitted, n = dim(Y)[1])
      ret.df[year.ix, 13:15] = colQuantiles(Y_res[ens.ix,], probs = c(0.025, 0.5, 0.975))
      ret.df[year.ix, 16:18] = cbind(run.pass.filt(x = Y_res, years = years, W = 20, type = "low", center = center, ens.ix = ens.ix))
      ret.df[year.ix, 19:21] = cbind(run.pass.filt(x = Y_res, years = years, W = 20, type = "high", center = center, ens.ix = ens.ix))
    }
    
  return(ret.df)
}






