
# Sebastian Sippel
# 11.12.2023

# Code for post-processing of simulations
# ------------------------------------------------------

# Define functions for post-processing of simulations: 
# ------------------------------------------------------

get.diff.window.cor <- function(x, y, years = 1850:2020, w.width = 51, W = 20, center = T, ens.means = NA) {
  
  x.low = pass.filt(y = x, W = W, type = "low", method = "Butterworth")
  x.high = pass.filt(y = x, W = W, type = "high", method = "Butterworth")
  y.low = pass.filt(y = y, W = W, type = "low", method = "Butterworth")
  y.high = pass.filt(y = y, W = W, type = "high", method = "Butterworth")
  if (center == T & !is.numeric(ens.means)) {
    x.low = x.low - mean(x.low[match(x = 1961:1990, table = years)])
    x.high = x.high - mean(x.high[match(x = 1961:1990, table = years)])
    y.low = y.low - mean(y.low[match(x = 1961:1990, table = years)])
    y.high = y.high - mean(y.high[match(x = 1961:1990, table = years)])
  } else if (center == T & is.numeric(ens.means)) {
    x.low = x.low - mean(ens.means$x.low[match(x = 1961:1990, table = years)])
    x.high = x.high - mean(ens.means$x.high[match(x = 1961:1990, table = years)])
    y.low = y.low - mean(ens.means$y.low[match(x = 1961:1990, table = years)])
    y.high = y.high - mean(ens.means$y.high[match(x = 1961:1990, table = years)])
  }
  # Get Differences & Correlations:
  # plot(y); lines(x-y); lines(x.low - y.low)
  ret.df = data.frame(x, x.high, x.low, y, y.high, y.low, 
                      x-y, x.high-y.high, x.low-y.low,  
                      rollapply(data = cbind(x, y), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA),
                      rollapply(data = cbind(x.high, y.high), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA), 
                      rollapply(data = cbind(x.low, y.low), width = w.width, FUN=function(x) { cor(x[,1], x[,2], use="complete.obs") }, by.column = F, fill = NA))
  names(ret.df) = c("x", "x.high", "x.low", "y", "y.high", "y.low", 
                    "diff_x_y", "diff_x.high_y.high", "diff_x.low_y.low", 
                    "cor_x_y", "cor_x.high_y.high", "cor_x.low_y.low")
  return(ret.df)
}

get.diff.window.cor_ens <- function(x, y, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix) {
  
  if (center.ens.means == T) {
    ens.means = get.diff.window.cor(colMedians(x[ens.ix,]), colMedians(y[ens.ix,]), years = years, w.width = w.width, W = W, center = T)
  } else {
    ens.means = NA
  }
  
  dat = list()
  for (en in 1:length(ens.ix)) {
    if (center.ens.means == T) {
      dat[[en]] = get.diff.window.cor(x[en,], y[en,], years = years, w.width = w.width, W = W, center = center, ens.means = ens.means) 
    } else {
      dat[[en]] = get.diff.window.cor(x[en,], y[en,], years = years, w.width = w.width, W = W, center = center, ens.means = ens.means) 
    }
  }
  
  # Get quantile estimates:
  diff_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x_y)), probs = c(0.025, 0.5, 0.975))
  diff_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.high_y.high)), probs = c(0.025, 0.5, 0.975))
  diff_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.low_y.low)), probs = c(0.025, 0.5, 0.975))
  
  cor_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x_y)), probs = c(0.025, 0.5, 0.975))
  cor_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.high_y.high)), probs = c(0.025, 0.5, 0.975))
  cor_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.low_y.low)), probs = c(0.025, 0.5, 0.975))
  
  ret.df = data.frame(diff_x_y, diff_x.high_y.high, diff_x.low_y.low, cor_x_y, cor_x.high_y.high, cor_x.low_y.low)
  names(ret.df) = c("diff_x_y_2.5", "diff_x_y_50", "diff_x_y_97.5", "diff_x.high_y.high_2.5", "diff_x.high_y.high_50", "diff_x.high_y.high_97.5", "diff_x.low_y.low_2.5", "diff_x.low_y.low_50", "diff_x.low_y.low_97.5", 
                    "cor_x_y_2.5", "cor_x_y_50", "cor_x_y_97.5", "cor_x.high_y.high_2.5", "cor_x.high_y.high_50", "cor_x.high_y.high_97.5", "cor_x.low_y.low_2.5", "cor_x.low_y.low_50", "cor_x.low_y.low_97.5")
  return(ret.df)
}

get.trend <- function(x, trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014) {
  
  if (all(is.na(x))) return(rep(NA, length(trend.years)))
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", 1:length(trend.years), sep="")
  return(trends.out)
}



get.trend_perioddiff <- function(x, trend.years = list(1900:1939, 1900:1950, 1980:2014), period.diff.years = list(c(1901:1920), c(1871:1890)), years = 1850:2014) {
  if (all(is.na(x))) return(rep(NA, length(trend.years)+1))
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", 1:length(trend.years), sep="")
  # period.diff:
  late.ix = match(x = period.diff.years[[1]], table = years)
  early.ix = match(x = period.diff.years[[2]], table = years)
  
  period.diff.out = mean(x[late.ix]) - mean(x[early.ix])
  names(period.diff.out) = "period.diff"
  trends.out = c(trends.out, period.diff.out)
  
  return(trends.out)
}




get.period.mean <- function(x, period.years = list(1901:1920), years = 1850:2020) {
  #if (all(is.na(x))) return(rep(NA, length(trend.years)+1))
  
  period.mean.out = sapply(X = 1:length(period.years), FUN=function(i) {
    ix = match(x = period.years[[i]], table = years)
    return( mean(x[ix]) )
  })
  names(period.mean.out) = paste("pm", 1:length(period.years), sep="")
  return(period.mean.out)
}







# get.trend(x = OBS.tos_$GMSST$ann$mod_p1_min[1,], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2020)
# -> further changes to do: Good name for trends, or meta-file...  
get.trend_ <- function(x, trend.length = 50, years = 1850:2014) {
  
  trend.years = lapply(X = years[1]:(tail(years, 1)-trend.length+1), FUN = function(cur.year) cur.year:(cur.year+trend.length-1))
  
  if (all(is.na(x))) return(rep(NA, length(trend.years)))
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", seq(1850, 1850+length(trend.years)-1), sep="")
  return(trends.out)
}







get.RE_score <- function(cur.region ) {
  cur.year = c(cur.region [[5]])
  year.ix = which(cur.year %in% 1870:2000)
  RE = sapply(X = 1:dim(cur.region[[2]])[2], FUN=function(i) {
    sum(cur.region[[2]][year.ix,i], na.rm = T)
  })
  return(RE) # which(RE > 0)
}







## 13. Define linear constraints based on ensemble or individual point estimate:
## ----------------------------------------

# Function to get linear model + constraint based on model relationship:
# y = CMIP6.trends$GMSST3_true
# x = CMIP6.trends$Tropics3_true
# x_new = ocean2k_trends[3]



get.linear.model.constraint <- function(y, x, x_new, plot.constraint = T) {
  dat = data.frame(y = y, x = x)
  new.dat = data.frame(x = x_new)
  reg.mod = lm(y ~ x, data = dat)
  
  n=length(reg.mod$residuals)
  RSS = sum((reg.mod$residuals)^2)
  s2 = 1/(n - 2) * RSS  # MSE
  sigma_x = sqrt( var(dat$x, na.rm=T) )
  new.dat.pred = predict(reg.mod, newdata = new.dat)
  sigma_f_x = sqrt(s2) * sqrt(1 + 1/n + (median(new.dat.pred) - mean(dat$x, na.rm = T)) / (n * sigma_x^2) )
  
  # combine variances in quadrature:
  mean_out =  median(new.dat.pred)
  sd_out =  sqrt(sigma_f_x^2) # + var(new.dat.pred))
  
  
  # plot emergent constraint for an ensemble:
  if(plot.constraint == T) {
    plot(x, y, col="red")
    abline(reg.mod, col ="red")
    
    plotCI(x = median(x_new), y = 0, li = quantile(x_new, 0.025), ui = quantile(x_new, 0.975), 
           err = "x", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
    plotCI(x = median(x_new), y = mean_out, li = mean_out - 2 * sd_out, ui = mean_out + 2 * sd_out, 
           err = "y", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
  }
  
  return(data.frame(mean_out=mean_out, sd_out = sd_out))
}


get.linear.model.constraint_ens <- function(y, x, x_new, plot.constraint = T) {
  
  dat = data.frame(y = y, x = x)
  new.dat = data.frame(x = x_new)
  reg.mod = lm(y ~ x, data = dat)
  
  n=length(reg.mod$residuals)
  RSS = sum((reg.mod$residuals)^2)
  s2 = 1/(n - 2) * RSS  # MSE
  sigma_x = sqrt( var(dat$x, na.rm = T) )
  new.dat.pred = predict(reg.mod, newdata = new.dat)
  sigma_f_x = sqrt(s2) * sqrt(1 + 1/n + (median(new.dat.pred) - mean(dat$x, na.rm = T)) / (n * sigma_x^2) )
  
  # combine variances in quadrature:
  mean_out =  median(new.dat.pred)
  sd_out =  sqrt(sigma_f_x^2 + var(new.dat.pred))
  
  
  # plot emergent constraint for an ensemble:
  if(plot.constraint == T) {
    plot(x, y, col="red")
    abline(reg.mod, col ="red")
    
    plotCI(x = median(x_new), y = 0, li = quantile(x_new, 0.025), ui = quantile(x_new, 0.975), 
           err = "x", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
    plotCI(x = median(x_new), y = mean_out, li = mean_out - 2 * sd_out, ui = mean_out + 2 * sd_out, 
           err = "y", add = T, col = "grey40", pch = 16, lwd = 2, lty = 2)
  }
  
  return(data.frame(mean_out=mean_out, sd_out = sd_out))
}




## Function to produce "matrix" of CMIP6 reconstructions:
get.CMIP6.recon.matrix <- function(CMIP6.tas_land_all.df, CMIP6.tos_all.df) {
  
    # select CMIP6 historical members:
    CMIP6.tas_land.all = data.frame(cbind(all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
    CMIP6.tos.all = data.frame(cbind(all = paste(CMIP6.tos_all.df$M$mod, "_", CMIP6.tos_all.df$M$scen, "_", CMIP6.tos_all.df$M$ens.mem, sep="")))
    ens.mem = data.frame(cbind(mod=CMIP6.tas_land_all.df$M$mod, scen = CMIP6.tas_land_all.df$M$scen, ens.mem = CMIP6.tas_land_all.df$M$ens.mem, 
                               all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
    ens.mem.un = unique(ens.mem)
    ens.mem.un = ens.mem.un[which(ens.mem.un$scen == "historical"),]
    # remove all ensemble members that contain NA's:
    na.mems = unique(CMIP6.tos.all$all[which(is.na(CMIP6.tos_all.df$ann$Yhat$GMST_FM))])
    omit.ix=na.omit(match(x = na.mems, table = ens.mem.un$all))
    if (length(omit.ix) > 0)  ens.mem.un = ens.mem.un[-omit.ix,]
    
    CMIP6.tas_land_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tos_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tas_land_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tos_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tas_land_hist_mod_p0 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tos_hist_mod_p0 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tas_land_hist_mod_p0_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.tos_hist_mod_p0_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.GMLSAT_NI = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    CMIP6.GMSST = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
    

    for (en in 1:dim(ens.mem.un)[1]) {
      print(en)
      
      ix_tas_land = which(CMIP6.tas_land.all$all == ens.mem.un$all[en] & CMIP6.tas_land_all.df$M$year %in% 1850:2014)
      ix_tos = which(CMIP6.tos.all$all == ens.mem.un$all[en] & CMIP6.tos_all.df$M$year %in% 1850:2014)
      
      CMIP6.tas_land_hist_mod_p1[en,] = CMIP6.tas_land_all.df$ann$Yhat$GMST_FM[ix_tas_land]
      CMIP6.tos_hist_mod_p1[en,] = CMIP6.tos_all.df$ann$Yhat$GMST_FM[ix_tos]
      CMIP6.tas_land_hist_mod_p1_pt1[en,] = CMIP6.tas_land_all.df$ann$Yhat_pt1$GMST_FM[ix_tas_land]
      CMIP6.tos_hist_mod_p1_pt1[en,] = CMIP6.tos_all.df$ann$Yhat_pt1$GMST_FM[ix_tos]
      
      CMIP6.tas_land_hist_mod_p0[en,] = CMIP6.tas_land_all.df$ann$Yhat_mod_p0$GMST_FM[ix_tas_land]
      CMIP6.tos_hist_mod_p0[en,] = CMIP6.tos_all.df$ann$Yhat_mod_p0$GMST_FM[ix_tos]
      CMIP6.tas_land_hist_mod_p0_pt1[en,] = CMIP6.tas_land_all.df$ann$Yhat_mod_p0_pt1$GMST_FM[ix_tas_land]
      CMIP6.tos_hist_mod_p0_pt1[en,] = CMIP6.tos_all.df$ann$Yhat_mod_p0_pt1$GMST_FM[ix_tos]
      
      # get additional variables for reviewer:
      CMIP6.GMLSAT_NI[en,] = CMIP6.tas_land_all.df$ann$Y$GMLSAT_NI[ix_tas_land]
      CMIP6.GMSST[en,] = CMIP6.tas_land_all.df$ann$Y$GMSST[ix_tas_land]
    }
    
    ret.list = list(CMIP6.tas_land_hist_mod_p1, CMIP6.tos_hist_mod_p1, CMIP6.tas_land_hist_mod_p1_pt1, CMIP6.tos_hist_mod_p1_pt1, 
                  CMIP6.tas_land_hist_mod_p0, CMIP6.tos_hist_mod_p0, CMIP6.tas_land_hist_mod_p0_pt1, CMIP6.tos_hist_mod_p0_pt1,
                  CMIP6.GMLSAT_NI, CMIP6.GMSST, ens.mem.un)
    names(ret.list) = c("tas_land_mod_p1", "tos_mod_p1", "tas_land_mod_p1_pt1", "tos_mod_p1_pt1",
                        "tas_land_mod_p0", "tos_mod_p0", "tas_land_mod_p0_pt1", "tos_mod_p0_pt1",
                        "CMIP6.GMLSAT_NI", "CMIP6.GMSST", "ens.mem.un")
    
    return(ret.list)
}








