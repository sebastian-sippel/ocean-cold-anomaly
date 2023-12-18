
# ------------------------------------------------------------------------------------
# Evaluate tos reconstruction based on ocean cells.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 25.10.2021
library(matrixStats)
library(zoo)
library(dplR)

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//code/_functions_CMIP6.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//code/_attribution_hildreth-lu.R")

# Load global observations / reconstruction:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03a_load_global_observations.R")
# load additional SST datasets:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/SST_datasets.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/TLand_datasets.RData")

# 00.(c) Load *new* reconstructions:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/03b_load_global_obs_reconstruction.R")

# 00.(d) Load forced response estimates from DAMIP simulations:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_ALL_ct_f.RData")
CMIP6.MMM = rowMeans((CMIP6.tas_ann_ALL_ct_f$AGMT_f_hist))
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HISTssp245_ct_f.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/CMIP6_DAMIP/_processed/CMIP6.tas_ann_HIST_ct_f.RData")
CMIP6.MMM = rowMeans(CMIP6.tas_ann_HISTssp245_ct_f$AGMT_f_hist)[1:171]


# 00.(e) CMIP6-reconstructions:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tas_land_all.df_v4.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedCMIP6_reconstr/CMIP6.tos_all.df_v4.RData")
CMIP6.tos_all.df$ann$Yhat$GSAT[which(CMIP6.tos_all.df$ann$Yhat$GSAT > 20)] = NA


# 00.(f) Define functions:
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

get.diff.window.cor_ens <- function(x, y, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix, probs = c(0.025, 0.5, 0.975)) {
    
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
  diff_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x_y)), probs = probs)
  diff_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.high_y.high)), probs = probs)
  diff_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$diff_x.low_y.low)), probs = probs)
  
  cor_x_y = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x_y)), probs = probs)
  cor_x.high_y.high = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.high_y.high)), probs = probs)
  cor_x.low_y.low = colQuantiles(x = t(sapply(X = dat, FUN=function(x) x$cor_x.low_y.low)), probs = probs)
  
  ret.df = data.frame(diff_x_y, diff_x.high_y.high, diff_x.low_y.low, cor_x_y, cor_x.high_y.high, cor_x.low_y.low)
  names(ret.df) = c("diff_x_y_2.5", "diff_x_y_50", "diff_x_y_97.5", "diff_x.high_y.high_2.5", "diff_x.high_y.high_50", "diff_x.high_y.high_97.5", "diff_x.low_y.low_2.5", "diff_x.low_y.low_50", "diff_x.low_y.low_97.5", 
                    "cor_x_y_2.5", "cor_x_y_50", "cor_x_y_97.5", "cor_x.high_y.high_2.5", "cor_x.high_y.high_50", "cor_x.high_y.high_97.5", "cor_x.low_y.low_2.5", "cor_x.low_y.low_50", "cor_x.low_y.low_97.5")
  return(ret.df)
}

get.trend <- function(x, trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014) {
  
  trends.out = sapply(X = 1:length(trend.years), FUN=function(i) {
    ix = match(x = trend.years[[i]], table = years)
    return(lm(x[ix] ~ years[ix])$coefficients[2] * length(ix))
  })
  names(trends.out) = paste("trend", 1:length(trend.years), sep="")
  return(trends.out)
}


## 0.(g) get mean-removed and PSL reconstruction:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processedOBS_reconstr/_old/GSAT.HadSLP2_psl_v2.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/data/03_processed_small_data/AGMT.tos_MR_.RData")


# GSAT Prepare data for plotting:
GSAT.tas_land = get.df(Y = OBS.tas_land_$GSAT$ann$mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2020)
GSAT.tos = get.df(Y = OBS.tos_$GSAT$ann$mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = 94:200, years.DA = 1850:2020)
GSAT.tos_ua = get.df(Y = HadSST_ua.annual$GSAT_mod_p1_min, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = NULL, years.DA = 1850:2020)
GSAT.hybrid36 = get.df(Y = OBS_hybrid36.tos_$GSAT$ann$mod_p1, f = CMIP6.MMM, years = 1850:2016, center = T, ens.ix = 1:200, years.DA = 1850:2020)
GSAT.tos_MR = get.df(Y = AGMT.tos_MR_$ann$mod_p1, f = CMIP6.MMM, years = 1850:2020, center = T, ens.ix = 1:200, years.DA = 1850:2020)



# plot(x = 1850:2020, y = cor_vec_hadsst4, type='l', ylim = c(-1, 1))
# lines(x = 1850:2020, y = cor_vec_hadsst4 + cor_vec_cmip6, col = "darkblue")
# lines(x = 1850:2016, y = cor_vec_hadsst4[1:167] + cor_vec_hybrid36, col = "darkgoldenrod4")

## 01. (a) Process observations
## ----------------------------------------
OBS_mod_p1_min = get.diff.window.cor_ens(x = OBS.tos_$GSAT$ann$mod_p1_min, y = OBS.tas_land_$GSAT$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 94:200)
OBS_mod_p1_min_PROB7765 = get.diff.window.cor_ens(x = OBS.tos_$GSAT$ann$mod_p1_min, y = OBS.tas_land_$GSAT$ann$mod_p1_min, years = 1850:2020, w.width = 51, W = 20, center = T, center.ens.means = T, ens.ix = 94:200, 
                                                  probs = c(0.2236, 0.5, 0.7765))

# Trends in observations:
alpha = 0.67
OBS.trends = data.frame(matrix(NA, nrow = 200, ncol = 15))
names(OBS.trends) = c("TMLSAT_tas_land1", "TMLSAT_tas_land2", "TMLSAT_tas_land3", 
                        "TMMSAT_tos1", "TMMSAT_tos2", "TMMSAT_tos3", 
                        "GSAT_tas_land1", "GSAT_tas_land2", "GSAT_tas_land3", 
                        "GSAT_tos1", "GSAT_tos2", "GSAT_tos3", 
                      "GSAT1", "GSAT2", "GSAT3")

for (i in 1:200) {
  print(i)
  # Get trends: 
  OBS.trends[i,1:3] = get.trend(x = OBS.tas_land_$TMLSAT_40S_40N$ann$mod_p1_min[i,], trend.years = list(1900:1939, 1900:1950, 1980:2020), years = 1850:2020)
  OBS.trends[i,4:6] = get.trend(x = OBS.tos_$TMMSAT_40S_40N$ann$mod_p1_min[i,], trend.years = list(1900:1939, 1900:1950, 1980:2020), years = 1850:2020)
  OBS.trends[i,7:9] = get.trend(x = OBS.tas_land_$GSAT$ann$mod_p1_min[i,], trend.years = list(1900:1939, 1900:1950, 1980:2020), years = 1850:2020)
  OBS.trends[i,10:12] = get.trend(x = OBS.tos_$GSAT$ann$mod_p1_min[i,], trend.years = list(1900:1939, 1900:1950, 1980:2020), years = 1850:2020)
  # GSAT average:
  OBS.trends$GSAT1[i] = (1 - alpha) * OBS.trends$GSAT_tas_land1[i] + alpha * OBS.trends$GSAT_tos1[i]
  OBS.trends$GSAT2[i] = (1 - alpha) * OBS.trends$GSAT_tas_land2[i] + alpha * OBS.trends$GSAT_tos2[i]
  OBS.trends$GSAT3[i] = (1 - alpha) * OBS.trends$GSAT_tas_land3[i] + alpha * OBS.trends$GSAT_tos3[i]
}


## 01. (b) Scale observations:
## ----------------------------------------
wr = c(mean(CRUTEM5.global.annual$Anomaly[142:171]) - mean(CRUTEM5.global.annual$Anomaly[8:51])) / c(mean(HadSST4.global.annual$Anomaly[142:171]) - mean(HadSST4.global.annual$Anomaly[8:51]))
alpha = 0.67
scale.land = alpha / wr + 1 - alpha   # see calculation!!
scale.sst = alpha + (1-alpha) * wr


GSAT_blended_ann_mod_p1 = colMedians(OBS.tos_$GSAT$ann$mod_p1_min) * sea_fraction + colMedians(OBS.tas_land_$GSAT$ann$mod_p1_min) * land_ice_fraction
ocean.raw.scale = HadSST4.global.annual$Anomaly[1:171] * scale.sst
land.raw.scale = CRUTEM5.global.annual$Anomaly[1:171] * scale.land


# regress CRUTEM5 and HadSST4:
raw.data.scale = get.diff.window.cor(x = ocean.raw.scale[8:171], y = land.raw.scale[8:171], years = 1857:2020, w.width = 51, W = 20, center = T)




## 01. (c) Process CMIP6
## ----------------------------------------
{
  # select CMIP6 historical members:
  CMIP6.tas_land.all = data.frame(cbind(all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  CMIP6.tos.all = data.frame(cbind(all = paste(CMIP6.tos_all.df$M$mod, "_", CMIP6.tos_all.df$M$scen, "_", CMIP6.tos_all.df$M$ens.mem, sep="")))
  ens.mem = data.frame(cbind(mod=CMIP6.tas_land_all.df$M$mod, scen = CMIP6.tas_land_all.df$M$scen, ens.mem = CMIP6.tas_land_all.df$M$ens.mem, 
                       all = paste(CMIP6.tas_land_all.df$M$mod, "_", CMIP6.tas_land_all.df$M$scen, "_", CMIP6.tas_land_all.df$M$ens.mem, sep="")))
  ens.mem.un = unique(ens.mem)
  ens.mem.un = ens.mem.un[which(ens.mem.un$scen == "historical"),]
  # remove all ensemble members that contain NA's:
  na.mems = unique(CMIP6.tos.all$all[which(is.na(CMIP6.tos_all.df$ann$Yhat$GSAT))])
  omit.ix=na.omit(match(x = na.mems, table = ens.mem.un$all))
  ens.mem.un = ens.mem.un[-omit.ix,]
  
  CMIP6.tas_land_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tos_hist_mod_p1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tas_land_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.tos_hist_mod_p1_pt1 = matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 165)
  CMIP6.DIFF = NA
  
  # Get trends for each ensemble member:
  CMIP6.trends = data.frame(matrix(NA, nrow = dim(ens.mem.un)[1], ncol = 15))
  names(CMIP6.trends) = c("TMLSAT_tas_land1", "TMLSAT_tas_land2", "TMLSAT_tas_land3", 
                        "TMMSAT_tos1", "TMMSAT_tos2", "TMMSAT_tos3", 
                        "GSAT_tas_land1", "GSAT_tas_land2", "GSAT_tas_land3", 
                        "GSAT_tos1", "GSAT_tos2", "GSAT_tos3", 
                        "GSAT1", "GSAT2", "GSAT3")
  
  
  for (en in 1:dim(ens.mem.un)[1]) {
    print(en)
    
    ix_tas_land = which(CMIP6.tas_land.all$all == ens.mem.un$all[en] & CMIP6.tas_land_all.df$M$year %in% 1850:2014)
    ix_tos = which(CMIP6.tos.all$all == ens.mem.un$all[en] & CMIP6.tos_all.df$M$year %in% 1850:2014)
    
    CMIP6.tas_land_hist_mod_p1[en,] = CMIP6.tas_land_all.df$ann$Yhat$GSAT[ix_tas_land]
    CMIP6.tos_hist_mod_p1[en,] = CMIP6.tos_all.df$ann$Yhat$GSAT[ix_tos]
    
    CMIP6.tas_land_hist_mod_p1_pt1[en,] = CMIP6.tas_land_all.df$ann$Yhat_pt1$GSAT[ix_tas_land]
    CMIP6.tos_hist_mod_p1_pt1[en,] = CMIP6.tos_all.df$ann$Yhat_pt1$GSAT[ix_tos]
    
    # Get trends: 
    CMIP6.trends[en,1:3] = get.trend(x = CMIP6.tas_land_all.df$ann$Yhat$TMLSAT_40S_40N[ix_tas_land], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014)
    CMIP6.trends[en,4:6] = get.trend(x = CMIP6.tos_all.df$ann$Yhat$TMMSAT_40S_40N[ix_tos], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014)
    CMIP6.trends[en,7:9] = get.trend(x = CMIP6.tas_land_all.df$ann$Yhat$GSAT[ix_tas_land], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014)
    CMIP6.trends[en,10:12] = get.trend(x = CMIP6.tos_all.df$ann$Yhat$GSAT[ix_tos], trend.years = list(1900:1939, 1900:1950, 1980:2014), years = 1850:2014)
    
    CMIP6.trends$GSAT1[en] = (1 - alpha) * CMIP6.trends$GSAT_tas_land1[en] + alpha * CMIP6.trends$GSAT_tos1[en]
    CMIP6.trends$GSAT2[en] = (1 - alpha) * CMIP6.trends$GSAT_tas_land2[en] + alpha * CMIP6.trends$GSAT_tos2[en]
    CMIP6.trends$GSAT3[en] = (1 - alpha) * CMIP6.trends$GSAT_tas_land3[en] + alpha * CMIP6.trends$GSAT_tos3[en]
    
    CMIP6.DIFF[en] = mean(CMIP6.tos_all.df$ann$Yhat$GSAT[ix_tos][51:90] - CMIP6.tas_land_all.df$ann$Yhat$GSAT[ix_tas_land][51:90])
  }
  
  CMIP6_mod_p1_min = get.diff.window.cor_ens(x = CMIP6.tos_hist_mod_p1, y = CMIP6.tas_land_hist_mod_p1, years = 1850:2014, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix = 1:599)
  CMIP6_mod_p1_min_pt1 = get.diff.window.cor_ens(x = CMIP6.tos_hist_mod_p1_pt1, y = CMIP6.tas_land_hist_mod_p1_pt1, years = 1850:2014, w.width = 51, W = 20, center = F, center.ens.means = F, ens.ix = 1:599)
  
  CMIP6_mod_p1_min_PROB7765 = get.diff.window.cor_ens(x = CMIP6.tos_hist_mod_p1, y = CMIP6.tas_land_hist_mod_p1, years = 1850:2014, w.width = 51, W = 20, center = T, center.ens.means = F, ens.ix = 1:599, 
                                             probs = c(0.2236, 0.5, 0.7765))
}



## Get correction vector(s):

# get correction vector from Cowtan
cor_vec_hybrid36 = GSAT.hybrid36$mod_p1_min_50[1:167] - GSAT.tos$mod_p1_min_50[1:167]
# plot(GSAT.hybrid36$mod_p1_min_50 - GSAT.tos$mod_p1_min_50)
cor_vec_hadsst4 = GSAT.tos$mod_p1_min_50 - GSAT.tos_ua$mod_p1_min_50

mean(c(cor_vec_hybrid36+cor_vec_hadsst4[1:167])[51:90])
mean(c(cor_vec_hadsst4[1:167])[51:90])

colMedians(OBS.tos_$GSAT$ann$mod_p1_min) - GSAT.tos_ua$mod_p1_min_50


# str(OBS.tos_$GSAT$ann$mod_p1)

cor_vec_cmip6 = rep(0, 171)
cor_vec_cmip6_ = rep(0, 171)
cor.ix = which((CMIP6_mod_p1_min$diff_x.low_y.low_2.5 - OBS_mod_p1_min$diff_x.low_y.low_97.5[1:165]) > 0)
cor_vec_cmip6[cor.ix] = (CMIP6_mod_p1_min$diff_x.low_y.low_2.5 - OBS_mod_p1_min$diff_x.low_y.low_97.5[1:165])[cor.ix]
cor.ix = which((CMIP6_mod_p1_min_PROB7765$diff_x.low_y.low_2.5 - OBS_mod_p1_min_PROB7765$diff_x.low_y.low_97.5[1:165]) > 0)
cor_vec_cmip6_[cor.ix] = (CMIP6_mod_p1_min_PROB7765$diff_x.low_y.low_2.5 - OBS_mod_p1_min_PROB7765$diff_x.low_y.low_97.5[1:165])[cor.ix]


# SST Adjustment:
{
  plot(x = 1850:2020, y = cor_vec_hadsst4, type='l', ylim = c(-1, 1), col = "blue")
  lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6), col = "blue", lty = 2)
  lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36), col = "blue", lty = 3)
  lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6_), col = "blue", lty = 2)
  
  corix = which((cor_vec_hadsst4[1:167] + cor_vec_hybrid36) > 0.5)  # 1900->1939
}
## ask nicolai about idea for correction...



DIFF_OBS_avg = sapply(X = 94:200, FUN=function(i) mean(OBS.tos_$GSAT$ann$mod_p1_min[i,corix] - OBS.tas_land_$GSAT$ann$mod_p1_min[i,corix]))
str(CMIP6.tas_land_all.df$ann$Yhat$GSAT)

str(OBS.tos_$GSAT$ann$mod_p1_min[94:200,])


str(CMIP6.tas_land_all.df)




# ------------------------------------------------------------------------------------
# 02. Plot cold anomaly plausibility in CMIP6
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/03_plausibility_cmip6/")

# mod_p1_min:
pdf(file = "01_plausibility_CMIP6_diff_cor.pdf", width = 8, height=10)
{
  par(mfrow=c(1, 1), mar=c(3,5,1,5))
  ylim = c(-3.6-1.25, 6.3); xlim = c(1848,2022)
  
  plot(x = 1850:2020, y = 1850:2020, type="n", 
       ylab = "", xlab = "", main = "", ylim = ylim, xlim = xlim, las=1, yaxt = "n", xaxt = "n", bty = "n", yaxs="i", xaxs="i")
  mtext(text = "Temperature difference, GSAT(ocean) minus GSAT(land) [°C]", side = 2, line = 3, at = 3.5, adj = 0.5)
  mtext(text = "Pearson Correlation, ocean/land (in 51-year window)", side = 2, line = 3, at = -2.75, adj = 0.5)
  
  {
    lines(x = c(1850, 1850), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1875, 1875), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1900, 1900), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1925, 1925), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1950, 1950), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(1975, 1975), y = c(-10, 10), col = "lightgrey", lty = 3)
    lines(x = c(2000, 2000), y = c(-10, 10), col = "lightgrey", lty = 3)
    
    axis(side = 1, at = c(1850, 1900, 1950, 2000, 2020))
    axis(side = 1, at = seq(1850, 2020, 5), tcl=0.2, labels=F)
    
    axis(side = 2, at = seq(4, 6, 0.5), labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(4, 6, 0.1), labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(5, 5), col = "darkgray", lty = 3)
    text(x = 1850, y = 6, labels = "a. Original Reconstruction", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(2, 4, 0.5) + 0.5, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 4, at = seq(2, 4, 0.1) + 0.5, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(3, 3) + 0.5, col = "darkgray", lty = 3)
    text(x = 1850, y = 4.1, labels = "b. Low-pass filtered \n (>20yr)", col = "grey40", pos = 4)
    
    axis(side = 2, at = seq(0, 2, 0.5) + 0.5*2, labels=seq(-1, 1, 0.5) / 2, las = 1)  # +0.9
    axis(side = 2, at = seq(0, 2, 0.1) + 0.5*2, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(1, 1) + 0.5*2, col = "darkgray", lty = 3)
    text(x = 1850, y = 1.7 + 0.5*2, labels = "c. High-pass filtered \n (<20yr)", col = "grey40", pos = 4)
    
    axis(side = 4, at = seq(-2, 0, 0.5) + 0.5*3, labels=seq(-1, 1, 0.5), las = 1, col = "blue", col.axis = "blue", col.ticks = "blue")   # +0.9
    axis(side = 4, at = seq(-2, 0, 0.1) + 0.5*3, labels=F, tcl=0.2)
    lines(x = c(1850, 2020), y = c(-1, -1) + 0.5*3, col = "darkgray", lty = 3)
    mtext(text = "Global SST Adjustment [°C]", side = 4, line = 3, at = 0.5, col = "blue")
    text(x = 1850, y = -0.25 + 0.5*3, labels = "d. SST Adjustments", col = "darkblue", pos = 4)
    
    lines(x = c(1850, 2020), y = c(-0.5, -0.5), col = "black", lty = 1)
  
    text(x = 1855, y = 0.5-1.25, labels = "e. Original Reconstruction, Correlation", col = "grey40", pos = 4)
    axis(side = 2, at = seq(-2, 0, 0.5) + 0.5 - 1.25, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 2, at = seq(-2, 0, 0.1) + 0.5 - 1.25, labels=F, tcl=0.2)
    
    axis(side = 4, at = seq(-4, -2, 0.5) + 0.5 - 1.25, labels=seq(0, 1, 0.25), las = 1)  # +0.9
    axis(side = 4, at = seq(-4, -2, 0.1) + 0.5 - 1.25, labels=F, tcl=0.2)
    text(x = 1855, y = -2.75, labels = "f. High-pass filtered (<20yr), Correlation", col = "grey40", pos = 4)
    
    # axis(side = 4, at = seq(-6, -4, 0.5) + 0.5*5, labels=seq(-1, 1, 0.5), las = 1)  # +0.9
    # axis(side = 4, at = seq(-6, -4, 0.1) + 0.5*5, labels=F, tcl=0.2)
    # lines(x = c(1850, 2020), y = c(-5, -5) + 0.5*5, col = "darkgray", lty = 3)
    # text(x = 1855, y = -4.5 + 0.5*5, labels = "Unforced, high-pass filtered (<20yr)", col = "grey40", pos = 4)
  }
  
  # Full reconstruction:
  {
    c = 5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x_y_2.5, rev(CMIP6_mod_p1_min$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x_y_2.5, rev(OBS_mod_p1_min$diff_x_y_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # Low-pass filtered:
  {
    c = 3 + 0.5; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x.low_y.low_2.5, rev(CMIP6_mod_p1_min$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x.low_y.low_2.5, rev(OBS_mod_p1_min$diff_x.low_y.low_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.low_y.low * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.low_y.low_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered:
  {
    c = 1+0.5*2; sc = 2
    # CMIP6:
    polygon(x = c(1850:2014, 2014:1850), y = c(CMIP6_mod_p1_min$diff_x.high_y.high_2.5, rev(CMIP6_mod_p1_min$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(1850:2020, 2020:1850), y = c(OBS_mod_p1_min$diff_x.high_y.high_2.5, rev(OBS_mod_p1_min$diff_x.high_y.high_97.5)) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$diff_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$diff_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # SST Adjustment:
  {
    c = -1+0.5*3; sc = 1
    
    lines(x = 1850:2020, y = cor_vec_hadsst4 * sc + c, type='l', ylim = c(-1, 1), col = "blue")
    #lines(x = 1850:2020, y = (cor_vec_hadsst4 + cor_vec_cmip6) * sc + c, col = "blue", lty = 2)
    lines(x = 1850:2016, y = (cor_vec_hadsst4[1:167] + cor_vec_hybrid36) * sc + c, col = "blue", lty = 3)
  }
  
  
  
  
  # Original Reconstruction, Correlation:
  {
    c = -2 + 0.5 - 1.25; sc = 2

    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p1_min_pt1$cor_x_y_2.5[26:140], rev(CMIP6_mod_p1_min_pt1$cor_x_y_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x_y_2.5[26:146], rev(OBS_mod_p1_min$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    # Original data:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x_y_2.5[26:146], rev(OBS_mod_p1_min$cor_x_y_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x_y_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x_y * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x_y_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  # High-pass filtered Reconstruction, Correlation:
  {
    c = -4 + 0.5 - 1.25; sc = 2
    
    # CMIP6:
    polygon(x = c(c(1850:2014)[26:140], c(2014:1850)[26:140]), y = c(CMIP6_mod_p1_min_pt1$cor_x.high_y.high_2.5[26:140], rev(CMIP6_mod_p1_min_pt1$cor_x.high_y.high_97.5[26:140])) * sc + c,
            col = make.transparent.color("red", alpha = 60), border = make.transparent.color("red", 150))
    
    # OBS:
    polygon(x = c(c(1850:2020)[26:146], c(2020:1850)[26:146]), y = c(OBS_mod_p1_min$cor_x.high_y.high_2.5[26:146], rev(OBS_mod_p1_min$cor_x.high_y.high_97.5[26:146])) * sc + c,
            col = make.transparent.color("grey40", alpha = 60), border = make.transparent.color("grey40", 150))
    lines(x = 1850:2014, y = CMIP6_mod_p1_min_pt1$cor_x.high_y.high_50 * sc + c, col = "red", lwd = 2)
    lines(x = 1857:2020, y = raw.data.scale$cor_x.high_y.high * sc + c, col = "black", lwd = 2, lty = 2)
    lines(x = 1850:2020, y = OBS_mod_p1_min$cor_x.high_y.high_50 * sc + c, col = "grey40", lwd = 2)
  }
  
  legend("bottom", c("Observational reconstruction", "Raw observations, scaled", "CMIP6 reconstruction",
                          "HadSST4 adjustment", "Implied Cowtan-hybrid36 adjustment"), 
                          col = c("grey40", "black", "red", "blue", "blue"),  
         cex = 0.8, ncol = 2, inset = 0.02, lwd = 2, lty = c(1, 2, 1, 1, 3))
}
dev.off()



# ------------------------------------------------------------------------------------
# Get the difference for correction:
# ------------------------------------------------------------------------------------
plot(CMIP6_mod_p1_min$diff_x.low_y.low_2.5 - OBS_mod_p1_min$diff_x.low_y.low_97.5[1:165])

cor_vec = rep(0, 171)
cor.ix = which((CMIP6_mod_p1_min$diff_x.low_y.low_2.5 - OBS_mod_p1_min$diff_x.low_y.low_97.5[1:165]) > 0)
cor_vec[cor.ix] = (CMIP6_mod_p1_min$diff_x.low_y.low_2.5 - OBS_mod_p1_min$diff_x.low_y.low_97.5[1:165])[cor.ix]
plot(cor_vec)


get.trend(x = colMedians(OBS.tas_land_$GSAT$ann$mod_p1_min), trend.years = list(1900:1939), years = 1850:2020)
get.trend(x = colMedians(OBS.tos_$GSAT$ann$mod_p1_min), trend.years = list(1900:1939), years = 1850:2020)
test.cor=get.trend(x = colMedians(OBS.tos_$GSAT$ann$mod_p1_min) + cor_vec, trend.years = list(1900:1939), years = 1850:2020)
test.cor2=get.trend(x = colMedians(OBS.tos_$TMMSAT_40S_40N$ann$mod_p1_min) + cor_vec, trend.years = list(1900:1939), years = 1850:2020)
# test.cor2=get.trend(x = colMedians(OBS.tos_$TMMSAT_40S_40N$ann$mod_p1_min), trend.years = list(1900:1939), years = 1850:2020)



### Plot the different trend distributions for ETCW:

plot(c(1,1), type='n', xlim = c(0, 1), ylim = c(0, 10))
lines(density(OBS.trends$GSAT_tas_land1), col = "darkorange")
lines(density(OBS.trends$GSAT_tos1), col = "blue")
lines(density(CMIP6.trends$GSAT1), col = "red")
lines(x = c(test.cor, test.cor), y = c(0, 10))





# ------------------------------------------------------------------------------------
# 03. Compare the patterns of the MR and Mean-included data projected on piControl simulations:
# ------------------------------------------------------------------------------------

load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_CMIP6_processed/CMIP6.tas_ann_piControl_ct.RData")
names(CMIP6.tas_ann_piControl_ct$Y)[1] = "GSAT"
trends_piControl = extract.nyear.trend_(XAX = CMIP6.tas_ann_piControl_ct, nyears = 40, trend.sep = 10, 
                                        var.names = c("GSAT_", "GMSST", "GMLSAT_NI", "TMMSAT_40S_40N", "TMLSAT_40S_40N"))
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/scripts/02c_process_piControl_v4.R")
TMLSAT_40S_40N.tas_land_piControl = extract.nyear.trend_from_reconstruction(XAX = CMIP6.tas_ann_piControl_ct, nyears = 40, trend.sep = 10, 
                                                                            var.names = c("TMLSAT_40S_40N"), 
                                                                            cur.file.name = "tas_land_predTMLSAT_v3", start.year = NULL, start.year.mask = 1900)

GSAT.tos_piControl = extract.nyear.sequence_from_reconstruction(XAX = CMIP6.tas_ann_piControl_ct, nyears = 70, trend.sep = 20, 
                                            var.names = c("GSAT"), 
                                            cur.file.name = "tos_predGSAT_v3", start.year = NULL, start.year.mask = 1885)
## try to extract nyear.sequence from MR-reconstruction:
GSAT.tos_MR_piControl = extract.nyear.sequence_from_reconstruction_4MR(XAX = CMIP6.tas_ann_piControl_ct, nyears = 70, trend.sep = 20, 
                                                                var.names = c("GSAT"), 
                                                                cur.file.name = "tos", start.year = NULL, start.year.mask = 1885)

str(GSAT.tos_MR_piControl)


plot(GSAT.tos_piControl$Y, GSAT.tos_piControl$Yhat)
plot(GSAT.tos_MR_piControl$Y, GSAT.tos_MR_piControl$Yhat)


exclude.ix = unique(c(which(GSAT.tos_piControl$M$mod == "E3SM-1-1"), unique(which(GSAT.tos_MR_piControl$Yhat < -1.2, arr.ind=T)[,1])))


plot(GSAT.tos_piControl$Yhat[-exclude.ix,], GSAT.tos_piControl$Y[-exclude.ix,])
cor(c(GSAT.tos_piControl$Yhat[-exclude.ix,]), c(GSAT.tos_piControl$Y[-exclude.ix,]), use = "complete.obs")  # R = 0.94

plot(GSAT.tos_MR_piControl$Yhat[-exclude.ix,], GSAT.tos_MR_piControl$Y[-exclude.ix,])
cor(c(GSAT.tos_MR_piControl$Yhat[-exclude.ix,]), c(GSAT.tos_MR_piControl$Y[-exclude.ix,]), use = "complete.obs")  # R = 0.55


plot(GSAT.tos_piControl$Yhat[-exclude.ix,], GSAT.tos_MR_piControl$Yhat[-exclude.ix,])
abline(0, 1, col = "red")
cor(c(GSAT.tos_piControl$Yhat[-exclude.ix,]), c(GSAT.tos_MR_piControl$Yhat[-exclude.ix,]), use = "complete.obs")   # correlation of 0.6 -> so is there decadal variability?


GSAT.tos_Y_r10 = t(sapply(X = 1:1913, FUN=function(i) rollmean(x = GSAT.tos_piControl$Y[i,], k = 20, fill = NA)))
GSAT.tos_Yhat_r10 = t(sapply(X = 1:1913, FUN=function(i) rollmean(x = GSAT.tos_piControl$Yhat[i,], k = 20, fill = NA)))
GSAT.tos_MR_Yhat_r10 = t(sapply(X = 1:1913, FUN=function(i) rollmean(x = GSAT.tos_MR_piControl$Yhat[i,], k = 20, fill = NA)))
       
plot(GSAT.tos_Yhat_r10[-exclude.ix,], GSAT.tos_MR_Yhat_r10[-exclude.ix,])
cor(c(GSAT.tos_Yhat_r10[-exclude.ix,]), c(GSAT.tos_MR_Yhat_r10[-exclude.ix,]), use = "complete.obs")  # R = 0.54 correlation... 


## draw all 1900 lines and see whether decoupling is possible in the extent seen in obs:
plot(c(1,1), type="n", xlim = c(1,70), ylim = c(-0.7, 0.7))
sapply(X = c(1:1913)[-exclude.ix], FUN = function(i) lines(x = 1:70, y = GSAT.tos_Yhat_r10[i,]-GSAT.tos_MR_Yhat_r10[i,]))

lines(x = 1:70, y = c(GSAT.tos$res_lp_50 - GSAT.tos_MR$res_lp_50)[36:105], type='l', col = "red", lwd = 2)   # it still is far outside from anything... 


plot(c(1,1), type="n", xlim = c(1,70), ylim = c(-0.7, 0.7))
sapply(X = c(1:1913)[-exclude.ix], FUN = function(i) lines(x = 1:70, y = GSAT.tos_Yhat_r10[i,]))
lines(x = 1:70, y = c(GSAT.tos$res_lp_50)[36:105], type='l', col = "red", lwd = 2)   # it still is far outside from anything... 


### something useful??  -> is the ocean pattern itself unusual?





sd(rollmean(x = GSAT.tos_piControl$Y[1,], k = 10, fill = NA)[seq(5, 65, 5)])
sd(rollmean(x = GSAT.tos_piControl$Yhat[1,], k = 10, fill = NA)[seq(5, 65, 5)])



GSAT.tos_piControl$M[unique(which(GSAT.tos_piControl$Yhat < -0.5 & GSAT.tos_MR_piControl$Yhat > 0.5, arr.ind = T)[,1]),]






plot(x = 1850:2020, y = GSAT.tos$res_50, type='l')
lines(x = 1850:2020, y = GSAT.tos_MR$res_50, type='l', col = "blue")
lines(x = c(1885, 1885), y = c(-1, 1))


plot(x = 1850:2020, y = GSAT.tos$res_hp_50, type='l')
lines(x = 1850:2020, y = GSAT.tos_MR$res_hp_50, type='l', col = "blue")
cor(GSAT.tos$res_hp_50, GSAT.tos_MR$res_hp_50)


# ------------------------------------------------------------------------------------
# 03. Land-ocean warming ratio:
# ------------------------------------------------------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/03_plausibility_cmip6/")


par(mfrow=c(1, 1), mar=c(5,5,1,1))
ylim = c(0, 4); xlim = c(0,1)

plot(x = 1850:2020, y = 1850:2020, type="n", 
     ylab = "Land-ocean warming ratio [-]", xlab = "GSAT trend [°C]", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxs="i")

axis(side = 1, at = seq(0, 1, 0.02), tcl=0.2, labels=F)
axis(side = 2, at = seq(0, 4, 0.1), tcl=0.2, labels=F)
lines(x = c(0, 1), y = c(wr, wr))  # Observed long-term warming ratio

# ?? Byrne line for theory??

# Byrne line for theory:
gamma = 0.72 # ERA-Interim
dq_O = 0.11 / 1000 # in kg/kg per decade
dT_O = 0.12 # in K pre decade
Lv = 2257  # kJ / kg latent heat of vaporization
cp =  1.0 # kJ/kg/K specific heat capactiy of air at constant pressure
WR_byrne = 1 + (1 - gamma) * Lv / cp * dq_O / dT_O # for given ocean warming and ocean rel. humidity change.
# lines(x = c(0, 1), y = c(WR_byrne, WR_byrne), col = "darkblue")  # when derived from ocean data

# Observational estimate, early warming:
points(x = median(OBS.trends$GSAT1), y = median(OBS.trends$TMLSAT_tas_land1) / median(OBS.trends$TMMSAT_tos1), col = "grey40", pch = 16)
points(x = OBS.trends$GSAT1[94:200], y = OBS.trends$TMLSAT_tas_land1[94:200] / OBS.trends$TMMSAT_tos1[94:200], col = "grey40", pch = 16, cex = 0.5)

# Observational estimate, late warming:
points(x = median(OBS.trends$GSAT3), y = median(OBS.trends$TMLSAT_tas_land3) / median(OBS.trends$TMMSAT_tos3), col = "grey40", pch = 4)
points(x = OBS.trends$GSAT3[94:200], y = OBS.trends$TMLSAT_tas_land3[94:200] / OBS.trends$TMMSAT_tos3[94:200], col = "grey40", cex = 0.5, pch = 4)


# CMIP6 estimates:
points(x = CMIP6.trends$GSAT1, y = CMIP6.trends$TMLSAT_tas_land1 / CMIP6.trends$TMMSAT_tos1, col = "red", pch = 16)

# Late warming:
points(x = CMIP6.trends$GSAT3, y = CMIP6.trends$TMLSAT_tas_land3 / CMIP6.trends$TMMSAT_tos3, col = "red", pch = 4)



### hypothetical correction:
# x = test.cor * 2/3 + 1/3 * 0.271
# y = median(OBS.trends$TMLSAT_tas_land1) / test.cor2
# points(x = x, y = y, col = "blue", pch = 16)
# points(x = median(OBS.trends$GSAT3), y = median(OBS.trends$TMLSAT_tas_land3) / median(OBS.trends$TMMSAT_tos3), col = "grey40", pch = 4)


#### How does land-warming alone constrain the ocean warming?

xlim = c(-1, 1)
ylim = c(-1, 1)

plot(x = 1850:2020, y = 1850:2020, type="n", 
     ylab = "Ocean warming [°C]", xlab = "Land warming [°C]", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxs="i")
abline(0,1)

axis(side = 1, at = seq(-1, 1, 0.05), tcl=0.2, labels=F)
axis(side = 2, at = seq(-1, 1, 0.05), tcl=0.2, labels=F)

points(x = median(OBS.trends$TMLSAT_tas_land1), y = median(OBS.trends$TMMSAT_tos1), col = "grey40", pch = 4)
points(x = OBS.trends$TMLSAT_tas_land1[94:200], y = OBS.trends$TMMSAT_tos1[94:200], col = "grey40", cex = 0.5, pch = 4)

points(x = CMIP6.trends$TMLSAT_tas_land1, y = CMIP6.trends$TMMSAT_tos1, col = "red", pch = 16)

points(x = trends_piControl$Y$TMLSAT_40S_40N * 40, y = trends_piControl$Y$TMMSAT_40S_40N * 40, col = "orange")

trends_piControl$Y$TMLSAT_40S_40N
trends_piControl$Y$TMMSAT_40S_40N



plot(x = 1850:2020, y = 1850:2020, type="n", 
     ylab = "Ocean warming [°C]", xlab = "Land warming [°C]", main = "", ylim = ylim, xlim = xlim, las=1, bty = "n", yaxs="i", xaxs="i")
abline(0,1)

axis(side = 1, at = seq(-1, 1, 0.05), tcl=0.2, labels=F)
axis(side = 2, at = seq(-1, 1, 0.05), tcl=0.2, labels=F)

points(x = median(OBS.trends$GSAT_tas_land1), y = median(OBS.trends$GSAT_tos1), col = "grey40", pch = 4)
points(x = OBS.trends$GSAT_tas_land1[94:200], y = OBS.trends$GSAT_tos1[94:200], col = "grey40", cex = 0.5, pch = 4)

points(x = CMIP6.trends$GSAT_tas_land1, y = CMIP6.trends$GSAT_tos1, col = "red", pch = 16)












##### THE END ########
  

# Cold ocean anomaly shows up clearly. Reasons why this cold ocean anomaly may be unrealistic:
# (1) literature: biases in ocean are prevalent; comparison to coastal stations appears to show that ocean anomaly is unrealistic.
# (2) "decoupling" of low-frequency variability (poor match between land- and ocean reconstruction) from high-frequency variability (very good match).
# (3) none of the cmip6 models shows a temporal "decoupling pattern" between low- and high-frequency variability.
# (4) none of the cmip6 models supports an ocean-land temperature difference close to that observed.





