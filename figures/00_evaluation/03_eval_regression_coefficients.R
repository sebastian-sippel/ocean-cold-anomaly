
# ------------------------------------------------------------------------------------
# Evaluate tos and tas reconstruction(s) based on CMIP6.
# ------------------------------------------------------------------------------------

# Sebastian Sippel
# 23.03.2022


# screen -S make_plots
# R

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/code/_functions_CMIP6.R")

# get plotting for fingerprints:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_plot_projected_worldmap_v2.R")


# 00.(c) Load GSAT/GMST reconstructions:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tas_land.df_v5.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tos.df_v5.RData")



# 00.(d) Define functions for evaluation:
# cmip6_mod = CMIP6.tas_land.df
# mon.ix = 6
# date.ix = 45
{
  # Get the beta coefficients:
  get.beta_coefs <- function(cmip6_mod, cur.var = "GMST_FM", mon.ix, date.ix) {
    
    ## mod_gta
    beta_land_gta = rep(NA, 72*36); beta_land_gta_ = raster.template
    beta_land_gta[cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$grid.ix] = cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$mod_gta
    values(beta_land_gta_) = c(matrix(beta_land_gta, 72, 36)[,36:1])
    # plot(beta_land_gta_, zlim = c(-0.04, 0.04), col = col)
    
    ## mod_p0
    beta_land_p0 = rep(NA, 72*36); beta_land_p0_ = raster.template
    beta_land_p0[cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$grid.ix] = cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$mod_p0[,2]
    values(beta_land_p0_) = c(matrix(beta_land_p0, 72, 36)[,36:1])
    # plot(beta_land_p0_, zlim = c(-0.01, 0.01), col = col)
    
    ## mod_p1
    beta_land_p1 = rep(NA, 72*36); beta_land_p1_ = raster.template
    beta_land_p1[cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$grid.ix] = cmip6_mod$mon[[mon.ix]]$beta[[cur.var]][[date.ix]]$mod_p1[,2]
    values(beta_land_p1_) = c(matrix(beta_land_p1, 72, 36)[,36:1])
    # plot(beta_land_p1_, zlim = c(-0.01, 0.01), col = col)
    
    ## mod_p3
    beta_land_p2 = rep(NA, 72*36); beta_land_p2_ = raster.template
    # beta_land_p2[cmip6_mod[[mon.ix]][[date.ix]]$grid.ix] = cmip6_mod[[mon.ix]][[date.ix]]$mod_p2$beta[,2]
    # values(beta_land_p2_) = c(matrix(beta_land_p2, 72, 36)[,36:1])
    # plot(beta_land_p3_, zlim = c(-0.04, 0.04), col = col)
    
    ret.list = list(beta_land_gta_ = beta_land_gta_, beta_land_p0_ = beta_land_p0_, beta_land_p1_ = beta_land_p1_, beta_land_p2_ = beta_land_p2_)
    
    return(ret.list)
  }
}



### 1. REGRESSION COEFFICIENTS:
### -----------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/")

# how are regression coefficients changing based on different training?

library(RColorBrewer)
col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(99)
map.col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(99)


# 1A. MAP of regression coefficients (tas):
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/03_tas_eval_regression_coefficients/")
{
  # for (date.ix in seq(6, 171, 10)) {
  date.ix = 46
    for (mon.ix in 1:12) {
      ## Define all respective grid cells and plot as a map:
      cur.betas = get.beta_coefs(cmip6_mod = CMIP6.tas_land.df, mon.ix = mon.ix, date.ix = date.ix)
      cur.title = paste(month.name[mon.ix], " ", c(1850:2020)[date.ix], sep="")
      cur.file.name = paste("_", c(1850:2020)[date.ix], "-", formatC(mon.ix, width = 2, format = "d", flag = "0"), ".png", sep="")
      
      # 1A. Compare map of regression coefficients:      
      png(file = cur.file.name, width = 4.5, height = 9, units = "in", res = 200)
      par(mar=c(1.5,1,2,1), mfrow=c(3,1), oma = c(0.5, 0,0,0))
      
      zlim. = round(max(c(quantile(cur.betas$beta_land_p0_, probs=0.98), abs(quantile(cur.betas$beta_land_p0_, probs=0.02)))), 2)
      zlim = c(-zlim., zlim.)
      
      plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_gta_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (GTA)", main = cur.title)
      plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p0_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training no unc.)")
      plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p1_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training CMIP6+unc+bias)")
      # plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p3_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training CMIP6+3x unc+ 3x bias)")
      
      dev.off()
    }
}





### Make all plots for tas_regression coefficients:
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/03_tos_eval_regression_coefficients/")
{
  # for (date.ix in seq(6, 171, 10)) {
  date.ix = 46
  for (mon.ix in 1:12) {
    ## Define all respective grid cells and plot as a map:
    cur.betas = get.beta_coefs(cmip6_mod = CMIP6.tos.df, mon.ix = mon.ix, date.ix = date.ix)
    cur.title = paste(month.name[mon.ix], " ", c(1850:2020)[date.ix], sep="")
    cur.file.name = paste("_", c(1850:2020)[date.ix], "-", formatC(mon.ix, width = 2, format = "d", flag = "0"), ".png", sep="")
    
    # 1A. Compare map of regression coefficients:      
    png(file = cur.file.name, width = 4.5, height = 9, units = "in", res = 200)
    par(mar=c(1.5,1,2,1), mfrow=c(3,1), oma = c(0.5, 0,0,0))
    
    zlim. = round(max(c(quantile(cur.betas$beta_land_p0_, probs=0.98), abs(quantile(cur.betas$beta_land_p0_, probs=0.02)))), 2)
    zlim = c(-zlim., zlim.)
    
    plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_gta_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (GTA)", main = cur.title)
    plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p0_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training no unc.)")
    plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p1_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training CMIP6+unc+bias)")
    # plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p3_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training CMIP6+3x unc+ 3x bias)")
    
    dev.off()
  }
}


#### END !!! ####



##### Illustrative example for:
# (1) fingerprint(s)
# (2) Oberservations
# (3) bias
# (4) uncertainties


#### Illustative plots for one time step:
# Predicted vs. Observed
# Bias-variance trade-off illustration
# Zonal statistics


## load required metrics:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/HadCRUT5_mon1.RData")
# load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.GSAT.tas_land_4eval_v2.RData")

setwd("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3//figures/00_cmip_evaluation/02_illustration_time-step/")

pdf(file = cur.file.name, width = 6, height = 9)
par(mar=c(4,4,2,1), mfrow=c(1,1))

## CONTINUE HERE FOR TOS COEFS.





### Bias-Variance trade-off  (evaluated *with* uncertainties):
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/1895-06.RData")
pred_early = mod_p1$Yhat$pt1[,2]
obs = CMIP6.tas_land.df$mon[[6]]$Y$GMST_FM



pdf("02_bias_variance_zonal.pdf", width = 7, height = 9)
{
  par(mar=c(4,5,2,1), mfrow=c(3,2))
  
  ### June 1895 Predicted vs. observed:
  load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/1895-06.RData")
  pred = mod_p1$Yhat$pt1[,2]
  obs = CMIP6.tas_land.df$mon[[6]]$Y$GMST_FM
  xlim <- ylim <- c(-1, 1.5)

  plot(obs, pred, xlim = xlim, ylim = ylim, main = "June 1895 Coverage Mask+Unc./Bias",
       xlab = "CMIP6 GMST [°C]", ylab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16)   # Monthly predicted vs. observed for 1895-06 time step.
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs)^2, 2)), sep=""), inset = 0.02)
  
  
  ### June 1995 Predicted vs. observed:
  load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/1995-06.RData")
  pred = mod_p1$Yhat$pt1[,2]
  obs = CMIP6.tas_land.df$mon[[6]]$Y$GMST_FM
  
  plot(obs, pred, xlim = c(-1, 1.5), ylim = c(-1, 1.5), main = "June 1995 Coverage Mask+Unc./Bias",
       xlab = "CMIP6 GMST [°C]", ylab = "Predicted GMST [°C]", 
       col = make.transparent.color("grey50", alpha = 40), pch = 16)   # Monthly predicted vs. observed for 1895-06 time step.
  axis(side = 1, at = seq(min(xlim), max(xlim), 0.1), tcl=0.2, labels=F)
  axis(side = 2, at = seq(min(ylim), max(ylim), 0.1), tcl=0.2, labels=F)
  abline(0, 1, col = "red")
  legend("topleft", c(paste("MSE = ", round(mse(pred, obs), 3), sep=""), paste("R2 = ", round(cor(pred, obs)^2, 2)), sep=""), inset = 0.02)
  
  
  # June 1895
    load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/1895-06.RData")
    xlim = c(0, 0.04); ylim = c(0, 0.02)
    
    plot(c(1,1), type="n", las = 1, xlim = xlim, ylim = ylim,
         xlab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), ylab = expression("Mean Squared Bias per bias ensemble realization [" ~ K^2 ~ "]"),
         main = "")
    axis(side = 1, at = seq(min(xlim), max(xlim), 0.002), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.001), tcl=0.2, labels=F)
    
    
    # Without uncertainties+biases:
    lines(x = mod_p1$MSE$pt0, y = (colMeans((mod_p1$me.by.ensreal$pt0)^2)), col = "darkblue", lty = 3, lwd = 2)
    lines(x = mod_p0$MSE$pt0, y = (colMeans((mod_p0$me.by.ensreal$pt0)^2)), col = "darkorange", lty = 3, lwd = 2)
    points(x = mod_p1$MSE$pt0[mod_p1$lambda.1_05], y = colMeans((mod_p1$me.by.ensreal$pt0)^2)[mod_p1$lambda.1_05], col = "darkblue", pch = 15, cex = 1.5)
    points(x = mod_p0$MSE$pt0[mod_p0$lambda.1_05], y = colMeans((mod_p0$me.by.ensreal$pt0)^2)[mod_p0$lambda.1_05], col = "darkorange", pch = 15, cex = 1.5)
    points(x = mod_gta$MSE$pt0, y = (mean((mod_gta$me.by.ensreal$pt0)^2)), col = "grey50", pch = 15, cex = 1.5)
    
    # With uncertainties+biases:
    lines(x = mod_p1$MSE$pt1, y = (colMeans((mod_p1$me.by.ensreal$pt1)^2)), col = "darkblue", lwd = 2)
    lines(x = mod_p0$MSE$pt1, y = (colMeans((mod_p0$me.by.ensreal$pt1)^2)), col = "darkorange", lwd = 2)
    points(x = mod_p1$MSE$pt1[mod_p1$lambda.1_05], y = colMeans((mod_p1$me.by.ensreal$pt1)^2)[mod_p1$lambda.1_05], col = "darkblue", pch = 16, cex = 1.5)
    points(x = mod_p0$MSE$pt1[mod_p0$lambda.1_05], y = colMeans((mod_p0$me.by.ensreal$pt1)^2)[mod_p0$lambda.1_05], col = "darkorange", pch = 16, cex = 1.5)
    points(x = mod_gta$MSE$pt1, y = (mean((mod_gta$me.by.ensreal$pt1)^2)), col = "grey50", pch = 16, cex = 1.5)
    
    ## add arrow for adding uncertainties:
    arrows(x0 = mod_p1$MSE$pt0[mod_p1$lambda.1_05], y0 = colMeans((mod_p1$me.by.ensreal$pt0)^2)[mod_p1$lambda.1_05], 
           x1 = mod_p1$MSE$pt1[mod_p1$lambda.1_05], y1 = colMeans((mod_p1$me.by.ensreal$pt1)^2)[mod_p1$lambda.1_05], col = "darkblue")
    arrows(x0 = mod_p0$MSE$pt0[mod_p0$lambda.1_05], y0 = colMeans((mod_p0$me.by.ensreal$pt0)^2)[mod_p0$lambda.1_05], 
           x1 = mod_p0$MSE$pt1[mod_p0$lambda.1_05], y1 = colMeans((mod_p0$me.by.ensreal$pt1)^2)[mod_p0$lambda.1_05], col = "darkorange")
    arrows(x0 = mod_gta$MSE$pt0, y0 = (mean((mod_gta$me.by.ensreal$pt0)^2)), 
           x1 = mod_gta$MSE$pt1, y1 = (mean((mod_gta$me.by.ensreal$pt1)^2)), col = "grey50")
    
    # legend("topleft", c("Evaluation without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
    #       col = c("black", "grey50", "darkorange", "darkblue", "black", "black", "grey50", "darkorange", "darkblue"), 
     #      lty = c(NA, NA, 3, 3, NA, NA, NA, 1, 1), lwd = c(NA, NA,2,2,NA, NA, NA,2,2), pch = c(NA, 15, 15, 15, NA, NA, 16, 16, 19), inset = 0.02, cex = 0.8)
    
    
    # June 1995
    load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGMST_v3/1995-06.RData")
    plot(c(1,1), type="n", las = 1, xlim = c(0, 0.04), ylim = c(0, 0.02),
         xlab = expression("Monthly temperature reconstruction MSE [" ~ K^2 ~ "]"), ylab = expression("Mean Squared Bias per bias ensemble realization [" ~ K^2 ~ "]"),
         main = "")
    axis(side = 1, at = seq(min(xlim), max(xlim), 0.002), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 0.001), tcl=0.2, labels=F)
    
    # Without uncertainties+biases:
    lines(x = mod_p1$MSE$pt0, y = (colMeans((mod_p1$me.by.ensreal$pt0)^2)), col = "darkblue", lty = 3, lwd = 2)
    lines(x = mod_p0$MSE$pt0, y = (colMeans((mod_p0$me.by.ensreal$pt0)^2)), col = "darkorange", lty = 3, lwd = 2)
    points(x = mod_p1$MSE$pt0[mod_p1$lambda.1_05], y = colMeans((mod_p1$me.by.ensreal$pt0)^2)[mod_p1$lambda.1_05], col = "darkblue", pch = 15, cex = 1.5)
    points(x = mod_p0$MSE$pt0[mod_p0$lambda.1_05], y = colMeans((mod_p0$me.by.ensreal$pt0)^2)[mod_p0$lambda.1_05], col = "darkorange", pch = 15, cex = 1.5)
    points(x = mod_gta$MSE$pt0, y = (mean((mod_gta$me.by.ensreal$pt0)^2)), col = "grey50", pch = 15, cex = 1.5)
    
    # With uncertainties+biases:
    lines(x = mod_p1$MSE$pt1, y = (colMeans((mod_p1$me.by.ensreal$pt1)^2)), col = "darkblue", lwd = 2)
    lines(x = mod_p0$MSE$pt1, y = (colMeans((mod_p0$me.by.ensreal$pt1)^2)), col = "darkorange", lwd = 2)
    points(x = mod_p1$MSE$pt1[mod_p1$lambda.1_05], y = colMeans((mod_p1$me.by.ensreal$pt1)^2)[mod_p1$lambda.1_05], col = "darkblue", pch = 16, cex = 1.5)
    points(x = mod_p0$MSE$pt1[mod_p0$lambda.1_05], y = colMeans((mod_p0$me.by.ensreal$pt1)^2)[mod_p0$lambda.1_05], col = "darkorange", pch = 16, cex = 1.5)
    points(x = mod_gta$MSE$pt1, y = (mean((mod_gta$me.by.ensreal$pt1)^2)), col = "grey50", pch = 16, cex = 1.5)
    
    legend("topleft", c("Evaluation, without unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias", "", "Evaluation, incl. unc./biases:", "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), 
           col = c("black", "grey50", "darkorange", "darkblue", "black", "black", "grey50", "darkorange", "darkblue"), 
           lty = c(NA, NA, 3, 3, NA, NA, NA, 1, 1), lwd = c(NA, NA,2,2,NA, NA, NA,2,2), pch = c(NA, 15, 15, 15, NA, NA, 16, 16, 19), inset = 0.02, cex = 0.8)
    

    # Zonal L2 norm:
      xlim = c(0,round(max((L2_zonal_p0[,546]))*1.1, 4)); ylim = c(-90, 90)
      
      plot(c(1,1), type='n', xlim = xlim, ylim = ylim, las = 1,
           ylab = "Latitude [°N]", xlab = "Sum of squared coefficients in each lat. band")
      axis(side = 1, at = seq(min(xlim), max(xlim), 0.0002), tcl=0.2, labels=F)
      axis(side = 2, at = seq(min(ylim), max(ylim), 10), tcl=0.2, labels=F)
      
      lines(x = (L2_zonal_gta[,546]), y = z[-(nzone+1)] + 7.5, type='l', col = "grey50", lwd = 2)
      lines(x = (L2_zonal_p0[,546]), y = z[-(nzone+1)] + 7.5, type='l', col = "darkorange", lwd = 2)
      lines(x = (L2_zonal_p1[,546]), y = z[-(nzone+1)] + 7.5, type='l', col = "darkblue", lwd = 2)
      
      legend("bottomright", c( "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), col = c("grey50", "darkorange", "darkblue"), lwd = 2, cex = 0.8, inset = 0.02)
    
    
    # 1995 Zonal L2 norm:
    xlim. = round(max((L2_zonal_p0[,546]))*1.1, 4)
    xlim = c(0,round(max((L2_zonal_p0[,546]))*1.1, 4)); ylim = c(-90, 90)
    
    plot(c(1,1), type='n', xlim = xlim, ylim = ylim, las = 1,
         ylab = "Latitude [°N]", xlab = "Sum of squared coefficients in each lat. band")
    axis(side = 1, at = seq(min(xlim), max(xlim), 0.0002), tcl=0.2, labels=F)
    axis(side = 2, at = seq(min(ylim), max(ylim), 10), tcl=0.2, labels=F)
    
    lines(x = (L2_zonal_gta[,1746]), y = z[-(nzone+1)] + 7.5, type='l', col = "grey50", lwd = 2)
    lines(x = (L2_zonal_p0[,1746]), y = z[-(nzone+1)] + 7.5, type='l', col = "darkorange", lwd = 2)
    lines(x = (L2_zonal_p1[,1746]), y = z[-(nzone+1)] + 7.5, type='l', col = "darkblue", lwd = 2)
    
    # legend("bottomright", c( "GTA", "CMIP-Train, no unc./bias", "CMIP-Train, incl. unc./bias"), col = c("grey50", "darkorange", "darkblue"), lwd = 2, cex = 0.8, inset = 0.02)
}
    dev.off()



    
    
    


