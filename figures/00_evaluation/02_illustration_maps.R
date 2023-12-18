
### Illustrate time step reconstruction conceptually:
# Sebastian Sippel
# 02.09.2022

# 00.(a) load  respective functions & code:
source("/net/h2o/climphys1/sippels/_code/tools/frenchcolormap.R")
# source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_functions_CMIP5_extr.R")
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_functions_CMIP6.R")

# get plotting for fingerprints:
source("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v2/code/_plot_projected_worldmap_v2.R")




# 00.(c) Load AGMT reconstructions:
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/03_processed_CMIP6_4evaluation/CMIP6.tas_land.df_v3.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/_processed_CRU/CRUTEM5_mon6.RData")
load("/net/h2o/climphys1/sippels/_projects/global_mean_reconstr_v3/data/02_trained_models/tas_land_predGSAT_v3/1895-06.RData")



library(RColorBrewer)
error.col = colorRampPalette((brewer.pal(n = 9, name = "Blues")))(99)
map.col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(99)
Oranges = colorRampPalette((brewer.pal(n = 9, name = "Oranges")))(99)

obs <- unc <- bias <- raster.template
values(obs) <- c(matrix(CRUTEM5_[46,], 72, 36)[,36:1])
values(unc) = c(matrix(sqrt(CRUTEM5_sampling_unc_[46,]^2), 72, 36)[,36:1])
values(bias) = c(matrix(rowSds(sapply(X = CRUTEM5_ENS_anom_, FUN=function(x) x[46,])), 72, 36)[,36:1])

beta_land_gta = rep(NA, 72*36); beta_land_gta_ = raster.template
beta_land_gta[CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$grid.ix] = CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$mod_gta
values(beta_land_gta_) = c(matrix(beta_land_gta, 72, 36)[,36:1])


### CONTINUE HERE WITH MAKING FULL PLOT:
### ------------------------------------
setwd("/net/h2o/climphys1/sippels/_projects/ocean_cold_anomaly/figures/00_evaluation/")

# 1A. Compare map of regression coefficients:      
png(file = "02_illustration_maps_1895-06.png", width = 8, height = 5, units = "in", res = 300)
par(mar=c(1.5,1,2,1), mfrow=c(2,2), oma = c(0.5, 0,0,0))

# Betas
zlim. = round(max(c(quantile(values(beta_land_gta_), probs=0.98, na.rm=T), abs(quantile(values(beta_land_gta_), probs=0.02, na.rm=T)))), 2)
zlim = c(0, zlim.)
plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(beta_land_gta_, zlim = zlim), disagg = 1, zlim = zlim, legend.text = "", main = "Regression Coefficients", 
                        col = Oranges)

# Obs 06/1895 
zlim = c(-4, 4)
plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(obs, zlim = zlim), disagg = 1, zlim = zlim, legend.text = "", main = "", col = map.col)

# Uncertainties:
zlim = c(0, 1)
plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(unc, zlim = zlim), zlim = zlim, legend.text = "", main = "", col = error.col)

# Bias ensemble realization:
zlim = c(0, 1)
plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(bias, zlim = zlim), zlim = zlim, legend.text = "", main = "", col = error.col)

# plot_projected_worldmap(file.name = NULL, beta = adjust_raster_values(cur.betas$beta_land_p3_, zlim = zlim), zlim = zlim, legend.text = "Temperature coefficients (training CMIP6+3x unc+ 3x bias)")

dev.off()





cur.betas = get.beta_coefs(cmip6_mod = CMIP6.tas_land.GSAT.beta, mon.ix = mon.ix, date.ix = date.ix)
cur.title = paste(month.name[mon.ix], " ", c(1850:2020)[date.ix], sep="")
cur.file.name = paste("_", c(1850:2020)[date.ix], "-", formatC(mon.ix, width = 2, format = "d", flag = "0"), ".png", sep="")




