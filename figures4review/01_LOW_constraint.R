
### what do I need to do?

# is there a correlation between decadal decoupling and each climate model's land-ocean warming ratio?

# solution:
# -> make plot of: x-axis: CMIP6 warming ratio, y-axis: decadal decoupling across the 600 simulations.

get.trend(x = CMIP6.tos_all.df$ann$Y$GMST_FM[ix_tos], trend.years = trend.years, years = 1850:2014)

# get land-ocean warming ratios:
CMIP6.OW = apply(X = CMIP6_matrix$CMIP6.GMSST, MARGIN = 1, FUN=get.trend, trend.years = list(1980:2014), years = 1850:2020)
CMIP6.LW = apply(X = CMIP6_matrix$CMIP6.GMLSAT_NI, MARGIN = 1, FUN=get.trend, trend.years = list(1980:2014), years = 1850:2020)

CMIP6.LOW = CMIP6.LW / CMIP6.OW
hist(CMIP6.LOW)


## get the metric for decoupling:

decoupl = rowMeans(CMIP6_matrix$tos_mod_p1[,51:80]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,51:80])
decoupl_large = rbind(
  rowMeans(CMIP6_matrix$tos_mod_p1[,1:30]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,51:30]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,11:40]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,11:40]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,21:50]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,21:50]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,31:60]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,31:60]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,41:70]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,41:70]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,51:80]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,51:80]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,61:90]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,61:90]),
  rowMeans(CMIP6_matrix$tos_mod_p1[,71:100]) - rowMeans(CMIP6_matrix$tas_land_mod_p1[,71:100]))



# stratify by model: 
mod.un = unique(substr(CMIP6_matrix$ens.mem.un$mod, 1,3))

CMIP6.LOW_ = rep(NA, length(mod.un))
decoupl_ = rep(NA, length(mod.un))
decoupl_large_ = rep(NA, length(mod.un))
ens.mem.un_ = rep(NA, length(mod.un))

for (m in 1:length(mod.un)) {
  ix = which(substr(CMIP6_matrix$ens.mem.un$mod, 1, 3) == mod.un[m]) 
  # if (length(ix) < 3) next;
  
  CMIP6.LOW_[m] = mean(CMIP6.LOW[ix], na.rm=T)
  decoupl_[m] = sd(decoupl[ix], na.rm = T)
  decoupl_large_[m] = sd(decoupl_large[,ix], na.rm = T)
  ens.mem.un_
}





library(RColorBrewer)
# Define colors using RColorBrewer
num_colors <- 27
palette <- colorRampPalette(brewer.pal(9, "Set1"))(num_colors)
label_colors <- setNames(palette, unique(substr(ens.mem.un$mod, 1, 3)))



setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures4review/")


# Plotting the points with colors based on labels
png(filename = "01_LOW_constraint.png", width = 4, height = 8, units = "in", res = 600)
par(mfrow=c(3, 1), mar = c(5,5,1,1))

plot(data$CMIP6.LW, data$CMIP6.OW, col = label_colors[data$label], pch = 19,
     xlab = "Land warming [°C]", ylab = "SST Warming [°C]", main = "")

# Add a legend
legend("topleft", legend = names(label_colors), col = label_colors, pch = 19, ncol = 3, cex = 0.7)

plot(data$CMIP6.LOW, data$decoupl, col = label_colors[data$label], pch = 19,
     xlab = "Land-ocean warming ratio", ylab = "Ocean - Land Difference [°C], \n 1901-1930 average ", main = "")

plot(CMIP6.LOW_, decoupl_large_, col = label_colors, pch = 19, cex = 1.5,
     xlab = "Land-ocean warming ratio", ylab = "SD of (Ocean - Land) Difference", main = "Aggregated by model")

dev.off()





plot(CMIP6.LOW, (decoupl))
cor(CMIP6.LOW, abs(decoupl), use = "complete.obs")




