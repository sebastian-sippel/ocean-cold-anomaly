
## Sebastian Sippel
# 02.07.2024

# This plot is to look at monthly reconstructions
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures4review/")

# run: 04a_master_load_reconstructions.R



# decadal differences in the annual cycle in required corrections:
dec_diff_ann_cycle = matrix(NA, nrow = 17, ncol = 12)
dec_diff_ann_cycle_ua = matrix(NA, nrow = 17, ncol = 12)

for (i in 1:17) {
  ix = ((i-1) * 10 * 12 + 1):(i * 10 * 12)
  dec_diff_ann_cycle[i,] = rowMeans(matrix(OBS.tos$GMST_FM$mon$mod_p1_min[ix] - OBS.tas_land$GMST_FM$mon$mod_p1_min[ix], 
                  nrow = 12))
  
  #dec_diff_ann_cycle_ua[i,] = rowMeans(matrix(HadSST_ua.monthly$GMST_FM_mod_p1_min[ix] - OBS.tas_land$GMST_FM$mon$mod_p1_min[ix], 
  #                                         nrow = 12))
  dec_diff_ann_cycle_ua[i,] = rowMeans(matrix(OBS.tos$GMST_FM$mon$mod_p1_min[ix] - HadSST_ua.monthly$GMST_FM_mod_p1_min[ix], 
                                              nrow = 12))
}




# 02a: monthly temperature reconstructions:
pdf(file="02a_seasonality_GMST_reconstruction.pdf", width = 7, height = 5.5)
  par(mfrow=c(2,1), mar = c(2,5,1,1))
  
  ix = 241:(1440) # 
  years = (seq(1870, 1970, 1/12))[-1201]
    
  plot(x = years, y = OBS.tas_land$GMST_FM$mon$mod_p1_min[ix], type='l', col = "darkorange", 
       ylab = "Monthly Temperature \n Anomaly [°C]", xlab = "", ylim = c(-1, 0.5))
  lines(x = years, y = OBS.tos$GMST_FM$mon$mod_p1_min[ix], col = "darkblue")
  lines(x = years, y = HadSST_ua.monthly$GMST_FM_mod_p1_min[ix], col = "darkorchid")
  legend("bottomright", c("HadSST4-rec", "HadSST4-unadj-rec", "CRUTEM5-rec"), col = c("darkblue", "darkorchid", "darkorange"), cex = 0.8, lty = 1, inset = 0.02)
  
  # differences:
  plot(x = years, y = OBS.tos$GMST_FM$mon$mod_p1_min[ix] - OBS.tas_land$GMST_FM$mon$mod_p1_min[ix], 
       type='l', col = "black", 
       ylab = "Temperature Difference [°C], \n HadSST4-rec - CRUTEM5-rec", xlab = "", ylim = c(-1.2, 0.6))
  lines(x = years, y = HadSST_ua.monthly$GMST_FM_mod_p1_min[ix] - OBS.tas_land$GMST_FM$mon$mod_p1_min[ix], col = "darkorchid")
  
  lines(x = c(1,2000), y = c(0,0), col = "darkgray")
dev.off()




# 02b: monthly temperature reconstructions:
pdf(file="02b_seasonality_reconstruction.pdf", width = 7, height = 3.5)
  par(mfrow=c(1,2), mar = c(4, 4, 2, 0.5))
  
  plot(1:12, dec_diff_ann_cycle[1,1:12], type='n', ylim = c(-0.5, 0.5),
       ylab = "Difference by month [°C]", xlab ="Month", 
       main = c("HadSST4-rec - CRUTEM5-rec."))
  for (i in 3:5) {
    lines(1:12, dec_diff_ann_cycle[i,1:12], type='l')
  }
    lines(1:12, colMeans(dec_diff_ann_cycle[3:5,1:12]), type='l', lwd = 2)
    
  lines(1:12, dec_diff_ann_cycle[6,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle[7,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle[8,1:12], type='l', col = "red")
    lines(1:12, colMeans(dec_diff_ann_cycle[6:8,1:12]), type='l', lwd = 2, col = "red")
  
  for (i in 9:11) {
    lines(1:12, dec_diff_ann_cycle[i,1:12], type='l', col = "darkgray")
  }
    lines(1:12, colMeans(dec_diff_ann_cycle[9:11,1:12]), type='l', lwd = 2, col = "darkgray")

    legend("top", c("1870s, 1880s, 1890s", "1900s, 1910s, 1920s", "1930s, 1940s, 1950s"), 
           col = c("black", "red", "darkgray"), lty = 1, cex = 0.8, inset = 0.02)
    
    # Add differences as well:
    plot(1:12, dec_diff_ann_cycle[1,1:12], type='n', ylim = c(-0.5, 0.5),
         ylab = "Difference by month [°C], \n compared to prev. period", xlab ="Month")
    lines(1:12, colMeans(dec_diff_ann_cycle[6:8,1:12]) - colMeans(dec_diff_ann_cycle[3:5,1:12]), type='l', lwd = 2, col = "darkred")
    lines(1:12, colMeans(dec_diff_ann_cycle[9:11,1:12]) - colMeans(dec_diff_ann_cycle[6:8,1:12]), type='l', lwd = 2, col = "grey40")
    
    legend("top", c("1900-1930 vs. 1870-1900", "1930-1960 vs. 1900-1930"), 
           col = c("darkred", "grey40"), lty = 1, cex = 0.8, inset = 0.02)
    
dev.off()





# 02c: monthly temperature reconstructions:
pdf(file="02c_seasonality_reconstruction_HadSST4-ua.pdf", width = 7, height = 3.5)
  par(mfrow=c(1,2), mar = c(4, 4, 2, 0.5))
  
  plot(1:12, dec_diff_ann_cycle[1,1:12], type='n', ylim = c(-0.5, 0.5),
       ylab = "Difference by month [°C]", xlab ="Month",
       main = c("HadSST4-rec - CRUTEM5-rec."))
  for (i in 3:5) {
    lines(1:12, dec_diff_ann_cycle[i,1:12], type='l')
  }
  lines(1:12, colMeans(dec_diff_ann_cycle[3:5,1:12]), type='l', lwd = 2)
  
  lines(1:12, dec_diff_ann_cycle[6,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle[7,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle[8,1:12], type='l', col = "red")
  lines(1:12, colMeans(dec_diff_ann_cycle[6:8,1:12]), type='l', lwd = 2, col = "red")
  
  for (i in 9:11) {
    lines(1:12, dec_diff_ann_cycle[i,1:12], type='l', col = "darkgray")
  }
  lines(1:12, colMeans(dec_diff_ann_cycle[9:11,1:12]), type='l', lwd = 2, col = "darkgray")
  
  legend("top", c("1870s, 1880s, 1890s", "1900s, 1910s, 1920s", "1930s, 1940s, 1950s"), 
         col = c("black", "red", "darkgray"), lty = 1, cex = 0.8, inset = 0.02)
  
  # Add differences as well:
  plot(1:12, dec_diff_ann_cycle_ua[1,1:12], type='n', ylim = c(-0.5, 0.5),
       ylab = "Difference by month [°C]", xlab ="Month",
       main = "HadSST4-rec - HadSST4-unadj-rec")
  for (i in 3:5) {
    lines(1:12, dec_diff_ann_cycle_ua[i,1:12], type='l')
  }
  lines(1:12, colMeans(dec_diff_ann_cycle_ua[3:5,1:12]), type='l', lwd = 2)
  
  lines(1:12, dec_diff_ann_cycle_ua[6,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle_ua[7,1:12], type='l', col = "red")
  lines(1:12, dec_diff_ann_cycle_ua[8,1:12], type='l', col = "red")
  lines(1:12, colMeans(dec_diff_ann_cycle_ua[6:8,1:12]), type='l', lwd = 2, col = "red")
  
  for (i in 9:11) {
    lines(1:12, dec_diff_ann_cycle_ua[i,1:12], type='l', col = "darkgray")
  }
  lines(1:12, colMeans(dec_diff_ann_cycle_ua[9:11,1:12]), type='l', lwd = 2, col = "darkgray")
  
  legend("bottom", c("1870s, 1880s, 1890s", "1900s, 1910s, 1920s", "1930s, 1940s, 1950s"), 
         col = c("black", "red", "darkgray"), lty = 1, cex = 0.8, inset = 0.02)
  
  

dev.off()








