

# ---------------------------------------------------------------
# Address Reviewer 1, Comment 8
# ---------------------------------------------------------------

# Sebastian Sippel
# 02.07.2024


# setwd to project folder:
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/")

# 00.(a) load  respective functions & code:
source("code/_convenience/frenchcolormap.R")
source("code/_functions_CMIP6.R")

# get plotting for fingerprints:
source("code/_convenience/_plot_projected_worldmap_v2.R")


# 00.(c) Load GMST reconstructions:
load("data/03_processed_CMIP6_4evaluation/CMIP6.tas_land.df_v5.RData")
load("data/01_processed4train_CRU/CRUTEM5_mon6.RData")
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
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly/figures4review/")

library(ggplot2)
library(ggsci)
library(gridExtra)



# investigate for time step June 1895:
beta_p0 = rep(NA, 72*36); beta_p0_ = raster.template
beta_p0[CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$grid.ix] = CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$mod_p0[,2]
values(beta_p0_) = c(matrix(beta_p0, 72, 36)[,36:1])

beta_p1 = rep(NA, 72*36); beta_p1_ = raster.template
beta_p1[CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$grid.ix] = CMIP6.tas_land.df$mon[[6]]$beta$GMST_FM[[46]]$mod_p1[,1]
values(beta_p1_) = c(matrix(beta_p1, 72, 36)[,36:1])

# check bias values:
hist(values(bias))


# Define vector of bias quantiles:
quant_vec = quantile(values(bias), probs = c(0.25, 0.5, 0.75), na.rm=T)

ix1 = which(values(bias) < quant_vec[1])
ix2 = which(values(bias) > quant_vec[1] & values(bias) < quant_vec[2])
ix3 = which(values(bias) > quant_vec[2] & values(bias) < quant_vec[3])
ix4 = which(values(bias) > quant_vec[3])


ratio_L2 = rep(NA, 4)
ratio_L2[1] = sum(values(beta_p1_)[ix1]^2) / sum(values(beta_p0_)[ix1]^2)
ratio_L2[2] = sum(values(beta_p1_)[ix2]^2) / sum(values(beta_p0_)[ix2]^2)
ratio_L2[3] = sum(values(beta_p1_)[ix3]^2) / sum(values(beta_p0_)[ix3]^2)
ratio_L2[4] = sum(values(beta_p1_)[ix4]^2) / sum(values(beta_p0_)[ix4]^2)

L2_biasincl = c(sum(values(beta_p1_)[ix1]^2), sum(values(beta_p1_)[ix2]^2), sum(values(beta_p1_)[ix3]^2), sum(values(beta_p1_)[ix4]^2))
L2_nobias = c(sum(values(beta_p0_)[ix1]^2), sum(values(beta_p0_)[ix2]^2), sum(values(beta_p0_)[ix3]^2), sum(values(beta_p0_)[ix4]^2))


# Prepare data for both barplots:
data1 <- data.frame(
  Category = rep(c("L2 norm (no bias)", "L2 norm (bias incl.)"), each = 4),
  Subcategory = rep(c("Bias Q1", "Bias Q2", "Bias Q3", "Bias Q4"), times = 2),
  Values = c(L2_nobias, L2_biasincl)
)

data2 <- data.frame(
  Category = c("Q1", "Q2", "Q3", "Q4"),
  Values = c(ratio_L2)
)


# Generate plots:
plot1 <- ggplot(data1, aes(x = Category, y = Values, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "L2 norm for land-based recon. June 1895", x = "", y = "L2 norm") +
  theme_minimal() +
  scale_fill_manual(values = pal_npg("nrc")(4))

# Save the plot as a PDF file
ggsave("03a_L2norm.pdf", plot = plot1, device = "pdf", width = 5, height = 5)



plot2 <- ggplot(data2, aes(x = Category, y = Values)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "blue") +
  labs(title = "", x = "Bias Quartiles", y = "Relative L2 norm change") +
  theme_minimal()

# Save the plot as a PDF file
ggsave("03b_L2norm-ratio.pdf", plot = plot2, device = "pdf", width = 5, height = 5)



