## Proxy plots early 20th century for Sebastian
setwd("/net/h2o/climphys1/sippels/_projects/ocean-cold-anomaly")

source("figures/Extended_Data_Figures/Extended_Data_Figure6/R-Functions_early20.R")
library(RNetCDF)
library(RColorBrewer)
library(rnaturalearth)
library(sp)

#definitions ---------
ref.period<-c(1871,1890)
early20.period<-c(1901,1920)


# load the proxies and metadata------------------
# PAGES2k 2017 temperature proxies, highres-data averaged over Apr-March window, as used in Neukom et al. 2019
proxies<-read.ts("data/00_DATASET/paleoclimate/neukom_indiv_proxies/proxy_ama_2.0.0.txt",sep="\t",header=T)

nproxies<-dim(proxies)[2]

# #interpolate gappy records
# proxies.int<-apply(proxies,2,function(x) approx(time(proxies),x,xout=min(which(!is.na(x))):max(which(!is.na(x)))))
# proxies.int<-lapply(proxies.int,function(x) ts(x$y,start=x$x[1]))
# proxies.int<-cbind.list(proxies.int)


#get metadata
meta<-as.matrix(read.table("data/00_DATASET/paleoclimate/neukom_indiv_proxies/metadata_2.0.0.txt",sep="\t"))
meta.names<-meta[,1]
meta<-meta[,-1]

lons.proxies<-as.numeric(meta[3,])
lats.proxies<-as.numeric(meta[2,])
coords.proxies<-cbind(lons.proxies,lats.proxies)

#select marine records
sel.marine<-which(meta[4,] %in% c("marine sediment","coral","sclerosponge","bivalve"))

##screened network from Neukom et al. 2019 based on regional+FDR screening (PAGES2k Consortium 2017)
proxies.fdr<-read.ts("data/00_DATASET/paleoclimate/neukom_indiv_proxies/proxy_ama_2.0.0_calib-selection_1881_1916_1995_0.67_infilled_DINEOF_PAGES-crit-regional+FDR.txt",sep="\t",header=T)

sel.fdr<-match(colnames(proxies),colnames(proxies.fdr))
length(which(!is.na(sel.fdr)))

selind.fdr<-which(!is.na(sel.fdr))

#remove Ocean2kHR_007 as it has no values in the ref.period
selind.fdr<-selind.fdr[-177]
proxies.fdr<-proxies.fdr[,-177]

selind.fdr.marine<-which(!is.na(match(selind.fdr,sel.marine)))

# instrumental data (Target from Neukom et al. 2019; Apr-March annual window)-------------------
instrfile<-"HadCRUT4.3_GraphEM_SP80_18502014_Apr-Mar_corr.nc"
nc<-open.nc(instrfile)  
lats<-var.get.nc(nc,1)
lons<-var.get.nc(nc,0)
instr<-var.get.nc(nc,3)
close.nc(nc)
lons.2d<-as.vector(lons)
lats.2d<-as.vector(lats)
lons.grid<-lons
lats.grid<-lats
lons<-unique(lons.2d)
lats<-unique(lats.2d)

nlon<-length(lons)
nlat<-length(lats)

ncells<-length(lons.2d)

lonsx<-lons
lonsx[which(lons>180)]<-lons[which(lons>180)]-360


lons.instr.proxies<-unlist(lapply(lons.proxies,function(x) which.min(abs(x-lonsx))))
lats.instr.proxies<-unlist(lapply(lats.proxies,function(x) which.min(abs(x-lats))))

# Proxy sign adjustment (so that correlation with local temp is positive for each proxy)------------
#import sign adjusment from matlab db (From PAGES2k 2017 database paper)
proxy_sign<-as.matrix(read.table("TSID_proxy_interpDirection_matlab_export.csv",sep=";",header=T))
smatch<-match(meta[6,],proxy_sign[,1])
reverse<-which(proxy_sign[smatch,4]=="negative")

##corrections:
#laguna chepical is already correct
reverse<-reverse[-which(reverse==686)]

proxies.sign<-proxies
for(pr in reverse){
  proxies.sign[,pr]<-proxies[,pr]*(-1)
}

#select adjusted fdr proxies
proxies.sign.fdr<-proxies.sign[,selind.fdr]


# Make the early 20 calculations --------------
#scale proxies to reference period
proxies.sign.fdr.ref<-scaletoperiod(proxies.sign.fdr,ref.period[1],ref.period[2])

#get anomalies during early 20th century
early20.mean<-apply(window(proxies.sign.fdr.ref,early20.period[1],early20.period[2]),2,mean,na.rm=T)

plot(early20.mean)

#check distributions
plot(density(early20.mean[selind.fdr.marine]))
densityline(density(early20.mean[-selind.fdr.marine]),col=2)

#write data----------
landmarine.fdr<-rep("terrestrial",209)
landmarine.fdr[selind.fdr.marine]<-"marine"
outdata<-rbind(meta[,selind.fdr],
               landmarine.fdr,
              round(early20.mean,3))
outdata.names<-c(meta.names,"land/marine","Anomaly 1901-1920 wrt 1871-1890")
write.table(cbind(outdata.names,outdata),file = "Proxies_Early20_anomalies.txt",quote = F,sep = ";",row.names = F,col.names = F)


#Plot definitions -----------
range.cols<-c(seq(-2,2,by=0.01))
bwr<-colorRampPalette(c('dark blue','white','dark red'), space = "Lab")
allcols.early20.mean<-bwr(length(range.cols)-1)

cols.early20.mean<-unlist(lapply(early20.mean,function(x) which.min(abs(x-range.cols))))

neganom.early20<-which(early20.mean<0)
neganom.early20.marine<-which(early20.mean[selind.fdr.marine]<0)
neganom.early20.terrestrial<-which(early20.mean[-selind.fdr.marine]<0)


#Plots ------------------
# sp::plot(ne_countries())
# points(lons.proxies[selind.fdr],lats.proxies[selind.fdr],pch=21,bg=allcols.early20.mean[cols.early20.mean],cex=2)
# points(lons.proxies[selind.fdr][selind.fdr.marine],lats.proxies[selind.fdr][selind.fdr.marine],pch=21,bg=allcols.early20.mean[cols.early20.mean[selind.fdr.marine]],cex=2)
# 
# points(lons.proxies[selind.fdr][neganom.early20],lats.proxies[selind.fdr][neganom.early20],pch=25,bg=allcols.early20.mean[cols.early20.mean[neganom.early20]],cex=2)
# points(lons.proxies[selind.fdr][-neganom.early20],lats.proxies[selind.fdr][-neganom.early20],pch=24,bg=allcols.early20.mean[cols.early20.mean[-neganom.early20]],cex=2)

png("Marine_proxies_negative_anomalies.png",width=300,height=200,res=200,units="mm")
sp::plot(ne_countries(returnclass = "sv"),axes=F)
points(lons.proxies[selind.fdr][selind.fdr.marine][neganom.early20.marine],lats.proxies[selind.fdr][selind.fdr.marine][neganom.early20.marine],pch=25,bg=allcols.early20.mean[cols.early20.mean[selind.fdr.marine][neganom.early20.marine]],cex=2)
dev.off()

png("Marine_proxies_positive_anomalies.png",width=300,height=200,res=200,units="mm")
sp::plot(ne_countries(returnclass = "sv"),axes=F)
points(lons.proxies[selind.fdr][selind.fdr.marine][-neganom.early20.marine],lats.proxies[selind.fdr][selind.fdr.marine][-neganom.early20.marine],pch=24,bg=allcols.early20.mean[cols.early20.mean[selind.fdr.marine][-neganom.early20.marine]],cex=2)
dev.off()

png("Terrestrial_proxies_negative_anomalies.png",width=300,height=200,res=200,units="mm")
sp::plot(ne_countries(returnclass = "sv"),axes=F)
points(lons.proxies[selind.fdr][-selind.fdr.marine][neganom.early20.terrestrial],lats.proxies[selind.fdr][-selind.fdr.marine][neganom.early20.terrestrial],pch=25,bg=allcols.early20.mean[cols.early20.mean[-selind.fdr.marine][neganom.early20.terrestrial]],cex=2)
dev.off()

png("Terrestrial_proxies_positive_anomalies.png",width=300,height=200,res=200,units="mm")
sp::plot(ne_countries(returnclass = "sv"),axes=F)
points(lons.proxies[selind.fdr][-selind.fdr.marine][-neganom.early20.terrestrial],lats.proxies[selind.fdr][-selind.fdr.marine][-neganom.early20.terrestrial],pch=24,bg=allcols.early20.mean[cols.early20.mean[-selind.fdr.marine][-neganom.early20.terrestrial]],cex=2)
dev.off()

# some additional tests -------------------
#check percentages
length(neganom.early20.marine);length(early20.mean[selind.fdr.marine])
length(neganom.early20.terrestrial);length(early20.mean[-selind.fdr.marine])

# some timeseries checks
order(early20.mean)
i<-62
plot(window(proxies.sign.fdr.ref[,i],1800,2020))
abline(h=0,lty=3)
abline(v=c(1871,1890))
abline(v=c(1901,1920),lty=2)
early20.mean[i]
