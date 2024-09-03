
# R functions to plot spatial patterns of fingerprints
# --------------------------------------------------------

# Sebastian Sippel
# 05.09.2019


require(rgdal)
require(fields)
require(rgeos)
#require(hydroGOF)
require(RColorBrewer)
require(psych) # for weighted correlation
library(gdalUtils)
library(sf)
library(terra)


# source("../../code/_convenience/convert.to.eurocentric.R")
# setwd("/net/h2o/climphys1/sippels/_projects/low_freq_anchor_v2/code/")
land.polygon <- readOGR("code/spatial_utils/shp_global110/", "110m_land_0_360")
land.polygon_eucentric <- readOGR("code/spatial_utils/shp_global110/", "110m_land")


{
  # Adjust raster values to zlim:
  adjust_raster_values <- function(cur.raster, zlim) {
    
    # check for zlim[1]:
    if (any(which(values(cur.raster) < zlim[1]))) {
      print(paste("Minimum raster value: ", min(values(cur.raster)), sep=""))
      values(cur.raster)[which(values(cur.raster) < zlim[1])] = zlim[1] - 10 ^ (-8)
    } 
    # check for zlim[2]:
    if (any(which(values(cur.raster) > zlim[2]))) {
      print(paste("Maximum raster value: ", max(values(cur.raster)), sep=""))
      values(cur.raster)[which(values(cur.raster) > zlim[2])] = zlim[2] + 10 ^ -8
    } 
    return(cur.raster)
  }
  
  
  plot_projected_worldmap <- function(file.name = F, beta, disagg = 1, col = map.col, zlim, useRaster = T, 
                                      legend.text = "Temperature coefficients", main = "", nlab = 5, legend = T, asp = 1, 
                                      points.to.add=NULL) {
    
    # 01. (Disaggregate) and Project Raster:
    # PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs +pm=180de +ellps=WGS84 +datum=WGS84"
    # PROJ <- "+proj=robin"
    PROJ_deshift <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs +pm=180de")
    PROJ = CRS("+proj=robin +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    PROJ1 <- CRS("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    beta = disaggregate(beta, fact = disagg)
    
    beta_ = beta
    extent(beta_) = extent(-180, 180, -90, 90)
    
    beta_proj = projectRaster(from=beta_, to=projectExtent(object=beta_, crs=PROJ))
    # beta_proj = projectRaster(from=beta, to=projectExtent(object=beta, crs=PROJ))
    # plot(beta_proj)
    
    # 02. Define surrounding polygon:
    temp = beta_
    values(temp) = 1
    pp <- rasterToPolygons(temp, na.rm=T, dissolve=T)
    pp_proj = spTransform(x = pp, CRSobj = PROJ)
    
    # 03. Cut raster to surrounding polygon:
    temp_proj = beta_proj
    values(temp_proj) = 1
    temp = rasterize(x = gDifference(spgeom1 =  rasterToPolygons(temp_proj, na.rm=T, dissolve=T),
                                     spgeom2 = pp_proj), y = beta_proj)
    values(beta_proj)[which(values(temp) == 1)] = NA
    
    # 04. Run plotting routine:
    leg.col.diff = abs(zlim[1]-zlim[2])/(length(col)-2)
    col.axis.lab = seq(zlim[1], zlim[2], length.out=nlab)
    col.axis.lab[1] = paste("< ", col.axis.lab[1], sep="")
    col.axis.lab[length(col.axis.lab)] = paste("> ", tail(col.axis.lab, 1), sep="")
    leg.const = 2  # Factor to add
    
    # par(mar=c(1,1,1,1), oma=c(0,0,0,0), mfrow=c(1,1))
    
    if (!is.null(file.name)) {
      png(file = file.name, width = 9, height = 6, units = "in", res = 150)
      par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    } 
    
    plot(beta_proj, axes=F, box=FALSE, add = F, horizontal=T, legend=F,
         zlim = c(zlim[1]-leg.col.diff, zlim[2]+leg.col.diff), interpolate=FALSE,
         col = col, colNA = "gray92", useRaster = T, main = main, breaks =  c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), asp = asp)
    # axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab))) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    mtext(legend.text, side=1, line=-0.5)
    plot(temp, add =T, horizontal=F, legend=F, useRaster=T, col="white", axes=F, box=FALSE, interpolate = F)  # temp simply color grid cells white that lie ouside of rectangle, 20190616

    
    # get the shape around the plot:
    # lines(spTransform(gridlines(x = beta, easts = seq(-180, 180, 60), norths = seq(-90, 90, 45)), CRSobj = PROJ), col = "darkgrey")
    options(warn=-1)
    land.polygon1=gBuffer(land.polygon, byid=TRUE, width=0)
    options(warn=0)
    plot(spTransform(crop(land.polygon1, extent(land.polygon)-0.0001 ), CRSobj = PROJ_deshift), add = T)
    
    # get the horizontal lines:
    hori.lines = (spLines(rbind(c(180, 0), c(-180, 0)), rbind(c(180, 45), c(-180, 45)), rbind(c(180, -45), c(-180, -45)), crs = proj4string(beta)))
    plot(spTransform(hori.lines, CRSobj = PROJ), add=T, col = "darkgrey")
    vert.lines =  spLines(cbind(-120, seq(-90, 90, by = 0.1)), cbind(-60, seq(-90, 90, by = 0.1)), cbind(0, seq(-90, 90, by = 0.1)), cbind(60, seq(-90, 90, by = 0.1)), cbind(120, seq(-90, 90, by = 0.1)), crs = proj4string(beta))
    lines(spTransform(vert.lines, CRSobj = PROJ), col = "darkgrey")
    plot(pp_proj, add=T, lwd = 2)
    
    if (!is.null(points.to.add)) plot(spTransform(SpatialPoints(points.to.add, proj4string = longlat), CRSobj = PROJ), pch = 16, cex = 0.2, add=T)
    
    # Add legend from image.plot to avoid skewed labels in raster::plot
    if (legend == T) {
    image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
               zlim = zlim, interpolate=FALSE,
               col = col, colNA = "gray95", useRaster = T, main = main, breaks = c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
               axis.args=list(labels = col.axis.lab, at = seq(zlim[1], zlim[2], length.out=nlab)),
               legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    }
    
    if (!is.null(file.name)) dev.off() 
    
    return(NULL)
  }
  
  
  
  plot_projected_worldmap_eucentric <- function(file.name = F, beta, disagg = 1, col = map.col, zlim, useRaster = T, 
                                      legend.text = "Temperature coefficients ", main = "", nlab = 5, legend = T, asp = 1, 
                                      points.to.add=NULL) {
    
    # 01. (Disaggregate) and Project Raster:
    # PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs +pm=180de +ellps=WGS84 +datum=WGS84"
    # PROJ <- "+proj=robin"
    PROJ = ("+proj=robin +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    PROJ1 <- ("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    beta = disaggregate(beta, fact = disagg)
    
    beta_ = beta
    extent(beta_) = extent(-180, 180, -90, 90)
    
    # beta_proj = projectRaster(from=beta_, to=projectExtent(object=beta_, crs=PROJ))
    beta_proj = raster(project(x=rast(beta_), y=PROJ))
    
    # plot(beta_proj)
    
    # 02. Define surrounding polygon:
    temp = beta_
    values(temp) = 1
    pp <- rasterToPolygons(temp, na.rm=T, dissolve=T)
    pp_proj = spTransform(x = pp, CRSobj = PROJ)
    
    # 03. Cut raster to surrounding polygon:
    temp_proj = beta_proj
    values(temp_proj) = 1
    temp = rasterize(x = gDifference(spgeom1 =  rasterToPolygons(temp_proj, na.rm=T, dissolve=T),
                                     spgeom2 = pp_proj), y = beta_proj)
    values(beta_proj)[which(values(temp) == 1)] = NA
    
    # 04. Run plotting routine:
    leg.col.diff = abs(zlim[1]-zlim[2])/(length(col)-2)
    col.axis.lab = seq(zlim[1], zlim[2], length.out=nlab)
    col.axis.lab[1] = paste("< ", col.axis.lab[1], sep="")
    col.axis.lab[length(col.axis.lab)] = paste("> ", tail(col.axis.lab, 1), sep="")
    leg.const = 2  # Factor to add
    
    # par(mar=c(1,1,1,1), oma=c(0,0,0,0), mfrow=c(1,1))
    
    if (!is.null(file.name)) {
      png(file = file.name, width = 3.5, height = 2.333, units = "in", res = 150, pointsize = 5)
      par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    } 
    
    # Create a mask for NA values
    beta_proj_NA <- beta_proj  # Create a logical mask where NA values are TRUE
    values(beta_proj_NA)[values(is.na(beta_proj))] <- -9999  # Create a logical mask where NA values are TRUE
    
    # Create a custom color palette including NA color
    custom_colors <- c("gray93", col)  # 'gray92' for NA values, followed by your regular color palette
    
    # Modify the breaks to include an extra range for NA values
    extended_breaks <- c(zlim[1] - leg.col.diff * leg.const, 
                         seq(zlim[1], zlim[2], length.out = length(col) - 1), 
                         zlim[2] + leg.col.diff * leg.const)
    
    # Plot the raster with the custom palette
    plot(beta_proj, 
         axes = FALSE, 
         box = FALSE, 
         add = FALSE, 
         horizontal = TRUE, 
         legend = FALSE,
         zlim = c(zlim[1] - leg.col.diff, zlim[2] + leg.col.diff), 
         interpolate = FALSE,
         col = col,  # Use the custom color palette
         colNA = "gray93",  # Set the NA color directly
         useRaster = T, 
         main = main, 
         breaks = extended_breaks,  # Use the extended breaks
         asp = asp)
    
    # axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab))) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    mtext(legend.text, side=1, line=-0.3)
    plot(temp, add =T, horizontal=F, legend=F, useRaster=T, col="white", axes=F, box=FALSE, interpolate = F)  # temp simply color grid cells white that lie ouside of rectangle, 20190616
    

    lines(spTransform(land.polygon_eucentric, CRSobj = PROJ), add=T, col = "grey", lwd = 0.6)
    # plot(spTransform(land.polygon, CRSobj = PROJ), add = T)
    
    hori.lines = (spLines(rbind(c(180, 0), c(-180, 0)), rbind(c(180, 45), c(-180, 45)), rbind(c(180, -45), c(-180, -45)), crs = longlat))
    plot(spTransform(hori.lines, CRSobj = PROJ), add=T, col = "darkgrey", lwd = 0.6)
    vert.lines =  spLines(cbind(-120, seq(-90, 90, by = 0.1)), cbind(-60, seq(-90, 90, by = 0.1)), cbind(0, seq(-90, 90, by = 0.1)), cbind(60, seq(-90, 90, by = 0.1)), cbind(120, seq(-90, 90, by = 0.1)), crs = longlat)
    lines(spTransform(vert.lines, CRSobj = PROJ), col = "darkgrey", lwd = 0.6)
    plot(pp_proj, add=T, lwd = 1.5)
    
    if (!is.null(points.to.add)) plot(spTransform(SpatialPoints(points.to.add, proj4string = longlat), CRSobj = PROJ), pch = 16, cex = 0.2, add=T)
    
    # Add legend from image.plot to avoid skewed labels in raster::plot
    if (legend == T) {
      image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
                 zlim = zlim, interpolate=FALSE,
                 col = col, colNA = "gray95", useRaster = F, main = main, breaks = c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
                 axis.args=list(labels = col.axis.lab, at = seq(zlim[1], zlim[2], length.out=nlab)),
                 legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    }
    
    if (!is.null(file.name)) dev.off() 
    
    return(NULL)
  }
  

  plot_projected_worldmap_eucentric_DIFF <- function(file.name = F, beta, disagg = 1, col = map.col, zlim, useRaster = T, 
                                                legend.text = "Temperature coefficients ", main = "", nlab = 5, legend = T, asp = 1, 
                                                raster.to.add) {
    
    # 01. (Disaggregate) and Project Raster:
    # PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs +pm=180de +ellps=WGS84 +datum=WGS84"
    # PROJ <- "+proj=robin"
    PROJ = ("+proj=robin +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    PROJ1 <- ("+proj=robin +lon_0=180 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    beta = disaggregate(beta, fact = disagg)
    
    beta_ = beta
    extent(beta_) = extent(-180, 180, -90, 90)
    
    #beta_proj = projectRaster(from=beta_, to=projectExtent(object=beta_, crs=PROJ))
    beta_proj = raster(project(x=rast(beta_), y=PROJ))
    # beta_proj = projectRaster(from=beta, to=projectExtent(object=beta, crs=PROJ))
    # plot(beta_proj)
    
    # 02. Define surrounding polygon:
    temp = beta_
    values(temp) = 1
    pp <- rasterToPolygons(temp, na.rm=T, dissolve=T)
    pp_proj = spTransform(x = pp, CRSobj = PROJ)
    
    # 03. Cut raster to surrounding polygon:
    temp_proj = beta_proj
    values(temp_proj) = 1
    temp = rasterize(x = gDifference(spgeom1 =  rasterToPolygons(temp_proj, na.rm=T, dissolve=T),
                                     spgeom2 = pp_proj), y = beta_proj)
    values(beta_proj)[which(values(temp) == 1)] = NA
    
    # 04. Run plotting routine:
    leg.col.diff = abs(zlim[1]-zlim[2])/(length(col)-2)
    col.axis.lab = seq(zlim[1], zlim[2], length.out=nlab)
    col.axis.lab[1] = paste("< ", col.axis.lab[1], sep="")
    col.axis.lab[length(col.axis.lab)] = paste("> ", tail(col.axis.lab, 1), sep="")
    leg.const = 2  # Factor to add
    
    # par(mar=c(1,1,1,1), oma=c(0,0,0,0), mfrow=c(1,1))
    
    if (!is.null(file.name)) {
      png(file = file.name, width = 9, height = 6, units = "in", res = 150)
      par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    } 
    
    # Create a mask for NA values
    beta_proj_NA <- beta_proj  # Create a logical mask where NA values are TRUE
    values(beta_proj_NA)[values(is.na(beta_proj))] <- -9999  # Create a logical mask where NA values are TRUE

    # Create a custom color palette including NA color
    custom_colors <- c("gray93", col)  # 'gray92' for NA values, followed by your regular color palette
    
    # Modify the breaks to include an extra range for NA values
    extended_breaks <- c(-10000, zlim[1] - leg.col.diff * leg.const, 
                         seq(zlim[1], zlim[2], length.out = length(col) - 1), 
                         zlim[2] + leg.col.diff * leg.const)
    
    # Plot the raster with the custom palette
    plot(beta_proj_NA, 
         axes = FALSE, 
         box = FALSE, 
         add = FALSE, 
         horizontal = TRUE, 
         legend = FALSE,
         zlim = c(zlim[1] - leg.col.diff, zlim[2] + leg.col.diff), 
         interpolate = FALSE,
         col = custom_colors,  # Use the custom color palette
         colNA = "gray96",  # Set the NA color directly
         useRaster = FALSE, 
         main = main, 
         breaks = extended_breaks,  # Use the extended breaks
         asp = asp)
    
    
    # plot(beta_proj, axes=F, box=FALSE, add = F, horizontal=T, legend=F,
    #     zlim = c(zlim[1]-leg.col.diff, zlim[2]+leg.col.diff), interpolate=FALSE,
    #     col = col, colNA = "gray92", useRaster = F, main = main, breaks =  c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), asp = asp)
  
    # axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab))) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    mtext(legend.text, side=1, line=-0.3)
    plot(temp, add =T, horizontal=F, legend=F, useRaster=T, col="white", axes=F, box=FALSE, interpolate = F)  # temp simply color grid cells white that lie ouside of rectangle, 20190616
    
    lines(spTransform(land.polygon_eucentric, CRSobj = PROJ), add=T, col = "grey", lwd = 0.6)
    # plot(spTransform(land.polygon, CRSobj = PROJ), add = T)
    
    hori.lines = (spLines(rbind(c(180, 0), c(-180, 0)), rbind(c(180, 45), c(-180, 45)), rbind(c(180, -45), c(-180, -45)), crs = longlat))
    plot(spTransform(hori.lines, CRSobj = PROJ), add=T, col = "darkgrey", lwd = 0.6)
    vert.lines =  spLines(cbind(-120, seq(-90, 90, by = 0.1)), cbind(-60, seq(-90, 90, by = 0.1)), cbind(0, seq(-90, 90, by = 0.1)), cbind(60, seq(-90, 90, by = 0.1)), cbind(120, seq(-90, 90, by = 0.1)), crs = longlat)
    lines(spTransform(vert.lines, CRSobj = PROJ), col = "darkgrey", lwd = 0.6)
    plot(pp_proj, add=T, lwd = 1.5)

    if (!is.null(raster.to.add)) {
      
      ix1 = which(values(!is.na(raster.to.add)))
      points.to.add.coord = coordinates(raster.to.add)[ix1,]
      points.to.add.val = values(raster.to.add)[ix1]
      
      for (i in 1:length(points.to.add.val)) {
        col.seq = seq(zlim[1], zlim[2], length.out = length(col))
        col.ix = which.min(abs(points.to.add.val[i] - col.seq))
        point.to.add = SpatialPoints(data.frame(x = points.to.add.coord[i,1], y = points.to.add.coord[i,2]))
        points(spTransform(SpatialPoints(point.to.add, proj4string = longlat), CRSobj = PROJ), bg =  col[col.ix], col = "grey40", pch = 21, cex = 0.8, lwd = 0.2)
      }
    }
    
    
    # Add legend from image.plot to avoid skewed labels in raster::plot
    if (legend == T) {
      image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
                 zlim = zlim, interpolate=FALSE,
                 col = col, colNA = "gray95", useRaster = T, main = main, breaks = c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
                 axis.args=list(labels = col.axis.lab, at = seq(zlim[1], zlim[2], length.out=nlab)),
                 legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    }
    
    if (!is.null(file.name)) dev.off() 
    
    return(NULL)
  }
  
  
    
  
  
  
  plot_projected_fingerprint <- function(beta, disagg = 2, col = map.col, zlim, useRaster = T, 
                                         legend.text = "Temperature coefficients ", main = "", file.name, png = T, nlab = 5,
                                         points.to.add=NULL) {
    
    # 01. (Disaggregate) and Project Raster:
    PROJ_deshift <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs +pm=180de"
    PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs"
    beta = disaggregate(beta, fact = disagg)
    extent(beta)=extent(c(-180,180,-90,90)) # shift beta to standard extent for robinson projection
    
    beta_proj = projectRaster(from=beta, to=projectExtent(object=beta, crs=PROJ))
    
    # 02. Define surrounding polygon:
    temp = beta
    values(temp) = 1
    pp <- rasterToPolygons(temp, na.rm=T, dissolve=T)
    pp_proj = spTransform(x = pp, CRSobj = PROJ)
    
    # 03. Cut raster to surrounding polygon:
    temp_proj = beta_proj
    values(temp_proj) = 1
    temp = rasterize(x = gDifference(spgeom1 =  rasterToPolygons(temp_proj, na.rm=T, dissolve=T),
                                     spgeom2 = pp_proj), y = beta_proj)
    values(beta_proj)[which(values(temp) == 1)] = NA
    
    # 04. Run plotting routine:
    if (png == T) {
      png(file = file.name, width = 9, height = 6, units = "in", res = 150)
    } else {
      pdf(file = file.name, width = 9, height = 6)
    }
    leg.col.diff = abs(zlim[1]-zlim[2])/(length(col)-2)
    col.axis.lab = seq(zlim[1], zlim[2], length.out=nlab)
    col.axis.lab[1] = paste("< ", col.axis.lab[1], sep="")
    col.axis.lab[length(col.axis.lab)] = paste("> ", tail(col.axis.lab, 1), sep="")
    leg.const = 2  # Factor to add
    
    par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    plot(beta_proj, axes=F, box=FALSE, add = F, horizontal=T, legend=F,
         zlim = c(zlim[1]-leg.col.diff, zlim[2]+leg.col.diff), interpolate=FALSE,
         col = col, colNA = "gray92", useRaster = T, main = main, breaks =  c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
         axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab))) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    mtext(legend.text, side=1, line=-0.5)
    plot(temp, add =T, horizontal=T, legend=F, useRaster=T, col="white", axes=F, box=FALSE, interpolate = F)  # temp simply color grid cells white that lie ouside of rectangle, 20190616
    
    # add gridlines:
    lines(spTransform(gridlines(x = beta, easts = seq(-180, 180, 60), norths = seq(-90, 90, 45)), CRSobj = PROJ), col = "darkgrey")
    options(warn=-1)
    land.polygon1=gBuffer(land.polygon, byid=TRUE, width=0)
    options(warn=0)
    plot(spTransform(crop(land.polygon1, extent(land.polygon)-0.0001 ), CRSobj = PROJ_deshift), add = T)
    plot(pp_proj, add=T, lwd = 2)
    
    if (!is.null(points.to.add)) plot(spTransform(SpatialPoints(points.to.add, proj4string = longlat), CRSobj = PROJ), pch = 16, cex = 0.2, add=T)
    
    # Add legend from image.plot to avoid skewed labels in raster::plot
    #image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
    #           zlim = zlim, interpolate=FALSE,
    #           col = col, colNA = "gray92", useRaster = T, main = main, breaks = seq(zlim[1], zlim[2], length.out = length(col)+1), 
    #           axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab)),
    #           legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    
    image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
               zlim = zlim, interpolate=FALSE,
               col = col, colNA = "gray92", useRaster = T, main = main, breaks = c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
               axis.args=list(labels = col.axis.lab, at = seq(zlim[1], zlim[2], length.out=nlab)),
               legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    
    dev.off()
    
    return(NULL)
  }
    
  

  
  
  
  
  plot_fingerprint_noproj <- function(beta, disagg = 1, col = map.col, zlim, useRaster = T, 
                                              legend.text = "Temperature coefficients ", main = "", file.name, png = T, nlab = 5) {
    
    # 01. (Disaggregate) and Project Raster:
    beta = disaggregate(beta, fact = disagg)
    beta_proj = beta
    
    temp = beta; values(temp) = 1; pp <- rasterToPolygons(temp, na.rm=T, dissolve=T)
    
    # 04. Run plotting routine:
    leg.col.diff = abs(zlim[1]-zlim[2])/(length(col)-2)
    col.axis.lab = seq(zlim[1], zlim[2], length.out=nlab)
    col.axis.lab[1] = paste("< ", col.axis.lab[1], sep="")
    col.axis.lab[length(col.axis.lab)] = paste("> ", tail(col.axis.lab, 1), sep="")
    leg.const = 2  # Factor to add
    
    # par(mar=c(1,1,1,1), oma=c(0,0,0,0))
    plot(beta_proj, axes=F, box=FALSE, add = F, horizontal=T, legend=F,
         zlim = c(zlim[1]-leg.col.diff, zlim[2]+leg.col.diff), interpolate=FALSE,
         col = col, colNA = "gray95", useRaster = T, main = main, breaks =  c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const),
         asp = 1.3)
    
    # axis.args=list(labels = seq(zlim[1], zlim[2], length.out=nlab), at = seq(zlim[1], zlim[2], length.out=nlab))) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    mtext(legend.text, side=1, line=1)
    
    lines(land.polygon)
    lines(pp, lwd = 2)
    
    # Define different gridlines:
    hori.lines = list(st_linestring(rbind(c(360, 0), c(0, 0))), st_linestring(rbind(c(360, 45), c(0, 45))), st_linestring(rbind(c(0, -45), c(360, -45))))
    sapply(hori.lines, FUN=function(x) lines(x, col = "darkgrey"))
    
    vert.lines = list(st_linestring(rbind(c(60, -90), c(60, 90))), rbind(c(120, -90), c(120, 90)), rbind(c(180, -90), c(180, 90)), rbind(c(240, -90), c(240, 90)), rbind(c(300, -90), c(300, 90)) )
    sapply(vert.lines, FUN=function(x) lines(x, col = "darkgrey"))

    # Add legend from image.plot to avoid skewed labels in raster::plot
    image.plot(axes=F, box=F, add = F, horizontal=T, legend.only = T,
               zlim = zlim, interpolate=FALSE,
               col = col, colNA = "gray95", useRaster = T, main = main, breaks = c(zlim[1]-leg.col.diff*leg.const, seq(zlim[1], zlim[2], length.out = length(col)-1), zlim[2]+leg.col.diff*leg.const), 
               axis.args=list(labels = col.axis.lab, at = seq(zlim[1], zlim[2], length.out=nlab)),
               legend.width = 0.7, legend.shrink = 0.6) # round(seq(zlim[1], zlim[2], length.out = length(col)+1), 4))   # colNA argument can be switched off
    
    return(NULL)
  }
  
  
  
  
  # Adjust raster values to zlim:
  adjust_raster_values <- function(cur.raster, zlim) {
    
    # check for zlim[1]:
    if (any(which(values(cur.raster) < zlim[1]))) {
      print(paste("Minimum raster value: ", min(values(cur.raster)), sep=""))
      values(cur.raster)[which(values(cur.raster) < zlim[1])] = zlim[1] - 10 ^ (-8)
    } 
    # check for zlim[2]:
    if (any(which(values(cur.raster) > zlim[2]))) {
      print(paste("Maximum raster value: ", max(values(cur.raster)), sep=""))
      values(cur.raster)[which(values(cur.raster) > zlim[2])] = zlim[2] + 10 ^ -8
    } 
    return(cur.raster)
  }
  
}

