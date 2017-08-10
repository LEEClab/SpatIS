##########################################################################
## PLOT GOOGLE 2.1									##
## Script for Plotting Animal Movement Data over GoogleMaps Imagery	## 
## Google Earth Image									##
## by Carlos Andre Zucco - September/2014 - cazucco14@gmail.com		##
## creatively based on previous work of Luiz Gustavo Oliveira-Santos	##
## No rights reserved. Feel free to modify and share				##
##########################################################################

################################################################################
##  USAGE												##
## plot.google(sp,zone,points=TRUE,transpp=0.3,cex=0.3,pcol="yellow",		##
##		   track=TRUE,lcol="white",transpl=0.3,lwd=0.3,add=F)				##
## plot.pol.google(spatpol,zone,,zone,col="red",transp,border=NA,			##
##			   lwd=0.3,scale=2,add=F)							##
##													##
## FUNCTION ARGUMENTS										##
##  sp - set of data in SpatialPoints or SpatialPointsDataFrame 		 	##
##	   with UTM projection									##
##  spatpol - set of data in SpatialPolygons or SpatialPolygonsDataFrame 	##
##	   with UTM projection									##
##  zone - UTM zone of data 									##
##  points - logical, whether to plot points of fixes					##
##  tracks - logical, whether to plot lines connecting points			##
##  transpp - transparence of plotted points						##
##  transpl - transparence of plotted lines				 		##
##  transp - transparence of polygon fill							##
##  cex - size of plotted points								##
##  lwd - width of plotted lines								##													##
##  pcol - color of points									##
##  lcol - color of lines									##
##  border - color of polygon border; 'NA' sets polygon with no border line	##
##  add 0 logical, indicating if plot should be over opened window 
##  VALUES												##
##  A plot of data over Google Maps satelite imagery 					##
##													##
################################################################################
## Version 2.0 Updates
# Plot in UTM CRS
# Plot of points, tracks or both
# Resolution of png image of 1280x1280 pixels
## Version 2.1 Updates
# 'add' argument for overploting points/tracks and polygons
### ATTENTION
# Be sure of full compatibility between 
# R version and packages RgoogleMaps' and 'png'


plot.google<- function(sp,points=TRUE,transpp=0.3,cex=0.3,pcol="yellow",
				track=TRUE,lcol="white",transpl=0.3,lwd=0.3,scale=2,add=F){
 require(sp)
 require(RgoogleMaps)
 require(rgdal)
 if(inherits(sp,"SpatialPoints")==FALSE){stop("sp argument should be a SpatialPoints or SpatialPointsDataFrame object")}
 if(is.na(proj4string(sp))){stop("CRS of sp argument should be assigned")}
 if(add==F){
## coordinate range of polygons 
 bboxUTM<-SpatialPoints(t(bbox(sp)),CRS(proj4string(sp)))
 bboxLL<-spTransform(bboxUTM,CRSobj=CRS("+proj=longlat +datum=WGS84"))
 satmap<-GetMap.bbox(latR=range(coordinates(bboxLL)[,2]),lonR=range(coordinates(bboxLL)[,1]),
 			   destfile ="MyTyle.png",maptype='satellite',SCALE=scale)
 tyle=readPNG("MyTyle.png")
 bbox=SpatialPoints(matrix(unlist(satmap$BBOX)[c(2,4,1,3)],2,2)
		,proj4string=CRS("+proj=longlat +datum=WGS84"))
 bbox=data.frame(spTransform(bbox,CRSobj=CRS(proj4string(sp))))
## Return data to UTM CRS format  ##
 lat<-coordinates(sp)[,2]
 long<-coordinates(sp)[,1]
## Plot Tyle and Data ##
par(mar=c(3,3,1,1))
plot(bbox[,1],bbox[,2],xaxs="i",yaxs="i",type="n",xlab="Longitude",ylab="Latitude",cex.axis=0.8)
tyle2<-rasterImage(tyle,bbox[1,1],bbox[1,2],bbox[2,1],bbox[2,2])
if(points==TRUE){points(long,lat,pch=16,col=adjustcolor(pcol,transpp),cex=cex)}
if(track==TRUE){points(long,lat,col=adjustcolor(lcol,transpl),type='l',lwd=lwd)}
 } # end add=T
 if(add==T){
lat<-coordinates(sp)[,2]
long<-coordinates(sp)[,1]
if(points==TRUE){points(long,lat,pch=16,col=adjustcolor(pcol,transpp),cex=cex)}
if(track==TRUE){points(long,lat,col=adjustcolor(lcol,transpl),type='l',lwd=lwd)}
} # end of add=F
}


plot.pol.google<- function(spatpol,col="red",transp=0.3,border=NA,lwd=0.3,scale=2,add=F){
 require(RgoogleMaps)
 require(rgdal)
 require(png)
 if(inherits(spatpol,"SpatialPolygons")==FALSE){stop("spatpol should be a SpatialPolygon or SpatialPolygonDataFrame object")}
 if(is.na(proj4string(spatpol))){stop("CRS of sp argument should be assigned")}
 if(add==F){
## coordinate range of polygons 
 bboxUTM<-SpatialPoints(t(bbox(spatpol)),CRS(proj4string(spatpol)))
 bboxLL<-spTransform(bboxUTM,CRSobj=CRS("+proj=longlat +datum=WGS84"))

## get PNG tyle from Google Maps and extract UTM bbox of tyle
 satmap<-GetMap.bbox(latR=range(coordinates(bboxLL)[,2]),lonR=range(coordinates(bboxLL)[,1]),
 			   destfile ="MyTyle.png",maptype='satellite',SCALE=scale)
 tyle=readPNG("MyTyle.png")
 bbox=SpatialPoints(matrix(unlist(satmap$BBOX)[c(2,4,1,3)],2,2)
		,proj4string=CRS("+proj=longlat +datum=WGS84"))
 bbox=data.frame(spTransform(bbox,CRSobj=CRS(proj4string(spatpol))))

## Plot Tyle and Data
 par(mar=c(3,3,1,1))
 plot(bbox[,1],bbox[,2],xaxs="i",yaxs="i",type="n",xlab="Longitude",ylab="Latitude",cex.axis=0.8)
 tyle2<-rasterImage(tyle,bbox[1,1],bbox[1,2],bbox[2,1],bbox[2,2])
 pols = unlist(lapply(spatpol@polygons, function(i)slot(i, 'Polygons')))
 pols = lapply(pols, function(pols)slot(pols, 'coords'))
 for(i in 1:length(pols)){
 polygon(pols[[i]],density = NULL,col=adjustcolor('red',transp),border=border,lwd=lwd)}
 } # end add=T
 if(add==T){
 pols = unlist(lapply(spatpol@polygons, function(i)slot(i, 'Polygons')))
 pols = lapply(pols, function(pols)slot(pols, 'coords'))
 for(i in 1:length(pols)){
 polygon(pols[[i]],density = NULL,col=adjustcolor('red',transp),border=border,lwd=lwd)}
 } #end of add=F
}


plot.google2<- function(sp,points=TRUE,transpp=0.3,cex=0.3,pcol="yellow",
                       track=TRUE,lcol="white",transpl=0.3,lwd=0.3,scale=2,add=F,
                       rangex=0, rangey=0){
  require(sp)
  require(RgoogleMaps)
  require(rgdal)
  if(inherits(sp,"SpatialPoints")==FALSE){stop("sp argument should be a SpatialPoints or SpatialPointsDataFrame object")}
  if(is.na(proj4string(sp))){stop("CRS of sp argument should be assigned")}
  if(add==F){
    ## coordinate range of polygons 
    bboxUTM<-SpatialPoints(t(bbox(sp)),CRS(proj4string(sp)))
    bboxLL<-spTransform(bboxUTM,CRSobj=CRS("+proj=longlat +datum=WGS84"))
    if(sum(rangey) == 0){
      satmap<-GetMap.bbox(latR=range(coordinates(bboxLL)[,2]),lonR=range(coordinates(bboxLL)[,1]),
                          destfile ="MyTyle2.png",maptype='satellite',SCALE=scale)
    } else {
      satmap<-GetMap.bbox(latR=rangey,lonR=range(coordinates(bboxLL)[,1]),
                          destfile ="MyTyle2.png",maptype='satellite',SCALE=scale)
    }
    tyle=readPNG("MyTyle2.png")
    bbox=SpatialPoints(matrix(unlist(satmap$BBOX)[c(2,4,1,3)],2,2)
                       ,proj4string=CRS("+proj=longlat +datum=WGS84"))
    bbox=data.frame(spTransform(bbox,CRSobj=CRS(proj4string(sp))))
    ## Return data to UTM CRS format  ##
    lat<-coordinates(sp)[,2]
    long<-coordinates(sp)[,1]
    ## Plot Tyle and Data ##
    par(mar=c(3,3,1,1))
    plot(bbox[,1],bbox[,2],xaxs="i",yaxs="i",type="n",xlab="Longitude",ylab="Latitude",cex.axis=0.8)
    tyle2<-rasterImage(tyle,bbox[1,1],bbox[1,2],bbox[2,1],bbox[2,2])
    if(points==TRUE){points(long,lat,pch=16,col=adjustcolor(pcol,transpp),cex=cex)}
    if(track==TRUE){points(long,lat,col=adjustcolor(lcol,transpl),type='l',lwd=lwd)}
  } # end add=T
  if(add==T){
    lat<-coordinates(sp)[,2]
    long<-coordinates(sp)[,1]
    if(points==TRUE){points(long,lat,pch=16,col=adjustcolor(pcol,transpp),cex=cex)}
    if(track==TRUE){points(long,lat,col=adjustcolor(lcol,transpl),type='l',lwd=lwd)}
  } # end of add=F
}