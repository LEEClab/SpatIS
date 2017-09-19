######################################################
#
# Linking individual specialization to space use
# for Sturnira lilium bats
#
# Patricia Rogeri - pa_bio04 at yahoo.com.br
# Bernardo Niebuhr - bernardo_brandaum at yahoo.com.br
# Renata Muylaert - remuylaert at gmail.com
#
# July 2017
# No copyrights - feel free to use, modify, and share
#######################################################

##########################################
# 0)  Preparing data
##########################################

# Loading pacakges
if(!require(adehabitatLT)) install.packages("adehabitatLT", dep=T); library(adehabitatLT)
if(!require(date)) install.packages("date", dep=T); library(date)
if(!require(adehabitatHR)) install.packages("adehabitatHR", dep=T); library(adehabitatHR)
if(!require(sp)) install.packages("sp", dep=T); library(sp)
if(!require(rgdal)) install.packages("rgdal", dep=T); library(rgdal)
if(!require(proj4)) install.packages("proj4", dep=T); library(proj4)
if(!require(RgoogleMaps)) install.packages("RgoogleMaps", dep=T); library(RgoogleMaps)
if(!require(png)) install.packages("png", dep=T); library(png)
if(!require(tcltk)) install.packages("tcltk", dep=T); library(tcltk)
if(!require(bipartite)) install.packages("bipartite", dep=T); library(bipartite)
if(!require(colortools)) install.packages("colortools", dep=T); library(colortools)

# Path to code folder
codedir <- "/home/leecb/Github/SpatIS/code"

# Path to data folder
datadir <- "/home/leecb/Github/SpatIS/data"

# Path to output folder
outputdir <- "/home/leecb/Github/SpatIS/output"

# Loading function to load Google Maps images
# Functions built by Luiz Gustavo Oliveira-Santos and Carlos Andre Zucco
# Change to the folder where he code is located
setwd(codedir)
source("plot_google_2_1_source_code.R")
source("spatis_source_code.R")

# Loading data
# Change to the data folder
setwd(datadir)

# Load land use map
# landuse.map <- readOGR("map_area/", layer = "map_saocarlos_wgs84_utm23S")
landuse.map <- readOGR("map_area/", layer = "map_cut_wgs84_utm23S")
landuse.map <- spTransform(landuse.map, CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))

# Load foraging areas map
resource.areas <- readOGR("map_area/", layer = "resource_areas_wgs84_utm23S")
resource.areas <- spTransform(resource.areas, CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))

# Load bat locations
data_ini <- read.table("bat_locations_ufscar_VHF.csv", header = T, sep = "\t", dec = ",")
data <- subset(data_ini) # We are not going to consider individuals 33 and 37
colnames(data) <- c("ID", "when_day", "when_hour", "x", "y", "position")
head(data,20)
str(data)

# Change from the VHF tag ID to the individual ID
n <- length(unique(data$ID))
data$tag_ID <- data$ID
indivs <- 1:n
data$ID <- ifelse(indivs < 10, paste(0, indivs, sep = ""), indivs)[as.factor(data$ID)]

#######################################
# Description of the data
# 
# ID: ID of bat individual
# when: day when the location was sampled, in the format ddmmYY
# hour: time of the day of the location, in the format hhmm
# x: longitude position, in UTM projection, zone 23S, datum WGS 84
# y: latitude position, in UTM projection, zone 23S, datum WGS 84
# position: whether the position of the animal was horizontal (H, flying), 
#           vertical (V, perched), and A (beginning of the movement)
# tag_ID: ID of the VHF tag used in the individual
#
#######################################

##########################
# Organizing data

# Coordinates
data$x <- round(data$x, 0)
data$y <- round(data$y, 0)

# Date
data$day <- ifelse(nchar(as.character(data$when_day)) == 5, paste("0", substr(as.character(data$when_day), 1,1), sep = ""), 
                   substr(as.character(data$when_day), 1,2))
data$month <- ifelse(nchar(as.character(data$when_day)) == 5, substr(as.character(data$when_day), 2,3), 
                     substr(as.character(data$when_day), 3,4))
data$year <- as.integer(ifelse(nchar(as.character(data$when_day)) == 5, substr(as.character(data$when_day), 4,5), 
                     substr(as.character(data$when_day), 5,6)))
data$year <- ifelse(data$year < 2000, 2000 + data$year, data$year)

# Date reorganized
data$date <- paste(data$year, data$month, data$day, sep="-")
# Julian day
data$julian <- mdy.date(as.numeric(data$month), as.numeric(data$day), data$year)

# Time of the day
data$hour <- substr(as.character(data$when_hour), 1,2)
data$minute <- substr(as.character(data$when_hour), 3,4)
data$second <- ifelse(nchar(as.character(data$when_hour)) == 4, "00", substr(as.character(data$when_hour), 5,6))

data$time <- ifelse(is.na(data$hour), NA, paste(data$hour, data$minute, data$second, sep=":"))

# We are going to work with the following data
data.f <- subset(data, select = c("ID", "x", "y", "date", "julian", "time", "tag_ID"))
head(data.f)

# Transforming data into Spatial Points
time <- as.POSIXct(paste(data.f$date, data.f$time), format="%Y-%m-%d %H:%M:%S")
spdados <- SpatialPointsDataFrame(data.f[,c(2,3)], data = data.f[,c(1,4:6)],
                                  proj4string=CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))
spdados <- spdados[order(spdados$ID, time),]  ### order data chronologically

##########################
# Extracting map information for each location

mapvalues <- over(spdados, landuse.map)
spdados$polygon_ID <- mapvalues$OBJECTID
spdados$land_use_class <- mapvalues$Class

areavalues <- over(spdados, resource.areas)
spdados$resource_area <- areavalues$ID

##########################
# Visualizing data

# Checking location of data
ids <- unique(data.f$ID)
(n <- length(unique(data.f$ID))) # number of individuals
cor <- rainbow(n)

# No background
plot(data.f$x, data.f$y, type = "n")
points(spdados, pch = 20, col = cor[as.factor(spdados$ID)])

# One at a time
for(i in 1:n){
  plot(data.f$x, data.f$y, type = "n")
  points(spdados[spdados$ID == ids[i],], pch = 20, col = cor[i])
  mtext(paste("ID = ", ids[i], sep = ""), side = 3, line = -1)
  Sys.sleep(2)
}

# With land use map as background
# nclass <- length(levels(landuse.map$Class))
# plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
levels(landuse.map$Class)
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols)
points(spdados, pch = 20, col = cor[as.factor(spdados$ID)])

# One at a time
for(i in 1:n){
  # plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
  plot(landuse.map, col = cols)
  points(spdados[spdados$ID == ids[i],], pch = 20, col = cor[i])
  mtext(paste("ID = ", ids[i], sep = ""), side = 3, line = 1)
  Sys.sleep(2)
}

# With resource areas as background
plot(resource.areas, border = "red", lwd = 2)
points(spdados, pch = 20, col = cor[as.factor(spdados$ID)])

# With land use map and foraging areas as background
# With land use map as background
# nclass <- length(levels(landuse.map$Class))
# plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
levels(landuse.map$Class)
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols, border = cols, lwd = 0.1)
points(spdados, pch = 20, col = cor[as.factor(spdados$ID)])
plot(resource.areas, border = "red", lwd = 2, add = T)

# Google maps image as background
plot.google(spdados, track=F, points = T, transpp = 0.5, pcol = cor[as.factor(spdados$ID)], cex=0.5) # google map image of the area

##########################################
# 1)  Assessing individual specialization using polygons of different
#     land use classes as resources
##########################################

# Output folder
setwd(outputdir)

##########################
# Calculating points por individual and polygon and organizing it as a network

# Aggregating number of points per individual per polygon
points.polygons.long <- aggregate(rep(1, nrow(spdados)), by = list(ID = spdados$ID, polygon = spdados$polygon_ID), FUN = sum)

# Transforming it into a matrix (network)
points.polygons.wide <- reshape(points.polygons.long, timevar = "polygon", idvar = "ID", direction = "wide")
points.polygons.wide <- points.polygons.wide[order(points.polygons.wide$ID),]
rownames(points.polygons.wide) <- points.polygons.wide$ID
points.polygons.wide <- points.polygons.wide[-1]

colnames(points.polygons.wide) <- substr(colnames(points.polygons.wide), 3,5)

# Filling NA's with zeros in the matrix
for(i in 1:nrow(points.polygons.wide)) {
  for(j in 1:ncol(points.polygons.wide)) {
    if(is.na(points.polygons.wide[i,j])) points.polygons.wide[i,j] <- 0
  }
}

freq.polygons.wide <- points.polygons.wide

# Frequency of points per individual in each polygon
tot.points.indiv <- rowSums(points.polygons.wide)
for(i in 1:nrow(freq.polygons.wide)) {
    freq.polygons.wide[i,] <- points.polygons.wide[i,]/tot.points.indiv[i]
}

# Exporting it to calculate individual specialization indices on DIETA
write.table(points.polygons.wide, "network_matrix_npoints_ids_polygons_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(points.polygons.wide, "network_matrix_npoints_ids_polygons_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(freq.polygons.wide, "network_matrix_freqpoints_ids_polygons_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(freq.polygons.wide, "network_matrix_freqpoints_ids_polygons_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)

###########################
# Plotting the network

points.polygons.decrease <- points.polygons.wide[order(rowSums(points.polygons.wide), decreasing = T), order(colSums(points.polygons.wide), decreasing = T)]

# As a bipartite network
plotweb(web = points.polygons.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Sites", 3)

# As a binary network
binary.points.polygons <- ifelse(points.polygons.decrease > 1, 1, 0)
bin.points.polygons.decrease <- binary.points.polygons[order(rowSums(binary.points.polygons), decreasing = T), order(colSums(binary.points.polygons), decreasing = T)]
plotweb(web = bin.points.polygons.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Sites", 3)

# As a matrix
cplace <- 1:ncol(points.polygons.decrease)-0.5
cnames <- colnames(points.polygons.decrease)
rplace <- 1:nrow(points.polygons.decrease)-0.5
rnames <- rownames(points.polygons.decrease)
visweb(as.matrix(points.polygons.decrease), type="none", square="interaction", prednames=F, preynames=F, boxes=F, clear = F, circles=TRUE, circle.col=1, circle.max=2.5, xlab="Plants")
text(cplace, rep(nrow(points.polygons.decrease)+1,ncol(points.polygons.decrease)), cnames, cex=0.5, )
text(rep(-1,nrow(points.polygons.decrease)), rplace, rev(rnames), cex=0.6)
mtext("Individuals", 2, line = 0.3)
mtext("Sites", 3, line = -3)

# Export
setwd(outputdir)
png("network_indiv_polygons.png", width = 15, height = 10, units = "cm", res = 300)
plotweb(web = points.polygons.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Sites", 3)
dev.off()


##########################################
# 2)  Assessing individual specialization using foraging sites
#     (considering kernel 50% of individuals and land use class edges)
##########################################

# Removing individual 06 (since the its home range did not reach stability - see below)
spdados.ud <- spdados[spdados$ID != "06",]
# Re-numbering individuals
spdados.ud$ID <- c(1:length(unique(spdados.ud$ID)))[as.factor(spdados.ud$ID)]
spdados.ud$ID <- ifelse(spdados.ud$ID < 10, paste(0, spdados.ud$ID, sep = ""), spdados.ud$ID)

# Output folder
setwd(outputdir)

##########################
# Calculating points por individual and polygon and organizing it as a network

# Aggregating number of points per individual per polygon
points.areas.long <- aggregate(rep(1, nrow(spdados.ud)), by = list(ID = spdados.ud$ID, area = spdados.ud$resource_area), FUN = sum)

# Transforming it into a matrix (network)
points.areas.wide <- reshape(points.areas.long, timevar = "area", idvar = "ID", direction = "wide")
points.areas.wide <- points.areas.wide[order(points.areas.wide$ID),]
rownames(points.areas.wide) <- points.areas.wide$ID
points.areas.wide <- points.areas.wide[-1]

colnames(points.areas.wide) <- substr(colnames(points.areas.wide), 3,5)

# Filling NA's with zeros in the matrix
for(i in 1:nrow(points.areas.wide)) {
  for(j in 1:ncol(points.areas.wide)) {
    if(is.na(points.areas.wide[i,j])) points.areas.wide[i,j] <- 0
  }
}

freq.areas.wide <- points.areas.wide

# Frequency of points per individual in each polygon
tot.points.area.indiv <- rowSums(points.areas.wide)
for(i in 1:nrow(freq.areas.wide)) {
  freq.areas.wide[i,] <- freq.areas.wide[i,]/tot.points.area.indiv[i]
}

# Exporting it to calculate individual specialization indices on DIETA
write.table(points.areas.wide, "network_matrix_npoints_ids_areas_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(points.areas.wide, "network_matrix_npoints_ids_areas_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(freq.areas.wide, "network_matrix_freqpoints_ids_areas_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(freq.areas.wide, "network_matrix_freqpoints_ids_areas_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)

###########################
# Plotting the network

points.areas.decrease <- points.areas.wide[order(rowSums(points.areas.wide), decreasing = T), order(colSums(points.areas.wide), decreasing = T)]

# As a bipartite network
plotweb(web = points.areas.decrease, method="cca", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)

# As a binary network
binary.points.areas <- ifelse(points.areas.decrease > 1, 1, 0)
bin.points.areas.decrease <- binary.points.areas[order(rowSums(binary.points.areas), decreasing = T), order(colSums(binary.points.areas), decreasing = T)]
plotweb(web = bin.points.areas.decrease, method="cca", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)

# As a matrix
cplace <- 1:ncol(points.areas.decrease)-0.5
cnames <- colnames(points.areas.decrease)
rplace <- 1:nrow(points.areas.decrease)-0.5
rnames <- rownames(points.areas.decrease)
visweb(as.matrix(points.areas.decrease), type="none", square="interaction", prednames=F, preynames=F, boxes=F, clear = F, circles=TRUE, circle.col=1, circle.max=2.5, xlab="Plants")
text(cplace, rep(nrow(points.areas.decrease)+1,ncol(points.areas.decrease)), cnames, cex=0.5, )
text(rep(-1,nrow(points.areas.decrease)), rplace, rev(rnames), cex=0.6)
mtext("Individuals", 2, line = -2)
mtext("Sites", 3, line = -1)

# Export
setwd(outputdir)
png("network_indiv_areas.png", width = 15, height = 10, units = "cm", res = 300)
plotweb(web = points.areas.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)
dev.off()

setwd(outputdir)
png("network_indiv_areas_binary.png", width = 15, height = 10, units = "cm", res = 300)
bin.points.areas.decrease <- binary.points.areas[order(rowSums(binary.points.areas), decreasing = T), order(colSums(binary.points.areas), decreasing = T)]
plotweb(web = bin.points.areas.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)
dev.off()

##########################################
# 3)  Assessing individual specialization using overlap between areas of use
#     (uilization distributions)
##########################################

##########################
# Accumulation curves - to check for sampling sufficiency

setwd(outputdir)

ids <- unique(spdados$ID)

### MCP 100%
### Parameters to control
cumHRmcp <- list() ## List with cumulative sample size for all individuals
for(i in 1:length(ids)){  ## loop for individuals
  temp <- spdados[which(spdados$ID == ids[i]),]
  temp <- SpatialPoints(coordinates(temp), CRS(proj4string(spdados)))
  cumulative <- vector()
  for(k in 5:length(temp)){  ##loop for sample size from 5 locations to all locations
    cumulative[k] <- mcp.area(temp[1:k,], percent=100, plotit=F)
  }  
  cumHRmcp[[i]] <- data.frame(hr=unlist(cumulative),ssize=5:length(temp))
}

names(cumHRmcp) <- ids
cumHRmcp

par(mar = c(5, 4, 4, 2) + 0.1)
# Seeing cummulative MCP area plots
for(i in 1:length(ids)) { ## plot all curves with 2 seconds interval
  plot(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, cex=0.5, pch=16, main=ids[i],
       xlab="Number of locations",ylab="MCP 100% area (ha)" )
  points(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, type="l", lwd=0.7, lty=2)
  Sys.sleep(2)
}

# To save all plots
# for(i in 1:length(ids)){
#   jpeg(paste("Cumulative_", ids[i], ".jpg", sep=""), width=20, height=20, units="cm", res = 300)
#   plot(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, cex=0.5, pch=16, main=ids[i],
#        xlab="Number of locations",ylab="MCP 100% area (ha)" )
#   points(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, type="l", lwd=0.7, lty=2)
#   dev.off() 
# }

### KERNEL 100%
cumHRkde <- list() ## List with cumulative sample size for all individuous
pb<-tkProgressBar(max=length(ids))
for(i in 1:length(ids)){  ## loop for individuals
  setTkProgressBar(pb, value=i, title="Kernel Estimation Progress", label=paste("Estimating",ids[i]))
  temp <- spdados[which(spdados$ID == ids[i]),]
  temp <- SpatialPoints(coordinates(temp), CRS(proj4string(spdados)))
  cumulative <- vector()
  for(k in 5:length(temp)){    			## loop for sample size from 5 locations to all locations
    href <- kernelUD(temp[1:k,], grid=200)@h$h  		## get reference smoother
    kde <- kernelUD(temp[1:k,], grid=200, h=href*0.8)  	## run kernel with 0.8*href
    cumulative[k] <- kernel.area(kde, percent=100)
  }  
  cumHRkde[[i]] <- data.frame(hr=na.omit(cumulative), ssize=5:length(temp))
}

names(cumHRkde) <- ids
close(pb)
cumHRkde

## Plot
for(i in 1:length(ids)){ ## plot all curves with 2 seconds interval
  plot(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,cex=0.5,pch=16,main=ids[i],
       xlab="Numero de Localizoes",ylab="Area do kernel 95%(ha)" )
  points(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,type="l",lwd=0.7,lty=2)
  Sys.sleep(2)
}

# INDIVIDUAL ID = 06 (tag_ID = 49) DID NOT REACH STABILITY ON ITS AREA OF USE - WE ARE GOING TO EXCLUDE IT FROM THE SPACE USE ANALYSES

# Removing individual 06 (since the its home range did not reach stability - see below)
spdados.ud <- spdados[spdados$ID != "06",]
# Re-numbering individuals
spdados.ud$ID <- c(1:length(unique(spdados.ud$ID)))[as.factor(spdados.ud$ID)]
spdados.ud$ID <- ifelse(spdados.ud$ID < 10, paste(0, spdados.ud$ID, sep = ""), spdados.ud$ID)

##########################
# Home range analysis

# Gathering all points to represent the population
pop <- spdados.ud
pop$ID <- "all" # ID all represents all population points
spdados.ud <- rbind(spdados.ud, pop)

ids <- unique(spdados.ud$ID)

# MCP

# MCP values
cpi <- mcp(spdados.ud[,1], percent=95)
areas <- mcp.area(spdados.ud[,1], percent=95, plotit = F)
cpi
summary(cpi)

plot(1:length(areas), areas, axes = F)
axis(2)
mtext("Area of use (MPC 95%, ha)", side = 2, line = 2)
axis(1, at = 1:11, names(areas))
mtext("Individual", side = 1, line = 2)

# Plotting only MCPs
cor <- c(rainbow(length(areas)-1), "black") # the last one is for the whole population
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols)
for(i in 1:length(areas)) {
  cp_id <- cpi[as.data.frame(cpi)[,1] == ids[i],]
  plot(cp_id, add = T, lwd = 2, border = cor[i])
}

# Plotting MCPs + points
cor <- c(rainbow(length(areas)-1), "black") # the last one is for the whole population
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols)
for(i in 1:length(areas)) {
  cp_id <- cpi[as.data.frame(cpi)[,1] == ids[i],]
  if(ids[i] != "all") points(spdados.ud[spdados.ud$ID == ids[i],], pch = 20, cex = 0.7, col = cor[i]) # points
  plot(cp_id, add = T, lwd = 2, border = cor[i]) # polygons
}

# Fixed Kernel
ids <- unique(spdados.ud$ID)
kernels <- list()
for(i in 1:length(ids)) {
  kudl <- adehabitatHR::kernelUD(spdados.ud[,1][spdados.ud$ID == ids[i],], h = "href")
  h = kudl[[1]]@h$h
  print(h)
  kudl <- adehabitatHR::kernelUD(spdados.ud[,1][spdados.ud$ID == ids[i],], h = 1*h, extent = 1.5, grid = 100)
  kernels[[i]] <- kudl
}

# Value of activity areas - kernel 95%
kern95 <- data.frame(-1, -1)[-1,]
colnames(kern95) <- c("id", "area")
homeranges <- list()
cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 95)
  homeranges[[i]] <- homerange
  kern95 <- rbind(kern95, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern95

# Value of core activity areas - kernel 50%
kern50 <- data.frame(-1, -1)[-1,]
colnames(kern50) <- c("id", "area")
homeranges <- list()
cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 50)
  homeranges[[i]] <- homerange
  kern50 <- rbind(kern50, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern50

# Export kernel 50%
setwd(outputdir)
tiff("Fig_S1_kernel50.tif", width = 15, height = 15, units = "cm", res = 300)
par(mar = c(0, 0, 0, 0))
kern50 <- data.frame(-1, -1)[-1,]
colnames(kern50) <- c("id", "area")
homeranges <- list()
cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 50)
  homeranges[[i]] <- homerange
  kern50 <- rbind(kern50, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern50
dev.off()

# Exporting values
area.use.vals <- cbind(kern95, kern50[2], t(areas))
colnames(area.use.vals) <- c("id", 'kern95', "kern50", "MCP95")

setwd(outputdir)
write.table(area.use.vals, "area_of_use_vals.txt", sep = "\t", row.names = F, col.names = T)

##########################
# Calculating the spatial use overlap
# SpatIS - Spatial individual specialization

calculated.spatis <- SpatIS(spdados.ud, individuals.col = "ID", population.ID = "all", grid = 200, extent = 1.5)

# Individual SpatIS
SpatIS.ind <- calculated.spatis$SpatIS.individual
# Population SpatIS
SpatIS.pop <- calculated.spatis$SpatIS.population

# Export individual SpatIS
# setwd(outputdir)
# write.table(SpatpIS.ind, "individual_spatis.csv", sep = "\t", row.names = F, col.names = T)

# Bootstrap - mixing points to assess if SpatIS is significantly greater than zero
randomized <- SpatIS_randomize(calculated.spatis, iterations = 999)

# Doing it manually
# naleat <- 100
# SpatIS_aleat <- c()
# for(i in 1:naleat)
# {
#   print(i)
#   spdados2 <- spdados.ud[,1][spdados$ID != "all",]
#   indivs <- sample(spdados2$ID)
#   spdados2$ID <- indivs
#   
#   pop <- spdados2
#   pop$ID <- "all"
#   spdados2 <- rbind(spdados2, pop)
#   
#   over <- kerneloverlap(spdados2[,1], grid = 100, extent = 1.5, method = "VI")
#   # Overlap of each individual with the whole population utilization distribution
#   SpatIS.ind.aux <- over[-population.line,population.line]
#   # Individual SpatIS = 1 - overlap of the individual with the population
#   SpatIS.ind <- 1 - SpatIS.ind.aux
#   # Population SpatIS = average of individual SpatIS
#   SpatIS <- mean(SpatIS.ind)
#   SpatIS_aleat <- c(SpatIS_aleat, SpatIS)
# }

# P-value
(p <- randomized$p)

# Bootstrap - mixing points to assess if SpatIS is significantly greater than zero
#             while keeping the roofing points
naleat <- 999
SpatIS_aleat_center <- SpatIS.pop
population.line <- 
for(i in 1:naleat)
{
  print(i)
  spdados2 <- spdados.ud[,1][spdados$ID != "all",]
  
  home_size <- 60 # in meters
  lines <- which(sqrt((coordinates(spdados2)[,1] - 203441)**2) < home_size & sqrt((coordinates(spdados2)[,2] - 7567795)**2) < home_size)
  home <- spdados2[lines,]
  spdados2 <- spdados2[-lines,]
  
  indivs <- sample(spdados2$ID)
  spdados2$ID <- indivs
  spdados2 <- rbind(home, spdados2)

  pop <- spdados2
  pop$ID <- "all"
  spdados2 <- rbind(spdados2, pop)
  
  over <- kerneloverlap(spdados2[,1], grid = 100, extent = 1.5, method = "VI")
  population.line <- which(rownames(over) == "all")
  # Overlap of each individual with the whole population utilization distribution
  SpatIS.ind.aux <- over[-population.line,population.line]
  # Individual SpatIS = 1 - overlap of the individual with the population
  SpatIS.ind <- 1 - SpatIS.ind.aux
  # Population SpatIS = average of individual SpatIS
  SpatIS <- mean(SpatIS.ind)
  SpatIS_aleat_center <- c(SpatIS_aleat_center, SpatIS)
}

# P-value
(p <- sum(SpatIS.pop <= SpatIS_aleat_center)/(naleat+1)) # proportion of random values that are greater than the observed value
hist(SpatIS_aleat_center, main = "", xlab = "Spatial Individual Specialization", ylab = "Frequency",
     xlim = c(min(SpatIS_aleat_center), SpatIS.pop))
abline(v = SpatIS.pop, col = "red")

# Export
setwd(outputdir)
# png("p_value_SpatIS_random.png", width = 10, height = 10, units = "cm", res = 300)
tiff("p_value_SpatIS_random.tif", width = 10, height = 10, units = "cm", res = 300)
hist(SpatIS_aleat_center[-1], main = "", xlab = "Spatial Individual Specialization", ylab = "Frequency",
     xlim = c(min(SpatIS_aleat_center), SpatIS.pop), breaks = 20)
abline(v = SpatIS.pop, col = "red")
dev.off()


##########################################
# 4)  Creating figure with results
##########################################

# Output folder
setwd(outputdir)

# png("Fig_results.png", width = 20, height = 20, units = "cm", res = 300)
tiff("Fig_results.tif", width = 20, height = 20, units = "cm", res = 300)
par(mfrow = c(2,2), mar = c(0, 0, 0, 0))

plot(0, 0, xlim = c(0,1), ylim = c(0,1), type = "n", axes = F)
#plot(landuse.map, col = "white", border = "white", lwd = 0.1)
mtext("A", side = 3, at = 0.1, line = -3, cex = 1.5)

# Figure A - network
# As a bipartite network
cols <- complementary("purple4", plot = F)
plotweb(web = points.areas.decrease, method="cca", text.rot=0, 
        col.interaction="grey65",
        col.high = cols[1], col.low=cols[2],
        bor.col.interaction ="grey65", bor.col.high=cols[1],
        bor.col.low=cols[2], labsize = 1.5)
mtext("Individuals", 1, line = 0)
mtext("Foraging sites", 3, line = 0)
mtext("B", 3, at = 0, cex = 1.5)

# Figure B - map + space use
# Valor das HR
kern95 <- data.frame(-1, -1)[-1,]
colnames(kern95) <- c("id", "area")
homeranges <- list()
cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 95)
  homeranges[[i]] <- homerange
  kern95 <- rbind(kern95, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    # x <- c(196250, 209856)
    # y <- c(7561021, 7574393)
    # plot(x, y, type = "n", bty = "n", axes = F)
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern95

mtext("C", 3, at = 198500, cex = 1.5)

# Figure C - map + space use (kernel 50%)
# Valor das HR
kern95 <- data.frame(-1, -1)[-1,]
colnames(kern95) <- c("id", "area")
homeranges <- list()
cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 50)
  homeranges[[i]] <- homerange
  kern50 <- rbind(kern50, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    # x <- c(196250, 209856)
    # y <- c(7561021, 7574393)
    # plot(x, y, type = "n", bty = "n", axes = F)
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern95

mtext("D", 3, at = 198500, cex = 1.5)

# legend("topright", legend = c())

# # Figure D - significance of SpatIS
# # P-value
# (p <- sum(SpatIS_real <= SpatIS_aleat_center)/naleat) # proportion of random values that are greater than the observed value
# hist(SpatIS_aleat_center, main = "", xlab = "Spatial Individual Specialization", ylab = "Frequency",
#      xlim = c(min(SpatIS_aleat_center), SpatIS_real))
# abline(v = SpatIS_real, col = "red")
# mtext("Spatial Individual Specialization", 1, line = 1.8)
# mtext("", 2, line = 1.8)
# mtext("C", 3, at = 0.2, cex = 1.5)

par(mfrow = c(1,1))
dev.off()

