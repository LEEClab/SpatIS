######################################################
#
# Linking individual specialization to space use
# for Sturnira lilium bats
#
# Patricia Rogeri - pa_bio04 at yahoo.com.br
# Bernardo Niebuhr - bernardo_brandaum at yahoo.com.br
# Renata Muylaert - remuylaert at gmail.com
#
# Jun 2020
# No copyrights - feel free to use, modify, and share
#######################################################

##########################################
# 0)  Preparing data
##########################################

# Loading pacakges
if(!require(tidyverse)) install.packages("tidyverse", dep=T); library(tidyverse)
if(!require(date)) install.packages("date", dep=T); library(date)
if(!require(adehabitatHR)) install.packages("adehabitatHR", dep=T); library(adehabitatHR)
if(!require(sp)) install.packages("sp", dep=T); library(sp)
if(!require(rgdal)) install.packages("rgdal", dep=T); library(rgdal)
if(!require(proj4)) install.packages("proj4", dep=T); library(proj4)
if(!require(png)) install.packages("png", dep=T); library(png)
if(!require(tcltk)) install.packages("tcltk", dep=T); library(tcltk)
if(!require(bipartite)) install.packages("bipartite", dep=T); library(bipartite)
if(!require(colortools)) install.packages("colortools", dep=T); library(colortools)
if(!require(igraph)) install.packages("igraph", dep=T); library(igraph)

# Loading function to calculate SpatIS and SpatICS
# Change to the folder where he code is located
source("code/spatis_source_code_v1_0.R")

# Loading data

# Load land use map
landuse.map <- readOGR(dsn = "data/map_area", layer = "map_cut_wgs84_utm23S")
landuse.map <- spTransform(landuse.map, CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))

# Load foraging areas map
resource.areas <- readOGR("data/map_area", layer = "resource_areas_wgs84_utm23S")
resource.areas <- spTransform(resource.areas, CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))

# Load bat information
info <- read.table("data/spatis_paper_individual_bat_info.csv", header= T, sep = ",")

# Load bat locations
data.f <- read.table("data/spatis_paper_bat_positions.csv", header = T, sep = ",")

# Transforming data into Spatial Points
time <- as.POSIXct(paste(data.f$date, data.f$time), format="%Y-%m-%d %H:%M:%S")
spdados <- SpatialPointsDataFrame(data.f[,c(3,4)], data = data.f[,c(1, 2, 4:6)],
                                  proj4string=CRS("+proj=utm +datum=WGS84 +zone=23 +south +ellps=WGS84 +towgs84=0,0,0"))
spdados <- spdados[order(spdados$Animal_ID, time),]  ### order data chronologically

##########################
# Extracting map information for each location

mapvalues <- over(spdados, landuse.map)
spdados$polygon_ID <- mapvalues$OBJECTID
spdados$land_use_class <- mapvalues$Classe

areavalues <- over(spdados, resource.areas)
spdados$resource_area <- areavalues$ID

##########################
# Visualizing data

# Checking location of data
ids <- unique(spdados$Animal_ID)
(n <- length(unique(data.f$Animal_ID))) # number of individuals
cor <- rainbow(n)

# No background
plot(data.f$x, data.f$y, type = "n",
     ylab = 'y', xlab = 'x')
points(spdados, pch = 20, col = cor[as.factor(spdados$Animal_ID)])

# One at a time
for(i in 1:n){
  plot(data.f$x, data.f$y, type = "n",
       ylab = 'y', xlab = 'x')
  points(spdados[spdados$Animal_ID == ids[i],], pch = 20, col = cor[i])
  mtext(paste("ID = ", ids[i], sep = ""), side = 3, line = -1)
  Sys.sleep(1)
}

# With land use map as background
# nclass <- length(levels(landuse.map$Class))
# plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
levels(landuse.map$Class)
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols)
points(spdados, pch = 20, col = cor[as.factor(spdados$Animal_ID)])

# One at a time
for(i in 1:n){
  # plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
  plot(landuse.map, col = cols)
  points(spdados[spdados$Animal_ID == ids[i],], pch = 20, col = cor[i])
  mtext(paste("ID = ", ids[i], sep = ""), side = 3, line = 1)
  Sys.sleep(1)
}

# With resource areas as background
plot(resource.areas, border = "red", lwd = 2)
points(spdados, pch = 20, col = cor[as.factor(spdados$Animal_ID)])

# With land use map and foraging areas as background
# With land use map as background
# nclass <- length(levels(landuse.map$Class))
# plot(landuse.map, col = grey.colors(nclass)[landuse.map$Class])
levels(landuse.map$Class)
rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols, border = cols, lwd = 0.1)
points(spdados, pch = 20, col = cor[as.factor(spdados$Animal_ID)])
plot(resource.areas, border = "red", lwd = 2, add = T)


##########################################
# 1)  Assessing individual specialization using polygons of different
#     land use classes as resources
##########################################

# Removing individual 11 (Tag_ID == 49) 
# (since the its home range did not reach stability - see below)
spdados.ud <- spdados[spdados$Animal_ID != 11,]

# Compare old and new numeric codes
table(spdados$Animal_ID)
table(spdados.ud$Animal_ID)

##########################
# Calculating points por individual and polygon and organizing it as a network

# Aggregating number of points per individual per polygon
points.polygons.long <- aggregate(rep(1, nrow(spdados.ud)), 
                                  by = list(Animal_ID = spdados.ud$Animal_ID, polygon = spdados.ud$polygon_ID), 
                                  FUN = sum)

# Transforming it into a matrix (network)
points.polygons.wide <- reshape(points.polygons.long, timevar = "polygon", 
                                idvar = "Animal_ID", direction = "wide")
points.polygons.wide <- points.polygons.wide[order(points.polygons.wide$Animal_ID),]
rownames(points.polygons.wide) <- points.polygons.wide$Animal_ID
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
write.table(points.polygons.wide, "output/network_matrix_npoints_ids_polygons_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(points.polygons.wide, "output/network_matrix_npoints_ids_polygons_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(freq.polygons.wide, "output/network_matrix_freqpoints_ids_polygons_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(freq.polygons.wide, "output/network_matrix_freqpoints_ids_polygons_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)

# return to the long version to calculate proportion of use of each site by each individual,
# with no arcsin sqrt transformation
props <- points.polygons.long %>% 
  dplyr::rename(success = x) %>% 
  dplyr::arrange(Animal_ID) %>% 
  dplyr::group_by(Animal_ID) %>% 
  dplyr::mutate(sample.size = sum(success),
                failure = sample.size - success) %>% 
  dplyr::left_join(
    dplyr::select(info, Animal_ID, sex, body_mass_g)
  )
props

write.csv(props, file = "data/proportion_use_landcover_polygons.csv", row.names = F)

###########################
# Plotting the network

points.polygons.decrease <- points.polygons.wide[order(rowSums(points.polygons.wide), decreasing = T), order(colSums(points.polygons.wide), decreasing = T)]

# As a bipartite network
plotweb(web = points.polygons.decrease, method="cca", text.rot=90,
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
png("output/network_indiv_polygons.png", width = 15, height = 10, units = "cm", res = 300)
plotweb(web = points.polygons.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Sites", 3)
dev.off()

# unipartite graph
g <- graph.incidence(points.polygons.decrease, weighted = T)
g
V(g)$type
V(g)$name
E(g)$weight
plot(g, layout = layout_as_bipartite)
projected_g <- bipartite_projection(g, multiplicity = TRUE)
individuals.polygons <- projected_g$proj1

png('output/network_unipatite_individuals_polygons.png', width = 15, height = 15, units = "cm", res = 300)
plot(individuals.polygons, vertex.color = 'brown3', vertex.label.color = 'white',
     edge.color = 'grey50', edge.width = E(individuals.polygons)$weight**2/4, 
     layout = layout_with_graphopt(individuals.polygons, charge = 0.01))
dev.off()

# dynamic plot
#tkplot(individuals, edge.color = 'grey30', edge.width = E(individuals)$weight)

##########################################
# 2)  Assessing individual specialization using foraging sites
#     (considering kernel 50% of individuals and land use class edges)
##########################################

##########################
# Calculating points por individual and polygon and organizing it as a network

# Aggregating number of points per individual per polygon
points.areas.long <- aggregate(rep(1, nrow(spdados.ud)), by = list(Animal_ID = spdados.ud$Animal_ID, area = spdados.ud$resource_area), FUN = sum)

# Transforming it into a matrix (network)
points.areas.wide <- reshape(points.areas.long, timevar = "area", idvar = "Animal_ID", direction = "wide")
points.areas.wide <- points.areas.wide[order(points.areas.wide$Animal_ID),]
rownames(points.areas.wide) <- points.areas.wide$Animal_ID
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
write.table(points.areas.wide, "output/network_matrix_npoints_ids_areas_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(points.areas.wide, "output/network_matrix_npoints_ids_areas_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)
write.table(freq.areas.wide, "output/network_matrix_freqpoints_ids_areas_HEADER.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(freq.areas.wide, "output/network_matrix_freqpoints_ids_areas_NOHEADER.txt", sep = "\t", row.names = T, col.names = F, quote = F)

# return to the long version to calculate proportion of use of each site by each individual,
# with no arcsin sqrt transformation
props.areas <- points.areas.long %>% 
  dplyr::rename(success = x) %>% 
  dplyr::arrange(Animal_ID) %>% 
  dplyr::group_by(Animal_ID) %>% 
  dplyr::mutate(sample.size = sum(success),
                failure = sample.size - success) %>% 
  dplyr::left_join(
    dplyr::select(info, Animal_ID, sex, body_mass_g)
  )
props.areas

write.csv(props.areas, file = "data/proportion_use_foraging_sites.csv", row.names = F)

###########################
# Plotting the network

points.areas.decrease <- points.areas.wide[order(rowSums(points.areas.wide), decreasing = T), order(colSums(points.areas.wide), decreasing = T)]

# As a bipartite network
plotweb(web = points.areas.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)

# Transform the previous bipartite object into an igraph object

g <- graph.incidence(points.areas.decrease, weighted = T)
g
V(g)$type
V(g)$name
E(g)$weight
plot(g, layout = layout_as_bipartite)
projected_g <- bipartite_projection(g, multiplicity = TRUE)
individuals.areas <- projected_g$proj1

png('output/network_unipatite_individuals_areas.png', width = 15, height = 15, units = "cm", res = 300)
plot(individuals.areas, vertex.color = 'brown3', vertex.label.color = 'white',
     edge.color = 'grey50', edge.width = E(individuals.areas)$weight**2/2, 
     layout = layout_with_graphopt(individuals.areas, charge = 0.01))
dev.off()

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
text(cplace, rep(nrow(points.areas.decrease)+1,ncol(points.areas.decrease)), cnames, cex=0.5)
text(rep(-1,nrow(points.areas.decrease)), rplace, rev(rnames), cex=0.6)
mtext("Individuals", 2, line = -2)
mtext("Sites", 3, line = -1)

# Export
png("output/network_indiv_areas.png", width = 15, height = 10, units = "cm", res = 300)
plotweb(web = points.areas.decrease, method="normal", text.rot=90,
        col.interaction="grey50",
        col.high = "darkolivegreen3", col.low="brown3",
        bor.col.interaction ="grey50", bor.col.high="darkolivegreen3",
        bor.col.low="brown3")
mtext("Individuals", 1)
mtext("Foraging sites", 3)
dev.off()

png("output/network_indiv_areas_binary.png", width = 15, height = 10, units = "cm", res = 300)
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
#     (utilization distributions)
##########################################

##########################
# Accumulation curves - to check for sampling sufficiency

ids <- unique(spdados$Animal_ID)

### MCP 95%
### Parameters to control
cumHRmcp <- list() ## List with cumulative sample size for all individuals
for(i in 1:length(ids)){  ## loop for individuals
  temp <- spdados[which(spdados$Animal_ID == ids[i]),]
  temp <- SpatialPoints(coordinates(temp), CRS(proj4string(spdados)))
  cumulative <- vector()
  for(k in 5:length(temp)){  ##loop for sample size from 5 locations to all locations
    cumulative[k] <- mcp.area(temp[1:k,], percent=95, plotit=F)
  }  
  cumHRmcp[[i]] <- data.frame(hr=unlist(cumulative),ssize=5:length(temp))
}

names(cumHRmcp) <- ids
cumHRmcp

par(mar = c(5, 4, 4, 2) + 0.1)
# Seeing cummulative MCP area plots
for(i in 1:length(ids)) { ## plot all curves with 2 seconds interval
  plot(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, cex=0.5, pch=16, main=ids[i],
       xlab="Number of locations",ylab="MCP 95% area (ha)" )
  points(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, type="l", lwd=0.7, lty=2)
  Sys.sleep(1)
}

# To save each plot
# for(i in 1:length(ids)){
#   jpeg(paste("Cumulative_", ids[i], ".jpg", sep=""), width=20, height=20, units="cm", res = 300)
#   plot(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, cex=0.5, pch=16, main=ids[i],
#        xlab="Number of locations",ylab="MCP 95% area (ha)" )
#   points(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, type="l", lwd=0.7, lty=2)
#   dev.off()
# }

# Plot for each individual (removing the 11th which did not reach an asymptote)
png('output/accumulation_curves_MCP.png', width = 15, height = 15, units = 'cm', res = 600)
par(mfrow = c(4,3), mar = c(2,3,3,1) + 0.1, oma = c(3, 3, 0, 0))
for(i in 1:length(ids)){ ## plot all curves with 2 seconds interval
  if(i != 11) {
    ind_num <- ifelse(i < 6, i, i-1)
    plot(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, cex=0.5, pch=16, main=ids[i],
         xlab="Number of locations",ylab="MCP 95% area (ha)" )
    points(cumHRmcp[[i]]$hr ~ cumHRmcp[[i]]$ss, type="l", lwd=0.7, lty=2)
  }
}
mtext(text = "Number of locations", side = 1, outer = T, line = 1)
mtext(text = "Space use area (95% MCP, in ha)", side = 2, outer = T, line = 1)

dev.off()

### KERNEL 99%
# reorder data based on distance to roost
roost.coords <- c(203441, 7567795)
distance <- sqrt((coordinates(spdados)[,1] - roost.coords[1])**2 + (coordinates(spdados)[,2] - roost.coords[2])**2)
spdados <- spdados[order(distance),]

cumHRkde <- list() ## List with cumulative sample size for all individuous
pb<-tkProgressBar(max=length(ids))
for(i in 1:length(ids)){  ## loop for individuals
  setTkProgressBar(pb, value=i, title="Kernel Estimation Progress", label=paste("Estimating",ids[i]))
  temp <- spdados[which(spdados$Animal_ID == ids[i]),]
  temp <- SpatialPoints(coordinates(temp), CRS(proj4string(spdados)))
  cumulative <- vector()
  for(k in 5:length(temp)){    			## loop for sample size from 5 locations to all locations
    href <- kernelUD(temp[1:k,], grid=200)@h$h  		## get reference smoother
    kde <- kernelUD(temp[1:k,], grid=200, h=href)#*0.8)  	## run kernel with 0.8*href
    cumulative[k] <- kernel.area(kde, percent=99)
  }  
  cumHRkde[[i]] <- data.frame(hr=na.omit(cumulative), ssize=5:length(temp))
}

names(cumHRkde) <- ids
close(pb)
cumHRkde

## Plot
for(i in 1:length(ids)){ ## plot all curves with 2 seconds interval
  plot(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,cex=0.5,pch=16,main=ids[i],
       xlab="Numero de Localizoes",ylab="Area do kernel 99%(ha)" )
  points(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,type="l",lwd=0.7,lty=2)
  Sys.sleep(1)
}

# Plot for each individual
png('output/accumulation_curves.png', width = 15, height = 15, units = 'cm', res = 600)
par(mfrow = c(4,3), mar = c(2,3,3,1) + 0.1, oma = c(3, 3, 0, 0))
for(i in 1:length(ids)){ ## plot all curves with 2 seconds interval
  if(i != 6) {
    ind_num <- ifelse(i < 6, i, i-1)
    plot(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,cex=0.5,pch=16,main=ind_num,
         xlab="",ylab="")
    points(cumHRkde[[i]]$hr~cumHRkde[[i]]$ss,type="l",lwd=0.7,lty=2)
  }
}
mtext(text = "Number of locations", side = 1, outer = T, line = 1)
mtext(text = "Space use area (99% KDE, in ha)", side = 2, outer = T, line = 1)

dev.off()

# test with KDE
iterations <- 20

### KDE 95%
### Parameters to control
KDE_percentage <- 95
cumHRkde <- list() ## List with cumulative sample size for all individuals
for (i in 1:length(ids)) {  ## loop for individuals
  print(i)
  cumHRkde[[i]] <- list()
  for(j in 1:iterations) { # Loop for iterations
    temp <- spdados[which(spdados$Animal_ID == ids[i]),]
    # Here we do not use the points as they are, but randomize their order
    temp <- SpatialPoints(coordinates(temp)[sample(length(temp)),], CRS(proj4string(spdados)))
    cumulative <- vector()
    for(k in 5:length(temp)){  ##loop for sample size from 5 locations to all locations
      UD <- kernelUD(temp[1:k,])
      cumulative[k] <- kernel.area(UD, percent = KDE_percentage)
    }  
    cumulative <- cumulative[5:length(cumulative)]
    cumHRkde[[i]][[j]] <- data.frame(hr = unlist(cumulative), ssize = 5:length(temp))
  }
}

names(cumHRkde) <- ids
cumHRkde

# Plotting
par(mar = c(5, 4, 4, 2) + 0.1)
# Seeing cummulative kde area plots
for(i in 1:length(ids)) { ## plot all curves with 2 seconds interval
  cum.df <- as.data.frame(cumHRkde[[i]])
  cumHR.df.mult <- cum.df[,c(2,seq(1, ncol(cum.df), 2))]
  
  cumHR.df.mult2 <- data.frame(npoins = rep(cumHR.df.mult[,1], iterations), HR = unlist(c(cumHR.df.mult[,2:ncol(cumHR.df.mult)])))
  plot(cumHR.df.mult2$npoins, cumHR.df.mult2$HR, cex=0.5, pch=16, col = 'grey', main=ids[i],
       xlab = "Number of locations",ylab = paste0("KDE ", KDE_percentage, "% area (ha)"))
  points(unique(cumHR.df.mult2$npoins), apply(cumHR.df.mult[,2:ncol(cumHR.df.mult)], MARGIN = 1, FUN = median), type="l", lwd=3, lty=1)
  
  Sys.sleep(1)
}

# INDIVIDUAL ID = 11 (tag_ID = 49) DID NOT REACH STABILITY ON ITS AREA OF USE - WE ARE GOING TO EXCLUDE IT FROM THE SPACE USE ANALYSES

# Removing individual 11 (since the its home range did not reach stability - see below)
spdados.ud <- spdados[spdados$Animal_ID != 11,]
spdados.ud <- spdados.ud[order(spdados.ud$Animal_ID),]

# Save object
saveRDS(spdados.ud, file = "code/bat_spdados2.RDS")

##########################
# Home range analysis

# Gathering all points to represent the population
pop <- spdados.ud
pop$Animal_ID <- "all" # ID all represents all population points
spdados.ud <- rbind(spdados.ud, pop)

ids <- unique(spdados.ud$Animal_ID)

# MCP

# MCP values
cpi <- mcp(spdados.ud[, "Animal_ID"], percent = 95)
areas <- mcp.area(spdados.ud[, "Animal_ID"], percent = 95, plotit = F)
cpi
summary(cpi)

plot(1:length(areas), areas, axes = F, xlab = "", ylab = "")
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
  if(ids[i] != "all") points(spdados.ud[spdados.ud$Animal_ID == ids[i],], pch = 20, cex = 0.7, col = cor[i]) # points
  plot(cp_id, add = T, lwd = 2, border = cor[i]) # polygons
}

# Fixed Kernel
ids <- unique(spdados.ud$Animal_ID)
kernels <- list()
for(i in 1:length(ids)) {
  kudl <- adehabitatHR::kernelUD(spdados.ud[, "Animal_ID"][spdados.ud$Animal_ID == ids[i],], h = "href")
  h = kudl[[1]]@h$h
  print(h)
  kudl <- adehabitatHR::kernelUD(spdados.ud[, "Animal_ID"][spdados.ud$Animal_ID == ids[i],], h = 1*h, extent = 1.5, grid = 100)
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

# Export kde50 as a shapefile
for(i in 1:length(homeranges)) {
  rgdal::writeOGR(homeranges[[i]], "output", paste0("kde50_h10_", as.character(homeranges[[i]]$id)),
                   driver = "ESRI Shapefile")
  
}

# spdf <- raster::union(homeranges[[1]], homeranges[[2]])
# for(i in 3:(length(homeranges)-1)) {
#   spdf <- raster::union(spdf, homeranges[[i]])
# }


# Export activity points
rgdal::writeOGR(spdados.ud[spdados.ud$Animal_ID != "all",-3], dsn = "output", layer = 'bat_locations', driver = "ESRI Shapefile")

# Export kernel 50%
tiff("output/Fig_S1_kernel50.tif", width = 15, height = 15, units = "cm", res = 300)
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

write.table(area.use.vals, "output/area_of_use_vals.txt", sep = "\t", row.names = F, col.names = T)

##########################
# Calculating the spatial use overlap
# SpatIS - Spatial individual specialization

calculated.spatis <- SpatIS(spdados.ud, individuals.col = "Animal_ID", 
                            population.ID = "all", grid = 200, extent = 1.5)

# calculated.spatics <- SpatICS(spdados.ud, individuals.col = "ID", population.ID = "all", grid = 200, extent = 1.5)

# Individual SpatIS
SpatIS.ind <- calculated.spatis$SpatIS.individual
# Population SpatIS
SpatIS.pop <- calculated.spatis$SpatIS.population

# Individual SpatIS
SpatICS.ind <- calculated.spatis$SpatICS.individual
# Population SpatIS
SpatICS.pop <- calculated.spatis$SpatICS.population

# Export individual SpatIS
# setwd(outputdir)
# write.table(SpatpIS.ind, "individual_spatis.csv", sep = "\t", row.names = F, col.names = T)

# Bootstrap - mixing points to assess if SpatIS is significantly greater than zero
randomized <- SpatIS.randomize(calculated.spatis, iterations = 1000)

# Let's look at the results
str(randomized)
# Individual values
randomized$SpatIS.individual.random
randomized$SpatIS.individual.observed
# Population value
randomized$SpatIS.population.random
randomized$SpatIS.population.observed
# Significance
randomized$SpatIS.significance
# Power analysis
randomized$SpatIS.power
randomized$SpatIS.power.curve

# keeping the nest
spdados2 <- spdados.ud[, "Animal_ID"][spdados.ud$Animal_ID != "all",]

nest.coords <- c(203441, 7567795)
home_size <- 60 # in meters
lines <- which(sqrt((coordinates(spdados2)[,1] - nest.coords[1])**2) < home_size & sqrt((coordinates(spdados2)[,2] - nest.coords[2])**2) < home_size)
spdados2$nest <- 0
spdados2$nest[lines] <- 1

calculated.spatis2 <- SpatIS(spdados2, individuals.col = "Animal_ID", 
                             population.ID = NULL, grid = 200, extent = 1.5)
randomized.nest <- SpatIS.randomize(calculated.spatis2, iterations = 1000, 
                                    not.randomize.col = "nest", not.randomize.val = 1)

# SpatIS
randomized.nest[[1]]$SpatIS.population.random %>% mean

# SpatICS
randomized.nest[[2]]$SpatICS.population.random %>% mean
randomized.nest[[2]]$SpatICS.power

##########################################
# 4)  Creating figure with results
##########################################

# Output folder

# Fig. 3 from the manuscript
# png("Fig_results.png", width = 20, height = 20, units = "cm", res = 300)
tiff("output/Fig_results_base.tif", width = 20, height = 30, units = "cm", res = 300)
png("output/Fig_3_results_base.png", width = 20, height = 30, units = "cm", res = 300)
par(mfrow = c(3,2), mar = c(0, 0, 0, 0))

# Figure B - map + space use
# Valor das HR
kern95 <- data.frame(-1, -1)[-1,]
colnames(kern95) <- c("id", "area")
homeranges <- list()
cor <- c(RColorBrewer::brewer.pal(10, 'Set3'), 'black')#c(rainbow(length(ids)-1), "black") # the last one is for the whole population
for(i in (1:length(ids))) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 97)
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
    raster::scalebar(5000, xy = c(200300, 7562500), divs = 4, type = 'bar', label = c(0, '2.5', '5 km'))
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
#legend('topright', legend = 1:10, col = cor, lwd = 2)
kern95

mtext("A", 3, at = 199250, line = -3, cex = 1.5)

# legend

plot(0, 0, xlim = c(0,1), ylim = c(0,1), type = "n", axes = F)
legend(x = 0.6, y = 1, legend = 1:10, col = cor, lwd = 3, title = "Individuals", cex = 1.7)
legend(x = 0, y = 1, legend = levels(landuse.map$Classe), fill = cols, border = cols,
       title = "Land use class", cex = 1.7)
#plot(landuse.map, col = "white", border = "white", lwd = 0.1)
# mtext("A", side = 3, at = 0.1, line = -3, cex = 1.5)

# As a bipartite network - land use polygons
cols <- cols[landuse.map$Classe[landuse.map$OBJECTID %in% colnames(points.polygons.decrease)]]

plotweb(web = points.polygons.decrease, method="normal", text.rot=0, 
        col.interaction="grey65",
        col.high = cols, col.low=cor[as.numeric(row.names(points.polygons.decrease))],
        arrow = 'down.center',
        bor.col.interaction ="grey65", bor.col.high=cols[1],
        bor.col.low=cols[2], labsize = 2, 
        high.lablength = 0)
mtext("Individuals", 1, line = 0)
mtext("Habitat polygons", 3, line = 0)

# cols <- complementary("purple4", plot = F)
# plotweb(web = points.polygons.decrease, method="normal", text.rot=0, 
#         col.interaction="grey65",
#         col.high = cols[1], col.low=cols[2],
#         arrow = 'down.center',
#         bor.col.interaction ="grey65", bor.col.high=cols[1],
#         bor.col.low=cols[2], labsize = 1.5, 
#         high.lablength = 0)
# mtext("Individuals", 1, line = 0)
# mtext("Habitat polygons", 3, line = 0)

mtext("B", 3, at = 0, cex = 1.5)


# Figure D - unipartite network - land use polygons
cols <- complementary("purple4", plot = F)
plot(individuals.polygons, vertex.color = cor[as.numeric(row.names(points.polygons.decrease))], #vertex.color = cols[2], vertex.label.color = 'white',
     #vertex.label.
     edge.color = 'grey50', edge.width = E(individuals.polygons)$weight**2/5, 
     layout = layout_with_graphopt(individuals.polygons, charge = 0.01))

mtext("C", 3, at = -1.7, cex = 1.5)


# Figure E - as a bipartite network - foraging sites
cols <- complementary("purple4", plot = F)
plotweb(web = points.areas.decrease, method="normal", text.rot=0, 
        col.interaction="grey65",
        col.high = rgb(1, 0, 0, alpha = 0.3), col.low=cor[as.numeric(row.names(points.areas.decrease))],
        bor.col.interaction ="grey65", bor.col.high='red',
        arrow = 'down.center',
        bor.col.low=cor[as.numeric(row.names(points.areas.decrease))], labsize = 2)
mtext("Individuals", 1, line = 0)
mtext("Foraging sites", 3, line = 0)

mtext("D", 3, at = 0, cex = 1.5)


# Figure F - unipartite network - foraging polygons
cols <- complementary("purple4", plot = F)
plot(individuals.areas, vertex.color = cor[as.numeric(row.names(points.areas.decrease))], #cols[2], vertex.label.color = 'white',
     #vertex.label.
     edge.color = 'grey50', edge.width = E(individuals.areas)$weight**2/2, 
     layout = layout_with_graphopt(individuals.areas, charge = 0.01))

mtext("E", 3, at = -1.7, cex = 1.5)


# # map + space use (kernel 50%)
# # Valor das HR
# kern95 <- data.frame(-1, -1)[-1,]
# colnames(kern95) <- c("id", "area")
# homeranges <- list()
# cor <- c(rainbow(length(ids)-1), "black") # the last one is for the whole population
# for(i in (1:length(ids))) {
#   homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 50)
#   homeranges[[i]] <- homerange
#   kern50 <- rbind(kern50, homerange@data)
#   if(i == 1){
#     rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
#     cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
#               rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
#     # x <- c(196250, 209856)
#     # y <- c(7561021, 7574393)
#     # plot(x, y, type = "n", bty = "n", axes = F)
#     plot(landuse.map, col = cols, border = cols, lwd = 0.1)
#   }
#   plot(homerange, border=cor[i], lwd = 2, add = T)
# }
# kern95

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

# Fig. S5 from the manuscript - SpatIS significance
# png("output/spatis_significance.png", width = 15, height = 15, units = "cm", res = 600)
tiff("output/spatis_significance.tif", width = 15, height = 15, units = "cm", res = 600)

brk_plot <- seq(0, 1, 0.05)
hist(unlist(randomized.nest[[1]]$SpatIS.individual.random), col = 'grey', freq = F, breaks = brk_plot, xlim = c(0, 1), 
     main = '', xlab = 'Spatial Individual Specialization', ylab = 'Probability density')
hist(randomized.nest[[1]]$SpatIS.individual.observed, freq = F, breaks = brk_plot, col = rgb(1, 0, 0, alpha = 0.4), add = T)
abline(v = mean(unlist(randomized.nest[[1]]$SpatIS.individual.random)), col = 1)
abline(v = mean(unlist(randomized.nest[[1]]$SpatIS.individual.observed)), col = 2)

dev.off()

# Fig. S6 from the manuscript - SpatICS significance
png("output/spatics_significance.png", width = 15, height = 15, units = "cm", res = 600)
#tiff("output/spatics_significance.tif", width = 15, height = 15, units = "cm", res = 600)

brk_plot <- seq(0, 1, 0.05)
hist(unlist(randomized.nest[[2]]$SpatICS.individual.random), col = 'grey', freq = F, breaks = brk_plot, xlim = c(0, 1), 
     main = '', xlab = 'Spatial Individual Complementary Specialization', ylab = 'Probability density')
hist(randomized.nest[[2]]$SpatICS.individual.observed, freq = F, breaks = brk_plot, col = rgb(1, 0, 0, alpha = 0.4), add = T)
abline(v = mean(unlist(randomized.nest[[2]]$SpatICS.individual.random)), col = 1)
abline(v = mean(unlist(randomized.nest[[2]]$SpatICS.individual.observed)), col = 2)

dev.off()

###################
# Fig. S2

# Export kernel 50%
png("output/Fig_S1.png", width = 18, height = 18, units = "cm", res = 300)
par(mfrow = c(2,2), mar = c(0, 0, 0, 0))

rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
          rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
plot(landuse.map, col = cols, border = 'black', lwd = 1)

mtext("A", 3, at = 199200, line = -3, cex = 1.5)

par(mar = c(0, 0, 0, 0))
kern50 <- data.frame(-1, -1)[-1,]
colnames(kern50) <- c("id", "area")
homeranges <- list()
cor <- c(RColorBrewer::brewer.pal(10, 'Set3'), 'black')
for(i in 1:(length(ids)-1)) {
  homerange <- adehabitatHR::getverticeshr(kernels[[i]], percent = 50)
  homeranges[[i]] <- homerange
  kern50 <- rbind(kern50, homerange@data)
  if(i == 1){
    rgb255 <- function(r, g, b) rgb(r, g, b, maxColorValue = 255)
    cols <- c(rgb255(103, 100, 107), rgb255(73, 70, 77), rgb255(130, 130, 130), 
              rgb255(225, 225, 225), rgb255(178, 178, 178), rgb255(204, 204, 204))
    plot(landuse.map, col = cols, border = cols, lwd = 0.1)
    raster::scalebar(2000, xy = c(200500, 7563000), divs = 4, type = 'bar', label = c(0, '1', '2 km'))
    
  }
  plot(homerange, border=cor[i], lwd = 2, add = T)
}
kern50

mtext("B", 3, at = 199200, line = -3, cex = 1.5)

plot(landuse.map, col = cols, border = cols, lwd = 0.1)
points(spdados, pch = 20, col = cor[as.factor(spdados$Animal_ID)])
plot(resource.areas, border = "red", lwd = 2, add = T)

mtext("C", 3, at = 199200, line = -3, cex = 1.5)

plot(0, 0, xlim = c(0,1), ylim = c(0,1), type = "n", axes = F)
legend(x = 0.6, y = 1, legend = 1:10, col = cor, lwd = 3, title = "Individuals", cex = 1.2)
legend(x = 0, y = 1, legend = levels(landuse.map$Classe), fill = cols, border = cols,
       title = "Land use class", cex = 1.2)

dev.off()
