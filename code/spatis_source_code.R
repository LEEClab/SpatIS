################################################################################
#
# SpatIS: Spatial Individual Specialization index
# 
# Patricia K. Rogeri - pa_bio04@yahoo.com.br
# Bernardo B. S. Niebuhr - bernardo_brandaum@yahoo.com.br
# Renata L. Muylaert - renatamuy@gmail.com 
#
# Function SpatIS, to calculate Individual Specializations considering the use 
#   of space. The functional calculates utilization distributions (UD) for both 
#   the individuals and the whole population, based on telemetry location points.
#
#   The individual SpatIS is calculated as the volume of the individual UD 
#   (or a X% kernel area, X defined by the user) that do not overlap with the 
#   populational UD (or a correspondent kernel area).
#
# Call: SpatIS(data, individuals.col, population.ID = NULL, method, ...)
#
# Inputs:
# data:            a SpatialPointsDataFrame object containing the locations of 
#                  the individuals, their IDs, and other spatial or individual 
#                  information.
#
# individuals.col: name(string) or number of the column of the input 
#                  SpatialPointsDataFrame 'data' that represents the identity 
#                  of individuals.
#
# population.ID:   ID (string or number) that represents the population in the
#                  column 'individuals.col' of the input dataset 'data'. In case
#                  the locations of all individuals were not pooled and combined
#                  to the original data, population.ID is set to NULL (the default).
#                  In this case, the calculations for the whole population are
#                  made inside the SpatIS function.
#
# method:          method of overlap between utilization distributions (the same
#                  options as the function kerneloverlap from adehabitatHR 
#                  package. See ?kerneloverlap for more information). The default
#                  is "VI", which computes the volume of the intersection between 
#                  the individual and populational UD.
# 
# ...:             additional arguments passed to the functions kerneloverlap and
#                  kernelUD, from the package adehabitatHR.
#
# Output:
# A list of four elements:
# data:              the data used to calculate SpatIS (with locations corresponding
#                    to the whole population, independently of whether the population
#                    locations were already in the input data). 
#                    It is a SpatialPointsDataFrame.
#
# parms:             a list with parameters used as input to call SpatIS function.
#
# SpatIS.individual: a vector of Spatial individual specialization indices for
#                    each individual of the population (the overlap between
#                    each individual and the population UDs).
#
# SpatIS.population: a value of Spatial individual specialization for the population,
#                    calculated as the average of all SpatIS individual values.
#
# July 2017
################################################################################

# Load packages
if(!require(adehabitatHR)) install.packages("adehabitatHR", dep=T); library(adehabitatHR)
if(!require(sp)) install.packages("sp", dep=T); library(sp)
if(!require(rgdal)) install.packages("rgdal", dep=T); library(rgdal)

# Function SpatIS declaration
SpatIS <- function(data, individuals.col, population.ID = NULL, 
                   method = c("VI", "HR", "PHR", "BA", "UDOI", "HD"), ...)
{
  # Check if the data are of the class SpatialPoints
  if (!inherits(data, "SpatialPoints")) 
    stop("Data should inherit the class SpatialPoints.")
  
  # Copying data and transforming ID values in character variables
  data2 <- data
  data2[[individuals.col]] <- as.character(data2[[individuals.col]])
  
  # Method of overlap
  method = method[1]
  
  # If population.ID is NULL, the points representing all the population were not
  #   calculated yet. They will be calculated then.
  if(is.null(population.ID)) {
    population.ID <- "all"
    pop.dat <- data2 # Copying data
    pop.dat[[individuals.col]] <- population.ID # We create a new "individual" with the ID "all"
                                                # which represents all population points
    data.aux <- rbind(data2, pop.dat) # Combining data from individuals and the population
  } else {
    # If the population.ID was furnished by the user, check if it really exists in the input data
    if(!(population.ID %in% unique(data2[[individuals.col]]))) {
      # If the population.ID does not exist, show an error message
      stop(paste("The population ID \"",population.ID,"\" is absent from the input data. Plase set the parameter population.ID to NULL.", sep = ""))
    } else {
      # If it exists, consider the input data itself and go on
      data.aux <- data2
    }
  }
  
  # This is a matrix with the overlap of utilization distribution between each pair of individuals
  # and each individual and the whole population
  over <- kerneloverlap(data.aux[,1], method = method, ...)
  
  # Line in the overlap matrix that represents the population
  population.line <- which(rownames(over) == population.ID)
  
  # Overlap of each individual with the whole population utilization distribution
  SpatIS.ind.aux <- over[-population.line,population.line]
  # SpatIS for individuals = 1 - overlap of the individual with the population
  SpatIS.ind <- 1 - SpatIS.ind.aux
  # SpatIS = average of individual SpatIS
  SpatIS.pop <- mean(SpatIS.ind)
  
  # List of parameters
  parms <- list(individuals.col = individuals.col, population.ID = population.ID, method = method, ...)
  
  return( list(data = data.aux, parms = parms, SpatIS.individual = SpatIS.ind, SpatIS.population = SpatIS.pop) )
}


################################################################################
#
# SpatIS_randomize
# 
# Patricia K. Rogeri - pa_bio04@yahoo.com.br
# Bernardo B. S. Niebuhr - bernardo_brandaum@yahoo.com.br
# Renata L. Muylaert - renatamuyt@gmail.com
#
# Function SpatIS_randomize, to randomize the locations of individuals, calculate
#   their utilization distributions, calculate SpatIS for numerous times using a
#   bootstrap procedure and compare the randomized SpatIS with the observed one.
#   This frunction may be used to calculate whether the oserved SpatIS in 
#   different from a random SpatIS value (a comparison with the null model).
#
# Call: SpatIS_randomize(SpatIS.object, iterations = 99, bootstrap = TRUE, plot = TRUE)
#
# Inputs:
# data:       The output of SpatIS function, a list containing:
#             1) the original data points used, for individuals and
#             the whole population; 2) the list of parameters used as input
#             to SpatIS function; 3) the observed individual SpatIS values; 
#             and 4) the observed population SpatIS value.
#
# iterations: Number of iterations to be performed for assessing significance
#             of the population SpatIS.
#
# bootstrap:  logical. If TRUE, each individual location may be selected more than once
#             in the randomization process, and some of the real locations may not be
#             selected in a given iteration (bootstrap process). If FALSE (Default), 
#             there is no replacement in the randomization of locations, so individual
#             locations are permutated instead of bootstraped (permutation process).
#
# plot:       Whether or not to plot the comparison between observed and 
#             randomized population SpatIS values.
#
# Output:
# A list of three elements:
# SpatIS.population.random:   a vector with the observed population SpatIS as the 
#                             first element and the randomized population SpatIS
#                             values calculated during the bootstrap procedure.
#
# SpatIS.population.observed: the observed population SpatIS (used as input).
#
# p:                          the p-value which represents the one-tailed significance. It 
#                             answers the question: how many times randomized 
#                             individual locations present a population SpatIS
#                             value equal or greater than the observed value?
#
# July 2017
################################################################################

# Function SpatIS_randomize declaration
SpatIS_randomize <- function(SpatIS.object, iterations = 99, bootstrap = TRUE, plot = TRUE)
{
  # Check if the observed SpatIS was correctly calculated
  type <- class(SpatIS.object)
  class.input <- class(SpatIS.object[[1]])[1]
  nrow.input <- nrow(SpatIS.object[[1]])
  class.parms <- class(SpatIS.object[[2]])
  length.spatis.pop <- length(SpatIS.object[[4]])
  
  if(type != 'list' | class.input != 'SpatialPointsDataFrame' | nrow.input < 5 | 
     class.parms != 'list' | length.spatis.pop < 1) {
    stop('There was a problem in the calculation of the input observed SpatIS.')
  }
  
  data <- SpatIS.object[[1]]
  parms <- SpatIS.object[[2]]
  SpatIS.ind <- SpatIS.object[[3]]
  SpatIS.pop <- SpatIS.object[[4]]
  
  # Ramdomized SpatIS values - initialized with the observed value
  SpatIS.pop.random <- SpatIS.pop
  
  # Generate one randomized SpatIS value per permutation
  for(i in 1:iterations)
  {
    print(paste('Iteration =', i))
    
    # Remove population IDs and randomize individual IDs
    data2 <- data[,parms$individuals.col][data[[parms$individuals.col]] != parms$population.ID,]
    # If boostrap argument == TRUE, sample points with replacement - each location may be chosen more than once
    if(bootstrap) {
      data2[[parms$individuals.col]] <- sample(data2[[parms$individuals.col]], replace = TRUE) # should we bootstrap with replace = TRUE?
    } else { # If boostrap argument == FALSE, locations are permutated instead of bootstraped (no replacement)
      data2[[parms$individuals.col]] <- sample(data2[[parms$individuals.col]])
    }
    
    # Combine randomized individuals with population
    pop.dat <- data2 # Copying data
    pop.dat[[parms$individuals.col]] <- parms$population.ID # We create a new "individual" with the ID "all"
                                                            # which represents all population points
    data.aux <- rbind(data2, pop.dat) # Combining data from randomized individuals and the population
    
    # Parameters to call kerneloverlap - getting them from the input SpatIS data
    over.parms <- list(xy = data.aux[,1], method = parms$method)
    if(length(parms) > 3) {
      name <- names(parms[4:length(parms)])
      val <- parms[4:length(parms)]
      for(j in 1:length(val)) {
        over.parms[name[j]] <- val[j]
      }
    }
    
    # This is a matrix with the overlap of utilization distribution between each pair of  randomized individuals
    # and each individual and the whole population
    over <- do.call(kerneloverlap, over.parms)
    
    # Line in the overlap matrix that represents the population
    population.line <- which(rownames(over) == parms$population.ID)
    
    # Overlap of each randomized individual with the whole population utilization distribution
    SpatIS.ind.aux <- over[-population.line,population.line]
    # Randomzed individual SpatIS = 1 - overlap of the randomized individual with the population
    SpatIS.ind.random <- 1 - SpatIS.ind.aux
    # Randomzed population SpatIS = average of randomized individual SpatIS
    # Append that to the vector of randomized values
    SpatIS.pop.random <- c(SpatIS.pop.random, mean(SpatIS.ind.random))
  }
  
  # Calculating significance (p-value)
  (p <- sum(SpatIS.pop <= SpatIS.pop.random)/(iterations+1)) # proportion of random values that are greater than the observed value
  # Plot signigicance
  if(plot == TRUE) {
    hist(SpatIS.pop.random[-1], main = "", xlab = "Spatial Individual Specialization", ylab = "Frequency",
         xlim = c(min(SpatIS.pop.random), max(SpatIS.pop.random)))
    abline(v = SpatIS.pop, col = "red")
  }
  
  return( list(SpatIS.population.random = SpatIS.pop.random, SpatIS.population.observed = SpatIS.pop, p = p))
}