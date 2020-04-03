################################################################################
#
# SpatIS: Spatial Individual Specialization index
# 
# Patricia K. Rogeri - pa_bio04@yahoo.com.br
# Bernardo B. Niebuhr - bernardo_brandaum@yahoo.com.br
# Renata L. Muylaert - renatamuy@gmail.com 
#
# Function SpatIS, to calculate Individual Specializations considering the use 
#   of space. The function calculates utilization distributions (UD) for both 
#   the individuals and the whole population, based on telemetry location points.
#
#   The individual SpatIS is calculated as the volume of the individual UD 
#   (or a X% kernel area, X defined by the user) that do not overlap with the 
#   populational UD (or a correspondent population kernel area).
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
#
# License: GPLv2 (GNU General Public License v2.0)
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
# t.power
# 
# Bernardo B. Niebuhr - bernardo_brandaum@yahoo.com.br
# Modified from a function by Phil Spector (https://www.stat.berkeley.edu/~spector/).
#
# Function t.power, to calculate power of a t-test analysis by comparing each randomized 
# array of individual SpatIS values with the observed array of individual SpatIS values, 
# using a t-test. As a result, one has n_iterations results of the t-test. The function then 
# calculates power as the proportion of times we got the t.test p-value smaller than a given alpha. 
# This function was slightly modified from here: https://statistics.berkeley.edu/computing/r-t-tests.
# All the merits and authorship to Phil Spector, from the Department of Statistics, 
# University of California, Berkeley.
#
# Call: t.power(obs, exp, alpha = 0.05)
#
# Inputs:
# obs:        An array of observed individual SpatIS values, with as many
#             values as individuals sampled. SpatIS values are calculated by
#             the function SpatIS.
#
# exp:        List of n_iteration arrays of individual SpatIS values 
#             (one value for each individual), calculated by the SpatIS_randomize function.
#
# alpha:      Threshold for significance. Default is 0.05.
#
# Output:
# 
# power:      Power of the analysis, calculated as the proportion of significant p-values for
#             the t.tests comparing the observed and expected individuals SpatIS values
# 
# Aug 2018
################################################################################

# Function t.power declaration
t.power = function(obs, exp, alpha = 0.05){
  
  # Perform a t.test for comparing the observed and each of the expected arrays of 
  # individual SpatIS values
  ts = lapply(X = exp, FUN = t.test, y = obs, alternative = 'less', mu = 0, 
              var.equal = FALSE, conf.level = 1 - alpha)
  
  # Extract p-values
  pval <- function(x) x$p.value
  tps <- unlist(lapply(ts, pval))
  
  # Return the proportion of significant p-values as the power of the analysis
  sum(tps < alpha) / length(tps)
}

################################################################################
#
# t.power.vs.n
# 
# Bernardo B. Niebuhr - bernardo_brandaum@yahoo.com.br
# Modified from a function by Phil Spector (https://www.stat.berkeley.edu/~spector/).
#
# Function t.power.vs.n, to calculate power of a t-test analysis by comparing each randomized 
# array of individual SpatIS values with the observed array of individual SpatIS values, 
# using a t-test, but considering only n individuals sampled at a time.
# This function was slightly modified from here: https://statistics.berkeley.edu/computing/r-t-tests.
# All the merits and authorship to Phil Spector, from the Department of Statistics, 
# University of California, Berkeley.
#
# Call: function(n, observed, expected, alpha = 0.05, n.repeat = 100)
#
# Inputs:
# n:          The number of values to be sampled from observed and expected
#             distributions, before performing power analysis.
#
# observed:   An array of observed individual SpatIS values, from which n values
#             will be sampled. SpatIS values are calculated by the function SpatIS.
#
# expected:   List of arrays of individual SpatIS values (one value for each individual), 
#             calculated by the SpatIS_randomize function. Each of these arrays is used
#             for sampling n values for power analysis.
#
# alpha:      Threshold for significance. Default is 0.05.
#
# n.repeat:   Number of repetitions of the sampling procedire, for averaging purposes.
#
# Output:
# 
# power:      Average power of the analysis, calculated as the proportion of significant 
#             p-values for the t.tests comparing the observed and expected individuals SpatIS values
#             (but only the n values sampled from each distribution).
# 
# Aug 2018
################################################################################

# Function t.power.vs.n declaration
t.power.vs.n <- function(n, observed, expected, alpha = 0.05, n.repeat = 100) {
  
  # Declaration of array for values of power
  power.val <- c()
  
  # For each iteration
  for(i in 1:n.repeat) {
    # Sample n values from the observed distribution
    samp.obs <- sample(observed, n)
    # Sample n values from the expected distribution
    samp.exp <- lapply(expected, sample, n)
    # Calculate the power of a t.test comparing observed and expected distributions
    power.val <- c(power.val, t.power(samp.obs, samp.exp, alpha = alpha))
  }
  
  # list(power.val, mean(power.val))
  # Return the mean power for each sampling size n
  mean(power.val)
}

################################################################################
#
# SpatIS_randomize
# 
# Patricia K. Rogeri - pa_bio04@yahoo.com.br
# Bernardo B. Niebuhr - bernardo_brandaum@yahoo.com.br
# Renata L. Muylaert - renatamuyt@gmail.com
#
# Function SpatIS_randomize, to randomize the locations of individuals, calculate
#   their utilization distributions, calculate SpatIS numerous times using a
#   bootstrap procedure and compare the randomized SpatIS with the observed one.
#   This frunction may be used to calculate whether the observed SpatIS is 
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
# alpha:      significance threshold (default = 0.05).
#
# plot:       Whether or not to plot the comparison between observed and 
#             randomized population SpatIS values.
#
# brk_plot:   argument breaks for the histogram plot of the observed and expected (randomized)
#             distribution of individual SpatIS values (only if plot = TRUE).
#
# Output:
# A list of seven elements:
# SpatIS.individual.random:   a list with a number of array equals to 'iterations', i.e.,
#                             one array per randomization. Each array presents an individual
#                             SpatIS value for each individual sampled, but based on randomized
#                             individual locations.
#
# SpatIS.individual.observed: an array of observed individual SpatIS values, one for each sampled
#                             individual (this array is the same present in the input data for the
#                             SpatIS_randomize function).
#
# SpatIS.population.random:   an array with the expected population SpatIS
#                             values calculated during the bootstrap procedure, after the
#                             randomization of the individuals' locations.
#
# SpatIS.population.observed: the observed population SpatIS value (used as input).
#
# SpatIS.significance:        an object containing the t statistics, p-value and other
#                             results of a t-test, used for assessing SpatIS significance.
#                             The t-test compares the observed and expected (after randomization)
#                             distributions of SpatIS.
#
# SpatIS.power:               power of the t.test analysis.
#
# SpatIS.power.curve:         data frame with values for power considering only a subset of n individuals,
#                             to assess sampling sufficiency. n is varies between n and the 
#                             total number of individuals.
# 
# July 2017
# Updated Aug 2018
#
# License: GPLv2 (GNU General Public License v2.0)
################################################################################

# Function SpatIS_randomize declaration
SpatIS_randomize <- function(SpatIS.object, iterations = 99, bootstrap = FALSE, alpha = 0.95,
                             plot = TRUE, brk_plot = seq(0, 1, 0.05))
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
  # SpatIS.population.obs.rand is an array with population SpatIS values; the first value is the
  # observed one and the following are the expected ones (one value per iteration), after randomization
  # of individuals' locations
  SpatIS.population.obs.rand <- SpatIS.pop
  # SpatIS.individual.obs.rand is a list of arrays with individual SpatIS values; the first array
  # presents the observed individual SpatIS values, one for each individual, and the following ones
  # (one array per iteration) to the randomized SpatIS values, after randomization of 
  # individuals' locations 
  SpatIS.individual.obs.rand <- list()
  SpatIS.individual.obs.rand[[1]] <- SpatIS.ind
  
  # Generate one randomized SpatIS value per permutation
  for(i in 1:iterations)
  {
    print(paste('Iteration =', i))
    
    # Remove population IDs and randomize individual IDs
    data2 <- data[, parms$individuals.col][data[[parms$individuals.col]] != parms$population.ID,]
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
    # Randomized individual SpatIS = 1 - overlap of the randomized individual with the population
    SpatIS.ind.aleat <- 1 - SpatIS.ind.aux
    SpatIS.individual.obs.rand[[(i+1)]] <- SpatIS.ind.aleat
    # Randomized population SpatIS = average of randomized individual SpatIS
    # Append that to the vector of randomized values
    SpatIS.population.obs.rand <- c(SpatIS.population.obs.rand, mean(SpatIS.ind.aleat))
  }
  
  # Calculating significance and power
  observed <- SpatIS.individual.obs.rand[[1]] # The first array is the observed individual SpatIS values
  expected <- SpatIS.individual.obs.rand[2:(iterations+1)] # The following arrays are the randomized values
  expected.polled <- unlist(expected) # Here we pool everything in a single distribution
  
  # Calculating significance (p-value)
  # (p <- sum(SpatIS.pop <= SpatIS.population.obs.rand)/(iterations+1)) # proportion of random values that are greater than the observed value
  significance <- t.test(observed, expected.polled, alternative = 'greater', mu = 0,
                         var.equal = FALSE, conf.level = 1-alpha)
  
  # Plot signigicance
  if(plot == TRUE) {
    # hist(SpatIS.population.obs.rand[-1], main = "", xlab = "Spatial Individual Specialization", ylab = "Frequency",
    #      xlim = c(min(SpatIS.population.obs.rand), max(SpatIS.population.obs.rand)))
    # abline(v = SpatIS.pop, col = "red")
    hist(expected.polled, col = 'grey', freq = F, breaks = brk_plot, xlim = c(0, 1), 
         main = '', xlab = 'Spatial Individual Specialization', ylab = 'Frequency')
    hist(observed, freq = F, breaks = brk_plot, col = rgb(1, 0, 0, alpha = 0.4), add = T)
    abline(v = mean(SpatIS.population.obs.rand[2:iterations]), col = 1)
    abline(v = SpatIS.population.obs.rand[1], col = 2)
  }
  
  # Calculating power
  power.value <- t.power(obs = observed, exp = expected, alpha = alpha)
  
  # Calculating sampling sufficiency
  sampling.suff <- data.frame(n = 2:length(observed), power = unlist(lapply(2:length(observed), t.power.vs.n, observed = observed, expected = expected, alpha = alpha, n.repeat = 100)))
  
  # Return results
  return( list(SpatIS.individual.random = expected, SpatIS.individual.observed = observed, 
               SpatIS.population.random = SpatIS.population.obs.rand[2:iterations], SpatIS.population.observed = SpatIS.population.obs.rand[[1]], 
               SpatIS.significance = significance, SpatIS.power = power.value, SpatIS.power.curve = sampling.suff))
}