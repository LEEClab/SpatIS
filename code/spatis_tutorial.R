
# Loading packages
if(!require(adehabitatLT)) install.packages("adehabitatLT", dep=T); library(adehabitatLT)
if(!require(date)) install.packages("date", dep=T); library(date)
if(!require(adehabitatHR)) install.packages("adehabitatHR", dep=T); library(adehabitatHR)
if(!require(sp)) install.packages("sp", dep=T); library(sp)
if(!require(rgdal)) install.packages("rgdal", dep=T); library(rgdal)

# Loading packages
if(!require(plotrix)) install.packages("plotrix", dep=T); library(plotrix)


beta <- c(1.5, 2, 1.8, 2.2, 0.1) # Coefficient of attraction or bias - positive values correpond to attraction, negative values correspond to avoidance
rho <- 0.6 # Concentration parameter around the bias absolute angle
scale <- 1
shape <- 1

# Number of individuals
ntracks <- 5
# Number of steps per trajectory/individual
nsteps <- 100

# Matrices of locations
X <- matrix(NA, nsteps, ntracks)
Y <- matrix(NA, nsteps, ntracks)

# Coordinates of the point of attraction/repulsion for each individual
xh <- c(-40,30,40,-25,0) 
yh <- c(35,30,-42,-30,0)

for(i in 1:ntracks){
  x <- numeric(nsteps)
  y <- numeric(nsteps) 
  h <- numeric(nsteps)
  steps <- numeric(nsteps)
  
  h[1] <- runif(1,1,2*pi)
  x[1] <- rnorm(1,0,1)
  y[1] <- rnorm(1,0,1)
  
  for(t in 2:nsteps){  
    adj <- xh[i] - x[t-1]
    op  <- yh[i] - y[t-1]
    r   <- sqrt(adj^2 + op^2)
    ya <- sin(h[t-1]) + beta[i]*(op/r)
    xa <- cos(h[t-1]) + beta[i]*(adj/r)    
    m_t <- atan2(ya,xa)
    h[t] <- rwrappedcauchy(1,mu=circular(m_t),rho=rho)
    steps[t-1] <- rweibull(1,scale=scale, shape=shape)
    x[t] <- x[t-1] + cos(h[t])*steps[t-1]
    y[t] <- y[t-1] + sin(h[t])*steps[t-1]
  } 
  X[,i] <- x
  Y[,i] <- y
}

# Draw landscape with prefered resources


radius <- c(10, 5, 8, 10, 5)
cols <- grey.colors(length(xh), alpha = 0.5)
matplot(X,Y, type="n", xlim=c(1.1*min(X),1.1*max(X)),ylim=c(1.1*min(Y),1.1*max(Y)))
for(i in 1:length(xh)) {
  draw.circle(xh[i], yh[i], radius = radius[i], border = cols[i], col = cols[i])
}

matplot(X, Y, type="l", pch=16, col=1:ntracks, asp=1, 
        xlim=c(min(X),max(X)),ylim=c(min(Y),max(Y)), add = T)

