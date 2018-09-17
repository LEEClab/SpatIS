######################################################
#
# Linking individual specialization to space use
# for Sturnira lilium bats
#
# Testing for effects of body size and sex on
# interindividual variaion
#
# Patricia Rogeri - pa_bio04 at yahoo.com.br
# Bernardo Niebuhr - bernardo_brandaum at yahoo.com.br
# Renata Muylaert - remuylaert at gmail.com
#
# July 2017
# No copyrights - feel free to use, modify, and share
#######################################################

# Loading pacakges
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load('bbmle', 'lme4')

# Path to data folder
datadir <- "/home/leecb/Github/SpatIS/data"

# Loading data
# Change to the data folder
setwd(datadir)

######################################################################
# 1) Verifying the influence of body mass and sex on interindividual 
#    variation considering individual use of foraging sites (kernel 50%)
######################################################################

# Load data
use_vs_foraging_sites <- read.table("test_for_bodysize_sex_foraging_sites.csv", header = T, sep = "\t", dec = ",")

# Organize data
colnames(use_vs_foraging_sites)[c(1:3)] <- c("ID", "sex", "log.weight")
use_vs_foraging_sites$ID <- as.factor(use_vs_foraging_sites$ID)

############################################################
# Description of the data
# 
# ID: ID of bat individual
# sex: sex of the individual
# log.weight: logarithm of the body size (in grams)
# FS1-FS9: Arccossine of the square root of the 
#          proportion of use of each foraging site (1-9)
#          by each individual
#
#############################################################

# Transforming it into long format
use_vs_foraging_sites.long <- reshape(use_vs_foraging_sites, varying = 4:12, v.names = "frequency.use", direction = "long", idvar = "ID")
colnames(use_vs_foraging_sites.long)[4] <- "foraging.site"
head(use_vs_foraging_sites.long)

# General statistics
# Body size
range(10^(use_vs_foraging_sites.long$log.weight))
mean(10^(use_vs_foraging_sites.long$log.weight)); sd(10^(use_vs_foraging_sites.long$log.weight))

# Use of foraging sites
range(use_vs_foraging_sites.long$frequency.use)
mean(use_vs_foraging_sites.long$frequency.use); sd(use_vs_foraging_sites.long$frequency.use)

# Run linear models
# No effect model
m0 <- glm(frequency.use ~ 1, data = use_vs_foraging_sites.long)
# Sex + weight model
m1 <- glm(frequency.use ~ sex + log.weight, data = use_vs_foraging_sites.long)
# Sex + weight + ID as a fixed effect
m2 <- glm(frequency.use ~ sex + log.weight + ID, data = use_vs_foraging_sites.long)
# Sex + weight + ID as a random effect
m3 <- lmer(frequency.use ~ sex + log.weight + (1|ID), data = use_vs_foraging_sites.long)

# Comparing models
AICctab(m0, m1, m2, m3, delta = T, weights = T) # Null wins

# Summary of models
summary(m1)
summary(m2)
summary(m3)

# Checking for percentages explained
af <- anova(m2)
afss <- af$"Deviance"
print(cbind(af, PctExplained = afss/sum(afss, na.rm = T)*100))

######################################################################
# 2) Verifying the influence of body mass and sex on interindividual 
#    variation considering individual use of polygons (landscape elements)
######################################################################

# Load data
use_vs_poligons <- read.table("test_for_bodysize_sex_polygons.csv", header = T, sep = ";", dec = ",")

# Organize data
colnames(use_vs_poligons)
colnames(use_vs_poligons)[c(1:3)] <- c("ID", "sex", "log.weight")
use_vs_poligons$ID <- as.factor(use_vs_poligons$ID)

####################################################################
# Description of the data
# 
# ID: ID of bat individual
# sex: sex of the individual
# log.weight: logarithm of the body size (in grams)
# FS1-FS31: Arccossine of the square root of the 
#           proportion of use of each habitat type polygon (1-31)
#           by each individual
#
####################################################################

# Transforming it into long format
use_vs_poligons.long <- reshape(use_vs_poligons, varying = 4:ncol(use_vs_poligons), v.names = "frequency.use", direction = "long", idvar = "ID")
colnames(use_vs_poligons.long)
colnames(use_vs_poligons.long)[4] <- "polygon"
head(use_vs_poligons.long)
table(use_vs_poligons.long[4])

# General statistics
# Body size
range(10^(use_vs_poligons.long$log.weight))
mean(10^(use_vs_poligons.long$log.weight)); sd(10^(use_vs_foraging_sites.long$log.weight))

# Use of foraging sites
range(use_vs_poligons.long$frequency.use)
mean(use_vs_poligons.long$frequency.use); sd(use_vs_foraging_sites.long$frequency.use)

# Run linear models
# No effect model
m0p <- glm(frequency.use ~ 1, data = use_vs_poligons.long)
# Sex + weight model
m1p <- glm(frequency.use ~ sex + log.weight, data = use_vs_poligons.long)
# Sex + weight + ID as a fixed effect
m2p <- glm(frequency.use ~ sex + log.weight + ID, data = use_vs_poligons.long)
# Sex + weight + ID as a random effect
m3p <- lmer(frequency.use ~ sex + log.weight + (1|ID), data = use_vs_poligons.long)

# Comparing models
AICctab(m0p, m1p, m2p, m3p, delta = T, weights = T) # Null wins

# Summary of models
summary(m1p)
summary(m2p)
summary(m3p)

# Checking for percentages explained
af <- anova(m2p)
afss <- af$"Deviance"
print(cbind(af, PctExplained = afss/sum(afss, na.rm = T)*100))

