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
# Renata Muylaert - renatamuy at gmail.com
#
# July 2017
# No copyrights - feel free to use, modify, and share
#######################################################

# Loading pacakges
if(!require(install.load)) install.packages('install.load'); library(install.load)
install.load::install_load('bbmle', 'lme4')
install.load::install_load('sjPlot', 'glmmTMB', 'ggpubr')

######################################################################
# 1) Verifying the influence of body mass and sex on interindividual 
#    variation considering individual use of foraging sites (kernel 50%)
######################################################################

# Load data
use_vs_foraging_sites <- read.csv("data/proportion_use_foraging_sites.csv")

# Organize data
use_vs_foraging_sites$Animal_ID <- as.factor(use_vs_foraging_sites$Animal_ID)

# General statistics
# Body size
range(use_vs_foraging_sites$body_mass_g)
mean(use_vs_foraging_sites$body_mass_g); sd(use_vs_foraging_sites$body_mass_g)

# proportion
range(p.fs <- use_vs_foraging_sites$success/use_vs_foraging_sites$sample.size)
mean(p.fs); sd(p.fs)

# Run linear models
# No effect model
m0 <- glm(cbind(success, failure) ~ 1, family = binomial, data = use_vs_foraging_sites)
# Sex + weight model
m1 <- glm(cbind(success, failure) ~ sex + body_mass_g, family = binomial, data = use_vs_foraging_sites)
# Sex + weight + ID as a fixed effect
m2 <- glm(cbind(success, failure) ~ sex + body_mass_g + Animal_ID, family = binomial, data = use_vs_foraging_sites)
# Sex + weight + ID as a random intercept
m3 <- glmer(cbind(success, failure) ~ sex + body_mass_g + (1|Animal_ID), family = binomial, data = use_vs_foraging_sites)
# In the models above, residual deviance is much higher than the degrees of freedom
# what means the data is overdispersed; we should add the polygon as a random intercept
m4 <- glmer(cbind(success, failure) ~ sex + body_mass_g + (1|Animal_ID) + (1|area), family = binomial, data = use_vs_foraging_sites)
m5 <- glmer(cbind(success, failure) ~ 1 + (1|Animal_ID) + (1|area), family = binomial, data = use_vs_foraging_sites)

# Comparing models
AICctab(m0, m1, m2, m3, m4, m5, delta = T, weights = T) # m0 wins

# Summary of models
summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)

# Checking for percentages explained
af <- anova(m2)
afss <- af$"Deviance"
print(cbind(af, PctExplained = afss/sum(afss, na.rm = T)*100))

## PRODUCING DIAGNOSTIC PLOTS ##

## Logistic regression
plot(m2$fit, residuals(m2), pch = 19, las = 1, cex = 1.4)
abline(0,0,lwd = 1.5)
## check for no pattern - there is some pattern

plot(m1$fit, residuals(m1), pch = 19, las = 1, cex = 1.4)
abline(0,0,lwd = 1.5)
## check for no pattern - there is some pattern

## GLMM
plot(m2)
plot(m4)

# same plots using ggplot; however, for GLMM these plots do no make
# so sense
qqnorm(residuals(m4))
ggplot(data.frame(lev=hatvalues(m4),pearson=residuals(m4,type="pearson")),
       aes(x=lev,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(eta=predict(m4,type="link"),pearson=residuals(m4,type="pearson")),
       aes(x=eta,y=pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(x1=use_vs_foraging_sites$sex, pearson=residuals(m4,type="pearson")),
       aes(x=x1,y=pearson)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "glm")

ggplot(data.frame(x2=use_vs_foraging_sites$body_mass_g,pearson=residuals(m4,type="pearson")),
       aes(x=x2,y=pearson)) +
  geom_point() +
  theme_bw()

#############
# Prepare additional plots with ggplot for GLMM
# this is based on code by Ben Bolker: 
# https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html#diagnostics
# and on these plots by Raju Rimal:
# https://rpubs.com/therimalaya/43190

# std residuals vs fitted values
plot(m4, type=c("p", "smooth"))

p1 <- ggplot(m4, aes(.fitted, scale(.resid))) +
  geom_point() +
  stat_smooth(method = "loess") + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Fitted values") +
  ylab("Standardized residuals") +
  ggtitle("Residual vs Fitted")+
  theme_bw()
p1

# qqnorm residuals - we are not including this one because it does not make sense 
# for a binomial fit
qqnorm(residuals(m4))

p2 <- ggplot(m4, aes(qqnorm(scale(resid(m4)))[[1]], scale(resid(m4)))) + 
  geom_point(na.rm = TRUE) +
  geom_abline(aes(slope = 1, intercept = 0)) + 
  xlab("Theoretical Quantiles") +
  ylab("Standardized Residuals") + 
  ggtitle("Gaussian Q-Q") +
  theme_bw()
p2

# Scale-location plot
p3 <- ggplot(m4, aes(.fitted, sqrt(abs(scale(resid(m4))))))+
  geom_point(na.rm=TRUE) +
  stat_smooth(method="loess", na.rm = TRUE) +
  xlab("Fitted Value") +
  ylab(expression(sqrt("|Standardized residuals|"))) +
  ggtitle("Scale-Location") +
  theme_bw()
p3

# Cook's 
# not including those because of errors/warnings and difficulty to interpret them for binomial regression
p4 <- ggplot(m4, aes(seq_along(cooks.distance(m4)), cooks.distance(m4))) +
  geom_bar(stat="identity", position="identity") +
  xlab("Obs. Number") +
  ylab("Cook's distance") +
  ggtitle("Cook's distance") +
  theme_bw()
p4

# residuals vs leverage
# not including those because of errors/warnings and difficulty to interpret them for binomial regression
p5 <- ggplot(m4, aes(hatvalues(m4), scale(resid(m4)))) + 
  geom_point(aes(size = cooks.distance(m4)), na.rm=TRUE) +
  stat_smooth(method="loess", na.rm=TRUE) +
  xlab("Leverage") +
  ylab("Standardized Residuals") +
  ggtitle("Residual vs Leverage") +
  scale_size_continuous("Cook's Distance", range=c(1,5)) +
  theme_bw() +
  theme(legend.position="bottom")
p5

plot(resid(m4,type="pearson") ~ as.numeric(use_vs_foraging_sites$sex))

# Plot residuals vs explanatory variables
# just to understand, but we do not need to add them
p6 <- ggplot(m4, aes(as.numeric(use_vs_foraging_sites$sex), scale(resid(m4)))) +
  geom_point(na.rm=TRUE) +
  stat_smooth(method="loess", na.rm = TRUE) +
  xlab("Sex") +
  ylab("Standardized residuals") +
  ggtitle("Residuals vs. Sex") +
  theme_bw()
p6

p7 <- ggplot(m4, aes(as.numeric(use_vs_foraging_sites$body_mass_g), scale(resid(m4)))) +
  geom_point(na.rm=TRUE) +
  stat_smooth(method="loess", na.rm = TRUE) +
  xlab("Body mass (g)") +
  ylab("Standardized residuals") +
  ggtitle("Residuals vs. Body mass") +
  theme_bw()
p7

# random effects vs each of the fixed effect variables
(pred.sex <- sjPlot::plot_model(m4, type = "pred", terms = "sex", 
                                title = "Model prediction",
                                axis.title = c("Sex", "Predicted probability")))

(pred.bs <- sjPlot::plot_model(m4, type = "pred", terms = "body_mass_g", 
                               title = "Model prediction",
                               axis.title = c("Body mass (g)", "Predicted probability")))

## now check for approximate normality of random effects:
(re.qq <- sjPlot::plot_model(m4, type = "diag"))

re.qq.id <- re.qq$Animal_ID + 
  ggtitle("Q-Q plot: Individual") +
  theme(plot.title = element_text(size=11))
re.qq.area <- re.qq$area + 
  ggtitle("Q-Q plot: Foraging site") +
  theme(plot.title = element_text(size=11))

# final diagnostic figure
# grid: str res vs fitted, scale-location, predicted prob vs sex and bodysize,
# qqplot residuals area and ID
plot.fs <- ggpubr::ggarrange(pred.sex, pred.bs, p1, p3, re.qq.id, re.qq.area,
          labels = LETTERS[1:6],
          ncol = 2, nrow = 3, align = "v")

ggsave(plot.fs, filename = "output/plot_diagnostics_foraging_site.png", device = "png", 
       width = 14, height = 20, units = "cm", dpi = 300)

######################################################################
# 2) Verifying the influence of body mass and sex on interindividual 
#    variation considering individual use of polygons (landscape elements)
######################################################################

# Load data
use_vs_landuse_polygons <- read.csv("data/proportion_use_landcover_polygons.csv")

# Organize data
use_vs_landuse_polygons$Animal_ID <- as.factor(use_vs_landuse_polygons$Animal_ID)

# General statistics
# Body size
range(use_vs_landuse_polygons$body_mass_g)
mean(use_vs_landuse_polygons$body_mass_g); sd(use_vs_landuse_polygons$body_mass_g)

# proportion
range(p.pol <- use_vs_landuse_polygons$success/use_vs_landuse_polygons$sample.size)
mean(p.pol); sd(p.pol)

# Run linear models
# No effect model
m0.pol <- glm(cbind(success, failure) ~ 1, family = binomial, data = use_vs_landuse_polygons)
# Sex + weight model
m1.pol <- glm(cbind(success, failure) ~ sex + body_mass_g, family = binomial, data = use_vs_landuse_polygons)
# Sex + weight + ID as a fixed effect
m2.pol <- glm(cbind(success, failure) ~ sex + body_mass_g + Animal_ID, family = binomial, data = use_vs_landuse_polygons)
# Sex + weight + ID as a random intercept
m3.pol <- glmer(cbind(success, failure) ~ sex + body_mass_g + (1|Animal_ID), family = binomial, data = use_vs_landuse_polygons)
# In the models above, residual deviance is much higher than the degrees of freedom
# what means the data is overdispersed; we should add the polygon as a random intercept
m4.pol <- glmer(cbind(success, failure) ~ sex + body_mass_g + (1|Animal_ID) + (1|polygon), family = binomial, data = use_vs_landuse_polygons)
m5.pol <- glmer(cbind(success, failure) ~ 1 + (1|Animal_ID) + (1|polygon), family = binomial, data = use_vs_landuse_polygons)

# Comparing models
AICctab(m0.pol, m1.pol, m2.pol, m3.pol, m4.pol, m5.pol, delta = T, weights = T) # m0 wins

# Summary of models
summary(m1.pol)
summary(m2.pol)
summary(m3.pol)
summary(m4.pol)
summary(m5.pol)

# Checking for percentages explained
af <- anova(m2.pol)
afss <- af$"Deviance"
print(cbind(af, PctExplained = afss/sum(afss, na.rm = T)*100))

## PRODUCING DIAGNOSTIC PLOTS ##

# std residuals vs fitted values
plot(m4.pol, type=c("p", "smooth"))

p1.pol <- ggplot(m4.pol, aes(.fitted, scale(.resid))) +
  geom_point() +
  stat_smooth(method = "loess") + 
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Fitted values") +
  ylab("Standardized residuals") +
  ggtitle("Residual vs Fitted")+
  theme_bw()
p1.pol

# Scale-location plot
p3.pol <- ggplot(m4.pol, aes(.fitted, sqrt(abs(scale(resid(m4.pol))))))+
  geom_point(na.rm=TRUE) +
  stat_smooth(method="loess", na.rm = TRUE) +
  xlab("Fitted Value") +
  ylab(expression(sqrt("|Standardized residuals|"))) +
  ggtitle("Scale-Location") +
  theme_bw()
p3.pol

# random effects vs each of the fixed effect variables
(pred.sex.pol <- sjPlot::plot_model(m4.pol, type = "pred", terms = "sex", 
                                title = "Model prediction",
                                axis.title = c("Sex", "Predicted probability")))

(pred.bs.pol <- sjPlot::plot_model(m4.pol, type = "pred", terms = "body_mass_g", 
                               title = "Model prediction",
                               axis.title = c("Body mass (g)", "Predicted probability")))

## now check for approximate normality of random effects:
(re.qq.pol <- sjPlot::plot_model(m4.pol, type = "diag"))

re.qq.id.pol <- re.qq.pol$Animal_ID + 
  ggtitle("Q-Q plot: Individual") +
  theme(plot.title = element_text(size=11))
re.qq.pol <- re.qq.pol$polygon + 
  ggtitle("Q-Q plot: Habitat polygon") +
  theme(plot.title = element_text(size=11))

# final diagnostic figure
# grid: str res vs fitted, scale-location, predicted prob vs sex and bodysize,
# qqplot residuals area and ID
plot.pol <- ggpubr::ggarrange(pred.sex.pol, pred.bs.pol, p1.pol, p3.pol, 
                              re.qq.id.pol, re.qq.pol,
                              labels = LETTERS[1:6],
                              ncol = 2, nrow = 3, align = "v")

ggsave(plot.pol, filename = "output/plot_diagnostics_landuse_polygons.png", device = "png", 
       width = 14, height = 20, units = "cm", dpi = 300)


#######################################################################
# Old code with Gaussian fit and Arcsin Sqrt transformation - avoid using it!!

######################################################################
# 1) Verifying the influence of body mass and sex on interindividual 
#    variation considering individual use of foraging sites (kernel 50%)
######################################################################

# Load data
use_vs_foraging_sites <- read.table("data/test_for_bodysize_sex_foraging_sites.csv", header = T, sep = "\t", dec = ",")

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
AICctab(m0, m1, m2, m3, delta = T, weights = T) # m0 wins

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
use_vs_poligons <- read.table("data/test_for_bodysize_sex_polygons.csv", header = T, sep = ";", dec = ",")

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
# FS1-FS31: Arcsine of the square root of the 
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

