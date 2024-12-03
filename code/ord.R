
# After clustering and finding proper number of clusters of the data. We can start conducting ordination. 


##### Clean all things #########################################################
rm(list = ls())
cat("\014")  # This is equivalent to pressing Ctrl+L
if (!is.null(dev.list())) dev.off()
gc()
################################################################################

# Unconstrained Ordination ####################################################

# ~ PCA (principal component analysis) ----
# ** Assumes linear relation and uses Euclidean distance **
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)

# Source additional functions that will be used later in this
# chapter. Our scripts assume that files to be read are in
# the working directory
source("code/fun.R")


# Load the required packages

# Load the Doubs data
load("data/Doubs.RData")
# Remove empty site 8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]


# Load the oribatid mite data
load('data/mite.RData')

# A reminder of the content of the env dataset
summary(env)
dim(env)
# PCA based on a correlation matrix
# Argument scale=TRUE calls for a standardization of the variables

env.pca <- rda(env, scale = TRUE)
env.pca # 11 variables --> 11 components (产生11个新坐标/PC)

summary(env.pca) # Default scaling 2
summary(env.pca, scaling = 1)


# test for number of components equal to number of variables
env_sub <- env[ , c(2, 4, 5, 6, 7)]
env_sub.pca <- rda(env_sub, scale = TRUE)
env_sub.pca # 5 variables --> 5 components (产生5个新坐标/PC)
summary(env_sub.pca) # default scaling 2
summary(env_sub.pca, scaling = 1)

# plot 画图

# Eigenvalues
(ev <- env.pca$CA$eig)
# Scree plot and broken stick model
screeplot(env.pca, bstick = TRUE, npcs = length(env.pca$CA$eig))

# Plots using biplot.rda
par(mfrow = c(1, 2))
biplot(env.pca, scaling = 1, main = "PCA - scaling 1")
biplot(env.pca, main = "PCA - scaling 2") # Default scaling 2


# Plots using cleanplot.pca
# A rectangular graphic window is needed to draw the plots together
par(mfrow = c(1, 2))
cleanplot.pca(env.pca, scaling = 1, mar.percent = 0.08)
cleanplot.pca(env.pca, scaling = 2, mar.percent = 0.04)

dev.off()


# ~~~ PCA on Transformed Species Data ----
# Hellinger pre-transformation of the species data
spe.h <- decostand(spe, "hellinger")
(spe.h.pca <- rda(spe.h))
# Scree plot and broken stick model
screeplot(
  spe.h.pca,
  bstick = TRUE,
  npcs = length(spe.h.pca$CA$eig)
) 
# PCA biplots

spe.pca.sc1 <- vegan::scores(spe.h.pca, display = "species", scaling = 1)
spe.pca.sc2 <- vegan::scores(spe.h.pca, display = "species", scaling = 2)
par(mfrow = c(1, 2))
cleanplot.pca(spe.h.pca, scaling = 1, mar.percent = 0.06)
cleanplot.pca(spe.h.pca, scaling = 2, mar.percent = 0.06)

dev.off()
# A posteriori projection of environmental variables in a PCA
# A PCA scaling 2 plot is produced in a new graphic window.
biplot(spe.h.pca, main = "PCA fish abundances-scaling 2") 
# Scaling 2 is default
(spe.h.pca.env <- envfit(spe.h.pca, env, scaling = 2)) 
# Plot significant variables with a user-selected colour
plot(spe.h.pca.env, p.max = 0.05, col = 3)
# This has added the significant environmental variables to the
# last biplot drawn by R.
# BEWARE: envfit() must be given the same scaling as the plot to
# which its result is added



# ~ CA (Correspondence Analysis) (self learning) ----
# categorical data, often used with contingency tables
# use chi-squared distance

# ~ PCoA (Principal Coordinate Analysis) ----
# Does not assume linear relation, and allows many distance matrices
spe.bray <- vegdist(spe) # bray-curtis distance
spe.b.pcoa <- cmdscale(spe.bray, k = (nrow(spe)-1), eig = TRUE)
# Plot of the sites
ordiplot(vegan::scores(spe.b.pcoa, choices = c(1, 2)),
         type = "t",
         main = "PCoA with species weighted averages")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)



# Add weighted average projection of species
spe.wa <- wascores(spe.b.pcoa$points[, 1:2], spe)
text(spe.wa, rownames(spe.wa), cex = 0.7, col = "red")
# A posteriori projection of environmental variables
(spe.b.pcoa.env <- envfit(spe.b.pcoa, env))
# Plot significant variables with a user-selected colour
plot(spe.b.pcoa.env, p.max = 0.05, col = 3)



# using pcoa() from ape package
spe.h.pcoa <- pcoa(dist(spe.h)) # hellinger distance
# Biplots
par(mfrow = c(1, 2))
# First biplot: Hellinger-transformed species data
biplot.pcoa(spe.h.pcoa, spe.h, dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)
# Second biplot: standardized Hellinger-transformed species data
spe.std <- scale(spe.h)
biplot.pcoa(spe.h.pcoa, spe.std, dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

dev.off()



# Third biplot: standardized Hellinger-transformed species data;
# only four species in the plot (figure not shown in the book)
biplot.pcoa(spe.h.pcoa, spe.h[, c(2, 5, 11, 21)], dir.axis1 = -1)
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# PCoA on a Hellinger distance matrix
is.euclid(dist(spe.h))
summary(spe.h.pcoa)
spe.h.pcoa$values

# PCoA on a percentage difference dissimilarity matrix
is.euclid(spe.bray)
spe.bray.pcoa <- pcoa(spe.bray)
spe.bray.pcoa$values # Observe eigenvalues 18 and following

# PCoA on the square root of a percentage difference
# dissimilarity matrix (aka Bray-Curtis)
is.euclid(sqrt(spe.bray))
spe.braysq.pcoa <- pcoa(sqrt(spe.bray))
spe.braysq.pcoa$values  # Observe the eigenvalues

# PCoA on a percentage difference dissimilarity matrix with
# Lingoes correction
spe.brayl.pcoa <- pcoa(spe.bray, correction = "lingoes")
spe.brayl.pcoa$values      # Observe the eigenvalues, col. 1 and 2

# PCoA on a percentage difference dissimilarity matrix with
# Cailliez correction
spe.brayc.pcoa <- pcoa(spe.bray, correction = "cailliez")
spe.brayc.pcoa$values      # Observe the eigenvalues, col. 1 and 2


# ~ NMDS (Nonmetric Multidimensional Scaling) ----
# Does not assume linear relation, and allows many distance matrices
# For the case when exact distance is not important


# Canonical Ordination ########################################################

# ~ RDA (redundancy analysis) ----
library(ade4)
library(adegraphics)
library(adespatial)
library(vegan)
library(vegan3d)
library(MASS)
library(ellipse)
library(FactoMineR)
library(rrcov)


# Source additional functions that will be used later in this
# Chapter. Our scripts assume that files to be read are in
# the working directory.
source('code/fun.R')

# Load the Doubs data. The file Doubs.Rdata is assumed to be in 
# the working directory
load("data/Doubs.RData")  
# Remove empty site 8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]

# Set aside the variable 'dfs' (distance from the source) for 
# later use
dfs <- env[, 1]

# Remove the 'dfs' variable from the env data frame
env2 <- env[, -1]
# Recode the slope variable (slo) into a factor (qualitative) 
# variable to show how these are handled in the ordinations
slo2 <- rep(".very_steep", nrow(env))
slo2[env$slo <= quantile(env$slo)[4]] <- ".steep"
slo2[env$slo <= quantile(env$slo)[3]] <- ".moderate"
slo2[env$slo <= quantile(env$slo)[2]] <- ".low"
slo2 <- factor(slo2, 
               levels = c(".low", ".moderate", ".steep", ".very_steep"))
table(slo2)
# Create an env3 data frame with slope as a qualitative variable
env3 <- env2
env3$slo <- slo2

# Create two subsets of explanatory variables
# Subset 1: Physiography (upstream-downstream gradient)
envtopo <- env2[, c(1 : 3)]
names(envtopo)
# Subset 2: Water quality
envchem <- env2[, c(4 : 10)]
names(envchem)

# Hellinger-transform the species data set
spe.hel <- decostand(spe, "hellinger")

# simpleRDA <- rda(Y,X,W)
# 
# formulaRDA <- rda(Y ~ var1 + factorA + var2*var3 + Condition(var4),
#                   data = XWdata)

(spe.rda <- rda(spe.hel ~ ., env3))
summary(spe.rda) # Scaling 2 (default)
coef(spe.rda)

# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(spe.rda)$r.squared)
# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(spe.rda)$adj.r.squared)

# Scaling 1
plot(spe.rda, 
     scaling = 1, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ env3 - scaling 1 - lc scores"
)


spe.sc1 <-
  vegan::scores(spe.rda, 
         choices = 1:2, 
         scaling = 1, 
         display = "sp"
  )
arrows(0, 0, 
       spe.sc1[, 1] * 0.92, 
       spe.sc1[, 2] * 0.92, 
       length = 0, 
       lty = 1, 
       col = "red"
)

# Scaling 2
plot(spe.rda, 
     display = c("sp", "lc", "cn"), 
     main = "Triplot RDA spe.hel ~ env3 - scaling 2 - lc scores"
)
    spe.sc2 <-
  vegan::scores(spe.rda, 
         choices = 1:2, 
         display = "sp"
  )
arrows(0, 0, 
       spe.sc2[, 1] * 0.92, 
       spe.sc2[, 2] * 0.92,
       length = 0,
       lty = 1,
       col = "red"
)


# Select species with goodness-of-fit at least 0.6 in the 
# ordination plane formed by axes 1 and 2
spe.good <- goodness(spe.rda)
sel.sp <- which(spe.good[, 2] >= 0.6)
# Triplots with homemade function triplot.rda(), scalings 1 and 2
source('code/fun.R')
triplot.rda(spe.rda, 
            site.sc = "lc", 
            scaling = 1, 
            cex.char2 = 0.7, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05, 
            select.spe = sel.sp
)
triplot.rda(spe.rda,
            site.sc = "lc", 
            scaling = 2, 
            cex.char2 = 0.7, 
            pos.env = 3, 
            pos.centr = 1, 
            mult.arrow = 1.1, 
            mar.percent = 0.05, 
            select.spe = sel.sp
)

# permutation test for RDA results
## Global test of the RDA result
anova(spe.rda, permutations = how(nperm = 999))
## Tests of all canonical axes
anova(spe.rda, by = "axis", permutations = how(nperm = 999))


# Apply Kaiser-Guttman criterion to residual axes
spe.rda$CA$eig[spe.rda$CA$eig > mean(spe.rda$CA$eig)]


# ~~~ Conduct RDA using simple fabricated data set ----
# environmental data
env_data <- data.frame(
  Temperature = rnorm(10, mean = 20, sd = 5),
  pH = rnorm(10, mean = 7, sd = 0.5),
  Nutrients = rnorm(10, mean = 50, sd = 10)
)

# species data
species_data <- data.frame(
  Species1 = rpois(10, lambda = 5),
  Species2 = rpois(10, lambda = 10),
  Species3 = rpois(10, lambda = 15)
)

rda_result <- rda(species_data ~ ., data = env_data)
summary(rda_result)
plot(rda_result)

# choose RDA or CCA ----
dca_result <- decorana(species_data)
summary(dca_result) # DCA1, Axis length = 0.71978

# so we shall use RDA




# ~ CCA ( Canonical Correspondence Analysis) ----


# ~ LDA (Linear Discriminant Analysis) ----