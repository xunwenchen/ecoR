# Association measures and matrices ----
# Most methods of multivariate analysis, in particular ordination and clustering techniques, are explicitly or implicitly based on the comparison of all possible pairs of objects or descriptors. Since the ***subsequent analyses are done on association matrices***, the choice of an appropriate measure is crucial. 

##### clean all things ####################################
rm(list = ls())
cat("\014")  # This is equivalent to pressing Ctrl+L
if (!is.null(dev.list())) dev.off()
gc()
###########################################################

# Major categories of association matrices ----

# Q mode: comparing objects ---- 
# checking dissimilarity

# Load the required packages 
library(ade4) 
library(adespatial) 
library(vegan) 
library(gclus) 
library(cluster) 
library(FD) 

# Load  data 
load("data/Doubs.Rdata") 
spe # site 8 is empty
rowSums(spe) # check which site had no observation


# Remove empty site 8 
spe <- spe[-8, ] 
env <- env[-8, ] 
spa <- spa[-8, ] 


## Q-mode dissimilarity and distance measures for  
## (semi-)quantitative data 
# Percentage difference (aka Bray-Curtis) dissimilarity matrix 
# on raw species data 
spe.db <- vegdist(spe)  # method = "bray" (default) 
spe.db

source('code/fun.R')
coldiss(spe.db, byrank = FALSE, diag = TRUE) # Is there a site effect?


# Check how spe.db looks like
spe.db.mt <- as.matrix(spe.db)
head(spe.db.mt)
View(spe.db.mt)

# Percentage difference (aka Bray-Curtis) dissimilarity matrix 
# on log-transformed abundances. log1p(x) is log(x+1), for small x close to zero
spe.dbln <- vegdist(log1p(spe)) 
spe.dbln

log1p(87) == log(87+1)

coldiss(spe.dbln, byrank = FALSE, diag = TRUE) 

# Chord distance matrix 
spe.dc <- dist.ldc(spe, "chord") # package 'adespatial'

# Alternate, two-step computation in vegan for chord distance matrix: 
# spe.norm <- decostand(spe, "nor") 
# spe.dc <- dist(spe.norm)
# spe.dc

# Hellinger distance matrix
spe.dh <- dist.ldc(spe) # Hellinger is the default distance
spe.dh

# Alternate, two-step computation in vegan:
# spe.hel <- decostand(spe, "hel")
# spe.dh <- dist(spe.hel)
# spe.dh

# Log-chord distance matrix
spe.logchord <- dist.ldc(spe, "log.chord")
spe.logchord

# Alternate, three-step computation in vegan:
spe.ln <- log1p(spe)
spe.ln.norm <- decostand(spe.ln, "nor")
spe.logchord <- dist(spe.ln.norm)
spe.logchord

coldiss(spe.logchord, byrank = FALSE, diag = TRUE)


# ~~ Mixed Types ------------
# Including Categorical
# (Qualitative Multiclass) Variables
# Fictitious data for Gower (S15) index
# Random normal deviates with zero mean and unit standard deviation
var.g1 <- rnorm(30, 0, 1)
# Random uniform deviates from 0 to 5
var.g2 <- runif(30, 0, 5)
# Factor with 3 levels (10 objects each)
var.g3 <- gl(3, 10, labels = c("A", "B", "C"))
# Factor with 2 levels,  orthogonal to var.g3
var.g4 <- gl(2, 5, 30, labels = c("D", "E"))

dat2 <- data.frame(var.g1, var.g2, var.g3, var.g4)
summary(dat2)
head(dat2)

# Computation of a matrix of Gower dissimilarity using 
# function daisy()

# Complete data matrix (4 variables)
dat2.S15 <- daisy(dat2, "gower")
range(dat2.S15)
coldiss(dat2.S15, diag = TRUE)

# Data matrix with the two orthogonal factors only 只用两个因子
dat2[, 3:4]
dat2partial.S15 <- daisy(dat2[, 3:4], "gower")
coldiss(dat2partial.S15, diag = TRUE)
head(as.matrix(dat2partial.S15))

# What are the dissimilarity values in the dat2partial.S15 matrix?
levels(factor(dat2partial.S15))

# # Alternate, computation of a matrix of Gower dissimilarity using
# # function gowdis() of package FD
# ?gowdis
# dat2.S15.2 <- gowdis(dat2)
# range(dat2.S15.2)
# coldiss(dat2.S15.2, diag = TRUE)
# 
# # Data matrix with the two orthogonal factors only
# dat2partial.S15.2 <- gowdis(dat2[ , 3:4])
# coldiss(dat2partial.S15.2, diag = TRUE)
# head(as.matrix(dat2partial.S15.2))
# 
# # What are the dissimilarity values in the dat2partial.S15.2 matrix? 
# levels(factor(dat2partial.S15.2))


# Euclidean distance 欧几里得距离
# contruct matrix
mat <- matrix(c(11, 0, 7, 8, 0, 24, 37, 5, 18, 1, 20, 45, 19, 33, 21), nrow = 3, byrow = TRUE)
rownames(mat) <- c('SiteA', 'SiteB', 'SiteC')
colnames(mat) <- c('sp1', 'sp2', 'sp3', 'sp4', 'sp5')
head(mat)
# calculate euclidean distance
dist_euc <- vegdist(mat, method = "euclidean")
dist_euc

# verification using manual calculation 手动计算验证
# SiteA to SiteB
sqrt((11-24)^2 + (0-37)^2 + (7-5)^2 + (8-18)^2 + (0-1)^2)
# SiteA to SiteC
sqrt((11-20)^2 + (0-45)^2 + (7-19)^2 + (8-33)^2 + (0-21)^2)
# SiteB to SiteC
sqrt((20-24)^2 + (37-45)^2 + (5-19)^2 + (18-33)^2 + (1-21)^2)


# bray-curtis distance BC距离
dist_bc <- vegdist(mat, method = "bray")
dist_bc

# 不同距离的优劣介绍：
# check https://r.qcbs.ca/workshop09/book-en/types-of-distance-coefficients.html
# for pros and cons of different dist.


# R mode: comparing variables ---- 
# checking covariance or correlation coefficient

View(spe)
spe.t <- t(spe)
View(spe.t)


# ~~ Species abundance data ----
# Chi-square pre-transformation followed by Euclidean distance
# chi-square sum(Oi-Ei)^2/Ei # 两个值的差异
spe.t.chi <- decostand(spe.t, "chi.square")
class(spe.t.chi)
spe.t.D16 <- dist(spe.t.chi) # D16 means chi-square (others: Jaccard (S7), Sørensen (S8) and Ochiai (S14) )
coldiss(spe.t.D16, diag = TRUE) # Effects of one species/var. on another species/var.?


# ~~ Species 0/1 data ----

# Jaccard index on fish presence-absence
spe.t.S7 <- vegdist(spe.t, "jaccard", binary = TRUE)
coldiss(spe.t.S7, diag = TRUE)


# chi-square vs Jaccard matrix pattern 是否一致？


#########################################################
#########################################################

# Pearson r linear correlation among environmental variables
View(env)
env.pearson <- cor(env) # default method = "pearson"
round(env.pearson, 2)
# Reorder the variables prior to plotting
env.o <- order.single(env.pearson)

# pairs() is a function to plot a matrix of bivariate scatter 
# plots.panelutils.R is a set of functions that add useful 
# features topairs().
pairs(env[ ,env.o], 
      lower.panel = panel.smooth, 
      upper.panel = panel.cor,
      diag.panel = panel.hist, 
      main = "Pearson Correlation Matrix")


# Exercise: Use 'data/data.xlsx'(sheet = 2) to conduct pairwise pearson correlation analysis using the above codes
# 练习：试试用data.xlsx (sheet = 2)，模仿上述代码做参数之间的相关性分析
library(readxl)
ww_data <- read_xlsx('data/data.xlsx', sheet = 2)
View(ww_data)
colnames(ww_data)

# answer is in 'ex_ans.R'. Provide later.

# Take-home message ############################################################

# Q mode -- create dist among sites/treatments etc. 
# for later analyses like PCA, PCoA, NMDS... 

# R mode -- create correlations among variables such as 
# chemical concentration, soil organic matter...

# Select dist based on data characteristics and objectives
## 1/0, abundance, mixed-type, double-zero...

# Understand pros and cons of different dist






