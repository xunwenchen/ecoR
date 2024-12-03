# Data transformation of ecological data ----

##### clean all things ########################################################
rm(list = ls())
cat("\014")  # This is equivalent to pressing Ctrl+L
if (!is.null(dev.list())) dev.off()
gc()
###############################################################################


load('data/Doubs.RData') # load the fish data from: Borcard D, Gillet F, Legendre P (2018) Numerical Ecology with R, 2nd ed. 2018 edition. Springer, New York, NY

# check the data
head(spe)
head(env)
head(spa) # coordinates of sampling locations

## Simple transformations of spe data ----
# Partial view of the raw data (abundance codes)
spe[1:5, 2:4]

# Transform to presence-absence (1-0)
library(vegan)
?decostand() # check how to use decostand()

spe.pa <- decostand(spe, method = "pa") 
# 'pa' - convert abundance to presence/absence (1/0)
spe.pa[1:5, 2:4]

## Standardization by columns (species)
# Scale abundances by dividing them by the maximum value of 
# each species
# Note: MARGIN = 2 (column, default value) for argument "max"
spe.scal <- decostand(spe, "max") 

spe.scal[1:5, 2:4]

# Display the maximum in each transformed column
apply(spe.scal, 2, max)


# Scale abundances by dividing them by the species totals
# (relative abundance by species)
# Note: here, override the default MARGIN = 1 argument of "total"
spe.relsp <- decostand(spe, "total", MARGIN = 2)
spe.relsp[1:5, 2:4]
# Display the sum by column
# Classical: apply(spe.relsp, 2, sum)
colSums(spe.relsp)


## Environmental data ----

summary(env)
skimr::skim(env) # call skimr package without charing it into the Environment

## Exploring pairwise relationship ----
## Scatter plots for all pairs of environmental variables
# Bivariate plots with histograms on the diagonal and smooth 
# fitted curves
source('code/fun.R')
pairs(env, 
      panel = panel.smooth, 
      diag.panel = panel.hist,
      main = "Bivariate Plots with Histograms and Smooth Curves"
)


## Simple transformation of an environmental variable
range(env$slo)
# Log-transformation of the slope variable (y = ln(x))
# Compare histograms and boxplots of the raw and transformed values
par(mfrow = c(2, 2))
hist(env$slo, 
     col = "bisque", 
     right = FALSE
)
hist(log(env$slo), 
     col = "lightgreen", 
     right = FALSE, 
     main = "Histogram of ln(env$slo)"
)
boxplot(env$slo, 
        col = "bisque", 
        main = "Boxplot of env$slo", 
        ylab = "env$slo"
)
boxplot(log(env$slo), 
        col = "lightgreen", 
        main = "Boxplot of ln(env$slo)",
        ylab = "log(env$slo)"
)

dev.off()

## Standardization of all environmental variables
# Center and scale = standardize the variables (z-scores)
env.z <- decostand(env, "standardize")
apply(env.z, 2, mean)   # means = 0
apply(env.z, 2, sd) # standard deviations = 1

# Same standardization using the scale() function (which returns 
# a matrix)
env.z <- as.data.frame(scale(env))
env.z

# Conclusions ----
# Choosing methods of data transformation is data- and objective-dependent

# END ----
