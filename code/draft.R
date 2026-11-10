
save <- F



# 1. introduction ####

# draft codes here 
# Install and load the vegan package

library(vegan)

# Example data: rows are samples, columns are species counts
your_data <- matrix(c(5, 3, 0, 2, 1, 4, 0, 0, 3, 2, 1, 0, 4, 5, 2, 1), nrow = 4, byrow = TRUE)

# Set the rarefaction threshold
rarefaction_threshold <- 3

# Perform rarefaction
rarefied_data <- rarefy(your_data, sample = rarefaction_threshold)

# Print the rarefied data
print(rarefied_data)

# Plot rarefaction curves
rarecurve(your_data, step = 1, sample = rarefaction_threshold, col = "blue", cex = 0.6)


df1 <-data.frame(species = c('a', 'b', 'c'), value = c(1, 2, 3))
df2 <-data.frame(species = c('a', 'b', 'c'), value = c(4, 5, 6))
df <- merge(df1, df2, by = 'species')
df


# Construct a matrix, columns are Sp1 to Sp5, and rows are SiteA to SiteC, and fill the matrix with order of left to right, and top to bottom for the following numbers: first row: 11, 0, 7, 8, 0; second row: 24, 37, 5, 18, 11, 20; third row: 45, 19, 33, 21

mat <- matrix(c(11, 0, 7, 8, 0, 24, 37, 5, 18, 1, 20, 45, 19, 33, 21), nrow = 3, byrow = TRUE)
rownames(mat) <- c('SiteA', 'SiteB', 'SiteC')
colnames(mat) <- c('sp1', 'sp2', 'sp3', 'sp4', 'sp5')
head(mat)
# calculate euclidean distance
library(vegan)
dist_euc <- vegdist(mat, method = "euclidean")
dist_euc


# verification of the euclidean distance calculation using manual calculation
# SiteA to SiteB
sqrt((11-24)^2 + (0-37)^2 + (7-5)^2 + (8-18)^2 + (0-1)^2)
# SiteA to SiteC
sqrt((11-20)^2 + (0-45)^2 + (7-19)^2 + (8-33)^2 + (0-21)^2)
# SiteB to SiteC
sqrt((20-24)^2 + (37-45)^2 + (5-19)^2 + (18-33)^2 + (1-21)^2)


# bray-curtis distance
dist_bc <- vegdist(mat, method = "bray")
dist_bc

# ~~~ PCA Example borrowed from StatQuest by Joshua Starmer ----
## StatQuest by Joshua Starmer is Appreciated and well acknowledged.
## In this example, the data is in a matrix called
## data.matrix
## columns are individual samples (i.e. cells)
## rows are measurements taken for all the samples (i.e. genes)
## Just for the sake of the example, here's some made up data...
data.matrix <- matrix(nrow=100, ncol=10)
colnames(data.matrix) <- c(
  paste("wt", 1:5, sep=""),
  paste("ko", 1:5, sep=""))
rownames(data.matrix) <- paste("gene", 1:100, sep="")
for (i in 1:100) {
  wt.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  ko.values <- rpois(5, lambda=sample(x=10:1000, size=1))
  
  data.matrix[i,] <- c(wt.values, ko.values)
}
head(data.matrix)
dim(data.matrix)

pca <- prcomp(t(data.matrix), scale=TRUE) 

## plot pc1 and pc2
plot(pca$x[,1], pca$x[,2])

## make a scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

## now make a fancy looking plot that shows the PCs and the variation:
library(ggplot2)

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

pca$rotation[top_10_genes,1] ## show the scores (and +/- sign)

#######
##
## NOTE: Everything that follow is just bonus stuff.
## It simply demonstrates how to get the same
## results using "svd()" (Singular Value Decomposition) or using "eigen()"
## (Eigen Decomposition).
##
#######

############################################
##
## Now let's do the same thing with svd()
##
## svd() returns three things
## v = the "rotation" that prcomp() returns, this is a matrix of eigenvectors
##     in other words, a matrix of loading scores
## u = this is similar to the "x" that prcomp() returns. In other words,
##     sum(the rotation * the original data), but compressed to the unit vector
##     You can spread it out by multiplying by "d"
## d = this is similar to the "sdev" value that prcomp() returns (and thus
##     related to the eigen values), but not
##     scaled by sample size in an unbiased way (ie. 1/(n-1)).
##     For prcomp(), sdev = sqrt(var) = sqrt(ss(fit)/(n-1))
##     For svd(), d = sqrt(ss(fit))
##
############################################

svd.stuff <- svd(scale(t(data.matrix), center=TRUE))

## calculate the PCs
svd.data <- data.frame(Sample=colnames(data.matrix),
                       X=(svd.stuff$u[,1] * svd.stuff$d[1]),
                       Y=(svd.stuff$u[,2] * svd.stuff$d[2]))
svd.data

## alternatively, we could compute the PCs with the eigen vectors and the
## original data
svd.pcs <- t(t(svd.stuff$v) %*% t(scale(t(data.matrix), center=TRUE)))
svd.pcs[,1:2] ## the first to principal components

svd.df <- ncol(data.matrix) - 1
svd.var <- svd.stuff$d^2 / svd.df
svd.var.per <- round(svd.var/sum(svd.var)*100, 1)

ggplot(data=svd.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", svd.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", svd.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("svd(scale(t(data.matrix), center=TRUE)")

############################################





# 2. Load data  ####
# ~~ 2.1. Load fish data ####
##
## Now let's do the same thing with eigen()
##
## eigen() returns two things...
## vectors = eigen vectors (vectors of loading scores)
##           NOTE: pcs = sum(loading scores * values for sample)
## values = eigen values
##
############################################
cov.mat <- cov(scale(t(data.matrix), center=TRUE))
dim(cov.mat)

## since the covariance matrix is symmetric, we can tell eigen() to just
## work on the lower triangle with "symmetric=TRUE"
eigen.stuff <- eigen(cov.mat, symmetric=TRUE)
dim(eigen.stuff$vectors)
head(eigen.stuff$vectors[,1:2])

eigen.pcs <- t(t(eigen.stuff$vectors) %*% t(scale(t(data.matrix), center=TRUE)))
eigen.pcs[,1:2]

eigen.data <- data.frame(Sample=rownames(eigen.pcs),
                         X=(-1 * eigen.pcs[,1]), ## eigen() flips the X-axis in this case, so we flip it back
                         Y=eigen.pcs[,2]) ## X axis will be PC1, Y axis will be PC2
eigen.data

eigen.var.per <- round(eigen.stuff$values/sum(eigen.stuff$values)*100, 1)

ggplot(data=eigen.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", eigen.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", eigen.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("eigen on cov(t(data.matrix))")

if(save) {
  
  ggsave('xx.jpg')
  
}


# I have a df, in R, columns are N, P, K, and treatment, rows are 1-100, representing sample IDs. pls plot boxplots for N, P, K, and show ANOVA test results across treatments using ggplot2. Thanks
set.seed(123)
df <- data.frame(
  N = rnorm(100, mean=50, sd=10),
  P = rnorm(100, mean=30, sd=5),
  K = rnorm(100, mean=20, sd=3),
  treatment = sample(c('A', 'B', 'C'), 100, replace=TRUE)
)
library(ggplot2)
library(reshape2)
df_melt <- melt(df, id.vars = 'treatment')
ggplot(df_melt, aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = 'free_y') +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", position=position_dodge(width=0.75)) +
  stat_compare_means(aes(group=treatment), label = "p.format") +
  theme_minimal() +
  labs(title="Boxplots of N, P, K by Treatment", y="Value", x="Treatment")
if(save) {
  ggsave('boxplots_NPK_treatment.jpg')
}
# end of the script ####


# inport towngas_charge.xlsx
library(readxl)
towngas_data <- read_excel("data/towngas_charge.xlsx")
head(towngas_data)

# plot towngas_data/charge over month but order with 'order' column
library(ggplot2)
ggplot(towngas_data, aes(x = reorder(date, order), 
                         y = charge, group = 1)) + 
  geom_line() + 
  geom_point() + 
  theme_bw() + 
  labs( title = "Towngas Charge Over Time", x = "Date", y = "Charge (HKD)" )+
  #rotate date labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  # highlight the last data point in indianred
  geom_point(data = tail(towngas_data, 1), aes(x = date, y = charge), color = "indianred", size = 3)

