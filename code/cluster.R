# Cluster analysis 聚类分析 ----
##### Clean all things #########################################################
rm(list = ls())
cat("\014")  # This is equivalent to pressing Ctrl+L
if (!is.null(dev.list())) dev.off()
gc()
################################################################################
library(ade4)
library(adespatial)
library(vegan)
library(gclus)
library(cluster)
library(pvclust)
library(RColorBrewer)
library(labdsv)
library(rioja)
library(indicspecies)
library(mvpart) # removed from CRAN
library(MVPARTwrap) # removed from CRAN
library(dendextend)
library(vegclust)
library(colorspace)
library(agricolae)
library(picante)

source('code/fun.R')


# Load the data
# File Doubs.Rdata is assumed to be in the working directory
load("data/Doubs.Rdata")
# Remove empty site 8
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]
latlong <- latlong[-8, ]

# 4.3.1. Single linkage ----
# Compute matrix of chord distance among sites
spe.norm <- decostand(spe, "normalize")
spe.ch <- vegdist(spe.norm, "euc")

# Attach site names to object of class 'dist'
attr(spe.ch, "labels") <- rownames(spe)

# Compute single linkage agglomerative clustering
spe.ch.single <- hclust(spe.ch, method = "single")
# Plot a dendrogram using the default options
plot(spe.ch.single,
     labels = rownames(spe),
     main = "Chord - Single linkage")

# pros and cons
# merge group based on smallest distance
# simple, but can create snake-like pattern

# 4.3.2. Complete linkage ----
# Compute complete-linkage agglomerative clustering
spe.ch.complete <- hclust(spe.ch, method = "complete")
plot(spe.ch.complete,
     labels = rownames(spe),
     main = "Chord - Complete linkage")

# pros and cons
# merge group based on furthest distance
# create more compact clusters, but senstive to outliers
# compact many clusters into one big cluster

# 4.4. Average agglomerative clustering ----
# ~~ Compute UPGMA agglomerative clustering ----
spe.ch.UPGMA <- hclust(spe.ch, method = "average")
plot(spe.ch.UPGMA,
     labels = rownames(spe),
     main = "Chord - UPGMA")

# pros and cons
# merges clusters based on the average distance 
# more balanced between single and complete linkage


# ~~ Compute centroid clustering 矩心：“物体”的转动中心----
spe.ch.centroid <- hclust(spe.ch, method = "centroid")
plot(spe.ch.centroid,
     labels = rownames(spe),
     main = "Chord - Centroid")

# pros and cons
# Merge the pair of clusters whose centroids are closest.
# simple, but sensitive to outliers and initial conditions


# 4.5. Compute Ward's minimum variance clustering ----
spe.ch.ward <- hclust(spe.ch, method = "ward.D2")
plot(spe.ch.ward,
     main = "Chord - Ward")

# pros and cons
# reduce within-cluster variance, create roughly equal distance 
# effective for spherical clusters.
# but assumes clusters are spherical.

# 4.6. Flexible clustering ----
# Compute Ward's minimum variance clustering
# Compute beta-flexible clustering using cluster::agnes()
# beta = -0.25
spe.ch.beta2 <- agnes(spe.ch, method = "flexible", 
                      par.method = 0.625)
# Change the class of agnes object
class(spe.ch.beta2) # [1] "agnes" "twins"
spe.ch.beta2 <- as.hclust(spe.ch.beta2)
class(spe.ch.beta2) # [1] "hclust"
plot(spe.ch.beta2,
     labels = rownames(spe),
     main = "Chord - Beta-flexible (beta=-0.25)")

# pros and cons
# allows for adjusting how clusters are formed.
# Pros: Adaptable to different data structures and can be fine-tuned.
# Cons: Requires careful parameter tuning and understanding of the data.



###############################################################################
###############################################################################



# 4.7. Finding a proper clustering method #####################################
# 寻找合适的聚类方法

# 4.7.2. Using Cophenetic correlation to validate clustering ----
# ** original distance vs clustered distance **
# 根据共表象相关
# Single linkage clustering
spe.ch.single.coph <- cophenetic(spe.ch.single)
cor(spe.ch, spe.ch.single.coph)
# Complete linkage clustering
spe.ch.comp.coph <- cophenetic(spe.ch.complete)
cor(spe.ch, spe.ch.comp.coph)
# Average clustering
spe.ch.UPGMA.coph <- cophenetic(spe.ch.UPGMA)
cor(spe.ch, spe.ch.UPGMA.coph)
# Ward clustering
spe.ch.ward.coph <- cophenetic(spe.ch.ward)
cor(spe.ch, spe.ch.ward.coph)

# plot them: 
# Shepard-like diagrams
par(mfrow = c(2, 2))
plot(
  spe.ch,
  spe.ch.single.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("Single linkage", paste("Cophenetic correlation =",
                                   round(cor(spe.ch, spe.ch.single.coph), 3)))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.single.coph), col = "red")
plot(
  spe.ch,
  spe.ch.comp.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("Complete linkage", paste("Cophenetic correlation =",
                                     round(cor(spe.ch, spe.ch.comp.coph), 3)))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.comp.coph), col = "red")
plot(
  spe.ch,
  spe.ch.UPGMA.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, sqrt(2)),
  main = c("UPGMA", paste("Cophenetic correlation =",
                          round(cor(spe.ch, spe.ch.UPGMA.coph), 3)))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.UPGMA.coph), col = "red")
plot(
  spe.ch,
  spe.ch.ward.coph,
  xlab = "Chord distance",
  ylab = "Cophenetic distance",
  asp = 1,
  xlim = c(0, sqrt(2)),
  ylim = c(0, max(spe.ch.ward$height)),
  main = c("Ward", paste("Cophenetic correlation =",
                         round(cor(spe.ch, spe.ch.ward.coph), 3)))
)
abline(0, 1)
lines(lowess(spe.ch, spe.ch.ward.coph), col = "red")
dev.off()

# use Gower method to judge which clustering is better 
# Gower (1983) distance. The smallest the best
(gow.dist.single <- sum((spe.ch - spe.ch.single.coph) ^ 2))
(gow.dist.comp <- sum((spe.ch - spe.ch.comp.coph) ^ 2))
(gow.dist.UPGMA <- sum((spe.ch - spe.ch.UPGMA.coph) ^ 2))
(gow.dist.ward <- sum((spe.ch - spe.ch.ward.coph) ^ 2))


# 4.7.3.1. Fusion level ----

par(mfrow = c(2, 2))
# Plot the fusion level values of the complete linkage clustering
plot(
  spe.ch.complete$height, 
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Complete", 
  ylab = "k (number of clusters)", 
  xlab = "h (node height)", 
  col = "grey" 
)
text(spe.ch.complete$height,
     nrow(spe):2, 
     nrow(spe):2,
     col = "red",
     cex = 0.8)

# Plot the fusion level values of the UPGMA clustering
plot(
  spe.ch.UPGMA$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.UPGMA$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)

# Plot the fusion level values of the Ward clustering
plot(
  spe.ch.ward$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Ward",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)

text(spe.ch.ward$height,
     nrow(spe):2,
     nrow(spe):2,
     col = "red",
     cex = 0.8)


# Plot the fusion level values of the beta-flexible clustering
# (beta = -0.25)
plot(
  spe.ch.beta2$height,
  nrow(spe):2,
  type = "S",
  main = "Fusion levels - Chord - Beta-flexible",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
text(spe.ch.beta2$height,
nrow(spe):2,
nrow(spe):2,
col = "red",
cex = 0.8)

dev.off()

# 4.7.3.2. Comparing two dendrograms ----
# Objects of class "hclust" must be first converted into objects of 
# class "dendrogram"
class(spe.ch.ward)     # [1] "hclust"
dend1 <- as.dendrogram(spe.ch.ward)
class(dend1) # [1] "dendrogram"
dend2 <- as.dendrogram(spe.ch.complete)
dend12 <- dendlist(dend1, dend2)

# Plots a tanglegram plot of a side by side trees.
tanglegram(
  untangle(dend12),
  sort = TRUE,
  common_subtrees_color_branches = TRUE,
  main_left = "Ward method",
  main_right = "Complete linkage"
)

 

# 4.7.3.3. Muliscale bootstrap resampling ----
# A hypothesis testing approach

# Compute p-values for all clusters (edges) of the dendrogram
spech.pv <-
  pvclust(t(spe.norm),
          method.hclust = "ward.D2",
          method.dist = "euc",
          parallel = TRUE)

# Plot dendrogram with p-values
plot(spech.pv)
# Highlight clusters with high AU p-values
pvrect(spech.pv, alpha = 0.95, pv = "au")
lines(spech.pv)
pvrect(spech.pv, alpha = 0.91, border = 4)


#  Silhousette-based ----
# 4.7.3.4. Average Silhouette Widths (ASW) ----
# *** high ASW indicates well-separated clusters

# Choose and rename the dendrogram ("hclust" object)
hc <- spe.ch.ward

# Plot average silhouette widths (using Ward clustering) for all 
# partitions except for the trivial partitions (k = 1 or k = n)
Si <- numeric(nrow(spe))
for (k in 2:(nrow(spe) - 1))
{
  sil <- silhouette(cutree(hc, k = k), spe.ch)
  Si[k] <- summary(sil)$avg.width
}
k.best <- which.max(Si)
plot(
  1:nrow(spe),
  Si,
  type = "h",
  main = "Silhouette-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Average silhouette width"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(Si),
       pch = 16,
       col = "red",
       cex = 1.5)
# k = 2

#  - Matrix correlation-based  ----
# 4.7.3.5. Comparison Between the Dissimilarity Matrix and Binary Matrices Representing Partitions ----
# Optimal number of clusters according to matrix correlation 
# statistic (Pearson)
kt <- data.frame(k = 1:nrow(spe), r = 0)
for (i in 2:(nrow(spe) - 1)) {
  gr <- cutree(hc, i)
  distgr <- grpdist(gr)
  mt <- cor(spe.ch, distgr, method = "pearson")
  kt[i, 2] <- mt
}
k.best <- which.max(kt$r)
plot(
  kt$k,
  kt$r,
  type = "h",
  main = "Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(kt$r),
       pch = 16,
       col = "red",
       cex = 1.5)

# k = 6

# - IndVal-based ----
# 4.7.3.6. Species Fidelity Analysis ----
# Optimal number of clusters as per indicator species analysis
# (IndVal, Dufrene-Legendre; package: labdsv)
IndVal <- numeric(nrow(spe))
ng <- numeric(nrow(spe))
for (k in 2:(nrow(spe) - 1)) {
  iva <- indval(spe, cutree(hc, k = k), numitr = 1000)
  gr <- factor(iva$maxcls[iva$pval <= 0.05])
  ng[k] <- length(levels(gr)) / k
  iv <- iva$indcls[iva$pval <= 0.05]
  IndVal[k] <- sum(iv)
}
k.best <- which.max(IndVal[ng == 1]) + 1
col3 <- rep(1, nrow(spe))
col3[ng == 1] <- 3

par(mfrow = c(1, 2))
plot(
  1:nrow(spe),
  IndVal,
  type = "h",
  main = "IndVal-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "IndVal sum",
  col = col3
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(
  which.max(IndVal),
  max(IndVal),
  pch = 16,
  col = "red",
  cex = 1.5
)
plot(
  1:nrow(spe),
  ng,
  type = "h",
  xlab = "k (number of groups)",
  ylab = "Ratio",
  main = "Proportion of clusters with significant indicator
species",
  col = col3 
)
axis(1, 
     k.best, 
     paste("optimum", k.best, sep = "\n"), 
     col = "red", 
     col.axis = "red") 
points(k.best, 
       max(ng), 
       pch = 16, 
       col = "red", 
       cex = 1.5) 

dev.off()

# k = 2; k = 2

# silhouette-based k = 2
# matrix correlation-based k = 6
# IndVal-based k = 2
# criteria do not return the same solution
# good compromise seems to be k = 4



# *** making k = 4
# 4.7.3.7. Silhouette Plot of the Final Partition ----
# Choose the number of clusters
k <- 4
# Silhouette Plot of the Final Partition
spech.ward.g <- cutree(spe.ch.ward, k = k)
sil <- silhouette(spech.ward.g, spe.ch)
rownames(sil) <- row.names(spe)
plot(
  sil,
  main = "Silhouette plot - Chord - Ward",  
  cex.names = 0.8,
  col = 2:(k + 1),
  nmax = 100
)

# 4.7.3.8. Final Dendrogram with Graphical Options ----
# Reorder clusters
spe.chwo <- reorder.hclust(spe.ch.ward, spe.ch)

# Plot reordered dendrogram with group labels
plot(
  spe.chwo,
  hang = -1,
  xlab = "4 groups",
  sub = "",
  ylab = "Height",
  main = "Chord - Ward (reordered)",
  labels = cutree(spe.chwo, k = k)
)
rect.hclust(spe.chwo, k = k)

# Plot the final dendrogram with group colors (RGBCMY...)
# Fast method using the additional hcoplot() function:
hcoplot(spe.ch.ward, spe.ch, lab = rownames(spe), k = 4)




# Convert the "hclust" object into a "dendrogram" object
dend <- as.dendrogram(spe.chwo)

# Plot the dendrogram with coloured branches
dend %>% set("branches_k_color", k = k) %>% plot

# Use standard colours for clusters
clusters <- cutree(dend, k)[order.dendrogram(dend)]
dend %>%
  set("branches_k_color", k = k, value = unique(clusters) + 1) %>%
  plot
# Add a coloured bar
colored_bars(clusters + 1,
             y_shift = -0.5,
             rowLabels = paste(k, "clusters"))

# 4.7.3.9. Spatial Plot of the Clustering Result ----
# Plot the Ward clusters on a map of the Doubs River 
# (see Chapter 2)
drawmap(xy = spa,
        clusters = spech.ward.g,
        main = "Four Ward clusters along the Doubs River")

# 4.7.3.10. Heat Map and Ordered Community Table ----
# Heat map of the dissimilarity matrix ordered with the dendrogram
heatmap(
  as.matrix(spe.ch),
  Rowv = dend,
  symm = TRUE,
  margin = c(3, 3)
)


# Ordered community table
# Species are ordered by their weighted averages on site scores.
# Dots represent absences.
or <- vegemite(spe, spe.chwo)


# Heat map of the doubly ordered community table, with dendrogram
heatmap(
t(spe[rev(or$species)]),
Rowv = NA,
Colv = dend,
col = c("white", brewer.pal(5, "Greens")),
scale = "none",
margin = c(4, 4),
ylab = "Species (weighted averages of sites)",
xlab = "Sites"

)

# Conclusions #################################################################
# - many clustering method
# - based on data and research objectives
# - compare different clustering methods
# - understand principle, pros and cons of each method
# - No single truth 没有“标准答案”

# END ----
