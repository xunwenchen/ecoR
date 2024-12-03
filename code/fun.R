# Build own functions here ----
# or borrow others' functions with acknowledgement ----

# calculate chi-squared statistic using R base
# Oi is observed vector, and Ei is expected proportional value
cal_x2 <- function(Oi, Ei_prop){
  # Oi is the observed value
  # Ei is expected proportional value
  t <- sum(Oi) # total number of observation
  r <- Ei_prop/sum(Ei_prop) # ratio/probability
  exp <- r*t # expected value
  chi_squared <- sum((Oi - exp)^2 / exp)# Calculate chi-squared statistic
  
  # Print the result
  return(chi_squared)
}


# panelutils.R 
#
# License: GPL-2 
# Author:  Francois Gillet
#          23 August 2012

## Put Pearson, Spearman or Kendall correlations on the upper panel
panel.cor <- function(x, y, method="pearson", digits=3, cex.cor=1.2, no.col=FALSE)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method=method)
  ra <- cor.test(x, y, method=method)$p.value
  txt <- round(r, digits)
  prefix <- ""
  if(ra <= 0.1) prefix <- "."
  if(ra <= 0.05) prefix <- "*"
  if(ra <= 0.01) prefix <- "**"
  if(ra <= 0.001) prefix <- "***"
  if(no.col)
  {
    color <- 1
    if(r < 0) { if(ra <= 0.001) sig <- 4 else sig <- 3 }
    else { if(ra <= 0.001) sig <- 2 else sig <- 1 }
  }
  else
  {
    sig <- 1
    if(ra <= 0.001) sig <- 2
    color <- 2
    if(r < 0) color <- 4
  }
  txt <- paste(txt, prefix, sep="\n")
  text(0.5, 0.5, txt, cex = cex.cor, font=sig, col=color)
}


## Put histograms on the diagonal
panel.hist <- function(x, no.col=FALSE, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  his <- hist(x, plot=FALSE)
  breaks <- his$breaks; nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  if(no.col) rect(breaks[-nB], 0, breaks[-1], y, col="gray", ...)
  else rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


## Add black lowess curves to scatter plots
panel.smoothb <- function (x, y, col=par("col"), bg=NA, pch=par("pch"), 
                           cex=1, col.smooth="black", span=2/3, iter=3, ...) 
{
  points(x, y, pch=pch, col=col, bg=bg, cex=cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f=span, iter=iter), col=col.smooth, ...)
}


#Usage:
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
#pairs(num.mat, lower.panel=panel.smooth, upper.panel=panel.cor, method="kendall")

# coldiss()
# Color plots of a dissimilarity matrix, without and with ordering
#
# License: GPL-2
# Author:  Francois Gillet
#          23 August 2012 - rev. 07 June 2016

"coldiss" <- function(D,
                      nc = 4,
                      byrank = TRUE,
                      diag = FALSE) {
  require(gclus)
  
  D <- as.dist(as.matrix(D))
  
  if (max(D) > 1)
    D <- D / max(D)
  
  if (byrank) {
    spe.color <- dmat.color(1 - D, cm.colors(nc))
  }
  else {
    spe.color <- dmat.color(1 - D, byrank = FALSE, cm.colors(nc))
  }
  
  spe.o <- order.single(1 - D)
  speo.color <- spe.color[spe.o, spe.o]
  
  op <- par(mfrow = c(1, 2), pty = "s")
  
  if (diag) {
    plotcolors(
      spe.color,
      rlabels = attributes(D)$Labels,
      main = "Dissimilarity Matrix",
      dlabels = attributes(D)$Labels
    )
    plotcolors(
      speo.color,
      rlabels = attributes(D)$Labels[spe.o],
      main = "Ordered Dissimilarity Matrix",
      dlabels = attributes(D)$Labels[spe.o]
    )
  }
  else {
    plotcolors(spe.color, rlabels = attributes(D)$Labels,
               main = "Dissimilarity Matrix")
    plotcolors(speo.color,
               rlabels = attributes(D)$Labels[spe.o],
               main = "Ordered Dissimilarity Matrix")
  }
  
  par(op)
}

# Usage:
# coldiss(D = dissimilarity.matrix, nc = 4, byrank = TRUE, diag = FALSE)
# If D is not a dissimilarity matrix (max(D) > 1), then D is divided by max(D)
# nc 							number of colours (classes)
# byrank = TRUE		equal-sized classes
# byrank = FALSE	equal-length intervals
# diag = TRUE			print object labels also on the diagonal

# Example:
# coldiss(spe.dj, nc = 9, byrank = FALSE, diag = TRUE)


# Map of the Doubs river (see Chapter 2)
# License: GPL-2
# Author:  Francois Gillet

drawmap <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", tcol = "white") {
    
    # Draw the Doubs river
    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "x coordinate (km)",
      ylab = "y coordinate (km)"
    )
    lines(xy, col = "light blue")
    text(65, 20, "Upstream", cex = 1.2)
    text(15, 32, "Downstream", cex = 1.2)
    
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
        #        pch = i + 20,
        pch = i + 15,
        cex = 3,
        col = i + 1,
        bg = i + 1
      )
    }
    text(xy,
         row.names(xy),
         cex = 0.8,
         col = tcol,
         font = 2)
    legend(
      "bottomright",
      paste("Cluster", 1:k),
      #      pch = (1:k) + 20,
      pch = (1:k) + 15,
      col = 2:(k + 1),
      pt.bg = 2:(k + 1),
      pt.cex = 2,
      bty = "n"
    )
    
  }


# Map of the Doubs river (see Chapter 4)
# License: GPL-2
# Author:  Francois Gillet

drawmap3 <-
  function(xy = spa, clusters, main = "Clusters along the Doubs river", 
           colors = palette()[-1], pch = 21, tcol = "black") {
    
    # Draw the Doubs river
    plot(
      xy,
      asp = 1,
      type = "n",
      main = main,
      xlab = "x coordinate (km)",
      ylab = "y coordinate (km)"
    )
    lines(xy, col = "light blue")
    text(65, 20, "Upstream", cex = 1.2)
    text(15, 32, "Downstream", cex = 1.2)
    
    # Add the clusters
    k <- length(levels(factor(clusters)))
    for (i in 1:k)
    {
      points(
        xy[clusters == i, 1],
        xy[clusters == i, 2],
        pch = pch,
        cex = 3,
        col = "white",
        bg = colors[i]
      )
    }
    text(xy,
         row.names(xy),
         cex = 0.8,
         col = tcol,
         font = 2)
    legend(
      "bottomright",
      paste("Cluster", 1:k),
      pch = 22,
      col = "white",
      pt.bg = colors,
      pt.cex = 2,
      bty = "n"
    )
    
  }

# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author:  Francois Gillet, 23 August 2012
# Revised: Daniel Borcard, 31 August 2017

"hcoplot" <- function(tree, 
                      diss, 
                      lab = NULL,
                      k, 
                      title = paste("Reordered dendrogram from", 
                                    deparse(tree$call), 
                                    sep="\n"))
{
  require(gclus)
  gr <- cutree(tree, k=k)
  tor <- reorder.hclust(tree, diss)
  plot(tor, 
       labels = lab,
       hang=-1, 
       xlab=paste(length(gr),"sites"), 
       sub=paste(k,"clusters"), 
       main=title)
  so <- gr[tor$order]
  gro <- numeric(k)
  for (i in 1 : k)
  {
    gro[i] <- so[1]
    if (i<k) so <- so[so!=gro[i]]
  }
  rect.hclust(tor, 
              k = k, 
              border = gro + 1, 
              cluster = gr)
  legend("topright", 
         paste("Cluster", 1 :k ), 
         pch = 22, 
         col = 2 : (k + 1), 
         bty = "n")
}


# Function to compute and test the statistic 'a' by permutation.
# 'a' is the co-occurrence of two species across the sites.
#
# Parameters:
# mat = site-by-species matrix. Sites are rows, species are columns.
# nperm = number of permutations. Choose this number so that the smallest
#         p-values will remain significant after correction for multiple
#         testing, e.g. Holm correction
#
# License: GPL-2
# Author:  Pierre Legendre
#          2010

test.a <-
  function (mat, nperm = 999) {
    A <- system.time({
      mat <- as.matrix(mat)
      site.names <- rownames(mat)
      sp.names <- colnames(mat)
      
      # Transform the 'pa' or abundance data to presence-absence
      mat <- decostand(mat, "pa")
      
      n <- nrow(mat)
      p <- ncol(mat)
      a <- t(mat) %*% mat
      
      # Permutation tests for a
      p.a = matrix(1, p, p, dimnames = list(sp.names, sp.names))
      for (iperm in 1:nperm) {
        perm.mat = mat[sample(n), 1]
        for (j in 2:p) {
          vec <- mat[sample(n), j]
          perm.mat <- cbind(perm.mat, vec)
        }
        #
        a.perm <- t(perm.mat) %*% perm.mat
        for (j in 2:p) {
          for (jj in 1:(p - 1)) {
            if (a.perm[j, jj] >= a[j, jj])
              p.a[j, jj] <- p.a[j, jj] + 1
          }
        }
      }
      p.a <- p.a / (nperm + 1)
      
      for (j in 1:(p - 1)) {
        for (jj in (j + 1):p) {
          p.a[j, jj] <- NA
        }
      }
      diag(p.a) <- NA
      
    })
    A[3] <- sprintf("%2f", A[3])
    cat("Computation time =", A[3], " sec", '\n')
    
    out <- list(a = a,
                p.a = p.a,
                p.a.dist = as.dist(p.a))
    out
  }


bartlett.perm <- function(y, fact, centr="MEDIAN", nperm=999, alpha=0.05)
  
  # Computation of parametric, permutational and bootstrap versions of the
  # Bartlett test of homogeneity of variances.
  #
  # The data are centred to their within-group medians (default) or means 
  # before the tests.
  #
  # Prior to the computation of the test of homogeneity of variances,
  # a Shapiro-Wilk test of normality of residuals is computed. If the residuals
  # are not normally distributed, a warning is issued because this nonnormality
  # influences the type I error of the parametric test, which will very likely 
  # have an incorrect type I error.
  #
  # USAGE
  # bartlett.perm(y, fact, centr, nperm, alpha)
  #
  # ARGUMENTS
  # y       a numeric vector of data values
  # fact    a vector or factor object giving the group for the corresponding
  #         elements of y
  # centr   should the data, within groups, be centred on their medians ("MEDIAN")
  #         or on their means ("MEAN")?
  # nperm   number of permutations
  # alpha   level of rejection of the H0 hypothesis of normality of residuals
  #         in the Shapiro-Wilk test
  #
  # RESULT
  # Bartlett          Bartlett's K-squared test statistic
  # Param.prob        Parametric probability (P-value) of the test
  # Permut.prob       Permutational probability (P-value) of the test
  # Bootstrap.prob    Bootstrap probability (P-value) of the test
  #
  # DETAILS
  #
  # Centring the groups on their median or mean is very important for permutation
  # and bootstrap tests to be correct when the groups do not share the same 
  # position. Permuting groups with unequal mean or median artificially increases 
  # the within-group variances of the permuted data.
  #
  # License: GPL-2
  # Author:  Daniel Borcard
  #          1 February 2016
  
  # ------------------------------------------------------
# EXAMPLE
#
# Species abundance type data:
# y1 <- log1p(round(rlnorm(5,0.2,2)))
# y2 <- log1p(round(rlnorm(5,1,2)))
# y3 <- log1p(round(rlnorm(5,2,5)))
# yy <- c(y1,y2,y3)
#
# Factor
# fac <- gl(3,5, labels=c("groupe1","groupe2","groupe3"))
#
# Bartlett test with centring on the group medians
# bartlett.perm(yy, fac, centr="MEDIAN", nperm=999, alpha=0.05)
# ------------------------------------------------------


{
  
  fact <- as.factor(fact)
  
  normal <- shapiro.test(resid(lm(y~fact)))
  if(normal[2]<=alpha){
    cat("\n-------------------------------------------------------------------")
    cat("\nWARNING") 
    cat("\nThe residuals of the ANOVA model are not normally distributed.")
    cat("\nThis is likely to change the rate of type I error of the test")
    cat("\nof homogeneity of variances.")
    cat("\n-------------------------------------------------------------------")
    cat("\n")
  }
  
  # Trap for groups with 0 dispersion
  y.sd <- tapply(y, fact, sd)
  if(any(y.sd==0)) {
    cat("\n-------------------------------------------------------------------")
    cat("\nPROBLEM ENCOUNTERED") 
    cat("\nOne or more groups have zero variance. Please chek and correct.")
    cat("\nThe computations are meaningless if a group is made of observations")
    cat("\nthat all have the same value.")
    cat("\n-------------------------------------------------------------------")
    cat("\n")
    stop
  }
  
  
  CENTRE <- c("MEDIAN", "MEAN")
  centr <- match.arg(centr, CENTRE)
  
  
  # Within-group centring of data
  
  if(centr == "MEDIAN"){
    meds <- tapply(y, fact, median, na.rm=TRUE)
    y.c <- y - meds[fact]
  }
  else{
    means <- tapply(y, fact, mean, na.rm=TRUE)
    y.c <- y - means[fact]
  }
  
  
  bart <- bartlett.test(y.c,fact)
  
  # Permutation tests
  
  cat("\nPermutations running...")
  cat("\n")
  
  compt.perm <- 1
  
  for(i in 1:nperm) {
    
    yprime <- sample(y.c)
    
    bart.perm <- bartlett.test(yprime,fact)
    if(bart.perm[[1]] >= bart[[1]]){
      compt.perm=compt.perm+1}
    
  }
  
  # Bootstrap tests
  # Difference with permutation test: resampling is done with replacement
  
  cat("\nBootstrap running...")
  cat("\n")
  cat("\n")
  
  compt.boot <- 1
  
  for(i in 1:nperm) {
    
    yboot <- sample(y.c, replace=TRUE)
    
    bart.boot <- bartlett.test(yboot,fact)
    if(bart.boot[[1]] >= bart[[1]]){
      compt.boot=compt.boot+1}
    
  }
  
  
  Result <- matrix(0,1,4)
  colnames(Result) <- c("Statistic", "Param.prob", "Permut.prob", "Bootstrap.prob")
  rownames(Result) <- "Bartlett" 
  
  Result[1,1] <- round(bart[[1]],4)
  Result[1,2] <- round(bart[[3]],4)
  Result[1,3] <- compt.perm/(nperm+1)
  Result[1,4] <- compt.boot/(nperm+1)
  
  Result
  
}


### Function to perform (unpaired) Kruskal-Wallis or paired Wilcoxon test with 
### post-hoc multiple comparisons and boxplots with letters 
### (using agricolae::kruskal)

# Arguments:
# Y = a numeric response variable
# X = a factor (groups, treatments...) or a qualitative variable that can be 
# converted to a factor
# p.adj = correction of p-values for multiple comparisons 
# (default=none, bonferroni, holm...)
# paired = if TRUE and if 2 groups only, then Wilcoxon paired test is performed
# (objects are supposed to be ordered twice the same way in the data frame)

# Value:
# A summary table ($comparison) is printed and boxplots are drawn with result
# of the test ($p.value) and, if it is significant, of post-hoc tests as letters
# (in decreasing order)

# # Example:
# library(agricolae)
# source("boxplerk.R")
# library(stats)
# data(InsectSprays)
# boxplerk(
#   Y = InsectSprays$count,
#   X = InsectSprays$spray,
#   ylab = "count",
#   xlab = "spray",
#   bcol = "bisque",
#   p.adj = "holm"
# )

# License: GPL-2
# Author:  Francois Gillet
#          2021-08-22 (adapted to agricolae 1.3-5)


boxplerk <-
  function(Y,
           X,
           main = NULL,
           xlab = "X factor",
           ylab = "Y value",
           bcol = "bisque",
           p.adj = "none",
           cexy = 1,
           varwidth = TRUE,
           las = 1,
           paired = FALSE) {
    
    if (!is.factor(X)) {
      X <- as.factor(X)
    }
    aa <- levels(X)
    
    tt1 <- matrix(nrow = length(aa), ncol = 7)
    for (i in 1:length(aa)) {
      temp <- Y[X == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- median(temp, na.rm = TRUE)
      tt1[i, 7] <- length(temp)
    }
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
    
    boxplot(
      Y ~ X,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth)
    
    require(agricolae)
    comp <- kruskal(Y, X, p.adj = p.adj)
    gror <- comp$groups[aa, ]
    tt1$rank <- gror$Y
    tt1$group <- gror$groups
    sig <- "ns"
    
    if (paired == TRUE & length(aa) == 2) {
      coms <- wilcox.test(Y ~ X, paired = TRUE)
      pp <- coms$p.value
      tt1$group <- rep("a", 2)
      if (pp <= 0.05) 
        tt1$group[which(tt1$rank == min(tt1$rank))] <- "b"
    }
    else {
      pp <- comp$statistics$p.chisq
    }
    
    if (pp <= 0.1)
      sig <- "."
    if (pp <= 0.05)
      sig <- "*"
    if (pp <= 0.01)
      sig <- "**"
    if (pp <= 0.001)
      sig <- "***"
    if (pp <= 0.0001)
      sig <- "****"
    mtext(
      sig,
      side = 3,
      line = 0.5,
      adj = 0,
      cex = 2,
      font = 1
    )
    
    if (pp <= 0.1)
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 0.5,
        cex = 1,
        font = 4
      )
    
    list(comparison = tt1, p.value = pp, p.adjust = p.adj)
  }


### Function to perform (unpaired) ANOVA or paired t-test with post-hoc
### multiple comparisons and boxplots with letters (using agricolae::LSD.test)

# Arguments:
# Y = a numeric response variable
# X = a factor (groups, treatments...) or a qualitative variable that can be 
# converted to a factor
# p.adj = correction of p-values for multiple comparisons
# (default=none, bonferroni, holm...)
# paired = if TRUE and if 2 groups only, then paired t-test is performed
# (objects are supposed to be ordered twice the same way in the data frame)

# Value:
# A summary table ($comparison) is printed and boxplots are drawn with result
# of the test ($p.value) and, if it is significant, of post-hoc tests as letters
# (in decreasing order)

# # Example:
# library(agricolae)
# data(sweetpotato)
# # Check ANOVA assumptions
# shapiro.test(resid(aov(sweetpotato$yield ~ sweetpotato$virus)))
# bartlett.test(sweetpotato$yield, sweetpotato$virus)
# source("boxplert.R")
# boxplert(
#   Y = sweetpotato$yield,
#   X = sweetpotato$virus,
#   ylab = "yield",
#   xlab = "virus",
#   bcol = "orange",
#   p.adj = "holm"
# )

# License: GPL-2
# Author:  Francois Gillet
#          2021-08-22 (adapted to agricolae 1.3-5)


boxplert <-
  function(Y,
           X,
           main = NULL,
           xlab = "X factor",
           ylab = "Y value",
           bcol = "bisque",
           p.adj = "none",
           cexy = 1,
           varwidth = TRUE,
           las = 1,
           paired = FALSE) {
    
    if (!is.factor(X)) {
      X <- as.factor(X)
    }
    aa <- levels(X)
    
    tt1 <- matrix(nrow = length(aa), ncol = 6)
    for (i in 1:length(aa)) {
      temp <- Y[X == aa[i]]
      tt1[i, 1] <- mean(temp, na.rm = TRUE)
      tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
      tt1[i, 3] <- sd(temp, na.rm = TRUE)
      tt1[i, 4] <- min(temp, na.rm = TRUE)
      tt1[i, 5] <- max(temp, na.rm = TRUE)
      tt1[i, 6] <- length(temp)
    }
    tt1 <- as.data.frame(tt1)
    row.names(tt1) <- aa
    colnames(tt1) <- c("mean", "se", "sd", "min", "max", "n")
    
    boxplot(
      Y ~ X,
      main = main,
      xlab = xlab,
      ylab = ylab,
      las = las,
      col = bcol,
      cex.axis = cexy,
      cex.lab = cexy,
      varwidth = varwidth
    )
    
    require(agricolae)
    sig <- "ns"
    model <- aov(Y ~ X)
    
    if (paired == TRUE & length(aa) == 2) {
      coms <- t.test(Y ~ X, paired = TRUE)
      pp <- coms$p.value
      tt1$group <- rep("a", 2)
      if (pp <= 0.05) 
        tt1$group[which(tt1$mean == min(tt1$mean))] <- "b"
    }
    else {
      pp <- anova(model)$Pr[1]
      comp <- LSD.test(model,
                       "X",
                       alpha = 0.05,
                       p.adj = p.adj,
                       group = TRUE)
      gror <- comp$groups[aa, ]
      tt1$group <- gror$groups
    }
    
    if (pp <= 0.1)
      sig <- "."
    if (pp <= 0.05)
      sig <- "*"
    if (pp <= 0.01)
      sig <- "**"
    if (pp <= 0.001)
      sig <- "***"
    if (pp <= 0.0001)
      sig <- "****"
    mtext(
      sig,
      side = 3,
      line = 0.5,
      adj = 0,
      cex = 2,
      font = 1
    )
    
    if (pp <= 0.1) {
      mtext(
        tt1$group,
        side = 3,
        at = c(1:length(aa)),
        line = 0.5,
        cex = 1,
        font = 4
      )
    }
    
    list(comparison = tt1, p.value = pp, p.adjust = p.adj)
  }


# Function to compute a binary dissimilarity matrix from clusters
grpdist <- function(X){
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}


# cleanplot.pca() function
#
# A function to draw a PCA biplot (scaling 1 or scaling 2) from an object
# of class "rda" containing PCA results, computed with vegan's rda() function.
# A circle of equilibrium contribution is drawn on scaling type 1 biplots.
#
# ARGUMENTS
#
# ##### General parameters
# res.pca          An rda{vegan} object. This object contains a PCA result if a
#                  single matrix was used in the rda(). The function will still
#                  operate with aRDA result but a Warning will be printed.
# ax1, ax2         Canonical axes to be drawn as abscissa and ordinate.
#                  Defaults: 1 and 2.
# scaling          Scaling type: only 1 or 2 are supported. Default: 1.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be
#                  plotted.
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# cex.char1        Character size (for sites and response variables).
#
# ##### Label positions
# ## Positions: 1 = below the point, 2 = left, 3 = above, 4 = right. Default: 4.
# ## Note - Argument pos = NULL centres the label on the position of the object
# ##       (site point, species or environmental variable arrow, centroid) when
# ##       the object is not drawn.
# pos.sites        Position of site labels. 1 to 4, as above. Default: 2.
# pos.spe          Position of species labels. 1 to 4, as above. Default: 4.
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be
#                  drawn in the biplot, e.g. c(1, 2, 5, 8, 12). Draw all species
#                  if select.spe = NULL (default value). The species that are
#                  well represented in the RDA plot can be identified using
#                  goodness(pca.output.object, display = "species").
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all objects and
#                  labels. Positive values increase the margins around the plot,
#                  negative values reduce them.
# optimum          If TRUE, the longest species arrow is stretched to
#                  a length equal to the distance to the origin of the site
#                  farthest from the origin of the plot of (ax1, ax2). This is
#                  an optimal combined representation of the sites and species.
#                  The lengths of the species arrows can be further modified
#                  using the arguments mult.spe.
# move.origin      Move plot origin right-left and up-down. 
#                  Default: move.origin = c(0,0).
#                  Ex. move.origin = c(-1, 0.5) moves origin by 1 unit left and
#                  0.5 unit up.
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed.
#                  Default: TRUE. Interpretation: examine the code lines
#                  controlled by argument "silent".
#
# Reference
# Legendre, P. & L. Legendre. 2012. Numerical ecology, 3rd English edition.
#    Elsevier Science BV, Amsterdam.
#
# # Example 1 - Doubs river fish data provided with the NEwR book
# # Data also available in the R package ade4
#
# # Load the Doubs.RData file from the NEwR (2018) book material.
# load("Doubs.RData")
#
# # The fish community data at 30 sites are in object "spe" (30 sites x 27
# # species). Site 8 may be removed before PCA, where no fish had been caught.
# # We did not do it here.
# ### Not done ###   spe.29 <- spe[-8, ]
# # Examine the position of site 8 in the biplots produced by the example code
# # lines below.
#
# # Chord-transform the fish data
# library(vegan)
# spe.norm <- decostand(spe, "normalize")
# pca.out <- rda(spe.norm)
#
# # Scaling 1
# source("cleanplot.pca.R")
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.out, scaling = 1, optimum = FALSE)
# cleanplot.pca(pca.out, scaling = 1, optimum = TRUE)
#
# # With optimum = TRUE (right-hand graph), the length of the longest species
# # arrow is equal to the distance between the origin and the site farthest
# # from the origin. Minor disadvantage: the radius of the equilibrium
# # contribution circle is not equal to sqrt(2/27) = 0.2721655 in that plot.
# # It has this value in the left-hand graph.
#
# # Argument silent = FALSE provides additional information in the R console, 
# # including the new circle radius.
# cleanplot.pca(pca.out,
#               scaling = 1,
#               optimum = FALSE,
#               silent = FALSE)
# cleanplot.pca(pca.out,
#               scaling = 1,
#               optimum = TRUE,
#               silent = FALSE)
#
# # Compare scaling 1 and scaling 2 biplots
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.out, scaling = 1, optimum = TRUE)
# cleanplot.pca(pca.out, scaling = 2, optimum = TRUE)
#
#
# # Example 2, part 1
# # Cajo ter Braak's dune meadow vegetation data (20 sites x 30 species)
#
# library(vegan)
# data(dune)
# dune.hel <- decostand(dune, "hellinger")
# pca.dune.hel <- rda(dune.hel)
#
# # Compare scaling 1 and scaling 2 biplots
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.dune.hel, scaling = 1, optimum = TRUE)
# cleanplot.pca(pca.dune.hel, scaling = 2, optimum = TRUE)
#
# # Example 2, part 2
#
# # The species in the previous plots are very crowded
# # A RDA of a community matrix by itself is a PCA of that matrix.
# # A demonstration is found in Legendre & Legendre (2012), p. 637.
# # The goodness function can only be called on RDA or CCA results
# rda.dune.hel <- rda(dune.hel ~ ., data = dune.hel)
# tmp <- goodness(rda.dune.hel)
# ( sp.sel <- which(tmp[, 2] >= 0.4) )   # 14 species selected for plotting
#
# # Scaling 2
# # Make the plots less crowded by only printing the best represented species
# # Right-hand graph: selected species only, and also example of moving the 
# # centroid of the points in the biplot.
# dev.new(width = 12, height = 7, noRStudioGD = TRUE)
# par(mfrow = c(1, 2))
# cleanplot.pca(pca.dune.hel, scaling = 2)  # All species
# cleanplot.pca(
#   pca.dune.hel,
#   scaling = 2,
#   select.spe = sp.sel,
#   move.origin = c(-0.3, 0)
# )
#
#
# License: GPL-2
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre
#          2016–2020

'cleanplot.pca' <-
  function(res.pca,
           ax1 = 1,
           ax2 = 2,
           scaling = 1,
           plot.sites = TRUE,
           plot.spe = TRUE,
           label.sites = TRUE,
           label.spe = TRUE,
           cex.char1 = 0.7,
           pos.sites = 2,
           pos.spe = 4,
           mult.spe = 1,
           select.spe = NULL,
           mar.percent = 0.1,
           optimum = TRUE,
           move.origin = c(0, 0),
           silent = TRUE) {
    
    ### Internal functions
    'stretch' <-
      function(sites, mat, ax1, ax2, n, silent = silent) {
        # Compute stretching factor for the species arrows
        # First, compute the longest distance to centroid for the sites
        tmp1 <- rbind(c(0, 0), sites[, c(ax1, ax2)])
        D <- dist(tmp1)
        target <- max(D[1:n])
        # Then, compute the longest distance to centroid for the species arrows
        if("matrix" %in% class(mat)) {
          #        if (is.matrix(mat)) {
          p <- nrow(mat)   # Number of species to be drawn
          tmp2 <- rbind(c(0, 0), mat[, c(ax1, ax2)])
          D <- dist(tmp2)
          longest <- max(D[1:p])
        } else {
          tmp2 <- rbind(c(0, 0), mat[c(ax1, ax2)])
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if (!silent)
          cat("target =",
              target,
              " longest =",
              longest,
              " fact =",
              target / longest,
              "\n")
        fact <- target / longest
      }
    
    'larger.plot' <-
      function(sit.sc,
               spe.sc,
               percent,
               move.origin,
               ax1,
               ax2) {
        # Internal function to expand plot limits (adapted from code by Pierre
        # Legendre)
        mat <- rbind(sit.sc, spe.sc)
        range.mat <- apply(mat, 2, range)
        rownames(range.mat) <- c("Min", "Max")
        z <- apply(range.mat, 2, function(x)
          x[2] - x[1])
        range.mat[1,] <- range.mat[1,] - z * percent
        range.mat[2,] <- range.mat[2,] + z * percent
        if (move.origin[1] != 0)
          range.mat[, ax1] <- range.mat[, ax1] - move.origin[1]
        if (move.origin[2] != 0)
          range.mat[, ax2] <- range.mat[, ax2] - move.origin[2]
        range.mat
      }
    
    "pcacircle" <-
      function (pca, mult.spe, fact.spe, silent = silent) {
        # This function draws a circle of equilibrium contribution on a PCA plot
        # generated from the result file of a vegan rda() analysis.
        eigenv <- pca$CA$eig
        p <- length(eigenv)
        n <- nrow(pca$CA$u)
        tot <- sum(eigenv)
        radius <- (2 / p) ^ 0.5 * mult.spe * fact.spe
        symbols(
          0,
          0,
          circles = radius,
          inches = FALSE,
          add = TRUE,
          fg = 2
        )
        if (!silent) {
          cat(
            "\nSpecies arrows and the radius of the equilibrium circle are stretched ",
            "by a factor of",
            mult.spe * fact.spe
          )
          cat(
            "\nThe radius of the equilibrium circle is thus",
            (2 / p) ^ 0.5,
            "*",
            mult.spe,
            "*",
            fact.spe,
            "=",
            radius,
            "\n"
          )
        }
      }
    ### End internal functions
    
    if (!class(res.pca)[1] == "rda")
      stop("The input file is not a vegan output object of class 'rda'",
           call. = FALSE)
    if (!(is.null(res.pca$CCA)))
      stop(
        "The input file contains an RDA, not a PCA result. ",
        "Use function triplot.rda from the NEwR (2018) book to produce an RDA triplot."
      )
    if (scaling != 1 &
        scaling != 2)
      stop("Function only available for scaling 1 or 2", call. = FALSE)
    
    k <- length(res.pca$CA$eig)         # n. of PCA eigenvalues
    n.sp <- length(res.pca$colsum)      # n. of species
    ahead <- 0.05   # Length of arrow heads
    aangle <- 30    # Angle of arrow heads
    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1
    
    # This version reconstructs and uses the original RDA output of L&L 2012,
    # Section 11.1.3
    
    Tot.var = res.pca$tot.chi         # Total variance in response data Y
    eig.val = res.pca$CA$eig          # Eigenvalues of Y-hat
    Lambda = diag(eig.val)            # Diagonal matrix of eigenvalues
    eig.val.rel = eig.val / Tot.var   # Relative eigenvalues of Y-hat
    Diag = diag(sqrt(eig.val.rel))    # Diagonal matrix of sqrt(relative
    # eigenvalues)
    U.sc1 = res.pca$CA$v              # Species scores, scaling=1
    U.sc2 = U.sc1 %*% Lambda ^ (0.5)  # Species scores, scaling=2
    n = nrow(res.pca$CA$u)            # Number of observations
    Z.sc2 = res.pca$CA$u * sqrt(n - 1)# Site scores, scaling=2
    Z.sc1 = Z.sc2 %*% Lambda ^ (0.5)  # Site scores, scaling=1
    
    if (is.null(select.spe)) {
      vec <- 1:n.sp
    } else {
      vec <- select.spe
    }
    
    if (scaling == 1) {
      sit.sc <- Z.sc1
      spe.sc <- U.sc1[vec,]
    } else {
      # For scaling=2
      sit.sc <- Z.sc2
      spe.sc <- U.sc2[vec,]
    }
    if (is.null(rownames(sit.sc)))
      rownames(sit.sc) <- paste("Site", 1:n, sep = "")
    if (is.null(rownames(spe.sc)))
      rownames(spe.sc) <- paste("Sp", 1:n.sp, sep = "")
    
    fact.spe <- 1
    if (optimum) {
      fact.spe <-
        stretch(sit.sc[, 1:k], spe.sc[, 1:k], ax1, ax2, n, silent = silent)
    }
    if (!silent)
      cat("fact.spe =", fact.spe, "\n\n")
    spe.sc <- spe.sc * fact.spe * mult.spe
    
    lim <-
      larger.plot(
        sit.sc[, 1:k],
        spe.sc[, 1:k],
        percent = mar.percent,
        move.origin = move.origin,
        ax1 = ax1,
        ax2 = ax2
      )
    if (!silent)
      print(lim)
    
    # Draw the main plot
    mat <- rbind(sit.sc[, 1:k], spe.sc[, 1:k])
    plot(
      mat[, c(ax1, ax2)],
      type = "n",
      main = paste("PCA biplot - Scaling", scaling),
      xlim = c(lim[1, ax1], lim[2, ax1]),
      ylim = c(lim[1, ax2], lim[2, ax2]),
      xlab = paste("PCA ", ax1),
      ylab = paste("PCA ", ax2),
      asp = 1
    )
    abline(h = 0, v = 0, col = "grey60")
    
    # Draw the site scores
    if (plot.sites) {
      points(sit.sc[, ax1], sit.sc[, ax2], pch = 20)
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = pos.sites,
          cex = cex.char1
        )
    } else {
      if (label.sites)
        text(
          sit.sc[, ax1],
          sit.sc[, ax2],
          labels = rownames(sit.sc),
          col = "black",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # Draw the species scores
    if (plot.spe) {
      arrows(
        0,
        0,
        spe.sc[, ax1],
        spe.sc[, ax2],
        length = ahead,
        angle = aangle,
        col = "red"
      )
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = pos.spe,
          cex = cex.char1
        )
    } else {
      if (label.spe)
        text(
          spe.sc[, ax1],
          spe.sc[, ax2],
          labels = rownames(spe.sc),
          col = "red",
          pos = NULL,
          cex = cex.char1
        )
    }
    
    # If scaling = 1 draw circle of equilibrium contribution
    if (scaling == 1) {
      pcacircle(
        res.pca,
        mult.spe = mult.spe,
        fact.spe = fact.spe,
        silent = silent
      )
    }
  }


"PCA.newr" <- function(Y, stand = FALSE)
  
  # Principal component analysis (PCA) with option for variable standardization
  
  # stand = FALSE : center by columns only, do not divide by s.d.
  # stand = TRUE  : center and standardize (divide by s.d.) by columns
  #
  # License: GPL-2
  # Author:  Pierre Legendre
  #          May 2006
{
  Y <- as.matrix(Y)
  obj.names <- rownames(Y)
  var.names <- colnames(Y)
  size <- dim(Y)
  Y.cent <- apply(Y, 2, scale, center = TRUE, scale = stand)
  Y.cov <- cov(Y.cent)
  Y.eig <- eigen(Y.cov)
  k <- length(which(Y.eig$values > 1e-10))
  U <- Y.eig$vectors[, 1:k]
  F <- Y.cent %*% U
  U2 <- U %*% diag(Y.eig$value[1:k] ^ (0.5))
  G <- F %*% diag(Y.eig$value[1:k] ^ (-0.5))
  rownames(F) <- obj.names
  rownames(U) <- var.names
  rownames(G) <- obj.names
  rownames(U2) <- var.names
  axenames <- paste("Axis", 1:k, sep = " ")
  colnames(F) <- axenames
  colnames(U) <- axenames
  colnames(G) <- axenames
  colnames(U2) <- axenames
  
  # Fractions of variance
  varY <- sum(diag(Y.cov))
  eigval <- Y.eig$values[1:k]
  relative <- eigval / varY
  rel.cum <- vector(length = k)
  rel.cum[1] <- relative[1]
  for (kk in 2:k) {
    rel.cum[kk] <- rel.cum[kk - 1] + relative[kk]
  }
  
  out <-
    list(
      total.var = varY,
      eigenvalues = eigval,
      rel.eigen = relative,
      rel.cum.eigen = rel.cum,
      U = U,
      F = F,
      U2 = U2,
      G = G,
      stand = stand,
      obj.names = obj.names,
      var.names = var.names,
      call = match.call()
    )
  class(out) <- "PCA.newr"
  out
}



"print.PCA.newr" <- function(x, ...) {
  cat("\nPrincipal Component Analysis\n")
  cat("\nCall:\n")
  cat(deparse(x$call), '\n')
  if (x$stand)
    cat("\nThe data have been centred and standardized
        by column", '\n')
  cat("\nTotal variance in matrix Y: ", x$total.var, '\n')
  cat("\nEigenvalues", '\n')
  cat(x$eigenvalues, '\n')
  cat("\nRelative eigenvalues", '\n')
  cat(x$rel.eigen, '\n')
  cat("\nCumulative relative eigenvalues", '\n')
  cat(x$rel.cum.eigen, '\n')
  invisible(x)
}



"biplot.PCA.newr" <-
  function(x,
           scaling = 1,
           plot.axes = c(1, 2),
           color.obj = "black",
           color.var = "red",
           ...)
    # scaling = 1 : preserves Euclidean distances among the objects
    # scaling = 2 : preserves correlations among the variables
  {
    #### Internal function
    larger.frame <- function(mat, percent = 0.07)
      # Produce an object plot 10% larger than strictly necessary
    {
      range.mat <- apply(mat, 2, range)
      z <- apply(range.mat, 2, function(x)
        x[2] - x[1])
      range.mat[1, ] <- range.mat[1, ] - z * percent
      range.mat[2, ] <- range.mat[2, ] + z * percent
      range.mat
    }
    ####
    
    if (length(x$eigenvalues) < 2)
      stop("There is a single eigenvalue.
           No plot can be produced.")
    if (length(which(scaling == c(1, 2))) == 0)
      stop("Scaling must be 1 or 2")
    
    par(mai = c(1.0, 0.75, 1.0, 0.5))
    
    if (scaling == 1)
    {
      # Distance biplot, scaling type = 1: plot F for objects, U for variables
      # This projection preserves the Euclidean distances among the objects
      lf.F <- larger.frame(x$F[, plot.axes])
      biplot(
        x$F[, plot.axes],
        x$U[, plot.axes],
        col = c(color.obj, color.var),
        xlim = lf.F[, 1],
        ylim = lf.F[, 2],
        arrow.len = 0.05,
        asp = 1
      )
      title(main = c("PCA biplot", "scaling type 1"),
            line = 3)
    }
    else
    {
      # Correlation biplot, scaling type = 2: plot G for objects, U2 for variables
      # This projection preserves the correlation among the variables
      lf.G <- larger.frame(x$G[, plot.axes])
      biplot(
        x$G[, plot.axes],
        x$U2[, plot.axes],
        col = c(color.obj, color.var),
        xlim = lf.G[, 1],
        ylim = lf.G[, 2],
        arrow.len = 0.05,
        asp = 1
      )
      title(main = c("PCA biplot", "scaling type 2"),
            line = 3)
    }
    abline(h = 0, lty = 3)
    abline(v = 0, lty = 3)
    invisible()
  }


`CA.newr` <- 
  function(Y, use.svd=TRUE, cumfit.obj=TRUE, cumfit.var=TRUE)
    #
    # Compute correspondence analysis (CA).
    # Data table Y must contain frequencies or equivalent.
    #
    # use.svd: decomposition is done by svd (default). It can also be done by eigen.
    #          The signs of axes may differ between the two methods.
    # color  : color of the species symbols and labels in the biplots.
    #
    # For exercise, you may choose 'Table_9.11.txt' (small example from Chapter 9
    # of the manual) or 'Spiders_28x12_spe.txt' (larger data set, real data).
    #
    # License: GPL-2
    # Author:  Pierre Legendre
    #          January 2008
  {
    # Begin internal functions
    sq.length <- function(vec) sum(vec^2)
    #
    'cumul.fit.va' <- function(Fhat,n,p,k,var.names)
      # Compute the table of "Cumulative fit per variable" 
    {
      sp.var <- diag(var(Y))
      res <- matrix(NA,p,k)
      for(i in 1:p) {
        res[i,] <- cumsum(Fhat[i,]^2)/sq.length(Fhat[i,])
      }
      rownames(res) <- var.names
      colnames(res) <- paste("Cum.axis",1:k,sep=".")
      res
    }
    #
    'cumul.fit.ob' <- function(F,n,p,k,obj.names)
      # Compute the table of "Cumulative fit of the objects" 
    {
      res <- matrix(NA,n,k)
      for(i in 1:n) {
        res[i,] <- cumsum(F[i,]^2)/sq.length(F[i,])
      }
      rownames(res) <- obj.names
      colnames(res) <- paste("Cum.axis",1:k,sep=".")
      res
    }
    # End internal functions
    #
    Y = as.matrix(Y)
    if(min(Y) < 0) stop("Negative values not allowed in CA")
    #
    # Calculate basic parameters of Y
    n = nrow(Y)
    p = ncol(Y)
    #
    # Save the row and column names
    site.names = rownames(Y)
    sp.names = colnames(Y)
    #
    # Construct the Qbar matrix (contributions to chi-square)
    # Numerical ecology (1998), equations 9.31 and 9.32
    fi. = matrix(apply(Y,1,sum),n,1)
    f.j = matrix(apply(Y,2,sum),1,p)
    f.  = sum(fi.)
    pi. = as.vector(fi./f.)
    p.j = as.vector(f.j/f.)
    E = (fi. %*% f.j)/f.
    Qbar = (Y - E) * E^(-0.5) / sqrt(f.)
    inertia = sum(Qbar^2)
    #
    if(use.svd) {
      # Analyse Qbar by 'svd'
      svd.res = svd(Qbar)
      k = length(which(svd.res$d > 1e-8))
      values = svd.res$d[1:k]^2
      U = svd.res$v[,1:k]
      Uhat = svd.res$u[,1:k]
    } else {
      # Alternative analysis or Qbar by 'eigen'
      Qbar = as.matrix(Qbar)
      QprQ.eig = eigen( t(Qbar) %*% Qbar )
      k = length(which(QprQ.eig$values > 1e-16))
      values = QprQ.eig$values[1:k]
      U = QprQ.eig$vectors[,1:k]
      Uhat = Qbar %*% U %*% diag(values^(-0.5))
    }
    #
    rel.values = values/inertia
    cum.rel <- cumsum(rel.values)
    #
    # Construct matrices V, Vhat, F, and Fhat for biplots, scalings 1 and 2
    V = diag(p.j^(-0.5)) %*% U
    Vhat = diag(pi.^(-0.5)) %*% Uhat
    F = Vhat %*% diag(values^(0.5))
    Fhat = V %*% diag(values^(0.5))
    #
    # Matrices for biplot, scaling = 3 (Symmetric scaling in Canoco)
    spec3 = V %*% diag(values^(0.25))                # Species scores
    site3 = Vhat %*% diag(values^(0.25))             # Site scores
    #
    if(cumfit.var) {
      cfit.spe <- cumul.fit.va(Fhat,n,p,k,sp.names)
    } else {
      cfit.spe <- NULL
    }
    #
    if(cumfit.obj) {
      cfit.obj <- cumul.fit.ob(F,n,p,k,site.names)
    } else {
      cfit.obj <- NULL
    }
    #
    rownames(U) <- rownames(V) <- rownames(spec3) <- rownames(Fhat) <- sp.names
    rownames(Uhat) <- rownames(F) <- rownames(Vhat) <- rownames(site3) <- site.names
    ax.names <- paste("Axis",1:k,sep="")
    colnames(U) <- colnames(Uhat) <- colnames(V) <- colnames(spec3) <- colnames(Vhat) <- colnames(site3) <- colnames(F) <- colnames(Fhat) <- ax.names
    #
    general <- list(inertia=inertia, values=values, rel.values=rel.values, cum.rel=cum.rel)
    scaling1 <- list(species=V, sites=F)
    scaling2 <- list(species=Fhat, sites=Vhat)
    scaling3 <- list(species=spec3, sites=site3)
    scaling4 <- list(species=Fhat, sites=F)
    fit <- list(cumulfit.spe=cfit.spe, cumulfit.obj=cfit.obj)
    other <- list(U=U, Uhat=Uhat, F=F, Fhat=Fhat, site.names=site.names, sp.names=sp.names, Qbar=Qbar, call=match.call() )
    #
    out <- list(general=general, scaling1=scaling1, scaling2=scaling2, scaling3=scaling3, scaling4=scaling4, fit=fit, other=other)
    class(out) <- "CA.newr"
    out
  }

`print.CA.newr` <-
  function(x, kk=5, ...)
  {
    if (!inherits(x, "CA.newr")) stop("Object of class 'CA.newr' expected")
    cat("\nCorrespondence Analysis\n")
    cat("\nCall:\n")
    cat(deparse(x$other$call),'\n')
    cat("\nTotal inertia in matrix Qbar: ",x$general$inertia,'\n')
    cat("\nEigenvalues",'\n')
    cat(x$general$values,'\n')
    cat("\nRelative eigenvalues",'\n')
    cat(x$general$rel.values,'\n')
    cat("\nCumulative relative eigenvalues",'\n')
    cat(x$general$cum.rel,'\n')
    kk <- min(length(x$general$values), kk)
    if(!is.null(x$fit$cumulfit.spe)) {
      cat("\nCumulative fit per species (",kk,"axes)",'\n')
      print.default(x$fit$cumulfit.spe[,1:kk], digits=5)
    }
    if(!is.null(x$fit$cumulfit.obj)) {
      cat("\nCumulative fit of the objects (",kk,"axes)",'\n')
      print.default(x$fit$cumulfit.obj[,1:kk], digits=5)
    }
    cat('\n')
    invisible(x) 
  }

`biplot.CA.newr` <-
  function(x, xax=1, yax=2, scaling=1, aspect=1, cex=1, color.sites="black", color.sp="red",...)
    # xax and yax determine the axes that will be plotted.
    # Use aspect=NA to remove the effect of parameter 'asp' in the biplot.
  {
    if (!inherits(x, "CA.newr")) stop("Object of class 'CA.newr' expected")
    if(length(x$general$values) < 2) stop("There is a single eigenvalue. No plot can be produced.")
    #
    sp.names = x$other$sp.names
    si.names = x$other$site.names
    #
    
    if(scaling == 1) {
      
      # The sites are at the centroids (barycentres) of the species
      # This projection preserves the chi-square distance among the sites
      type = "scaling type 1"
      sp = x$scaling1$species
      si = x$scaling1$sites
      
    } else if(scaling == 2) {
      
      # The species are at the centroids (barycentres) of the sites
      # This projection preserves the chi-square distance among the species
      type = "scaling type 2"
      sp = x$scaling2$species
      si = x$scaling2$sites
      
    } else if(scaling == 3) {
      
      # Biplot, scaling = 3 (Symmetric scaling in Canoco)
      type = "scaling type 3"
      sp = x$scaling3$species
      si = x$scaling3$sites
      
    } else if(scaling == 4) {
      
      # For contingency tables --
      # Preserves the chi-square distance among the rows and among the columns
      type = "scaling type 4"
      sp = x$scaling2$species
      si = x$scaling1$sites
      
    } else { 
      
      stop("Program stopped: error in scaling type")
    }
    
    # Find the limits of the axes
    sp.range = apply(sp[,c(xax,yax)],2,range)
    si.range = apply(si[,c(xax,yax)],2,range)
    
    # Biplot: plot 'si' for sites, 'sp' for species
    ran.si = si.range[2,] - si.range[1,]
    ran.sp = sp.range[2,] - sp.range[1,]
    ran.x = max(ran.si[1], ran.sp[1])
    xmin = min(sp.range[1,1], si.range[1,1]) - ran.x/8
    xmax = max(sp.range[2,1], si.range[2,1]) + ran.x/3
    ymin = min(sp.range[1,2], si.range[1,2])
    ymax = max(sp.range[2,2], si.range[2,2])
    #
    plot(si[,c(xax,yax)], asp=aspect, pch=20, cex=cex, xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=paste("CA axis",xax), ylab=paste("CA axis",yax), col=color.sites)
    text(si[,c(xax,yax)], labels=si.names, cex=cex, pos=4, offset=0.5, col=color.sites)
    points(sp[,c(xax,yax)], pch=22, cex=cex, col=color.sp)
    text(sp[,c(xax,yax)], labels=sp.names, cex=cex, pos=4, offset=0.5, col=color.sp)
    title(main = c("CA biplot",type), family="serif")
    #
    invisible()
  }

# triplot.rda() function
#
# A function to draw a triplot (scaling 1 or scaling 2) from an object 
# of class "rda" containing RDA result from vegan's rda() function.
# 
# This version can handle RDA results with a single or more canonical axes. It can also
# draw any combination of canonical and residual axes, or pairs of canonical axes or pairs
# of residual axes.
# This function requires package vegan.
#
# ARGUMENTS
#
# ##### General parameters
# res.rda          An rda{vegan} object.
# ax1, ax2         Canonical axes to be drawn as abscissa and ordinate. Defaults: 1 and 2.
# site.sc          Can be set to "lc" (linear constraints or model scores, default) 
#                  or "wa" (weighted averages, default in vegan).
# scaling          Scaling type: only 1 or 2 are supported. Default: 1.
#
# ##### Items to be plotted
# plot.sites       If TRUE, the sites will be plotted as small circles.
# plot.spe         If TRUE, the species (or other response variables) will be plotted.
# plot.env         If TRUE, arrows for the explanatory variables will be plotted.
# plot.centr       If TRUE, symbols are plotted at the centroids of factor levels, if any. 
# arrows.only      if TRUE, plot arrows for quant. explanatory var. and factor classes
# label.sites      If TRUE, labels are added to the site symbols.
# label.spe        If TRUE, labels are added to the species arrows.
# label.env        If TRUE, labels are added to the environmental variable arrows.
# label.centr      If TRUE, labels are added to the centroids of factor levels.
# cex.char1        Character size (for sites and response variables).
# cex.char2        Character size (for explanatory variables).
# cex.point        Size of points representing centroids of factor levels.
# pca.ax           If there is a single canonical axis, use a PCA axis of the residuals
#                  as the ordinate of the triplot. Default: pca.ax=1
#
# ##### Label positions, colour centroids
# ## Positions: 1=below the point, 2=left, 3=above, 4=right. Default: 4.
# ## Note - Argument pos=NULL centres the label on the position of the object (site point,  
# ## species or environmental variable arrow, centroid) when the object is not drawn.
# pos.sites        Position of site labels. 1 to 4, as above. Default: 2.
# pos.spe          Position of species labels. 1 to 4, as above. Default: 4.
# pos.env          Position of env.variable labels. 1 to 4, as above. Default: 4.
# pos.centr        Position of centroid labels. 1 to 4, as above. Default: 4. 
# col.centr        Colour of centroid labels. Default: "blue". 
#                  Change to another colour as needed, for ex. "turquoise4"
#
# ##### Multipliers, selection of species to be plotted
# mult.spe         Multiplier for length of the species arrows. Default: 1.
# mult.arrow       Multiplier for length of the environmental arrows. Default: 1.
# select.spe       Vector containing a selection of the species numbers to be drawn in 
#                  the biplot, e.g. c(1,2,5,8,12). Draw all species if select.spe=NULL 
#                  (default value). The species that are well represented in the RDA plot 
#                  can be identified using goodness(RDA.output.object,display="species")
#
# ##### Position of the plot in frame, margins
# mar.percent      Factor to expand plot size to accomodate all items and labels. Positive 
#                  values, 0.g. 0.1, increase the margins around the plot, negative  
#                  values, e.g. -0.1, reduce them.
# optimum          If TRUE, the longest species and environmental arrows are stretched to 
#                  a length equal to the distance to the origin of the site farthest from 
#                  the origin in the plot of (ax1,ax2). This is an optimal combined 
#                  representation of the three elements. The lengths of the species and 
#                  environmental arrows can be further modified using the arguments 
#                  mult.spe and mult.arrow.
# move.origin      Move plot origin right-left and up-down. Default: move.origin=c(0,0).
#                  Ex. move.origin=c(-1,0.5) moves origin by 1 unit left and 0.5 unit up.
#
# ##### Varia
# silent           If FALSE, intermediate computation steps are printed. Default: TRUE.
#
# 
# # Example 1 - Table 11.3 of Legendre & Legendre (2012, p. 644), first 6 species only
#
# Y.mat = matrix(c(1,0,0,11,11,9,9,7,7,5,0,0,1,4,5,6,7,8,9,10,0,0,0,0,17,0,13,0,10,0,0, 
# 0,0,0,7,0,10,0,13,0,0,0,0,8,0,6,0,4,0,2,0,0,0,1,0,2,0,3,0,4),10,6)
# rownames(Y.mat) = paste("site",1:10,sep="")
# colnames(Y.mat) = paste("Sp",1:6,sep="")
# Depth = 1:10
# Sub. = as.factor(c(rep(1,3),3,2,3,2,3,2,3))     # 1=Sand, 2=Coral, 3=Other
# env = cbind(data.frame(Depth),data.frame(Sub.))
# 
# rda.out = rda(Y.mat~ Depth+Sub.,env)
# 
# # Scaling=1
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=1, mar.percent=0)
# triplot.rda(rda.out, scaling=1, move.origin=c(5,-5), mar.percent=-0.1)
#
# # Scaling=2
# par(mfrow=c(1,2))
# triplot.rda(rda.out, scaling=2, mar.percent=0.15, silent=FALSE)
# triplot.rda(rda.out, scaling=2, move.origin=c(0.4,-0.25), mar.percent=0.05,silent=FALSE)
# With silent=FALSE, additional information is provided on the axes positions and ranges
#
# # Example 2 - Same data
# 
# # A single quantitative explanatory variable; a single canonical axis is produced
# ( rda.out.1 = rda(Y.mat~ Depth, data=env) )
# 
# # Plot lc-scores. Compare vegan's plot.cca() and triplot.rda()
# par(mfrow = c(1, 2))
# plot(rda.out.1, scaling=1, display=c("sp","lc","cn"), col = "green4", cex = .7) 
# triplot.rda(rda.out.1, scaling=1)
# 
# # Any of the residual PCA axes can be used as the ordinate. Recommended: use ax2 = 2
#   i.e., the first non-canonical axis.
# par(mfrow = c(1, 2))
# triplot.rda(rda.out.1, scaling=1, ax1 = 1, ax2 = 2) # Ordinate: first residual axis
# triplot.rda(rda.out.1, scaling=2, ax1 = 1, ax2 = 3) # Ordinate: second residual axis
# 
# # Example 3 - Dune data from vegan
# 
# library(vegan)
# data(dune)
# data(dune.env)
# 
# rda.dune = rda(dune ~ .,dune.env)
# 
# Select the species well represented on RDA axes 1 and 2 of the triplot
# tmp = goodness(rda.dune)
# ( sp.sel = which(tmp[,2] >= 0.4) )
#
# # Scaling=2. Left: plot witl all species. Right: with selected species
# par(mfrow=c(1,2))
# triplot.rda(rda.dune, scaling=2, mar.percent=0)
# triplot.rda(rda.dune, scaling=2, 
#       select.spe=sp.sel, move.origin=c(-0.3,0), mar.percent=0.1)
#
# Use all species but plot the two first residual (non-canonical) axes, i.e., axes 13 
# and 14. Left plot with scaling = 1, right plot with scaling 2.
# dev.new(noRStudioGD = TRUE, width = 12, height = 6)
# par(mfrow = c(1,2))
# triplot.rda(rda.dune, ax1 = 13, ax2 = 14, scaling = 1)
# triplot.rda(rda.dune, ax1 = 13, ax2 = 14, scaling = 2)
#
######################################################################################
# Users interested in similar plots but in ggplot-style may try the sister function 
# ggtriplotRDA.R. This function requires the additional packages tidyverse and ggrepel
# with their dependencies.
######################################################################################
#
# #####
#
# License: GPL-2 
# Authors: Francois Gillet, Daniel Borcard & Pierre Legendre, 2016, 2021
#          This version is dev5.0.1 (4 June 2021)

'triplot.rda' <-
  function(res.rda, 
           ax1=1,
           ax2=2,
           site.sc="lc",
           scaling=1,
           plot.sites=TRUE, 
           plot.spe=TRUE,
           plot.env=TRUE,
           plot.centr=TRUE, 
           arrows.only=FALSE,
           label.sites=TRUE,
           label.spe=TRUE, 
           label.env=TRUE, 
           label.centr=TRUE,
           cex.char1=0.7, 
           cex.char2=1.0,
           cex.point=0.8,
           pos.sites=2,
           pos.spe=4,
           pos.env=4,  
           pos.centr=4,
           col.centr="blue",
           mult.spe=1, 
           mult.arrow=1,
           select.spe=NULL,
           mar.percent=0.15,  
           optimum=TRUE,
           move.origin=c(0,0),
           silent=TRUE) {
    
    ### Internal functions
    #
    'stretch' <- 
      function(sites, mat, ax1, ax2, n, silent=silent) {
        # Compute stretching factor for the species or environmental arrows
        # First, compute the longest distance to centroid for the sites
        tmp1 <- rbind(c(0,0), sites[,c(ax1,ax2)])
        D <- dist(tmp1)
        target <- max(D[1:n])
        # Then, compute the longest distance to centroid for the species or environmental arrows
        if (is.matrix(mat)) {
          p <- nrow(mat)   # Number of species or env. arrows to be drawn
          tmp2 <- rbind(c(0,0), mat[,c(ax1,ax2)])
          D <- dist(tmp2)
          longest <- max(D[1:p])
        } else {
          tmp2 <- rbind(c(0,0), mat[c(ax1,ax2)]) 
          longest <- dist(tmp2)
          # print(tmp2)
        }  # If a single row left in 'mat'
        #
        if(!silent)
          cat("target =", target,
              " longest =", longest,
              " fact =",target/longest,
              "\n")
        fact <- target/longest
      }
    #
    'larger.plot' <- 
      function(sit.sc, 
               spe.sc,
               BP.sc,
               percent,
               move.origin,
               ax1,
               ax2,
               k) {
        # Internal function to expand plot limits (adapted from code by Pierre Legendre)
        if(ax1 > k) {
          mat <- rbind(sit.sc, spe.sc)
        }
        else {
          mat <- rbind(sit.sc, spe.sc, BP.sc)
        }
        range.mat <- apply(mat, 2, range)
        rownames(range.mat) <- c("Min","Max")
        z <- apply(range.mat, 2, function(x) x[2]-x[1])
        range.mat[1,] <- range.mat[1,]-z*percent
        range.mat[2,] <- range.mat[2,]+z*percent
        if(move.origin[1] != 0) range.mat[,ax1] <- range.mat[,ax1] - move.origin[1]
        if(move.origin[2] != 0) range.mat[,ax2] <- range.mat[,ax2] - move.origin[2]
        range.mat
      }
    ### End of internal functions
    
    if(class(res.rda)[1]!="rda" & class(res.rda)[2]!="rda")
      stop("The input file is not a vegan rda output object")
    if(length(res.rda$colsum)==1)
      stop("Function triplot.rda is not compatible with results that contain no species scores")
    if(scaling!=1 & scaling!=2) 
      stop("Function only available for scaling = 1 or 2")
    if(ax1 > ax2)
      stop("Axes must be drawn in increasing order")
    
    k <- length(res.rda$CCA$eig)        # number of RDA eigenvalues
    if(ax1 <= k & site.sc=="lc") {
      cat(  "-----------------------------------------------------------------------")
      cat("\nSite constraints (lc) selected. To obtain site scores that are weighted") 
      cat("\nsums of species scores (default in vegan), argument site.sc must be set")
      cat("\nto wa.")
      cat("\n-----------------------------------------------------------------------\n\n")
    }
    n.sp <- length(res.rda$colsum)      # number of species
    # 'vec' will contain the selection of species to be drawn
    if (is.null(select.spe)) {
      vec <- 1:n.sp } else {
        vec <- select.spe
      }
    
    ahead <- 0.05   		# Length of arrow heads
    aangle <- 30    		# Angle of arrow heads
    
    # Properties of species and site scores in scalings 1 and 2 –
    # Scaling 1: the species scores have norms of 1
    # Scaling 1: the site scores are scaled to variances = can.eigenvalues
    # Scaling 2: the species scores have norms of sqrt(can.eigenvalues)
    # Scaling 2: the site scores are scaled to variances of 1
    # --------------------------------------------------------------------
    
    ### This version reconstructs the original RDA output in L&L 2012, Section 11.1.3
    
    Tot.var <- res.rda$tot.chi            # Total variance in response data Y
    eig.val <- c(res.rda$CCA$eig, res.rda$CA$eig) # all eigenvalues
    Lambda <- diag(eig.val)               # Diagonal matrix of eigenvalues
    eig.val.rel <- eig.val / Tot.var      # Relative eigenvalues of Y-hat
    Diag <- diag(sqrt(eig.val.rel))       # Diagonal matrix of sqrt(relative eigenvalues)
    U.sc1 <- cbind(res.rda$CCA$v, res.rda$CA$v) # All species scores, scaling=1
    U.sc2 <- U.sc1 %*% sqrt(Lambda)       # Species scores, scaling=2
    colnames(U.sc2) <- colnames(U.sc1)
    n <- nrow(res.rda$CCA$u)              # Number of observations
    Z.sc2 <- cbind(res.rda$CCA$u, res.rda$CA$u) * sqrt(n - 1)  # "lc" site scores, scaling=2
    Z.sc1 <- Z.sc2 %*% sqrt(Lambda)       # "lc" site scores, scaling=1
    colnames(Z.sc1) <- colnames(Z.sc2)
    F.sc2 <- cbind(res.rda$CCA$wa, res.rda$CA$u) * sqrt(n - 1)  # "wa" site scores, scaling=2
    F.sc1 <- F.sc2 %*% sqrt(Lambda)       # "wa" site scores, scaling=1
    colnames(F.sc1) <- colnames(F.sc2)
    BP.sc2 <- res.rda$CCA$biplot          # Biplot scores, scaling=2 ; cor(Z.sc1, X)
    BP.sc2 <- cbind(BP.sc2, matrix(0, nrow = nrow(BP.sc2), 
                                   ncol = length(eig.val) - k))
    colnames(BP.sc2) <- colnames(F.sc2)
    BP.sc1 <- BP.sc2 %*% Diag             # Biplot scores, scaling=1
    colnames(BP.sc1) <- colnames(BP.sc2)
    
    if (!is.null(res.rda$CCA$centroids)) {
      centroids.sc2 <- res.rda$CCA$centroids * sqrt(n - 1) # Centroids, scaling=2
      centroids.sc2 <- cbind(centroids.sc2, matrix(0, nrow = nrow(centroids.sc2), 
                                                   ncol = length(eig.val) - k))
      colnames(centroids.sc2) <- colnames(F.sc2)
      centroids.sc1 <- centroids.sc2 %*% sqrt(Lambda)      # Centroids, scaling=1
      colnames(centroids.sc1) <- colnames(centroids.sc2)
    }
    centroids.present <- TRUE
    if (is.null(res.rda$CCA$centroids)) {
      centroids.present <- FALSE
      if (plot.centr | label.centr) {
        cat("\nNo factor, hence levels cannot be plotted with symbols;")
        cat("\n'plot.centr' is set to FALSE\n")
        plot.centr  <- FALSE
        label.centr <- FALSE
      }
    }
    
    if (is.null(select.spe)) {vec <- 1:n.sp} else {vec <- select.spe}
    
    if (scaling == 1) {
      if (site.sc == "lc") {
        sit.sc <- Z.sc1
      } else {
        sit.sc <- F.sc1
      }
      spe.sc <- U.sc1[vec, ]
      BP.sc  <- BP.sc1
      if (centroids.present)
        centroids <- centroids.sc1
    } else {
      # For scaling 2
      if (site.sc == "lc") {
        sit.sc <- Z.sc2
      } else {
        sit.sc <- F.sc2
      }
      spe.sc <- U.sc2[vec, ]
      BP.sc  <- BP.sc2
      if (centroids.present)
        centroids <- centroids.sc2
    }
    
    fact.spe <- 1
    fact.env <- 1
    if (centroids.present & (plot.centr | label.centr)) {
      to.plot <- which(!(rownames(BP.sc) %in% rownames(centroids)))
    } else {
      to.plot <- 1:nrow(BP.sc)
    }
    
    if (optimum) {
      fact.spe <-
        stretch(sit.sc, spe.sc, ax1, ax2, n, silent = silent)
      if (arrows.only) {
        fact.env <-
          stretch(sit.sc, BP.sc, ax1, ax2, n, silent = silent)
      } else {
        # arrows only==FALSE
        quant.env.present <- FALSE
        if (length(to.plot) > 0) {
          quant.env.present <- TRUE
          fact.env <-
            stretch(sit.sc, BP.sc[to.plot, ], ax1, ax2, n, silent = silent)
        }
      }
    }
    
    if (!silent)
      cat("fac.spe =", fact.spe, "   fact.env =", fact.env, "\n")
    spe.sc <- spe.sc * fact.spe * mult.spe
    BP.sc <- BP.sc * fact.env * mult.arrow
    lim <-
      larger.plot(
        sit.sc[, ],
        spe.sc[, ],
        BP.sc[, ],
        percent = mar.percent,
        move.origin = move.origin,
        ax1 = ax1,
        ax2 = ax2,
        k = k
      )
    if (!silent)
      print(lim)
    
    
    ### Drawing the triplot begins here ###
    
    # Draw the main plot
    mat <- rbind.data.frame(sit.sc, spe.sc, BP.sc)
    #  mat <- rbind(sit.sc[,1:k], spe.sc[,1:k], BP.sc[,1:k])
    
    if(ax1 <= k){
      titre <- "RDA triplot - Scaling"
    }
    else{
      titre <- "Biplot of residuals of RDA - Scaling"
    }    
    
    plot(mat[,c(ax1,ax2)], 
         type="n", 
         main=paste(titre, scaling, "-", site.sc),
         xlim=c(lim[1,ax1], lim[2,ax1]),
         ylim=c(lim[1,ax2], lim[2,ax2]), 
         xlab=colnames(U.sc1)[ax1], 
         ylab=colnames(U.sc1)[ax2],
         asp=1)
    abline(h=0, 
           v=0, 
           col="grey60")
    
    # Draw the site scores ("lc" or "wa" for RDA, implicitly "wa" only for residuals (PCA))
    if(plot.sites) {
      points(sit.sc[,ax1], sit.sc[,ax2], pch=20)
      if(label.sites)
        text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=pos.sites, cex=cex.char1)
    } else {
      if(label.sites)
        text(sit.sc[,ax1], sit.sc[,ax2], labels = rownames(sit.sc), col="black", pos=NULL, cex=cex.char1)
    }
    
    # Draw the species scores
    if(plot.spe) {
      arrows(0, 0, spe.sc[,ax1], spe.sc[,ax2], length=ahead, angle=aangle, col="red")
      if(label.spe)
        text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=pos.spe, cex=cex.char1)
    } else {
      if(label.spe)
        text(spe.sc[,ax1], spe.sc[,ax2], labels = rownames(spe.sc), col="red", pos=NULL, cex=cex.char1)
    }
    
    # Draw the explanatory variables (only if at least one axis is canonical)
    #
    if(ax1 <= k) {
      if(!arrows.only) {
        # 1. Quantitative variables
        if(quant.env.present & plot.env) {   # Print arrows and labels for quantitative var.
          arrows(0, 0, BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, length=ahead, angle=aangle, col="blue")
          if(label.env)   # Print labels for the quantitative variables
            text(BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, labels = rownames(BP.sc)[to.plot], col="blue", pos=pos.env, cex=cex.char2)
        } else {
          if(quant.env.present & !plot.env & label.env)   # Only print labels for quant. var.
            text(BP.sc[to.plot,ax1]*mult.arrow, BP.sc[to.plot,ax2]*mult.arrow, labels = rownames(BP.sc)[to.plot], col="blue", pos=NULL, cex=cex.char2)
        }
        #
        # 2. Centroids and labels of factor levels
        if(centroids.present & plot.centr) {   # Print symbols and labels for factor classes
          points(centroids[,ax1], centroids[,ax2], pch=18, cex=cex.point, col=col.centr)
          if(label.centr)
            text(centroids[,ax1], centroids[,ax2], labels = rownames(centroids), col="blue", pos=pos.centr, cex=cex.char2)
        } else {
          #
          if(centroids.present & !plot.centr & label.centr)   # Only print labels for classes
            text(centroids[,ax1], centroids[,ax2], labels = rownames(centroids), col="blue", pos=NULL, cex=cex.char2)
        }
      }
      
      # 3. All env. var.: plot arrows and labels for all var. in 'BP.sc', quant. and factors
      if(arrows.only) {
        arrows(0, 0, BP.sc[,ax1]*mult.arrow, BP.sc[,ax2]*mult.arrow, length=ahead, angle=aangle, col="blue")
        if(label.env)  # Print labels for the quantitative variables
          text(BP.sc[,ax1]*mult.arrow, BP.sc[,ax2]*mult.arrow, labels = rownames(BP.sc), col="blue", pos=pos.env, cex=cex.char2)
      }
    }
    #
  }

# Function hcoplot()
# Reorder and plot dendrogram with colors for groups and legend
#
# Usage:
# hcoplot(tree = hclust.object, diss = dissimilarity.matrix, k = nb.clusters, 
#	title = paste("Reordered dendrogram from",deparse(tree$call),sep="\n"))
#
# License: GPL-2 
# Author:  Francois Gillet, 23 August 2012
# Revised: Daniel Borcard, 31 August 2017

"hcoplot" <- function(tree, 
                      diss, 
                      lab = NULL,
                      k, 
                      title = paste("Reordered dendrogram from", 
                                    deparse(tree$call), 
                                    sep="\n"))
{
  require(gclus)
  gr <- cutree(tree, k=k)
  tor <- reorder.hclust(tree, diss)
  plot(tor, 
       labels = lab,
       hang=-1, 
       xlab=paste(length(gr),"sites"), 
       sub=paste(k,"clusters"), 
       main=title)
  so <- gr[tor$order]
  gro <- numeric(k)
  for (i in 1 : k)
  {
    gro[i] <- so[1]
    if (i<k) so <- so[so!=gro[i]]
  }
  rect.hclust(tor, 
              k = k, 
              border = gro + 1, 
              cluster = gr)
  legend("topright", 
         paste("Cluster", 1 :k ), 
         pch = 22, 
         col = 2 : (k + 1), 
         bty = "n")
}

plot.lda = function(lda.out, groups, colour.vec=NULL, plot.sites=1, plot.centroids=0, xax=1, yax=2, plot.env=TRUE, plot.ell=TRUE, title="LDA predicted classes", mul.coef=2, pos.names=NULL, col.env="black", xlim=NULL, ylim=NULL) 
  ### 
  # lda.out : Output file of function lda() in {MASS}.
  # groups  : Vector listing the group number for each object (factor or numeric).
  # colour.vec : personal vector with colour names, to be used for the groups in
  #    the plot instead of the standard colours. Example vector with 7 colours:
  #    col.vec = c("gray60","bisque","brown4","red","blue","darkgreen","orange4")
  #    Function colors() makes 657 colours available to R users.
  # plot.sites: 0 = Do not plot the sites.
  #             1 = plot symbols for the sites
  #             2 = print the site names
  # plot.centroids: 0 = Do not plot the group centroids.
  #                 1 = Plot the group centroids (triangles) with assigned 
  #                     group colours.
  # xax : The axis to be used for abscissa of the plot.
  # yax : The axis to be used for ordinate of the plot.
  # plot.env=TRUE: plot the explanatory variables to the plot. FALSE: Do not 
  #                plot them.
  # plot.ell=TRUE: plot the 95% coverage ellipses of the groups. FALSE: Do not
  #                plot them.
  # title : Allows user to customize the title that will print above the plot.
  # mul.coef : Multiplication factor for the length of the variable arrows. Some
  #            trial and error is needed here.
  # pos.names : Offset the names of the binary variables: 1=bottom, 2=left, 
  #             3=top, 4=right.
  # col.env="black" : Colour for the environmental variable arrows and names.
  # xlim, ylim : Vectors describing the minimum and maximum values of the plotted region.
  #
  # License: GPL-2
  # Authors: Daniel Borcard and Pierre Legendre
  
{
  if(class(lda.out) != "lda") stop("File lda.out was not produced by function lda of MASS")
  if(min(summary(as.factor(groups))) < 2) stop("There is at least one group with less than 2 observations")
  library(ellipse)
  coef = lda.out$scaling   # Standardized discriminant function coefficients
  k <- ncol(coef)  # Number of canonical axes
  lev <- length(levels(as.factor(groups)))  # Number of groups
  # print(c(k,lev))
  Fp <- predict(lda.out)$x
  centre <- matrix(NA,lev,k)
  for(i in 1:lev) { # Compute the group centroids in LDA space
    #	centre[i,] <- apply(Fp[groups==i,], 2, mean) }  
    centre[i,] <- apply(Fp[groups==levels(as.factor(gr))[i],], 2, mean) }  
  # print(centre)
  class.num <- as.numeric(predict(lda.out)$class) # Assignment of sites to classes
  # print(class.num)
  if(xax > k) stop("Their are not enough canonical axes; change the xax value")
  if(yax > k) stop("Their are not enough canonical axes; change the yax value")
  xlab=paste("LDA axis",xax," ")
  ylab=paste("LDA axis",yax," ")
  plot(Fp[,xax], Fp[,yax], type="n", main=title, xlab=xlab, ylab=ylab, xlim, ylim, asp=1)
  if(plot.sites==1) {   # Plot symbols for the sites
    if(is.null(colour.vec)) {
      points(Fp[,xax], Fp[,yax], pch=21, bg=class.num+1)
    } else {
      colour.sel <- colour.vec[class.num]
      points(Fp[,xax], Fp[,yax], pch=21, bg=colour.sel)
    }
  } else if(plot.sites==2) {
    if(is.null(colour.vec)) {
      text(Fp[,xax], Fp[,yax], row.names(Fp), col=class.num+1)
    } else {
      colour.sel <- colour.vec[class.num]
      text(Fp[,xax], Fp[,yax], row.names(Fp), col=colour.sel)
    }
  }	
  if(plot.centroids) {
    if(is.null(colour.vec)) {
      points(centre[,xax], centre[,yax], pch=24, bg=(1:lev)+1, cex=1.5)
    } else {
      colour.sel <- colour.vec[1:lev]
      # print(colour.sel)
      points(centre[,xax], centre[,yax], pch=24, bg=colour.sel, cex=1.5)
    }
  }
  abline(v=0, lty="dotted")
  abline(h=0, lty="dotted")
  # Draw 95% ellipses around the groups
  if(plot.ell) {
    for(i in 1:length(levels(as.factor(groups)))) { 
      #		cov <- cov(Fp[groups==i,c(xax,yax)])
      cov <- cov(Fp[groups==levels(as.factor(gr))[i],
                    c(xax,yax)])
      #		centre <- apply(Fp[groups==i,c(xax,yax)], 2, mean)
      centre <- apply(Fp[groups==levels(as.factor(gr))[i],
                         c(xax,yax)], 2, mean)
      lines(ellipse(cov, centre=centre, level=0.95))
    }
  }
  if(plot.env) { 
    arrows(x0=0, y0=0, x1=coef[,xax]*mul.coef, y1=coef[,yax]*mul.coef, 
           col=col.env, code=2, lty=1, length=0.1, lwd=1) 	
    if(!is.null(rownames(coef))) {
      text(x=coef[,xax]*mul.coef, y=coef[,yax]*mul.coef, 
           rownames(coef), col=col.env, cex=1, pos=pos.names) }
  }
}


polyvars <- function(X, degr = 2, raw = FALSE) 
{
  
  # A function computing polynomials of vectors within a matrix.
  # Contrary to function poly() on which it is based, this function 
  # only computes polynomials separately for each vector of the provided matrix,
  # e.g. x, x^2, x^3, and not combinations such as xy, x^2y and so on.
  #
  # License: GPL-2
  # Author:  Daniel Borcard
  #          December 2014, March 2017
  
  #
  # Usage
  # -----
  # polymatrix(X = rawdatamatrix, degr = 3, raw = FALSE)
  #
  # Arguments
  # ---------
  #    X: a matrix or data frame containing quantitative variables
  #
  #    degr: the degree to which the variables must be raised. Default: 2
  #
  #    raw: logical; if TRUE raw polynomials are computed directly from 
  #         the raw variables. If FALSE (default), orthogonal polynomials 
  #         are computed.
  #
  # Value
  # -----
  # A data frame containing the polynomials. In the output matrix, each
  # variable appears in turn, followed by its polynomial terms, e.g.
  # v1, v2_square, v2, v2_square, and so on.
  #
  # Details
  # -------
  # When raw = FALSE, the function computes orthogonal polynomial terms 
  # of each variable separately. This means that in the resulting matrix
  # the polynomial terms of each variable are orthogonal, but that they
  # are not orthogonal to the terms of the other variables.
  
  class.verif <- apply(X, 2, class)
  if (any(class.verif == "factor") | any(class.verif == "character") == TRUE)
    stop("No factor or character variables allowed.", call. = FALSE)
  
  ## Store or assign variable names
  if(!is.null(colnames(X)))
  {
    var.names <- colnames(X)
  }
  else
  {
    var.names <- paste("v", 1 : ncol(X), sep = "")
  }
  
  ## Compute polynomial terms
  X.poly <- matrix(0, nrow(X), ncol(X) * degr)
  for(i in 0: (ncol(X) - 1)) 
  {
    toto <- poly(X[, (i + 1)], degr, raw=raw)
    X.poly[,(i * degr + 1) : ((i + 1) * degr)] <- toto
  }
  
  if((ncol(X) * degr) > (nrow(X) - 1) ) 
  {
    cat("\n------------------------------------------------------------------")
    cat("\nWARNING: the number of polynomial terms is equal to or larger than")
    cat("\nthe number of observations.")
    cat("\n------------------------------------------------------------------\n")
  }
  
  ## Create new column names
  indices <- rep(1 : degr, ncol(X))
  tmp <- vector(length = ncol(X.poly))
  for(j in 1 : ncol(X))
  {
    tmp[(j * degr - degr + 1) : (j * degr)] <- rep(var.names[j], degr)
  }
  var.poly <- paste(tmp, indices, sep = ".")
  colnames(X.poly) <- var.poly
  
  X.poly.df <- as.data.frame(X.poly)
  
  X.poly.df
  
}


## Examples

## Construction of a fictitious matrix of 5 observations and 4 variables:
# env <- matrix(1:20, 5)

## Computation of orthogonal polynomials of degree 3:
# env.ortho.deg3 <- polymatrix(env, degr = 3)

## Computation of a matrix of raw polynomials of degree 4:
# env.raw.deg4 <- polymatrix(env, degr = 4, raw = TRUE)

screestick <- function(ev, las = 1, positive = TRUE, 
                       evmean = FALSE, relative = FALSE) 
{
  
  # Description
  
  # A function to draw a screeplot of a vector of eigenvalues
  # with overlay of the values predicted by the broken stick model
  #
  # This function mimics the default output of function 
  # screeplot.cca {vegan} but it works on a vector of eigenvalues, 
  # which makes it independent of the function and package that has 
  # computed the ordination.
  
  # Arguments
  
  # ev       vector of eigenvalues
  # las      orientation of the x axis labels
  # ww       width of bars in barplot
  # ss       space between bars in barplot
  # evmean   logical: should the mean of the eigenvalues be overlaid?
  # relative logical: should the eigenvalues be divided by their
  #          sum to represent relative proportions of variation?
  
  # Licence: GPL-2
  # Author:  Daniel Borcard
  #          September 2017
  
  require(vegan)
  
  if(is.vector(ev) == FALSE) 
    stop("Please provide a *vector* of eigenvalues", call. = FALSE)
  
  evs <- sort(ev, decreasing = TRUE)
  
  if(sum(abs(ev - evs)) > 0)
  {
    cat("The values were not in decreasing order. They have been sorted.\n")
  }
  
  if(positive == TRUE)
  {
    ntemp <- length(evs)
    evs <- evs[evs > 0]
    n <- length(evs)
    
    if(ntemp > n) 
    {
      cat("Warning: presence of", ntemp - n, "negative eigenvalues. \n")
      cat("Only the", n, "positive eigenvalues are considered.\n")
    }
  }
  else
  {
    n <- length(evs)
  }
  
  if(relative == TRUE)
  {
    evs <- evs/sum(evs)
  }
  
  names(evs) <- paste("EV", 1 : n)
  
  # Broken stick model
  broken <- bstick(n) * sum(evs)
  
  if(relative == TRUE)
  {
    ylabel = "Relative eigenvalue"
  }
  else
  {
    ylabel = "Eigenvalue"
  }
  
  # Scree plot
  barplot(evs,
          ylim = c(0, max(c(evs, broken))),
          main = "Eigenvalues and broken stick model",
          ylab = ylabel,
          las = las)
  
  # Overlay the lines and points representing the broken stick model.
  # The bars of the bar plot have width ww = 1 and are preceded and
  # separated by an interval ss = 0.2. There are n bars. Argument 
  # 'absc' below is computed to ensure that the line breaks and the 
  # corresponding points of the broken stick model are aligned with 
  # the centers of the bars of the scree plot.
  # ww and ss might be added as arguments in a later version.
  ww <- 1
  ss <- 0.2
  ws <- ww + ss
  n05 <- ww/2
  absc <- seq((n05 + ss), ((n * ws) - n05), by = ws)
  lines(
    absc, 
    broken,
    type = "o",
    pch = 1,
    col = "red"
  )
  
  if(evmean == TRUE) 
  {
    abline (h = mean(evs), col = "darkgray")
    legend(
      "topright",
      legend = c("Broken Stick", "Mean of eigenvalues"),
      pch = c(1, NA_integer_),
      lty = 1, 
      col = c("red", "gray"),
      bty = "n")
  }
  else
  {
    legend(
      "topright", 
      "Broken Stick", 
      pch = 1, 
      lty = 1, 
      col = "red", 
      bty = "n")
  }
}

