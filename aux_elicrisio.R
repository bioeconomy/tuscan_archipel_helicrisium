#######################################################################################################
# to install packages

#######################################################################################################
# before to run code you have to install R packages:  

# please use pak::pkg_install function available in R packages pkg 
# install.packages("pak")
# pak::pkg_install(c("ordr","ordr.extra","devtools","tidyverse","knitr","cluster","factoextra",
#                   "flexclust","clustertend","ClusterR","Cluster","smacof",
#                    "gridExtra","MASS","paran","BBmisc","ggplot2","ggcorrplot","ggbiplot",
#                    "caret","heplots","spls","ggpubr","performance","FactoMineR",
#                    "CCA","flextable","kableExtra","klaR","hopkins","clusterSim",
#                    "matrixTests","multiColl","moments","psych","skimr","gtsummary"))
# pak::pkg_install(c("ordr","ordr.extra",
# or in ordinary way
# install.packages(c("devtools","tidyverse","knitr","cluster","factoextra",
#                   "flexclust","clustertend","ClusterR","Cluster","smacof",
#                    "gridExtra","MASS","paran","BBmisc","ggplot2","ggcorrplot","ggbiplot",
#                    "caret","heplots","spls","ggpubr","performance","FactoMineR",
#                    "CCA","flextable","kableExtra","klaR","hopkins","clusterSim",
#                    "matrixTests","multiColl","moments","psych","skimr","gtsummary"))
#
# install.packages(c("ordr","ordr.extra")
#################################################################################################################

##########################################################################
# load libraries

library(dplyr)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(openxlsx)
library(gridExtra)
library(cluster)
library(MASS)
library(FactoMineR)         # Multivariate Exploratory Data Analysis and Data Mining
library(factoextra)
library(ordr)
library(ordr.extra)
library(caret)
library(multiColl)
library(klaR)
library(folda)
########################################################################
library(knitr)
library(kableExtra)
library(matrixTests)
library(gtsummary)          # publication ready summary tables
library(ggpubr)             # publication ready data visualization in R
library(ggcorrplot)
library(dendextend)

############################################################################################
# functions

standardize <- function(x) {(x - mean(x))/sd(x)}

stat_groups=function(x,g,xlabv="Groups",ylabv="Compound") {
  myls=list();
  myls[[1]]=aov(x ~ g)
  myls[[2]]=shapiro.test(myls[[1]]$residuals)
  # not normal if 0.05 >
  myls[[3]]=bartlett.test(x ~ g) # not ok if 0.05 >
  # not homogeneous if 0.05 >
  myls[[4]]=kruskal.test(x ~ g)
  # differs 0.05 >
  myls[[5]]=boxplot(x ~ g, xlab=xlabv, ylab=ylabv)
  return(myls)
}

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}

#######################################################################################################
# useful function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}



#######################################################################################################
# useful function
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

##############################################################
# useful functions

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


discriminant.significance <- function(eigenvalues, p, k, N) {
  w <- N - 1 - .5 * (p + k)
  t <- sqrt((p^2 * (k - 1)^2 - 4) / (p^2 + (k - 1)^2 - 5))
  df1 <- p * (k - 1)
  df2 <- w * t - .5 * (p * (k - 1) - 2)
  lambda1 <- prod(1 / (1 + eigenvalues))
  f1 <- (1 - lambda1^(1/t)) / (lambda1^(1/t)) * df2 / df1
  p1 <- pf(f1, df1, df2, lower.tail = F)
  
  result <- NULL
  for (i in 2:p) {
    m <- i
    
    if (m == p) {
      t.i <- 1
    } else {
      t.i <- sqrt(((p - m + 1)^2 * (k - m)^2 - 4) / ((p - m + 1)^2 + (k - m)^2 - 5))
    }
    
    df1.i <- (p - m + 1) * (k - m)
    df2.i <- w * t.i - .5 * ((p - m + 1) * (k - m) - 2)
    lambda.i <- prod(1 / (1 + eigenvalues[i:p]))
    f.i <- (1 - lambda.i^(1/t.i)) / lambda.i^(1/t.i) * df2.i / df1.i
    p.i <- pf(f.i, df1.i, df2.i, lower.tail = F)    
    result <- rbind(result, data.frame(lambda.i, f.i, p.i))
  }
  res <- rbind(c(lambda1, f1, p1), result)
  colnames(res) <- c('Lambda', 'Approximate F', 'p-value')
  return(res)
}

# y           LDA score (`predict(iris.lda)$x`)
# a = 0       (no intercept since the data is centered)
# b1, b2, ... LDA coefficients (`$scaling`)
# x1, x2, ... centered data (at mean of group `$means`)

# LD model$scaling[,1:2] 
# eigvalues<-scale(model$scaling, F, sqrt(colSums(model$scaling^2)))
# N <- dim(data)[1] raws
# p <- dim(data)[2] - 1 column
# k <- length(unique(root$Tree.Number)) # gruppi


############################################################################################
# references


# [1] https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

