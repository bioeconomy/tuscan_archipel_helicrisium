#######################################################################################################
# to install packages



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
############################################################################################
# references


# [1] https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

