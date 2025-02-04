################################################################################
# load libraries & set working directory & set ramdom seeds

setwd("C:\\aaa_lavori\\lav_elicrisio_2025")

#######################################################################################################
                 
source("aux_elicrisio.R")

#######################################################################################################
# Set seeds 

set.seed(123)

#######################################################################################################

# dati_tot=read.csv("https://raw.githubusercontent.com/bioeconomy/tuscan_archipel_helicrisium/main/dati_wild_elicrisio.csv",dec = ",")

dati_tot=read.csv("dati_wild_elicrisio.csv")

cat("\014") # clean console

#######################################################################################################
# Purging data by EDA 

# N=10

soliton_vector<-c(122,123, #EII21 EII22
                  275,288, #GII1 e GII21
                  180,177, #EIV23 e EIV20
                  287,291, #GII20 e GII24
                  326,211) #GIII28 e EV28

dati_tot$sesquiterpene8=NULL # missing value presence

dati_sel=dati_tot[-c(which(dati_tot$miX==1),soliton_vector,337:340),] # purged


# dim(dati_sel)
# 316  46


#######################################################################################################

dati_sel_rel=100*dati_sel[,8:48]/dati_sel$TOT_mono.sesqui

dati_monosesqui=cbind(100*dati_sel[,8:27]/dati_sel$TOT_Mono,
                   100*dati_sel[,27:48]/(dati_sel$TOT_mono.sesqui-dati_sel$TOT_Mono))
      
dati_monosesquimiche=cbind(100*dati_sel[,8:27]/dati_sel$TOT_Mono,
                      100*dati_sel[,28:48]/(dati_sel$TOT_mono.sesqui))

#######################################################################################################
# imposing michelozzi's approach

dati_sel_rel=dati_monosesquimiche

# [1] "a.pinene"            "fenchene"            "camphene"           
# [4] "b.pinene"            "b.myrcene"           "a.phellandrene"     
# [7] "a.terpinene"         "d.limonene"          "cineolo"            
# [10] "b.ocimene"           "g.terpinene"         "p.cymene"           
# [13] "terpinolene"         "neryl.propionate"    "nerol"              
# [16] "neryl.isovaleriate"  "neryl.acetate"       "linalool"           
# [19] "terpinen.4.ol"       "a.terpineol"         "a.ylangene"         
# [22] "a.copaene"           "isoitalicene"        "italicene"          
# [25] "cis.a.bergamotene"   "trans.a.bergamotene" "b.caryophillene"    
# [28] "guaia.6.9.diene"     "b.farnesene"         "g.curcumene"        
# [31] "a.muurolene"         "b.bisabolene"        "d.cadinene"         
# [34] "b.selinene"          "b.curcumene"         "selinene"           
# [37] "a.curcumene"         "C15H26O2"            "guaiol"             
# [40] "sesquiterpene.4"     "rosifoliol"  

write.xlsx(list(aree=dati_sel,perc_relativi=dati_sel_rel),"matrice_aree_dati.xlsx") # to write data in local


# if scaled needed decomment

# dati_monosesquimiche_std <- dati_monosesquimiche %>% mutate_all(~(scale(.) %>% as.vector))

#######################################################################################################
# Correlation analisys

# correlation matrix and significance

mat_corr_sel <- cor(dati_sel[,8:48],method = "pearson" )

p.mat_wild <- cor.mtest(as.matrix(mat_corr_sel))

# palette

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

# plot

pdf(file="06.CORR_plot_pearson_tot.pdf",width=22,height=20)
corrplot(mat_corr_sel , method="color", col=col(200),  
         type="upper", 
         order="hclust",
         hclust.method = "average",
         addrect = 3, 
         tl.cex = 1,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", 
         tl.srt=45, #Text label color and rotation
         p.mat = p.mat_wild, # Combine with significance
         sig.level = 0.001, 
         insig = "blank", # hide correlation coefficient on the principal diagonal
         diag=F 
)
dev.off()

highcorrelated=c(13, #a_terpinolene
                 5,  #b_myrcene
                 6,  #a_phellandrene
                 34, #b_selinene
                 21, # a_ylangene
                 25, #cis_bergamotene
                 39, #guaiol
                 23  #iso_italicene
)
names(dati_sel_rel)[highcorrelated]



###################################################
# vif selection for LDA

multiColl::CN(dati_sel_rel) # 44.40231 col linearity!!!

a=multiCol(dati_sel_rel, graf = TRUE)


id_VIF=which(names(dati_sel_rel) %in% c("italicene","isoitalicene","a.ylangene","fenchene","d.limonene","b.pinene","a.muurolene")==T)

dati_sel_rel_LDA=dati_sel_rel[, -id_VIF]

a=multiCol(dati_sel_rel_LDA, graf = TRUE)

saveRDS(dati_sel_rel_LDA,"dati_sel_rel_LDA.rds")


###################################################################################
# Wild sampling classification use
#   WorkGgroup  (CI-CII-CIII  EI-EII-EIII-EIV-EV-EVI   GI-GII-GIII)
#   Island classification Isole
#   Capraia (CI-CII-CIII)
#   Elba  (EI-EII-EIII-EIV-EV-EVI)
#   Giglio  (GI-GII-GIII)
#   Novel classification   Work_prop2
#   Gruppo capraia  (CI-CII-CIII-EI) 
#   Gruppo Elba     (EIII)
#   Enfola          (EII-EIV)
#   Gruppo Giglio   (EV-EVI-GI-GII-GIII)
################################################################################

####################################################################################
# Clustering kmeans coupled by silhouette analisys by using all data

df=dati_sel_rel

silhouette_score <- function(k){
  km <- kmeans(df, centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(df))
  mean(ss[, 3])
}

k <- 2:10
avg_sil <- sapply(k, silhouette_score)

png("silhoutette_clustering_all_data.png")
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
dev.off()

res.hk <-hkmeans(dati_sel_rel, 4,hc.metric =  'manhattan')

dend_all=fviz_dend(res.hk,
          cex = 0.3,
          palette = c("orange", "green3", "red", "blue"), 
          rect = TRUE,
          rect_border = "jco", 
          rect_fill = TRUE)

dend_all

dati_sel$WorkGroup_dup[res.hk$cluster==4]  # Gruppo Elba in prevalenza ( Enfola "E_II" Elba E_IV)
dati_sel$WorkGroup_dup[res.hk$cluster==3]  # Gruppo Giglio in prevalenza  "E_V" "G_I" "G_II" "G_III"
dati_sel$WorkGroup_dup[res.hk$cluster==2]  # Gruppo Capraia in prevalenza 
dati_sel$WorkGroup_dup[res.hk$cluster==1]  # Gruppo Misto 1

dati_sel$Work_prop2[res.hk$cluster==4]  # 
dati_sel$Work_prop2[res.hk$cluster==3]  # Gruppo Giglio in prevalenza 
dati_sel$Work_prop2[res.hk$cluster==2]  # Gruppo Capraia in prevalenza Gruppo capraia (CI-CII-CIII-EI) 
dati_sel$Work_prop2[res.hk$cluster==1]  # Gruppo Misto 1


####################################################################################################
# Statistic summaries

  dati_sel_rel=100*dati_sel[,8:48]/dati_sel$TOT_mono.sesqui
  
  dati_sel_rel_group=data.frame(dati_sel_rel,WorkGroup=dati_sel$WorkGroup)
  dati_sel_rel_clus=data.frame(dati_sel_rel,cluster4=res.hk$cluster)
  dati_sel_rel_group_miche=data.frame(dati_monosesquimiche,WorkGroup=dati_sel$WorkGroup)
  
  summarise_group <- dati_sel_rel_group |> tbl_summary(by = WorkGroup,
                                              statistic = list(
                                                all_continuous() ~ "{mean} ± ({sd})"),digits = all_continuous() ~ 1)|>as_gt() |> 
                                                                                                                          gt::gtsave(filename = "summarise_group_data_single.docx")
  summarise_miche <- dati_sel_rel_group_miche |> tbl_summary(by = WorkGroup,
                                                       statistic = list(
                                                         all_continuous() ~ "{mean} ± ({sd})"),digits = all_continuous() ~ 1)|>as_gt() |> 
                                                                                                                                        gt::gtsave(filename = "summarise_group_data_miche.docx")
  
  summarise_clus <- dati_sel_rel_clus |> tbl_summary(by = cluster4,
                                                       statistic = list(
                                                         all_continuous() ~ "{mean} ± ({sd})"),digits = all_continuous() ~ 1)|>as_gt() |> gt::gtsave(filename = "summarise_cluster_data_single.docx")

################################################################################
# kruskal   
aa=col_kruskalwallis(dati_sel_rel_group_miche[1:41], dati_sel_rel_group_miche$WorkGroup)
aa  |> kable(format="html") |> cat(file = "ks_table.docx")  

# col_kruskalwallis(dati_sel_rel_group_miche[1:41], dati_sel_rel_clus$cluster4)
  
################################################################################
#  Clustering analysis on sampling group means # Lorenzo marini


dati_wild_2020_mean <- aggregate(dati_sel_rel, list(dati_sel$WorkGroup), mean,  na.rm=TRUE)
d_wild_2020 <- dist(dati_wild_2020_mean , method = "euclidean")
hr_wild_2020 <- hclust(d_wild_2020, method = "complete")
clustelist=as.numeric(cutree(hr_wild_2020, k = 4))
sing.COD <- levels(factor(dati_sel$WorkGroup)) 
hr_wild_2020$labels <- sing.COD
plot(hr_wild_2020, main="Wild populations 2020", hang = 0.5)

# Plot color branches based on cutting the tree into 4 clusters:

dend <- as.dendrogram(hr_wild_2020)

color_cluster_dendro <- c("#FF0000", #cluster 1
                          "#669900", #cluster 2
                          "#00CCCC", #cluster 3
                          "#9933FF") #cluster 4

dend4 <- dendextend::color_branches(dend, k = 4, col=color_cluster_dendro)

png("dendrogram_wild.png",width = 500, height = 500) #500 x 500

plot(dend4,main="HCA - Wild populations (k = 4)")

dev.off()

dati_sel_rel_clus=data.frame(dati_wild_2020_mean,cluster4=as.numeric(cutree(hr_wild_2020, k = 4)))
dati_sel_rel_clus=dati_sel_rel_clus[,-1]

summarise_clusters <-  dati_sel_rel_clus |> tbl_summary(by = cluster4,
                                                                   type=list(terpinolene~ "continuous",
                                                                             a.terpineol~ "continuous",
                                                                             b.bisabolene~ "continuous",
                                                                             d.cadinene~ "continuous",
                                                                             sesquiterpene.4~ "continuous"),
                                                        statistic = list(
                                                          all_continuous() ~ "{mean} ± ({sd})"),digits = everything() ~ 1)|>as_gt() |> 
                                                                                                                                    gt::gtsave(filename = "summarise_clus_data_mean.docx")







#######################################################################
# 12 rows and cluster obtained labels

sing.LOC <-  c("Cluster_4",
               "Cluster_4",
               "Cluster_2",
               "Cluster_2",
               "Cluster_1", # EII
               "Cluster_2",
               "Cluster_1", # EIV
               "Cluster_3",
               "Cluster_4",
               "Cluster_3", # 10
               "Cluster_3",
               "Cluster_3")


dati_wild_2020_mean_data=dati_wild_2020_mean[-1]
row.names(dati_wild_2020_mean_data)=dati_wild_2020_mean$Group.1
pca_prcomp_wild_2020 <- prcomp(dati_wild_2020_mean_data,scale = T) 


color_cluster_dendro <- c("#FF0000", #cluster 1
                          "#669900", #cluster 2
                          "#00CCCC", #cluster 3
                          "#9933FF") #cluster 4


type=c(rep("monoterpenes",21),rep("sesquiterpenes",20))

#######################################################################
# plotting PCA

# explained variance

scree_plot=factoextra::fviz_eig(pca_prcomp_wild_2020, main = "Wild population 2020 - Explained variance",addlabels = TRUE, ylim = c(0, 50))

scree_plot

ggsave(filename="scree_plot_PCA_wild_data.png",scree_plot,width=900,height=500,units = "px")

PCA_final_wild_2020 <- factoextra::fviz_pca_biplot(pca_prcomp_wild_2020, 
                                       # fill individuals by groups
                                       geom.ind = c("text","point"),
                                       col.ind = "black",
                                       fill.ind = sing.LOC,
                                       pointshape = 21, #i numeri definiscono  uno stile di punti!
                                       palette = "color_cluster_dendro", 
                                       addEllipses = T, 
                                       ellipse.level = 0.10,
                                       ellipse.type = "convex",
                                       geom.var = c("arrow", "text"), 
                                       arrowsize = 0.3,
                                       labelsize = 3,
                                       col.var = type,
                                       legend.title = list(fill = "Cluster", color = "Compounds"),
                                       title = "PCA - Wild population terpenes data",
                                       repel = T  
) +
  ggpubr::color_palette("color_cluster_dendro") +
  theme(legend.direction = "horizontal") +
  theme(legend.position = "bottom") +
  theme(legend.box = "vertical")

ggsave(filename="PCA_wild_data1.png")


#ggsave(filename="PCA_wild_data.png",PCA_final_wild_2020,width=900,height=600,units = "px")













####################################################################################








####################################################################################
# Supervised methods : LDA  
# LDA
# Prior probabilities of groups: 
# the proportion of training observations in each group. 
# Group means: group center of gravity. 
# Shows the mean of each variable in each group.
# Coefficients of linear discriminants: 
# Shows the linear combination of predictor variables that are used to form the LDA d
# 
#########################################################################
# define matrix of data

X=dati_sel_rel_LDA

#########################################################################

Y=dati_sel$WorkGroup # LDA per campionamenti


training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)

dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]

model_folda_campio <- folda(datX = X, response = Y, subsetMethod = "all")
model_folda_campio
model <- lda(Y~., data =dati_train)
predictions_train_lda <- model %>% predict(dati_train)
predictions_lda <- model %>% predict(dati_test)
predictions_lda_molda <- predict(model_folda_campio,dati_test)

model_accuracy_lda=mean(predictions_lda$class==dati_test$Y)
LD_1_2=model$scaling[,1:2] 
eigvalues<-scale(model$scaling, F, sqrt(colSums(model$scaling^2)))
manova1<-manova(as.matrix(X)~Y)
summary(manova1,test="Wilks")
aa=discriminant.significance(model$scaling,
                             dim(dati_train)[2]-1,
                             length(unique(dati_train$Y)),
                             dim(dati_train)[1])
gw_obj = greedy.wilks(Y ~ ., data=dati_train, niveau= 0.5) 
gw_obj

##################################################################################################################
# Accuracy statistics  classification

res_lda_campio_molda=caret::confusionMatrix(factor(predictions_lda_molda),factor(dati_test$Y))
res_lda_campio=caret::confusionMatrix(predictions_lda$class,factor(dati_test$Y))

kable(as.table(res_lda_campio_molda$table), "html") %>%
   kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_matrix_pops_sampling_sites.docx")

kable(as.table(round(res_lda_campio_molda$overall,2)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_accuracy_pops_sampling_sites.docx")

kable(as.table(format(res_lda_campio_molda$byClass,digits=1)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_acc_byClass_pops_sampling_sites.docx")

###########################################################################################################################################
# plots 

dataset = data.frame(popolazioni = Y[training.samples], lda = predictions_train_lda$x)
centroids <- aggregate(.~popolazioni, dataset, mean)[,1:3]

#############################################################################################################
# plot LDA sampling populations


eli_lda <- lda_ord(dati_train[,1:34], dati_train[,35], axes.scale = "standardized")

eli_lda %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  mutate_rows(
    species = grouping,
    data = ifelse(.element == "active", "centroid", "case")
  ) %>%
  ggbiplot() +
  theme_bw() +
  geom_rows_point(aes(
    color = grouping,
    size = data, 
    alpha = data
  )) +
  ggtitle('Linear discriminant analysis biplot of wild populations') +
  expand_limits(y = c(-3, 5))+
  labs(color = "Populations",subtitle="By sampling sites")

##################################################################################################################
# Isole LDA 

Y=dati_sel$Isole 


training.samples <- Y %>% createDataPartition(p = 0.8, list = FALSE)

dati_train=X[training.samples, ]
dati_test=X[-training.samples, ]
dati_test$Y=Y[-training.samples]
dati_train$Y=Y[training.samples]


model_folda_isole <- folda(datX = X, response = Y, subsetMethod = "all")
model_folda_isole
model <- lda(Y~., data =dati_train)
predictions_train_lda <- model %>% predict(dati_train)
predictions_lda <- model %>% predict(dati_test)
predictions_lda_molda <- predict(model_folda_isole,dati_test)
model_accuracy_lda=mean(predictions_lda$class==dati_test$Y)



LD_1_2=model$scaling[,1:2] 
eigvalues<-scale(model$scaling, F, sqrt(colSums(model$scaling^2)))
manova1<-manova(as.matrix(X)~Y)
summary(manova1,test="Wilks")
aa=discriminant.significance(model$scaling,
                             dim(dati_train)[2]-1,
                             length(unique(dati_train$Y)),
                             dim(dati_train)[1])
gw_obj = greedy.wilks(Y ~ ., data=dati_train, niveau= 0.5) 
gw_obj

##################################################################################################################
# Accuracy statistics  classification

res_lda_isole_molda=caret::confusionMatrix(factor(predictions_lda_molda),factor(dati_test$Y))
res_lda_isole=caret::confusionMatrix(predictions_lda$class,factor(dati_test$Y))

kable(as.table(res_lda_isole_molda$table), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_matrix_pops_isole.docx")

kable(as.table(round(res_lda_isole_molda$overall,2)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_accuracy_pops_isole.docx")

kable(as.table(format(res_lda_isole_molda$byClass,digits=1)), "html") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  cat(., file = "confusion_acc_byClass_pops_isole.docx")

###########################################################################################################################################
# LDA plots Isole

eli_lda <- lda_ord(dati_train[,1:34], dati_train[,35], axes.scale = "standardized")

eli_lda %>%
  as_tbl_ord() %>%
  augment_ord() %>%
  mutate_rows(
    species = grouping,
    data = ifelse(.element == "active", "centroid", "case")
  ) %>%
  ggbiplot() +
  theme_bw() +
  geom_rows_point(aes(
    color = grouping,
    size = data, 
    alpha = data
  )) +
  ggtitle('Linear discriminant analysis biplot of wild populations') +
  expand_limits(y = c(-3, 5))+
  labs(color = "Isole",,subtitle="By islands")

#############################################################################





##########################################################################
# Code References
# https://davetang.github.io/muse/pheatmap.html
# https://www.geeksforgeeks.org/draw-heatmap-with-clusters-using-pheatmap-in-r/
# https://stats.stackexchange.com/questions/592404/discriminant-analysis-when-assumptions-are-violated
# https://stats.stackexchange.com/questions/501367/does-lda-need-normality-assumption-for-its-dimensionality-reduction
# https://rpubs.com/sebnemer/discriminantR
# https://corybrunson.github.io/ordr/reference/lda-ord.html


# D. A. Belsley (1991). Conditioning diagnostics: collinearity and weak dara in regression. John Wiley & Sons, New York.


#' Kolde R (2019). _pheatmap: Pretty Heatmaps_. R package version 1.0.12,
#' <https://CRAN.R-project.org/package=pheatmap>.
#' 
#' 
#' @Manual{,
#'   title = {pheatmap: Pretty Heatmaps},
#'   author = {Raivo Kolde},
#'   year = {2019},
#'   note = {R package version 1.0.12},
#'   url = {https://CRAN.R-project.org/package=pheatmap},
#' }

#' Beck M (2022). _ggord: Ordination Plots with ggplot2_. R package version
#' 1.1.7.
#' 
#' @Manual{,
#'   title = {ggord: Ordination Plots with ggplot2},
#'   author = {Marcus W. Beck},
#'   year = {2022},
#'   note = {R package version 1.1.7},
#' }
#' 
#' Brunson J (2022). _ordr: A 'tidyverse' Extension for Ordinations and
#' Biplots_. R package version 0.1.1,
#' <https://CRAN.R-project.org/package=ordr>.
#' 
#' Una voce BibTeX per gli utenti LaTeX è
#' 
#' @Manual{,
#'   title = {ordr: A 'tidyverse' Extension for Ordinations and Biplots},
#'   author = {Jason Cory Brunson},
#'   year = {2022},
#'   note = {R package version 0.1.1},
#'   url = {https://CRAN.R-project.org/package=ordr},
#' }