library(pheatmap) ## for heatmap generation
#library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing 
library("cluster")
library("factoextra")
library("magrittr")
library(latticeExtra)
library(data.table)
library(Matrix)
library("data.table")
library("cluster")
library("factoextra")
library("magrittr")
library(NbClust)


Final_data=read.csv("Final_data_all.csv", row.names=1)

#  Select Optimal number of Cluster from Final-Data Using gap_stat method
# Elbow method

# gap_stat method
set.seed(23245)
gap_stat <- clusGap(Final_data, FUN = kmeans,
                    K.max = 25, B = 500)
print(gap_stat, method = "firstSEmax")
fviz_gap_stat(gap_stat)

km.res <- kmeans(Final_data, 8, nstart = 25)
# Visualize

fviz_cluster(km.res, data = my_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

#++++++++++++++++++++++++++++++++++++++++++++++
#+####################################
# GENE SUBSET

type2kclu = data.frame(
  my_data =substr(colnames(t(Final_data)),1,8),
  cluster=km.res$cluster)

Cluster_Tab=table(type2kclu)
colnames(Cluster_Tab) <- c("C1", "C2", "C3","C4", "C5", "C6","C7", "C8")
head(Cluster_Tab)

#write.csv(Cluster_Tab,"Cluster_Tab_8.csv")

##  Binary Matrix (M) generated from ( Final_data=Isolate*gene) (N1) and (Cluster_Tab_8=Gene *Clustes) (N2)
Final_data=read.csv("Final_data_all.csv", row.names=1)
N1=t(Final_data)

## Cluster and Genes
Cluster_gene=read.csv("Cluster_Tab_8.csv",row.names = 1)
#head(Cluster_gene[1:5])
#dim(Cluster_gene)
Data_2=as.matrix(Cluster_gene)
#dim(Data_2)
#head(Data_2[,1:6])
N2=(Data_2)
#head(N2[,1:6])
#dim(N2)

MM=matrix(NA,nrow =nrow(N1),ncol =ncol(N2))
dim(MM)
for(i in 1:nrow(N1)) {
  for(j in 1:ncol(N2)) {
    MM[i,j]=ifelse(any(N1[i,]&N2[,j]==1),1,0)
  }
}

#MM
Isolate_name<- rownames(N1)
Cluster_name<- colnames(N2)
rownames(MM) <- Isolate_name
colnames(MM) <- Cluster_name
Cross_sectional_data=MM
