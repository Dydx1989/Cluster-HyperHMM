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


#AMR_binary_data=read.csv("AMR_binary.csv", row.names=1)

AMR_binary_data=read.csv("AMR_binary.csv", row.names=1)
dim(AMR_binary_data)
#  Select Optimal number of Cluster from Final-Data Using gap_stat method
# Elbow method

# gap_stat method
set.seed(23245)
gap_stat <- clusGap(AMR_binary_data, FUN = kmeans,
                    K.max = 25, B = 500)
# Identify the optimal number of clusters KK

# Find the first maximum value
first_max_index <- which(diff(sign(diff(gap_stat$Tab[, "gap"]))) < 0) + 1
first_max_value <- gap_stat$Tab[, "gap"][first_max_index[1]]  # Selecting only the first maximum
KK=first_max_index[1]
#print(gap_stat, method = "firstSEmax")
fviz_gap_stat(gap_stat)

km.res <- kmeans(AMR_binary_data, KK, nstart = 25)
# Visualize

fviz_cluster(km.res, data = AMR_binary_data,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

#++++++++++++++++++++++++++++++++++++++++++++++
#+####################################
# GENE SUBSET
type2kclu = data.frame(
  my_data =substr(colnames(t(AMR_binary_data)),1,KK),
  cluster=km.res$cluster)

Cluster_Tab=table(type2kclu)
generate_column_names <- function(num_columns) {
  # Generate column names like C1, C2, C3, ...
  column_names <- paste0("C", 1:num_columns)
  return(column_names)
}

colnames(Cluster_Tab) <- generate_column_names(KK)
head(Cluster_Tab)

## Save new data as Gene_cluster: Genes and Cluster 
write.csv(Cluster_Tab,"Gene_cluster.csv")
#*******************************************************************************************###
#*** Use the save data (Gene_cluster) and AMR_binary_data to generate cross sectional data**###
#*******************************************************************************************###


