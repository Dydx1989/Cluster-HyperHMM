##  Binary Matrix (MM) generated from (AMR_binary_data=Isolate*gene) (N1) and ("Gene_cluster=Gene *Clustes) (N2)
Final_data=read.csv("AMR_binary.csv", row.names=1)
N1=t(Final_data)

## Cluster and Genes
Cluster_gene=read.csv("Gene_cluster.csv",row.names = 1)
#head(Cluster_gene[1:8])
#dim(Cluster_gene)
Cluster_gene=as.matrix(Cluster_gene)
#dim(Data_2)
#head(Data_2[,1:6])
N2=(Cluster_gene)
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
# Rename row names
row.names(MM) <- 1:nrow(MM)
Cross_sectional_data=MM
### Save as MM as cross-section as text file (for HyperHMM ) and csv for phylognetic estimation. 
write.csv(Cross_sectional_data,"Cross_sectional_data.csv")
# Write data frame to a text file without column and row names
data_matrix <- as.matrix(Cross_sectional_data)
#Concatenate the values of each row without gaps
row_values <- apply(data_matrix, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "Cross_sectional_data.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



