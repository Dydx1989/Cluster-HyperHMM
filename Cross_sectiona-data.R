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
