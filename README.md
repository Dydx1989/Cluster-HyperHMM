# Cluster-HyperHMM
##  HyperHMM Figure Illustration
Hypercubic Inference for Large-Scale Genomic Data.
Cluster-HyperHMM is an extension of HyperHMM [1], from https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098. 
![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/feb73be9-258f-4885-96b4-485dd57ce505)

## Cluster-HyperHMM Figure Illustration
### The extendent version of HyperHMM (Cluster-HyperHMM) is presented below:

![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/b0fabdb0-07e5-42fb-b2db-85b2f6ecfc26)

## Synthetic Illustration of Data Preprocessing in Cluster-HyperHMM 

```{r}
## Synthetic Data:  5 Isolates and  4 Genes
SYN1= matrix(c(1,1,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,0,1,0),nrow = 5)
SYN1
##      [,1] [,2] [,3] [,4]
## [1,]    1    0    0    1
## [2,]    1    0    0    0
## [3,]    0    1    1    0
## [4,]    0    0    0    1
## [5,]    1    0    1    0
## Synthetic Data: 4 Genes and 3 Clusters
SYN2=matrix(c(1,0,1,0,0,0,0,1,0,1,1,0),nrow = 4)
SYN2
##      [,1] [,2] [,3]
## [1,]    1    0    0
## [2,]    0    0    1
## [3,]    1    0    1
## [4,]    0    1    0
## Generate New Binary values with respect to Isolate and Cluster
##

MM=matrix(NA,nrow =nrow(SYN1),ncol =ncol(SYN2))
for(i in 1:nrow(SYN1)) {
  for(j in 1:ncol(SYN2)) {
    MM[i,j]=ifelse(any(SYN1[i,]&SYN2[,j]==1),1,0)
  }
}

MM
##      [,1] [,2] [,3]
## [1,]    1    1    0
## [2,]    1    0    0
## [3,]    1    0    1
## [4,]    0    1    0
## [5,]    1    0    1
```

## Klebsiella pneumnonia Data Preprocessing in Cluster-HyperHMM 

###   Data Preprocessing: Stage 1

```{r setup,echo=TRUE}
library(readxl)
raw_dff=read.csv("C://Users//HP//Documents//BIGSdb_3708197_1272842014_68324.csv")
```
### Subset out everything with a nonzero entry in the resistance info column


```{r,echo=TRUE}
data_dff_1= subset(raw_dff,resistance_info>0)
dim(data_dff_1)
head(data_dff_1[1:6])
#write.csv(data_dff_1,"data_dff_1_csv.csv")
```
### Subset out everything that is exactly K. pneumoniae (our particular species of interest)

```{r, echo=TRUE}
data_dff_2=subset(data_dff_1,taxonomic_designation=="K. pneumoniae")
dim(data_dff_2)
head(data_dff_2[1:6])
#write.csv(data_dff_2,"data_dff_2.csv")
```
### Assume (for now) that predictor genes are 0 (empty cell) or 1 (non-empty cell)

```{r}
library(dplyr)
data_dff_new <- data_dff_2 %>% replace(is.na(.), 0)
head(data_dff_new[1:5])
#write.csv(data_dff_new,"data_dff_new.csv")
data_dff_new[data_dff_new==""] <- 0
data_dff_new_1=data_dff_new[1:14]
dim(data_dff_new_1)
data_dff_new_2=data_dff_new[14:dim(data_dff_new)[2]]
dim(data_dff_new_2)
data_dff_new_3=data.frame(lapply(data_dff_new_2, function(x) { 
  if(is.factor(x)) {
    x <- as.character(x)
  }
  x[x!="0"]="1"
  x}))
head(data_dff_new_3[1:5])
#write.csv(data_dff_new_3,"data_dff_new_3.csv",row.names=FALSE)

Data_trans=data.frame(data_dff_new_1,data_dff_new_3)
head(Data_trans[1:5])

#write.csv(Data_trans,"Data_trans.csv")

```
####  Heatmap plot that reveal the structure and the pattern of the data after First preprocessing

```{r,warning=FALSE}
data_dff_3=read.csv("C://Users//HP//Documents//data_dff_new_3.csv")
head(data_dff_3[1:5])
class(data_dff_3$aac2p_Ia)
data_dff_3_1 = as.data.frame(sapply(data_dff_3, as.numeric))
class(data_dff_3_1$aac2p_Ia)
data_dff_3_2=as.matrix(data_dff_3_1)
#data_dff_3_2=scale(data_dff_3_2)
library(ggplot2)
library("gplots")
heatmap.2(data_matrix, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

```


![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/efa846cd-b20b-4570-b78e-2a07699dfb49)

###   Data Preprocessing: Stage 2
In this stage, we remove all the isolates with completely empty cells and retain only isolates with nonempty cells

```{r ,warning=FALSE,echo=FALSE,message=FALSE, fig.width=10,fig.height=7}

library(pheatmap) ## for heatmap generation
#AMR_gene_RR=Data_trans with only genes
AMR_gene_RR <- read.csv("AMR_gene_RR.csv", row.names=1)
#dim(AMR_gene_RR)
#head(AMR_gene_RR[1:5])

remove_zero_cols <- function(df) {
  rem_vec <- NULL
  for(i in 1:ncol(df)){
    this_sum <- summary(df[,i])
    zero_test <- length(which(this_sum == 0))
    if(zero_test == 6) {
      rem_vec[i] <- names(df)[i]
    }
  }
  features_to_remove <- rem_vec[!is.na(rem_vec)]
  rem_ind <- which(names(df) %in% features_to_remove)
  df <- df[,-rem_ind]
  return(df)
}

amrgene_1=remove_zero_cols(AMR_gene_RR)


```

![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/349a6f76-b7e4-4bb1-910d-ee356f6de2fc)



