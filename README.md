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

### Subset out everything with a nonzero entry in the resistance info column

### Subset out everything that is exactly K. pneumoniae (our particular species of interest)

### Assume (for now) that predictor genes are 0 (empty cell) or 1 (non-empty cell)
###  Heatmap plot that reveal the structure and the pattern of the data after First preprocessing



![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/efa846cd-b20b-4570-b78e-2a07699dfb49)

###   Data Preprocessing: Stage 2
In this stage, we remove all the isolates with completely empty cells and retain only isolates with nonempty cells


![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/349a6f76-b7e4-4bb1-910d-ee356f6de2fc)



