# Cluster-HyperHMM

Cluster large-scale data and perform evolutionary inference using HyperHMM https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098 (code https://github.com/StochasticBiology/hypercube-hmm ).

![image](https://github.com/user-attachments/assets/d14fadb2-9c8a-46bf-837f-1a12ed239f89)

Overview
---

A pipeline will generally involve acquiring and pre-processing data, clustering, then HyperHMM inference. The latter two steps are case-general and are performed in `cHHMM-expts_Kleborate.R`. 

Several case studies are included:

* _Kleborate_. Data on *Klebsiella pneumoniae* genomes from around the world. `Download_and_curate-kleborate_data.R` downloads and curates these data. Run `cHHMM-expts_Kleborate.R` with `expt = "Kleborate"` (default) for this.
* _Malaria_. Data on severe malaria progression in patients. Data is already included in `jallow_dataset_binary_with2s.csv`. Run `cHHMM-expts_Kleborate.R` with `expt = "malaria"` for this.

The outputs will include illustrations of the cluster selection and identities, summaries of the inferred dynamics under different evolutionary assumptions, and comparisons across assumptions and clustering protocols.

# Cluster-HyperHMM Steps
 ### Data preprocessing stage: See synthetic example and real-life data on K. Pneumania below. 
1. Read in binary matrix (**AMR_binary.csv: isolate and genes matrix**)
2. Perform clustering toget optimal number of cluster to generate  **Gene_cluster.csv:Cluster and Genes matrix**. The R codes is **Cluster_code.R**
3. Produce two "working" datasets in txt – cross-sectional (ie just cluster sets) and phylogenetic (estimating a relationship by similarity between bitstring): **Cross_sectional_data.txt** and **phylogeny_data.txt**. The R codes  are **Cross_sectiona-data.R** and **Phylogenetic-Estimation.R**
4. Construct transition sets for HyperHMM – for cross-sectional this is just 0000 -> (each observation); for phylogenetic it's (estimated) ancestor -> descendant pairs. 
5. Run HyperHMM on both cross-sectional and phylogenetic transition sets and plot output of each, and summary matrix folding both outputs together. Rcode for cross and phylogeny is **Run HyperHMM on both cross-sectional and phylogenetic.r** and for cross-continents comparison. 
##  HyperHMM Figure Illustration
Hypercubic Inference for Large-Scale Genomic Data.
Cluster-HyperHMM is an extension of HyperHMM [1], from https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098. 
![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/feb73be9-258f-4885-96b4-485dd57ce505)

## Cluster-HyperHMM Figure Illustration
### The extendent version of HyperHMM (Cluster-HyperHMM) is presented below:

![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/b0fabdb0-07e5-42fb-b2db-85b2f6ecfc26)
# Preprocessing steps for Raw data to form Binary Matrix 
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

## Klebsiella pneumnonia Data Preprocessing in Cluster-HyperHMM (Raw:BIGSdb_3708197_1272842014_68324.csv)

###   Data Preprocessing

### Subset out everything with a nonzero entry in the resistance info column

### Subset out everything that is exactly K. pneumoniae (our particular species of interest)

### Assume (for now) that predictor genes are 0 (empty cell) or 1 (non-empty cell)
###  Heatmap plot that reveal the structure and the pattern of the data after First preprocessing



![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/efa846cd-b20b-4570-b78e-2a07699dfb49)

###   Data Preprocessing 
In this stage, we remove all the isolates with completely empty cells and retain only isolates with nonempty cells


![image](https://github.com/Dydx1989/Cluster-HyperHMM/assets/53042175/349a6f76-b7e4-4bb1-910d-ee356f6de2fc)



