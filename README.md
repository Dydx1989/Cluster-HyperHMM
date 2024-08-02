# Cluster-HyperHMM

Cluster large-scale data and perform evolutionary inference using HyperHMM https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098 (code https://github.com/StochasticBiology/hypercube-hmm ).

![image](https://github.com/user-attachments/assets/1649970a-5b9e-42e8-8286-222d982cac92)

Overview
---

Accumulation models infer the dynamics by which features are acquired by a system over time. One flexible and efficient approach for this is HyperHMM. But it is limited in the number of features it can consider, with >20 or so features being intractable. The idea of Cluster-HyperHMM is to cluster datasets with large numbers of features, to reduce dimensionality to a smaller set of "effective features" or "feature clusters". HyperHMM can then be used to learn the dynamics by which the system acquires elements of these effective features.

Several points must be considered:
* Different clustering methods can be used to assign clusters based on which features appear in similar patterns across observations.
* Once clusters are assigned (i.e. which features correspond to which cluster), an occupancy protocol must be used. For example, if an observation has some but not all features of cluster X, does this count as having acquired cluster X or not?
* The relationship, if any, between observations must be considered. If they are fully independent, this is not an issue. But if they are, for example, linked by a shared ancestry, the possibility of similarity-by-descent should be addressed.

As such, Cluster-HyperHMM allows:
* A choice of clustering method: *gap statistic* or *k-means*.
* A choice of cluster occupancy rules for how many member features of a cluster are required for it to count as acquired: *at least one member*, *a majority of members*, *more than average number of members*
* A choice of accounting for relationships between observations: *all independent*, *estimated phylogeny based on clustering*, *independent phylogeny provided by user*

Pipeline outline
---

A pipeline will generally involve acquiring and pre-processing data, clustering, HyperHMM inference on the clustered data, then output and interpretation. All except pre-processing are case-general and are performed in `cHHMM-expts_Kleborate.R`. 

Several case studies are included:

* _Kleborate_. Data on *Klebsiella pneumoniae* genomes from around the world. `Download_and_curate-kleborate_data.R` downloads and curates these data. Run `cHHMM-expts_Kleborate.R` with `expt = "Kleborate"` (default) for this.
* _Malaria_. Data on severe malaria progression in patients. Data is already included in `jallow_dataset_binary_with2s.csv`. Run `cHHMM-expts_Kleborate.R` with `expt = "malaria"` for this.
* _Synthetic_. A simple, synthetically constructed dataset. Run `cHHMM-expts_Kleborate.R` with `expt = "synthetic"` for this.
* _Organelle DNA_. Data on evolutionary loss of genes from organelle DNA (mitochondrial and chloroplast DNA) -- with independently estimated phylogenies. Run `cHHMM-expts_Kleborate.R` with `expt = "mtDNA"` or `expt = "cpDNA"` for this.

The outputs will include illustrations of the cluster selection and identities, summaries of the inferred dynamics under different evolutionary assumptions, and comparisons across assumptions and clustering protocols.

Pipeline functions
---
| Function    | Description | Key arguments (* default)|
| -------- | ------- |------|
| `cHHMM.cluster.features`  | Takes binary data and returns an object describing the clustering that has been assigned    | `method` = `Gap`* (gap statistic) or `NbClust` (k-means) |
| `cHHMM.cross.sectional` | Takes a clustered object, interprets the data as cross-sectional (independent observtions) and returns a set of observations for use with HyperHMM    | `occupancy` = `any`* (at least one member), `majority` (majority of members), `relative` (above average member count) |
| `cHHMM.phylogenetic.estimation`    | Takes a clustered object, estimates a parsimonious phylogeny linking observations, and returns a set of observation transitions for use with HyperHMM    |`occupancy` (as above) |
| `cHHMM.phylogenetic.assignment`    | Takes a clustered object and a given phylogeny, and returns a set of observation transitions for use with HyperHMM    | `occupancy` (as above) |
| `HyperHMM` | Performs HyperHMM inference (see refs) | (see documentation) |
| `plot.cHHMM` | Takes a dataset and a HyperHMM fit and produces summary visualisations of the output | |
| `cHHMM.matrix.comparison` | Takes two HyperHMM fits (from different protocol choices) and returns a matrix summarising the similarities and differences in the inferred dynamics | |

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



