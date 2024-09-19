# Cluster-HyperHMM

Cluster large-scale data and perform evolutionary inference using HyperHMM https://academic.oup.com/bioinformatics/article/39/1/btac803/6895098 (code https://github.com/StochasticBiology/hypercube-hmm ).

![image](https://github.com/user-attachments/assets/f6fb08e1-d0c5-422d-81e6-20e0d706ac06)


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

A pipeline will generally involve acquiring and pre-processing data, clustering, HyperHMM inference on the clustered data, then output and interpretation. All except pre-processing are case-general and are performed in `cHHMM-expts.R`. 

Several case studies are included in `cHHMM-expts.R`:

* _Kleborate_. Data on *Klebsiella pneumoniae* genomes from around the world. `Download_and_curate-kleborate_data.R` downloads and curates these data. Run `cHHMM-expts.R` with `expt = "Kleborate"` for this.
* _Malaria_. Data on severe malaria progression in patients. Data is already included in `jallow_dataset_binary_with2s.csv`. Run `cHHMM-expts.R` with `expt = "malaria"` for this.
* _Synthetic_. A simple, synthetically constructed dataset. Run `cHHMM-expts.R` with `expt = "synthetic"` for this.
* _Organelle DNA_. Data on evolutionary loss of genes from organelle DNA (mitochondrial and chloroplast DNA) -- with independently estimated phylogenies. Run `cHHMM-expts.R` with `expt = "mtDNA"` or `expt = "cpDNA"` for this.

The outputs will include illustrations of the cluster selection and identities, summaries of the inferred dynamics under different evolutionary assumptions, and comparisons across assumptions and clustering protocols.

Pipeline functions
---
See the above figure for each key function's position in the workflow. These are all either in `cHHMM.R` or in the original HyperHMM code.

| Function    | Description | Key arguments (* default)|
| -------- | ------- |------|
| `cHHMM.cluster.features`  | Takes binary data and returns an object describing the clustering that has been assigned    | `method` = `Gap`* (gap statistic) or `NbClust` (k-means) |
| `cHHMM.cross.sectional` | Takes a clustered object, interprets the data as cross-sectional (independent observtions) and returns a set of observations for use with HyperHMM    | `occupancy` = `any`* (at least one member), `majority` (majority of members), `relative` (above average member count) |
| `cHHMM.phylogenetic.estimation`    | Takes a clustered object, estimates a parsimonious phylogeny linking observations, and returns a set of observation transitions for use with HyperHMM    |`occupancy` (as above) |
| `cHHMM.phylogenetic.assignment`    | Takes a clustered object and a given phylogeny, and returns a set of observation transitions for use with HyperHMM    | `occupancy` (as above) |
| `HyperHMM` | Performs HyperHMM inference (see refs) | (see documentation) |
| `plot.cHHMM` | Takes a dataset and a HyperHMM fit and produces summary visualisations of the output | |
| `cHHMM.matrix.comparison` | Takes two HyperHMM fits (from different protocol choices) and returns a matrix summarising the similarities and differences in the inferred dynamics | |
| `.Rdata file` |Consists of results generated through the following functions: `cHHMM.cross.sectional`, `cHHMM.cluster.features`, and `cHHMM.phylogenetic.estimation` for Klebsiella data, providing quick access and demonstration. | E.g. `load("cross.sectional.obs.RData");load("fit.cross.sectional.RData"); plot.cHHMM(cross.sectional.obs$cross_sectional_data, fit.cross.sectional, label="Independent\nMajority occupancy")` 
 |

