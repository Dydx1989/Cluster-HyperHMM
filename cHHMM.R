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
library("magrittr")
library(NbClust)
library(phytools)
library(phangorn)
library(ape)
### New packages
library(parameters)
library(easystats)
library(Rcpp)
library(RcppArmadillo)
library(ggpubr)
library(git2r)

# set to FALSE if you want to update HyperHMM from remote repo
use.local.hyperHMM = TRUE
if(use.local.hyperHMM == TRUE) {
  sourceCpp("hyperhmm-r.cpp")
  source("hypercube-plots.R")
} else {
  # pull and load most recent HyperHMM version
  git2r::clone("https://github.com/StochasticBiology/hypercube-hmm", "./tmp/")
  sourceCpp("tmp/hyperhmm-r.cpp")
  source("tmp/hypercube-plots.R")
}

# take data, return named list with data and cluster assignments
#method1: method="Gap"--Gap without Monte Carlo (Bootstrap) samples (Tibshirani et al.,2001)
#method2: method="NbClust"--NbClust method (Charrad et al., 2014)
cHHMM.cluster.features = function(binary_data,method="Gap",
                                  max.clusters = 25) {
  
  # initialise return structure
  r.list = list()
  
  if(is.null(colnames(t(binary_data)))) { rownames(binary_data) = 1:nrow(binary_data) }
  # read data from file
#  AMR_binary_data=read.csv(filename, row.names=1)
  #dim(AMR_binary_data)
  r.list[["data"]] = binary_data

  # cluster with gap statistic method, using maximum useful value of cluster and bootstrap params
  # New method without Monte Carlo (“bootstrap”) samples  Charrad et al., 2014
  
  if(method=="Gap") { 
    
    gap_stat <- n_clusters_gap(binary_data,gap_method = "firstSEmax", n_max = 25)
    
    #first_max_index <- which(diff(sign(diff(gap_stat$Gap))) < 0)+1
    #KK=first_max_index[1]
    
    first_max_index <- which(diff(sign(diff(gap_stat$Gap))) == -2)[1]+1
    
    KK=first_max_index
    
    #print(gap_stat, method = "firstSEmax")
    
    # store gap statistic analysis for interrogation
    r.list[["gap_stat"]] = gap_stat
    
    km.res <- kmeans(binary_data, KK, nstart = max.clusters)
    # Visualize
    
    r.list[["km_res"]] = km.res
    
    #++++++++++++++++++++++++++++++++++++++++++++++
    #+####################################
    # assign features to GENE SUBSET
    type2kclu = data.frame(
      my_data = colnames(t(binary_data)),#substr(colnames(t(binary_data)),1,KK),
      cluster=km.res$cluster)
    
    Cluster_Tab=table(type2kclu)
    generate_column_names <- function(num_columns) {
      # Generate column names like C1, C2, C3, ...
      column_names <- paste0("C", 1:num_columns)
      return(column_names)
    }
    
    colnames(Cluster_Tab) <- generate_column_names(KK)
    #head(Cluster_Tab)
    
    ## Save new data as Gene_cluster: Genes and Cluster 
    #write.csv(Cluster_Tab,"Gene_cluster.csv")
    r.list[["cluster_tab"]] = Cluster_Tab
    return(r.list)
    
  }
  
  
  
  if(method=="NbClust") 
  { Nb_stat<-NbClust(binary_data, min.nc=2, max.nc=15, 
                     method = "kmeans", index = "sdindex")
  
  KK=as.integer(Nb_stat$Best.nc[1])
  
  gap_stat=Nb_stat
  #print(gap_stat, method = "firstSEmax")
  
  # store gap statistic analysis for interrogation
  r.list[["gap_stat"]] = gap_stat
  
  km.res <- kmeans(binary_data, KK, nstart = max.clusters)
  # Visualize
  
  r.list[["km_res"]] = km.res
  
  #++++++++++++++++++++++++++++++++++++++++++++++
  #+####################################
  # assign features to GENE SUBSET
  type2kclu = data.frame(
    my_data =colnames(t(binary_data)), #substr(colnames(t(binary_data)),1,KK),
    cluster=km.res$cluster)
  
  Cluster_Tab=table(type2kclu)
  generate_column_names <- function(num_columns) {
    # Generate column names like C1, C2, C3, ...
    column_names <- paste0("C", 1:num_columns)
    return(column_names)
  }
  
  colnames(Cluster_Tab) <- generate_column_names(KK)
  #head(Cluster_Tab)
  
  ## Save new data as Gene_cluster: Genes and Cluster 
  #write.csv(Cluster_Tab,"Gene_cluster.csv")
  r.list[["cluster_tab"]] = Cluster_Tab
  return(r.list)}
  
  
  
}

# take clustering structure and return set of cross-sectional observations
cHHMM.cross.sectional = function(cluster.structure, occupancy="any") {
  r.list = list()
  ##  Binary Matrix (MM) generated from (AMR_binary_data=Isolate*gene) (N1) and ("Gene_cluster=Gene *Clustes) (N2)
  Final_data= cluster.structure[["data"]]   #read.csv("AMR_binary.csv", row.names=1)
  N1=t(Final_data)
  
  ## Cluster and Genes
  Cluster_gene= cluster.structure[["cluster_tab"]]  #read.csv("Gene_cluster.csv",row.names = 1)
  #head(Cluster_gene[1:8])
  #dim(Cluster_gene)
  Cluster_gene=as.matrix(Cluster_gene)
  #dim(Data_2)
  #head(Data_2[,1:6])
  N2=(Cluster_gene)
  #head(N2[,1:6])
  #dim(N2)
  
  MM=matrix(NA,nrow =nrow(N1),ncol =ncol(N2))
  #dim(MM)
  # loop through rows (isolates) in original data 
  for(i in 1:nrow(N1)) {
    # loop through columns (clusters) in cluster mapping
    for(j in 1:ncol(N2)) {
      # count number of members of this cluster
      num.members = sum(N2[,j])
      # impose occupancy requirements
      if(occupancy == "any") { 
        required = 1
      } else if(occupancy == "majority") {
        required = num.members/2
      } else if(occupancy == "relative") {
        this.cluster = which(N2[,j]==1)
        if(length(this.cluster) > 1) {
          occ.dist = rowSums(N1[,this.cluster])
        } else {
          occ.dist = N1[,this.cluster]
        }
        required = mean(occ.dist)
      }
      # count number of members in this record
      hits = length(which(N1[i,]==1 & N2[,j]==1))
      # assign occupancy accordingly
      MM[i,j] = ifelse(hits >= required, 1, 0)
      # if this isolate's row shares any 1s with this cluster's column, assign a 1
      #MM[i,j]=ifelse(any(N1[i,]&N2[,j]==1),1,0)
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
  #write.csv(Cross_sectional_data,"Cross_sectional_data1.csv")
  # Write data frame to a text file without column and row names
  data_matrix <- as.matrix(Cross_sectional_data)
  r.list[["cross_sectional_data"]] = data_matrix
  #Concatenate the values of each row without gaps
  row_values <- apply(data_matrix, 1, function(x) paste0(x, collapse = ""))
  #write.table(row_values , "Cross_sectional_data1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
  r.list[["bitstrings"]] = row_values
  return(r.list)
}

String_barcode <- function(x) {
  s <- str_split(x, "")[[1]]
  return(paste(s))
}

cHHMM.phylogenetic.assignment = function(cluster.structure, tree, ind.labels, occupancy="any") {
  r.list = list()
  cross.sectional = cHHMM.cross.sectional(cluster.structure, occupancy)
  dset = data.frame(label=ind.labels, cross.sectional$cross_sectional_data)
  ct = curate.tree(tree, dset)
  
  r.list[["full_tree"]] = ct$tree
  # Loop over all nodes in the tree
  finalpairs = data.frame()
  for (i in 1:nrow(ct$srcs)) {
    if (paste0(ct$srcs[i,], collapse="")==paste0(ct$dests[i,], collapse="")) {
      next
    }
    finalpairs = rbind(finalpairs, data.frame(Anc=paste0(ct$srcs[i,], collapse=""), 
                                              Desc=paste0(ct$dests[i,], collapse="")))
  }
 
  sorted_pairs <- finalpairs[order(finalpairs$Anc), ]

  r.list[["sorted_pairs"]] = sorted_pairs
  sp= sorted_pairs
  r.list[["srcs"]] = matrix(as.numeric(unlist(strsplit(sp$Anc, ""))), ncol=str_length(sp$Anc[1]), byrow = TRUE)
  r.list[["dests"]] = matrix(as.numeric(unlist(strsplit(sp$Desc, ""))), ncol=str_length(sp$Anc[1]), byrow= TRUE)
  
  return(r.list)
}
  
# take clustered structure and return estimated phylogeny and transitions from ancestral reconstruction
cHHMM.phylogenetic.estimation = function(cluster.structure, occupancy="any", true.tree = NA) {
  
  r.list = list()
  cross.sectional = cHHMM.cross.sectional(cluster.structure, occupancy)
  full_data = cross.sectional[["cross_sectional_data"]] # IGJ read.csv("Cross_sectional_data.csv",row.names = 1)
  #head(Realife_Kp)
  # Check if all values in the data frame are 0s or 1s
  
  # Function to check if a data frame contains a mixture of 0s and 1s
  is_mixture_of_01 <- function(df) {
    has_zero <- FALSE
    has_one <- FALSE
    
    for (i in 1:nrow(df)) {
      for (j in 1:ncol(df)) {
        if (df[i, j] == 0) {
          has_zero <- TRUE
        } else if (df[i, j] == 1) {
          has_one <- TRUE
        }
        
        # If both 0s and 1s are present, return TRUE
        if (has_zero && has_one) {
          return(TRUE)
        }
      }
    }
    
    # If only 0s or only 1s are present, return FALSE
    return(FALSE)
  }
  
  # Check if full_data contains all 0s or 1s
  if (is_mixture_of_01(t(full_data))) {
    out=pheatmap(t(full_data), show_rownames = TRUE,cluster_cols=T,cluster_rows=F,
                 cex=1,clustering_distance_rows = "manhattan", cex=1,
                 clustering_distance_cols = "manhattan",border_color = TRUE)
    
    
    
    my_tree <- as.phylo(out$tree_col) 
    r.list[["full_tree"]] = my_tree
    #  plot(my_tree)
    # read barcodes
    #barcodesfilename="Cross_sectional_data.csv"
    barcodes = as.data.frame(full_data) # IGJ read.csv(barcodesfilename, header=T)
    # barcodes$ID = 1:nrow(barcodes)
    
    # Remove rows with the same binary values
    # Subset Binary column names
    binary_cols <- colnames(barcodes)
    
    
    barcodes<- barcodes[!duplicated(barcodes[,binary_cols]), ]
    #write.table(barcodes,"barcodes_2.txt")
    # counts = data.frame(Species = barcodes$ID, bitwise = sapply(1, function(x)  
    #    do.call(paste0, as.data.frame(barcodes[-1][, x:(x+7)]))))
    
    pruned.out = pheatmap(t(barcodes), show_rownames = TRUE,cluster_cols=T,cluster_rows=F,
                          cex=1,clustering_distance_rows = "manhattan", cex=1,
                          clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE)
    
    #counts$bitwise
    
    # uncomment to dichotomise, but doesn't seem to improve performance
    # tree = multi2di(tree)
    tree <- as.phylo(pruned.out$tree_col) 
  
    tree$orig.label = tree$tip.label
    tree$tip.label = unlist(apply(barcodes[tree$tip.label,], 1, paste0, collapse=""))
    # Convert tree$tip.label to numeric indices
    #indices <- as.numeric(tree$tip.label)
    
    add_nodelable=c(1:tree$Nnode)
    tree$node.label = gsub("", "_", add_nodelable)
    r.list[["tree_2"]] = tree
    #plot(tree)
    #nodelabels(text=tree$node.label)
    
    tree.labels = c(tree$tip.label, tree$node.label)
    # set all nodes in the tree to negative value of trait by default
    #tree.vals = rep(0, length(tree.labels))
    tree.vals = c(tree$tip.label,rep(-1, length(tree$node.label)))
    cat("\n------- Painting ancestors...\n  ")
    
    change = T
    # while we're still not converged
    while(change == T) {
      change = F
      # loop through all nodes
      for(tree.ref in 1:length(tree.labels)) {
        if(tree.vals[tree.ref] == -1) {
          # if this node is null, check to see if its children are all characterised
          descendant.refs = Children(tree, tree.ref)
          if(all(tree.vals[descendant.refs] != -1)) {
            # if so, make the parent positive
            
            result <- sapply(seq_len(nchar(tree.vals[descendant.refs][[1]])), function(i) {
              ifelse(all(sapply(tree.vals[descendant.refs], function(x) substr(x, i, i) == "1")), "1", "0")
            })
            tree.vals[tree.ref] <- paste(result, collapse = "")
            #ifelse(String_barcode(tree.vals[descendant.refs][1])==String_barcode(tree.vals[descendant.refs][2]), 1,0) 
            
          }
          change = T
        }
      }
    }
    
    #print(tree.vals)
    
    
    #tiplab=tree$tip.label = gsub("_", " ",Bitcount)
    # Find strings in A not in B
    #nodelab<- tree.vals[tiplab != tree.vals]
    #nodelab<- tree.vals[(length(tiplab)+1): length(tree.vals)]
    #tree$tip.label = gsub("_", " ", tiplab)
    #tree$node.label = gsub("_", " ", nodelab)
    #plot(tree)
    #nodelabels(text=tree$node.label)
    
    # Convert the tree to Newick format
    #newick_tree <- write.tree(tree, file = "")
    
    
    # Read the Newick tree string
    #tree_string <- newick_tree
    #tree <- ape::read.tree(text = tree_string)
    
    # Get the tip labels and node labels
    labels <- tree.vals #c(as.character(tree$tip.label), as.character(tree$node.label))
    
    treePairs<-tree$edge
    
    # # Get the bitstrings of each node
    # bitstrings <- labels
    #
    # # Create an empty list to store the ancestor-descendant pairs
    # ancestor_descendant_pairs <- list()
    finalpairs=data.frame()
    
    # Loop over all nodes in the tree
    for (i in 1:nrow(treePairs)) {
      
      if (labels[treePairs[i,1]]==labels[treePairs[i,2]]) {
        next
      }
      finalpairs = rbind(finalpairs, data.frame(Anc=labels[treePairs[i,1]], Desc=labels[treePairs[i,2]]))
    }
    
    # Print the ancestor-descendant pairs with differences
    
    
    # Sort the pairs based on the "Anc" column
    sorted_pairs <- finalpairs[order(finalpairs$Anc), ]
    
    # Print the sorted pairs
    #  for (i in 1:nrow(sorted_pairs)) {
    #    cat(sorted_pairs$Anc[i], "-", sorted_pairs$Desc[i], "\n")
    #  }
    
    
    ## Save the data as phylogeny_data in txt for HyperHMM 
    # set your working directory or find your pc director getwd() and setwd("C:\\file\\")
    # for my PC use setwd("C:\\Users")
    #output_file <- "C:\\Users\\phylogeny_data.txt"
    
    # Write the sorted pairs to the text file
    #write.table(sorted_pairs, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
    r.list[["sorted_pairs"]] = sorted_pairs
    sp= sorted_pairs
    r.list[["srcs"]] = matrix(as.numeric(unlist(strsplit(sp$Anc, ""))), ncol=str_length(sp$Anc[1]), byrow = TRUE)
    r.list[["dests"]] = matrix(as.numeric(unlist(strsplit(sp$Desc, ""))), ncol=str_length(sp$Anc[1]), byrow= TRUE)
    
    return(r.list)
    
    
    
  } else {
    # If not all values are 0s or 1s, use pheatmap
    
    out=heatmap(t(full_data), 
                Rowv = NA,            # Do not cluster rows
                Colv = NA,            # Do not cluster columns
                col = heat.colors(256),  # Color palette
                scale = "none",       # Do not scale the data
                main = "Heatmap of Full Data",
                xlab = "Columns",
                ylab = "Rows")
    # Check if df contains a mixture of 0s and 1s
    if (!is_mixture_of_01(t(full_data))) {
      message("**WARNING!** Data frame does not contain a mixture of 0s and 1s. Hence, genes are present/absent across the isolates")
      message("Returning cross-sectional version to avoid breaking pipeline.")
      finalpairs = data.frame()
      df = cross.sectional$cross_sectional_data
      for(i in 1:nrow(df)) {
      finalpairs = rbind(finalpairs, data.frame(Anc=paste0(rep("0", ncol(df)), collapse=""), 
                                                Desc=paste0(df[i,], collapse="")))
      }
      sorted_pairs <- finalpairs[order(finalpairs$Anc), ]
      r.list[["sorted_pairs"]] = sorted_pairs
      sp= sorted_pairs
      r.list[["srcs"]] = matrix(as.numeric(unlist(strsplit(sp$Anc, ""))), ncol=str_length(sp$Anc[1]), byrow = TRUE)
      r.list[["dests"]] = matrix(as.numeric(unlist(strsplit(sp$Desc, ""))), ncol=str_length(sp$Anc[1]), byrow= TRUE)
      
      return(r.list)
    } else {
      cat("Genes present across the isolates\n")
      # Continue with other functions here
    }
  }
  
  
  
}


############ Maximum Parsimony #####################
####################################################

# take clustered structure and return estimated phylogeny and transitions from ancestral reconstruction
cHHMM.phylogenetic.estimation.parsimony = function(cluster.structure, occupancy="any", true.tree = NA) {
  
  r.list = list()
  cross.sectional = cHHMM.cross.sectional(cluster.structure, occupancy)
  full_data = cross.sectional[["cross_sectional_data"]] # IGJ read.csv("Cross_sectional_data.csv",row.names = 1)
  #head(Realife_Kp)
  # Check if all values in the data frame are 0s or 1s
  
  # Function to check if a data frame contains a mixture of 0s and 1s
  is_mixture_of_01 <- function(df) {
    has_zero <- FALSE
    has_one <- FALSE
    
    for (i in 1:nrow(df)) {
      for (j in 1:ncol(df)) {
        if (df[i, j] == 0) {
          has_zero <- TRUE
        } else if (df[i, j] == 1) {
          has_one <- TRUE
        }
        
        # If both 0s and 1s are present, return TRUE
        if (has_zero && has_one) {
          return(TRUE)
        }
      }
    }
    
    # If only 0s or only 1s are present, return FALSE
    return(FALSE)
  }
  
  # Check if full_data contains all 0s or 1s
  if (is_mixture_of_01(t(full_data))) {
    #out=pheatmap(t(full_data), show_rownames = TRUE,cluster_cols=T,cluster_rows=F,
    #           cex=1,clustering_distance_rows = "manhattan", cex=1,
    #            clustering_distance_cols = "manhattan",border_color = TRUE)
    
    
    
    #my_tree <- as.phylo(out$tree_col) 
    #r.list[["full_tree"]] = my_tree
    #  plot(my_tree)
    # read barcodes
    #barcodesfilename="Cross_sectional_data.csv"
    barcodes = as.data.frame(full_data) # IGJ read.csv(barcodesfilename, header=T)
    # barcodes$ID = 1:nrow(barcodes)
    
    # Remove rows with the same binary values
    # Subset Binary column names
    binary_cols <- colnames(barcodes)
    
    
    barcodes<- barcodes[!duplicated(barcodes[,binary_cols]), ]
    #write.table(barcodes,"barcodes_2.txt")
    # counts = data.frame(Species = barcodes$ID, bitwise = sapply(1, function(x)  
    #    do.call(paste0, as.data.frame(barcodes[-1][, x:(x+7)]))))
    
    ## Trees. Parsimony
    convert_to_named_list <- function(data) {
      # Transpose the data
      transposed_data <- t(data)
      
      # Convert the transposed matrix to a list
      named_list <- as.list(as.data.frame(transposed_data, stringsAsFactors = FALSE))
      
      # Convert each element in the list to a character vector
      named_list <- lapply(named_list, as.character)
      
      return(named_list)
    }
    result <- convert_to_named_list(barcodes)
    
    
    #"""""""""""""""""""""""
    #write.nexus.data(x, file="bitwise_11_C.nex")
    #nex <- read.nexus.data("bitwise_11_C.nex")
    #head(nex)
    # compute (Conversion among Sequence Formats)
    b.pd <- phyDat(result, type="USER", levels=c(0:1)) # convert nex to phyDat format
    b.pd #check what went in....
    ## 26 sequences with 8 character and 8 different site patterns.
    ## The states are 0 1
    ###### Parsimony tree
    # set.seed(1)
    #  C_tree <- random.addition(b.pd )
    #    primates_opt <- optim.parsimony(C_tree, b.pd)
    #    primates_opt <- acctran(primates_opt,  b.pd)
    treeRatchet  <- pratchet(b.pd, trace = 0, minit=100)
    
    #plot(primates_opt)
    ## Rooted  phylogenetic trees BKI
    tree <- root(treeRatchet, 1, r = TRUE)
    
    tree$orig.label = tree$tip.label
    tree$tip.label = unlist(apply(barcodes[tree$tip.label,], 1, paste0, collapse=""))
    # Convert tree$tip.label to numeric indices
    #indices <- as.numeric(tree$tip.label)
    
    add_nodelable=c(1:tree$Nnode)
    tree$node.label = gsub("", "_", add_nodelable)
    r.list[["tree_2"]] = tree
    #plot(tree)
    #nodelabels(text=tree$node.label)
    
    tree.labels = c(tree$tip.label, tree$node.label)
    # set all nodes in the tree to negative value of trait by default
    #tree.vals = rep(0, length(tree.labels))
    tree.vals = c(tree$tip.label,rep(-1, length(tree$node.label)))
    cat("\n------- Painting ancestors...\n  ")
    
    change = T
    # while we're still not converged
    while(change == T) {
      change = F
      # loop through all nodes
      for(tree.ref in 1:length(tree.labels)) {
        if(tree.vals[tree.ref] == -1) {
          # if this node is null, check to see if its children are all characterised
          descendant.refs = Children(tree, tree.ref)
          if(all(tree.vals[descendant.refs] != -1)) {
            # if so, make the parent positive
            
            result <- sapply(seq_len(nchar(tree.vals[descendant.refs][[1]])), function(i) {
              ifelse(all(sapply(tree.vals[descendant.refs], function(x) substr(x, i, i) == "1")), "1", "0")
            })
            tree.vals[tree.ref] <- paste(result, collapse = "")
            #ifelse(String_barcode(tree.vals[descendant.refs][1])==String_barcode(tree.vals[descendant.refs][2]), 1,0) 
            
          }
          change = T
        }
      }
    }
    
    #print(tree.vals)
    
    
    #tiplab=tree$tip.label = gsub("_", " ",Bitcount)
    # Find strings in A not in B
    #nodelab<- tree.vals[tiplab != tree.vals]
    #nodelab<- tree.vals[(length(tiplab)+1): length(tree.vals)]
    #tree$tip.label = gsub("_", " ", tiplab)
    #tree$node.label = gsub("_", " ", nodelab)
    #plot(tree)
    #nodelabels(text=tree$node.label)
    
    # Convert the tree to Newick format
    #newick_tree <- write.tree(tree, file = "")
    
    
    # Read the Newick tree string
    #tree_string <- newick_tree
    #tree <- ape::read.tree(text = tree_string)
    
    # Get the tip labels and node labels
    labels <- tree.vals #c(as.character(tree$tip.label), as.character(tree$node.label))
    
    treePairs<-tree$edge
    
    # # Get the bitstrings of each node
    # bitstrings <- labels
    #
    # # Create an empty list to store the ancestor-descendant pairs
    # ancestor_descendant_pairs <- list()
    finalpairs=data.frame()
    
    # Loop over all nodes in the tree
    for (i in 1:nrow(treePairs)) {
      
      if (labels[treePairs[i,1]]==labels[treePairs[i,2]]) {
        next
      }
      finalpairs = rbind(finalpairs, data.frame(Anc=labels[treePairs[i,1]], Desc=labels[treePairs[i,2]]))
    }
    
    # Print the ancestor-descendant pairs with differences
    
    
    # Sort the pairs based on the "Anc" column
    sorted_pairs <- finalpairs[order(finalpairs$Anc), ]
    
    # Print the sorted pairs
    #  for (i in 1:nrow(sorted_pairs)) {
    #    cat(sorted_pairs$Anc[i], "-", sorted_pairs$Desc[i], "\n")
    #  }
    
    
    ## Save the data as phylogeny_data in txt for HyperHMM 
    # set your working directory or find your pc director getwd() and setwd("C:\\file\\")
    # for my PC use setwd("C:\\Users")
    #output_file <- "C:\\Users\\phylogeny_data.txt"
    
    # Write the sorted pairs to the text file
    #write.table(sorted_pairs, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
    r.list[["sorted_pairs"]] = sorted_pairs
    sp= sorted_pairs
    r.list[["srcs"]] = matrix(as.numeric(unlist(strsplit(sp$Anc, ""))), ncol=str_length(sp$Anc[1]), byrow = TRUE)
    r.list[["dests"]] = matrix(as.numeric(unlist(strsplit(sp$Desc, ""))), ncol=str_length(sp$Anc[1]), byrow= TRUE)
    
    return(r.list)
    
    
    
  } else {
    # If not all values are 0s or 1s, use pheatmap
    
    out=heatmap(t(full_data), 
                Rowv = NA,            # Do not cluster rows
                Colv = NA,            # Do not cluster columns
                col = heat.colors(256),  # Color palette
                scale = "none",       # Do not scale the data
                main = "Heatmap of Full Data",
                xlab = "Columns",
                ylab = "Rows")
    # Check if df contains a mixture of 0s and 1s
    if (!is_mixture_of_01(t(full_data))) {
      message("**WARNING!** Data frame does not contain a mixture of 0s and 1s. Hence, genes are present/absent across the isolates")
      message("Returning cross-sectional version to avoid breaking pipeline.")
      finalpairs = data.frame()
      df = cross.sectional$cross_sectional_data
      for(i in 1:nrow(df)) {
        finalpairs = rbind(finalpairs, data.frame(Anc=paste0(rep("0", ncol(df)), collapse=""), 
                                                  Desc=paste0(df[i,], collapse="")))
      }
      sorted_pairs <- finalpairs[order(finalpairs$Anc), ]
      r.list[["sorted_pairs"]] = sorted_pairs
      sp= sorted_pairs
      r.list[["srcs"]] = matrix(as.numeric(unlist(strsplit(sp$Anc, ""))), ncol=str_length(sp$Anc[1]), byrow = TRUE)
      r.list[["dests"]] = matrix(as.numeric(unlist(strsplit(sp$Desc, ""))), ncol=str_length(sp$Anc[1]), byrow= TRUE)
      
      return(r.list)
    } else {
      cat("Genes present across the isolates\n")
      # Continue with other functions here
    }
  }
  
  
  
}


curate.tree = function(tree.src, data.src, losses = FALSE, data.header=TRUE) {
  if(is.character(tree.src)) {
    # read in Newick tree and root
    my.tree = read.tree(tree.src)
    my.rooted.tree = root(my.tree, 1, resolve.root = TRUE)
  } else {
    my.rooted.tree = tree.src
  }
  
  if(is.character(data.src)) {
    # read in barcode data
    my.data = read.csv(data.src, header=data.header)
    colnames(my.data)[1] = "label"
  } else {
    my.data = data.src
    colnames(my.data)[1] = "label"
  }
  
  match.set = match(my.data$label, my.rooted.tree$tip.label)
  if(any(is.na(match.set))) {
    message("Found observations that didn't correspond to tips of the tree!")
    my.data = my.data[-which(is.na(match.set)),]
    match.set = match(my.data$label, my.rooted.tree$tip.label)
  }
  
  if(any(duplicated(my.data$label))) {
    message("Duplicates in observation set!")
  }
  
  # prune tree to include only those tips in the barcode dataset
  tree = drop.tip(my.rooted.tree,
                  my.rooted.tree$tip.label[-match.set])
  
  tree$node.label = as.character(length(tree$tip.label) + 1:tree$Nnode)
  tree.labels = c(tree$tip.label, tree$node.label)
  
  if(any(duplicated(tree$tip.label))) {
    message("Duplicates in tree tips!")
  }
  
  cat("\n------- Painting ancestors...\n  ")
  
  # initialise "recursive" algorithm
  change = T
  new.row = my.data[1,]
  changes = data.frame()
  
  # while we're still making changes
  while(change == T) {
    change = F
    # loop through all nodes
    for(tree.ref in 1:length(tree.labels)) {
      #print(tree.ref)
      this.label = tree.labels[tree.ref]
      # see if this node exists in our barcode dataset
      if(!(this.label %in% my.data$label)) {
        # if not, check to see if its children are all characterised
        descendant.refs = Children(tree, tree.ref)
        if(all(tree.labels[descendant.refs] %in% my.data$label)) {
          #print(paste0("Doing ", this.label))
          ## ancestral state reconstruction
          # pull the rows in our barcode dataset corresponding to children of this node
          descendant.rows = which(my.data$label %in% tree.labels[descendant.refs])
          if(losses == FALSE) {
            # bitwise AND to construct the ancestral state
            new.barcode = apply(my.data[descendant.rows,2:ncol(my.data)], 2, prod)
          } else {
            # bitwise OR to construct the ancestral state
            new.barcode = apply(my.data[descendant.rows,2:ncol(my.data)], 2, max)
          }
          # add this ancestral state to the barcode dataset
          new.data = new.row
          new.data$label[1] = this.label
          new.data[1,2:ncol(new.data)] = new.barcode
          my.data = rbind(my.data, new.data)
          
          ## adding transitions to our observation set
          # loop through children
          for(d.ref in descendant.refs) {
            # pull the barcodes and branch lengths together and add to data frame
            d.row = which(my.data$label == tree.labels[d.ref])
            # get the reference for this edge
            e.ref = which(tree$edge[,2] == d.ref)
            changes = rbind(changes, 
                            data.frame(from=paste0(new.data[2:ncol(new.data)], collapse=""),
                                       to=paste0(my.data[d.row,2:ncol(new.data)], collapse=""),
                                       time=tree$edge.length[e.ref],
                                       from.node=this.label,
                                       to.node=tree.labels[d.ref]))
          }
          # we made a change, so keep the loop going
          change = T
        }
      }
    }
  }
  srcs = matrix(as.numeric(unlist(lapply(changes$from, strsplit, split=""))), byrow=TRUE, ncol=ncol(new.data)-1)
  dests = matrix(as.numeric(unlist(lapply(changes$to, strsplit, split=""))), byrow=TRUE, ncol=ncol(new.data)-1)
  
  rL = list("tree" = tree,
            "data" = my.data,
            "transitions" = changes,
            "srcs" = srcs,
            "dests" = dests,
            "times" = changes$time)
  return(rL)
}

# produce comparison matrix for different approaches
cHHMM.matrix.comparison = function(fitted.obj.1, fitted.obj.2) {
  
  r.list = list()
  
  stats.1 = fitted.obj.1$stats
  stats.1 = stats.1[order(stats.1$feature, stats.1$order),]
  A <- matrix(stats.1$mean, 
              nrow = max(stats.1$feature), 
              ncol = max(stats.1$order), byrow = TRUE)
  # Print the matrix
  #print(matrix_cross)
  
  PAij = matrix(0,nrow =nrow(A),ncol =ncol(A))
  for(i in 1:nrow(A)) {
    for(j in 1:nrow(A)) {
      for(ti in 1:nrow(A)) {
        for(tj in 1:nrow(A)) {
          elem = A[i,ti]*A[j,tj] #/(1-A[j,ti])
          if(ti < tj & i != j & !is.na(elem)) {
            #print(paste(c(i, "@", ti, "&", j, "@", tj, elem), collapse=" "))
            PAij[i,j] = PAij[i,j] + elem
          }
        }
      }
    }
  }
  
  num_rows <- nrow(PAij)
  num_cols <- ncol(PAij)
  rownames(PAij) <- as.character(1:num_rows)
  colnames(PAij) <- as.character(1:num_cols)
  #PAij
  ### Phylogeny: Klebsiella pneumonia dataset (PKp)
  
  stats.2 = fitted.obj.2$stats
  stats.2 = stats.2[order(stats.2$feature, stats.2$order),]
  B <- matrix(stats.2$mean, 
              nrow = max(stats.2$feature), 
              ncol = max(stats.2$order), byrow = TRUE)
  # Print the matrix
  #print(matrix_phylo)
  
  PBij =matrix(0,nrow =nrow(B),ncol =ncol(B))
  for(i in 1:nrow(B)) {
    for(j in 1:nrow(B)) {
      for(ti in 1:nrow(B)) {
        for(tj in 1:nrow(B)) {
          elem = B[i,ti]*B[j,tj] #/(1-A[j,ti])
          if(ti < tj & i != j & !is.na(elem)) {
            #print(paste(c(i, "@", ti, "&", j, "@", tj, elem), collapse=" "))
            PBij[i,j] = PBij[i,j] + elem
          }
        }
      }
    }
  }
  
  num_rows <- nrow(PBij)
  num_cols <- ncol(PBij)
  rownames(PBij) <- as.character(1:num_rows)
  colnames(PBij) <- as.character(1:num_cols)
  #PBij
  
  
  result_matrix=(PBij+PAij)/2
  num_rows <- nrow(result_matrix)
  num_cols <- ncol(result_matrix)
  rownames(result_matrix) <- as.character(1:num_rows)
  colnames(result_matrix) <- as.character(1:num_cols)
  # Print the result_matrix
  r.list[["results_matrix"]] = result_matrix
  
  return(r.list)
}






