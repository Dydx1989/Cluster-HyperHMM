library(pheatmap) ## for heatmap generation
#library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing 
library("cluster")
library("factoextra")
library("magrittr")
library(cluster)  # we'll use these packages
library(fpc)
library(mclust)
library(phytools)
library(phangorn)
library(ape)
## Read in cross-sectional data that contains isolate and clsuter
###########################################

Realife_Kp=read.csv("Cross_sectional_data.csv",row.names = 1)
#head(Realife_Kp)
out = pheatmap(t(Realife_Kp), show_rownames = TRUE,cluster_cols=T,cluster_rows=F,
               cex=1,clustering_distance_rows = "manhattan", cex=1,
               clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE)

my_tree <- as.phylo(out$tree_col) 
plot(my_tree)
# read barcodes
barcodesfilename="Cross_sectional_data.csv"
barcodes = read.csv(barcodesfilename, header=T)
names(barcodes)[names(barcodes) == "X"]<- "Species"
# Remove rows with the same binary values
# Subset Binary column names
#binary_cols <- colnames(barcodes[-1])
#barcodes<- barcodes[!duplicated(barcodes[,binary_cols]), ]
#write.table(barcodes,"barcodes_2.txt")
counts = data.frame(Species = barcodes$Species, bitwise = sapply(1, function(x)  
  do.call(paste0, as.data.frame(barcodes[-1][, x:(x+7)]))))

#counts$bitwise



# uncomment to dichotomise, but doesn't seem to improve performance
# tree = multi2di(tree)
tree <- as.phylo(out$tree_col) 
# Convert tree$tip.label to numeric indices
indices <- as.numeric(tree$tip.label)

Bitcount<- counts$bitwise[indices]

tree$tip.label = gsub("_", " ", Bitcount)
add_nodelable=c(1:tree$Nnode)
tree$node.label = gsub("_", " ", add_nodelable)
plot(tree)
nodelabels(text=tree$node.label)

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
        #print(paste(c("found node ", tree.ref, " which has children ", descendant.refs), collapse=""))      
        String_barcode <- function(x) {
          s <- str_split(x, "")[[1]]
          paste(s)
        }
        
        #parent_value <-paste(as.character(ifelse(String_barcode(tree.vals[descendant.refs][1])==1&String_barcode(tree.vals[descendant.refs][2])==1, 1,0) ), collapse="")
        
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


tiplab=tree$tip.label = gsub("_", " ",Bitcount)
# Find strings in A not in B
#nodelab<- tree.vals[tiplab != tree.vals]
nodelab<- tree.vals[(length(tiplab)+1): length(tree.vals)]
tree$tip.label = gsub("_", " ", tiplab)
tree$node.label = gsub("_", " ", nodelab)
#plot(tree)
nodelabels(text=tree$node.label)

# Convert the tree to Newick format
newick_tree <- write.tree(tree, file = "")


# Read the Newick tree string
tree_string <- newick_tree
tree <- ape::read.tree(text = tree_string)

# Get the tip labels and node labels
labels <- c(tree$tip.label, tree$node.label)

treePairs<-tree$edge

# # Get the bitstrings of each node
# bitstrings <- labels
#
# # Create an empty list to store the ancestor-descendant pairs
# ancestor_descendant_pairs <- list()
finalpairs=data.frame(Anc=character(),Desc=character())

# Loop over all nodes in the tree
for (i in 1:nrow(treePairs)) {
  
  if (labels[treePairs[i,1]]==labels[treePairs[i,2]]) {
    next
  }
  finalpairs[nrow(finalpairs)+1,]=c(labels[treePairs[i,1]],labels[treePairs[i,2]])
}

# Print the ancestor-descendant pairs with differences


# Sort the pairs based on the "Anc" column
sorted_pairs <- finalpairs[order(finalpairs$Anc), ]

# Print the sorted pairs
for (i in 1:nrow(sorted_pairs)) {
  cat(sorted_pairs$Anc[i], "-", sorted_pairs$Desc[i], "\n")
}


## Save the data as phylogeny_data in txt for HyperHMM 
# set your working directory or find your pc director getwd() and setwd("C:\\file\\")
# for my PC use setwd("C:\\Users")
output_file <- "C:\\Users\\phylogeny_data.txt"

# Write the sorted pairs to the text file
write.table(sorted_pairs, file = output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

