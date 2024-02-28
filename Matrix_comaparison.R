## Simple demos of R embedding of HyperHMM


### simple demos of R embedding of HyperHMM

# source code for inference
library(Rcpp)
library(RcppArmadillo)
sourceCpp("hyperhmm-r.cpp")

# source code for plots
library(ggpubr)
source("hypercube-plots.R")

# read in cross-sectional data and return a matrix
cube.read.crosssectional = function(fname) {
  data.raw = readLines(fname)
  data.mat = do.call(rbind, lapply(strsplit(data.raw, ""), as.numeric))
  return(data.mat)
}

# read in longitudinal data and return a list of two matrices
cube.read.longitudinal = function(fname) {
  data.list = list()
  data.raw = read.table(fname, header=FALSE, colClasses = "character")
  data.list$from = do.call(rbind, lapply(strsplit(data.raw[,1], ""), as.numeric))
  data.list$to = do.call(rbind, lapply(strsplit(data.raw[,2], ""), as.numeric))
  return(data.list)
}


### Cross-sectional: Klebsiella pneumonia dataset (CKp)
CKp.mat = cube.read.crosssectional("Data/Cross_sectional_data.txt")
fit.CKp = HyperHMM(CKp.mat)

matrix_cross <- matrix(fit.CKp$stats$mean, nrow = max(fit.CKp$stats$feature), ncol = max(fit.CKp$stats$order), byrow = TRUE)
# Print the matrix
print(matrix_cross)

A <- matrix_cross

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
PAij
### Phylogeny: Klebsiella pneumonia dataset (PKp)
PKp.list = cube.read.longitudinal("Data//phylogeny_data.txt")
fit.PKp = HyperHMM(PKp.list$to, initialstates=PKp.list$from, nboot=100)

matrix_phylo <- matrix(fit.PKp$stats$mean, nrow = max(fit.PKp$stats$feature), ncol = max(fit.PKp$stats$order), byrow = TRUE)
# Print the matrix
print(matrix_phylo)

B <- matrix_phylo

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
PBij

result_matrix=(PBij+PAij)/2
num_rows <- nrow(result_matrix)
num_cols <- ncol(result_matrix)
rownames(result_matrix) <- as.character(1:num_rows)
colnames(result_matrix) <- as.character(1:num_cols)
# Print the result_matrix
print(result_matrix)


library(pheatmap)
library(gridExtra)
library(grid)


## Create the heatmap:
#png("overalheatmap.png",width=1200, height=750)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(result_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
               cex=1,clustering_distance_rows = "manhattan", cex=1,
               clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T)

setHook("grid.newpage", NULL, "replace")
grid.text("Features/Clusters", y=-0.03, gp=gpar(fontsize=16))
#grid.text("ylabel example", x=-0.03, rot=90, gp=gpar(fontsize=16))
grid.text("Order", x=0.98, rot=270, gp=gpar(fontsize=16))



