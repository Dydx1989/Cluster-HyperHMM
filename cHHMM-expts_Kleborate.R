#Gap method without Monte Carlo (“bootstrap”) samples (B)
source("cHHMM.R")

# we must have source data with observations as rows and features as columns
# there was still specific code assuming 8 features (e.g. x+7) -- I think I have removed this
# applying the code to the other cases below throws errors
# we must allow flexibility in cluster assignment e.g. a single member feature gives me a 1 vs a majority of member features vs ...?
# with different seeds the clustering gives dramatically different results!

# this AMR cade study should reproduce, but the other examples below throw errors
# NB with different seeds the AMR case sometimes throws an error too
#Error in if (mm < k) stop("more cluster centers than distinct data points.") : 
#  missing value where TRUE/FALSE needed

#expt = "Kleborate"
expt = "malaria"
#expt = "PATHOGEN"
sf = 2


if(expt == "AMR") {
  other.data = t(read.csv("K.p.binary.csv"))
  other.data = other.data[2:nrow(other.data),1:ncol(other.data)]
  class(other.data) = "numeric"
  other.data[other.data==2] = 1
  src.data = as.data.frame(other.data)
  # Find rows that are not all zeros or all ones
  non_zero_one_rows <- which(!(rowSums(src.data == 0) == ncol(src.data) | rowSums(src.data == 1) == ncol(src.data)))
  
  # Find columns that are not all zeros or all ones
  non_zero_one_columns <- which(!(colSums(src.data == 0) == nrow(src.data) | colSums(src.data == 1) == nrow(src.data)))
  
  # Retain only rows and columns that are not all zeros or all ones
  src.data<- src.data[non_zero_one_rows, non_zero_one_columns]
  
} else if(expt == "synthetic") {
  set.seed(1)
  nr = 180
  nc = 180
  synth.data = matrix(0, nrow=nr, ncol=nc)
  for(i in 1:nrow(synth.data)) {
    if(i < nr/3) {
      synth.data[i,] = c(rbinom(nc/3, 1, 0.9), rbinom(nc/3, 1, 0.1), rbinom(nc/3, 1, 0.5))
    } else if(i < 2*nr/3) {
      synth.data[i,] = c(rbinom(nc/3, 1, 0.1), rbinom(nc/3, 1, 0.5), rbinom(nc/3, 1, 0.9))
    } else {
      synth.data[i,] = c(rbinom(nc/3, 1, 0.5), rbinom(nc/3, 1, 0.9), rbinom(nc/3, 1, 0.1))
    }
  }
  src.data = as.data.frame(t(synth.data))
} else if(expt == "malaria") {
  other.data = t(read.csv("jallow_dataset_binary_with2s.csv"))
  other.data = other.data[2:nrow(other.data),1:200]
  
  class(other.data) = "numeric"
  other.data[other.data==2] = 1
  src.data = as.data.frame(other.data)
} else if(expt == "PATHOGEN") {
  #Source: Germany 
  other.data = t(read.csv("pathogenwatch.csv"))
  other.data = other.data[2:nrow(other.data),1:ncol(other.data)]
  
  class(other.data) = "numeric"
  other.data[other.data==2] = 1
  src.data = as.data.frame(other.data)
  # Find rows that are not all zeros or all ones
  non_zero_one_rows <- which(!(rowSums(src.data == 0) == ncol(src.data) | rowSums(src.data == 1) == ncol(src.data)))
  
  # Find columns that are not all zeros or all ones
  non_zero_one_columns <- which(!(colSums(src.data == 0) == nrow(src.data) | colSums(src.data == 1) == nrow(src.data)))
  
  # Retain only rows and columns that are not all zeros or all ones
  src.data<- src.data[non_zero_one_rows, non_zero_one_columns]
  
}else if(expt == "Kleborate") {
  #Source: Germany 
  other.data = t(read.csv("kleborate-wide-binary.csv"))
  other.data = other.data[2:nrow(other.data),1:ncol(other.data)]
  
  class(other.data) = "numeric"
  other.data[other.data==2] = 1
  src.data = as.data.frame(other.data)
  # Find rows that are not all zeros or all ones
  non_zero_one_rows <- which(!(rowSums(src.data == 0) == ncol(src.data) | rowSums(src.data == 1) == ncol(src.data)))
  
  # Find columns that are not all zeros or all ones
  non_zero_one_columns <- which(!(colSums(src.data == 0) == nrow(src.data) | colSums(src.data == 1) == nrow(src.data)))
  
  # Retain only rows and columns that are not all zeros or all ones
  src.data<- src.data[non_zero_one_rows, non_zero_one_columns]
  
}











## Use plot for Gap method without Monte Carlo (“bootstrap”) samples (B)
#plot((clustered.structure[["gap_stat"]]))

png(paste0(expt, "-data.png"), width=600*sf, height=400*sf, res=72*sf)
print(pheatmap(src.data))
dev.off()

clustered.structure = cHHMM.cluster.features(src.data,method="Gap")

first_max_index <- which(diff(sign(diff(clustered.structure[["gap_stat"]]$Gap))) == -2)[1]+1
# Plot the gap statistic using ggplot
png(paste0(expt, "-gap-stat.png"), width=400*sf, height=300*sf, res=72*sf)
print(
  ggplot(clustered.structure[["gap_stat"]], aes(x = 1:length(Gap), y = Gap)) +
    geom_line(color = "blue") + 
    geom_point(color = "blue") +
    geom_vline(xintercept = first_max_index, color = "red", linetype = "dashed") +
    labs(x = "Number of Clusters", y = "Gap Statistic", title = "Gap Statistic vs. Number of Clusters") 
)
dev.off()

# plot cluster IDs
# convert any cell greater than 1 to 1
clustered.structure$cluster_tab[clustered.structure$cluster_tab > 1] <- 1

ids = clustered.structure$cluster_tab
df.plot = data.frame()
levels = rep(0, ncol(ids))
for(i in 1:nrow(ids)) {
  cid = which(ids[i,] == 1)
  df.plot = rbind(df.plot, data.frame(cluster=cid, level=levels[cid], name=rownames(ids)[i]))
  levels[cid] = levels[cid]+1
}
g.cids = ggplot(df.plot, aes(x=cluster, y=level, label=name)) + geom_text() +
  scale_x_continuous(breaks=1:ncol(ids)) + scale_y_continuous(breaks=NULL) +
  labs(x = "Cluster ID", y = "") + theme_minimal()
png(paste0(expt, "-cluster-IDs.png"), width=600*sf, height=300*sf, res=72*sf)
print(g.cids)
dev.off()

####
png(paste0(expt, "-clusters.png"), width=600*sf, height=400*sf, res=72*sf)
print(
  fviz_cluster(clustered.structure[["km_res"]], data = clustered.structure[["data"]],
               ellipse.type = "convex",
               palette = "viridis",
               ggtheme = theme_minimal())
)
dev.off()

cross.sectional.obs = cHHMM.cross.sectional(clustered.structure)
phylogenetic.obs = cHHMM.phylogenetic.estimation(clustered.structure)

cross.sectional.obs.majority = cHHMM.cross.sectional(clustered.structure, occupancy="majority")
phylogenetic.obs.majority = cHHMM.phylogenetic.estimation(clustered.structure, occupancy="majority")

cross.sectional.obs.relative = cHHMM.cross.sectional(clustered.structure, occupancy="relative")
phylogenetic.obs.relative = cHHMM.phylogenetic.estimation(clustered.structure, occupancy="relative")

fit.cross.sectional = HyperHMM(cross.sectional.obs$cross_sectional_data, nboot = 2)
fit.phylogenetic = HyperHMM(phylogenetic.obs$dests, initialstates = phylogenetic.obs$srcs, nboot = 2)
fit.cross.sectional.majority = HyperHMM(cross.sectional.obs.majority$cross_sectional_data, nboot = 2)
fit.phylogenetic.majority = HyperHMM(phylogenetic.obs.majority$dests, initialstates = phylogenetic.obs.majority$srcs, nboot = 2)
fit.cross.sectional.relative = HyperHMM(cross.sectional.obs.relative$cross_sectional_data, nboot = 2)
fit.phylogenetic.relative = HyperHMM(phylogenetic.obs.relative$dests, initialstates = phylogenetic.obs.relative$srcs, nboot = 2)

g.fluxes =  ggarrange(
  plot.hypercube.flux(fit.cross.sectional, thresh = 0.02),
  plot.hypercube.flux(fit.phylogenetic, thresh=0.02),
  labels=c("A. CS", "B. Phy"))
png(paste0(expt, "-out-fluxes.png"), width=600*sf, height=400*sf, res=72*sf)
print(g.fluxes)
dev.off()

plot.cHHMM(cross.sectional.obs$cross_sectional_data, fit.cross.sectional, label="Cross-sectional\nAny occupancy")
pheatmap(cross.sectional.obs$cross_sectional_data, cluster_cols = FALSE, cluster_rows = FALSE, color = c("white", "grey"), legend=FALSE)

#g.all = ggarrange( plot.standard(fit.cross.sectional, label="Cross-sectional\nAny occupancy"),
#           plot.standard(fit.phylogenetic, label="Phylogenetic\nAny occupancy"), 
#           plot.standard(fit.cross.sectional.majority, label="Cross-sectional\nMajority occupancy"),
#           plot.standard(fit.phylogenetic.majority, label="Phylogenetic\nMajority occupancy"),
#           plot.standard(fit.cross.sectional.relative, label="Cross-sectional\nRelative occupancy"),
#           plot.standard(fit.phylogenetic.relative, label="Phylogenetic\nRelative occupancy"),
    g.all = ggarrange( plot.cHHMM(cross.sectional.obs$cross_sectional_data, fit.cross.sectional, label="Cross-sectional\nAny occupancy"),
                       plot.cHHMM(phylogenetic.obs$dests, fit.phylogenetic, label="Phylogenetic\nAny occupancy"),
                       plot.cHHMM(cross.sectional.obs.majority$cross_sectional_data, fit.cross.sectional.majority, label="Cross-sectional\nMajority occupancy"),
                       plot.cHHMM(phylogenetic.obs.majority$dests, fit.phylogenetic.majority, label="Phylogenetic\nMajority occupancy"),
                       plot.cHHMM(cross.sectional.obs.relative$cross_sectional_data, fit.cross.sectional.relative, label="Cross-sectional\nRelative occupancy"),
                       plot.cHHMM(phylogenetic.obs.relative$dests, fit.phylogenetic.relative, label="Phylogenetic\nRelative occupancy"),
           nrow=3, ncol=2,
           labels=c("A.", "B.", "C.", "D.", "E.", "F."))
png(paste0(expt, "-out-all.png"), width=1600*sf, height=1200*sf, res=72*sf)
print(g.all)
dev.off()

comp.any.cs.phy = cHHMM.matrix.comparison(fit.cross.sectional, fit.phylogenetic)
comp.cs = cHHMM.matrix.comparison(fit.cross.sectional, fit.cross.sectional.majority)
comp.phy = cHHMM.matrix.comparison(fit.phylogenetic, fit.phylogenetic.majority)

g.any.cs.phy = as.ggplot(pheatmap(comp.any.cs.phy$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
         cex=1,clustering_distance_rows = "manhattan", cex=1,
         clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.cs = as.ggplot(pheatmap(comp.cs$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
                                  cex=1,clustering_distance_rows = "manhattan", cex=1,
                                  clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.phy = as.ggplot(pheatmap(comp.phy$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
                                  cex=1,clustering_distance_rows = "manhattan", cex=1,
                                  clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.matrices = ggarrange(g.any.cs.phy,
                       g.cs,
                       g.phy,
                       g.cids,
                       labels = c("A. Any cs v phy", "B. cs any vs maj", "C. phy any vs maj", ""))

png(paste0(expt, "-comparison-matrices.png"), width=1000*sf, height=600*sf, res=72*sf)
print(g.matrices)
dev.off()


############## Optimal cluster: NbClust

# NbClust method without Monte Carlo (“bootstrap”) samples  Charrad et al., 2014

clustered.structure.alt = cHHMM.cluster.features(src.data ,method="NbClust")

## NbClust method plot

sd_index <- data.frame((clustered.structure.alt[["gap_stat"]])$All.index)
# Cluster numbers
c_numbers <- 2:20
# Create data frame
cluster_df <- data.frame(sd_index,c_numbers )
colnames(cluster_df )[1]="sd_index"


# Find the index of the minimum cluster value
min_index <- which.min(cluster_df$sd_index)
# Plot
g.nbclust = ggplot(data = cluster_df, aes(x = c_numbers, y = sd_index)) +
  geom_point(size = 4) +  # Increase point size
  geom_line(size = 1.5) +  # Increase line size
  geom_point(data = cluster_df[min_index, ], aes(x = c_numbers, y = sd_index, color = "Optimal cluster"), size = 6) +  # Circle the minimum cluster value with blue color
  labs(title = "SD indices vs. Cluster Numbers",
       x = "Cluster Numbers",
       y = "SD validity index") +
  scale_x_continuous(breaks = c_numbers) +
  #scale_color_manual(values = c("Optimal cluster" = "blue")) +  # Manual legend for color
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white"),  # Set plot background color
    panel.background = element_rect(fill = "white"),     # Set panel background color
    panel.grid.major = element_line(color = "gray"),      # Set major grid line color
    panel.grid.minor = element_blank()                    # Remove minor grid lines
  )

png(paste0(expt, "-alt-nbclust.png"), width=400*sf, height=300*sf, res=72*sf)
print(g.nbclust)
dev.off()

####
png(paste0(expt, "-alt-clusters.png"), width=600*sf, height=400*sf, res=72*sf)
print(
fviz_cluster(clustered.structure.alt[["km_res"]], data = clustered.structure.alt[["data"]],
             ellipse.type = "convex",
             palette = "viridis",
             ggtheme = theme_minimal())
)
dev.off()

ids = clustered.structure.alt$cluster_tab
df.plot = data.frame()
levels = rep(0, ncol(ids))
for(i in 1:nrow(ids)) {
  cid = which(ids[i,] == 1)
  df.plot = rbind(df.plot, data.frame(cluster=cid, level=levels[cid], name=rownames(ids)[i]))
  levels[cid] = levels[cid]+1
}
g.cids.alt = ggplot(df.plot, aes(x=cluster, y=level, label=name)) + geom_text() +
  scale_x_continuous(breaks=1:ncol(ids)) + scale_y_continuous(breaks=NULL) +
  labs(x = "Cluster ID", y = "") + theme_minimal()
png(paste0(expt, "-alt-cluster-IDs.png"), width=600*sf, height=300*sf, res=72*sf)
print(g.cids.alt)
dev.off()

cross.sectional.obs.alt = cHHMM.cross.sectional(clustered.structure.alt)
phylogenetic.obs.alt = cHHMM.phylogenetic.estimation(clustered.structure.alt)

cross.sectional.obs.majority.alt = cHHMM.cross.sectional(clustered.structure.alt, occupancy="majority")
phylogenetic.obs.majority.alt = cHHMM.phylogenetic.estimation(clustered.structure.alt, occupancy="majority")

cross.sectional.obs.relative.alt = cHHMM.cross.sectional(clustered.structure.alt, occupancy="relative")
phylogenetic.obs.relative.alt = cHHMM.phylogenetic.estimation(clustered.structure.alt, occupancy="relative")

fit.cross.sectional.alt = HyperHMM(cross.sectional.obs.alt$cross_sectional_data, nboot = 2)
fit.phylogenetic.alt = HyperHMM(phylogenetic.obs.alt$dests, initialstates = phylogenetic.obs.alt$srcs, nboot = 2)
fit.cross.sectional.majority.alt = HyperHMM(cross.sectional.obs.majority.alt$cross_sectional_data, nboot = 2)
fit.phylogenetic.majority.alt = HyperHMM(phylogenetic.obs.majority.alt$dests, initialstates = phylogenetic.obs.majority.alt$srcs, nboot = 2)
fit.cross.sectional.relative.alt = HyperHMM(cross.sectional.obs.relative.alt$cross_sectional_data, nboot = 2)
fit.phylogenetic.relative.alt = HyperHMM(phylogenetic.obs.relative.alt$dests, initialstates = phylogenetic.obs.relative.alt$srcs, nboot = 2)

g.fluxes.alt =  ggarrange(
  plot.hypercube.flux(fit.cross.sectional.alt, thresh = 0.02),
  plot.hypercube.flux(fit.phylogenetic.alt, thresh=0.02),
  labels=c("A. CS", "B. Phy"))
png(paste0(expt, "-alt-out-fluxes.png"), width=600*sf, height=400*sf, res=72*sf)
print(g.fluxes.alt)
dev.off()

g.all.alt = ggarrange( plot.cHHMM(cross.sectional.obs.alt$cross_sectional_data, fit.cross.sectional.alt, label="Cross-sectional\nAny occupancy"),
                   plot.cHHMM(phylogenetic.obs.alt$dests, fit.phylogenetic.alt, label="Phylogenetic\nAny occupancy"),
                   plot.cHHMM(cross.sectional.obs.majority.alt$cross_sectional_data, fit.cross.sectional.majority.alt, label="Cross-sectional\nMajority occupancy"),
                   plot.cHHMM(phylogenetic.obs.majority.alt$dests, fit.phylogenetic.majority.alt, label="Phylogenetic\nMajority occupancy"),
                   plot.cHHMM(cross.sectional.obs.relative.alt$cross_sectional_data, fit.cross.sectional.relative.alt, label="Cross-sectional\nRelative occupancy"),
                   plot.cHHMM(phylogenetic.obs.relative.alt$dests, fit.phylogenetic.relative.alt, label="Phylogenetic\nRelative occupancy"),
                   nrow=3, ncol=2,
                   labels=c("A.", "B.", "C.", "D.", "E.", "F."))
png(paste0(expt, "-alt-out-all.png"), width=1600*sf, height=1200*sf, res=72*sf)
print(g.all)
dev.off()


comp.any.cs.phy.alt = cHHMM.matrix.comparison(fit.cross.sectional.alt, fit.phylogenetic.alt)
comp.cs.alt = cHHMM.matrix.comparison(fit.cross.sectional.alt, fit.cross.sectional.majority.alt)
comp.phy.alt = cHHMM.matrix.comparison(fit.phylogenetic.alt, fit.phylogenetic.majority.alt)

g.any.cs.phy.alt = as.ggplot(pheatmap(comp.any.cs.phy.alt$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
                                  cex=1,clustering_distance_rows = "manhattan", cex=1,
                                  clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.cs.alt = as.ggplot(pheatmap(comp.cs.alt$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
                          cex=1,clustering_distance_rows = "manhattan", cex=1,
                          clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.phy.alt = as.ggplot(pheatmap(comp.phy.alt$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
                           cex=1,clustering_distance_rows = "manhattan", cex=1,
                           clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T))
g.matrices.alt = ggarrange(g.any.cs.phy.alt,
                       g.cs.alt,
                       g.phy.alt,
                       g.cids.alt,
                       labels = c("A. Any cs v phy", "B. cs any vs maj", "C. phy any vs maj", ""))

png(paste0(expt, "-alt-comparison-matrices.png"), width=1000*sf, height=600*sf, res=72*sf)
print(g.matrices.alt)
dev.off()

