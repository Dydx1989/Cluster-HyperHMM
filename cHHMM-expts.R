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

expt = "malaria"
sf = 2

if(expt == "AMR") {
  src.data = read.csv("AMR_binary.csv", row.names=1)
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

fit.cross.sectional = HyperHMM(cross.sectional.obs$cross_sectional_data, nboot = 2)
fit.phylogenetic = HyperHMM(phylogenetic.obs$dests, initialstates = phylogenetic.obs$srcs, nboot = 2)
fit.cross.sectional.majority = HyperHMM(cross.sectional.obs.majority$cross_sectional_data, nboot = 2)
fit.phylogenetic.majority = HyperHMM(phylogenetic.obs.majority$dests, initialstates = phylogenetic.obs.majority$srcs, nboot = 2)

g.fluxes =  ggarrange(
  plot.hypercube.flux(fit.cross.sectional, thresh = 0.02),
  plot.hypercube.flux(fit.phylogenetic, thresh=0.02),
  labels=c("A. CS", "B. Phy"))
png(paste0(expt, "-out-fluxes.png"), width=600*sf, height=400*sf, res=72*sf)
print(g.fluxes)
dev.off()

g.all = ggarrange( plot.standard(fit.cross.sectional),
           plot.standard(fit.phylogenetic), 
           plot.standard(fit.cross.sectional.majority),
           plot.standard(fit.phylogenetic.majority),
           nrow=2, ncol=2,
           labels=c("A. CS any", "B. Phy any", "C. CS maj", "D. Phy maj"))
png(paste0(expt, "-out-all.png"), width=1600*sf, height=800*sf, res=72*sf)
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
g.matrices = ggarrange(g.any.cs.phys,
                       g.cs,
                       g.phy,
                       labels = c("A. Any cs v phy", "B. cs any vs maj", "C. phy any vs maj"))

png(paste0(expt, "-comparison-matrices.png"), width=1000*sf, height=600*sf, res=72*sf)
print(g.matrices)
dev.off()


############## Optimal cluster: NbClust

# NbClust method without Monte Carlo (“bootstrap”) samples  Charrad et al., 2014

# we must have source data with observations as rows and features as columns
# there was still specific code assuming 8 features (e.g. x+7) -- I think I have removed this
# applying the code to the other cases below throws errors
# we must allow flexibility in cluster assignment e.g. a single member feature gives me a 1 vs a majority of member features vs ...?
# with different seeds the clustering gives dramatically different results!

# this AMR cade study should reproduce, but the other examples below throw errors
# NB with different seeds the AMR case sometimes throws an error too
#Error in if (mm < k) stop("more cluster centers than distinct data points.") : 
#  missing value where TRUE/FALSE needed

clustered.structure = cHHMM.cluster.features(AMR.data ,method="NbClust")
clustered.structure = cHHMM.cluster.features(other.data ,method="NbClust")

## NbClust method plot

sd_index <- data.frame((clustered.structure[["gap_stat"]])$All.index)
# Cluster numbers
c_numbers <- 2:20
# Create data frame
cluster_df <- data.frame(sd_index,c_numbers )
colnames(cluster_df )[1]="sd_index"


# Find the index of the minimum cluster value
min_index <- which.min(cluster_df$sd_index)
# Plot
ggplot(data = cluster_df, aes(x = c_numbers, y = sd_index)) +
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
    plot.background = element_rect(fill = "lightblue"),  # Set plot background color
    panel.background = element_rect(fill = "white"),     # Set panel background color
    panel.grid.major = element_line(color = "gray"),      # Set major grid line color
    panel.grid.minor = element_blank()                    # Remove minor grid lines
  )




####
fviz_cluster(clustered.structure[["km_res"]], data = clustered.structure[["data"]],
             ellipse.type = "convex",
             palette = "viridis",
             ggtheme = theme_minimal())

cross.sectional.obs = cHHMM.cross.sectional(clustered.structure)
phylogenetic.obs = cHHMM.phylogenetic.estimation(clustered.structure)

fit.cross.sectional = HyperHMM(cross.sectional.obs$cross_sectional_data)

fit.phylogenetic = HyperHMM(phylogenetic.obs$dests, initialstates = phylogenetic.obs$srcs)

ggarrange( plot.standard(fit.cross.sectional),
           plot.standard(fit.phylogenetic), nrow=2 )

matrix.comp = cHHMM.matrix.comparison(fit.cross.sectional, fit.phylogenetic)

pheatmap(matrix.comp$results_matrix, show_rownames = TRUE,cluster_cols=F,cluster_rows=F,
         cex=1,clustering_distance_rows = "manhattan", cex=1,
         clustering_distance_cols = "manhattan", clustering_method = "complete",border_color = TRUE,display_numbers = T)





