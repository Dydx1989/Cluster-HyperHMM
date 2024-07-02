#Gap method without Monte Carlo (“bootstrap”) samples (B)
source("cHHMM_KD_methods.R")

# we must have source data with observations as rows and features as columns
# there was still specific code assuming 8 features (e.g. x+7) -- I think I have removed this
# applying the code to the other cases below throws errors
# we must allow flexibility in cluster assignment e.g. a single member feature gives me a 1 vs a majority of member features vs ...?
# with different seeds the clustering gives dramatically different results!

# this AMR cade study should reproduce, but the other examples below throw errors
# NB with different seeds the AMR case sometimes throws an error too
#Error in if (mm < k) stop("more cluster centers than distinct data points.") : 
#  missing value where TRUE/FALSE needed
AMR.data = read.csv("AMR_binary.csv", row.names=1)
clustered.structure = cHHMM.cluster.features(AMR.data,method="Gap")

### uncomment these for other test cases
### either
# other.data = t(read.csv("c4-curated.csv", header=FALSE))
### or
# other.data = t(read.csv("jallow_dataset_binary_with2s.csv"))
# other.data[other.data==2] = 1
# clustered.structure = cHHMM.cluster.features(other.data, n.boot = 10)


## Use plot for Gap method without Monte Carlo (“bootstrap”) samples (B)
#plot((clustered.structure[["gap_stat"]]))

first_max_index <- which(diff(sign(diff(clustered.structure[["gap_stat"]]$Gap))) == -2)[1]+1
# Plot the gap statistic using ggplot
ggplot(clustered.structure[["gap_stat"]], aes(x = 1:length(Gap), y = Gap)) +
  geom_line(color = "blue") + 
  geom_point(color = "blue") +
  geom_vline(xintercept = first_max_index, color = "red", linetype = "dashed") +
  labs(x = "Number of Clusters", y = "Gap Statistic", title = "Gap Statistic vs. Number of Clusters")

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



############## Optimal cluster: NbClust

# NbClust method without Monte Carlo (“bootstrap”) samples  Charrad et al., 2014
source("cHHMM_KD_methods.R")
# we must have source data with observations as rows and features as columns
# there was still specific code assuming 8 features (e.g. x+7) -- I think I have removed this
# applying the code to the other cases below throws errors
# we must allow flexibility in cluster assignment e.g. a single member feature gives me a 1 vs a majority of member features vs ...?
# with different seeds the clustering gives dramatically different results!

# this AMR cade study should reproduce, but the other examples below throw errors
# NB with different seeds the AMR case sometimes throws an error too
#Error in if (mm < k) stop("more cluster centers than distinct data points.") : 
#  missing value where TRUE/FALSE needed
AMR.data = read.csv("AMR_binary.csv", row.names=1)
clustered.structure = cHHMM.cluster.features(AMR.data ,method="NbClust")
### uncomment these for other test cases
### either
# other.data = t(read.csv("c4-curated.csv", header=FALSE))
### or
# other.data = t(read.csv("jallow_dataset_binary_with2s.csv"))
# other.data[other.data==2] = 1
# clustered.structure = cHHMM.cluster.features(other.data, n.boot = 10)

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





