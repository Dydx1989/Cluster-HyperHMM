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
AMR.data = read.csv("AMR_binary.csv", row.names=1)
clustered.structure = cHHMM.cluster.features(AMR.data)

### uncomment these for other test cases
### either
# other.data = t(read.csv("c4-curated.csv", header=FALSE))
### or
# other.data = t(read.csv("jallow_dataset_binary_with2s.csv"))
# other.data[other.data==2] = 1
# clustered.structure = cHHMM.cluster.features(other.data, n.boot = 10)

fviz_gap_stat(clustered.structure[["gap_stat"]])

fviz_cluster(clustered.structure[["km_res"]], data = clustered.structure[["data"]],
             ellipse.type = "convex",
             palette = "jco",
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



