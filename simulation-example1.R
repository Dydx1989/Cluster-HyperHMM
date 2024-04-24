# Set the parameters
n_samples <- 800  # Number of Genes
n_variables <- 200  # Number of Isolates
prob <- 0.5        # Probability of success (probability of 1)

# Simulate multivariate binary data
simulated_data <- matrix(rbinom(n_samples * n_variables, 1, prob), nrow = n_samples, ncol = n_variables)

# Name the rows of the matrix
rownames(simulated_data) <- paste0("row", 1:n_samples)

# Convert to data frame if needed
simulated_data <- as.data.frame(simulated_data)

# View the first few rows of the simulated data
head(simulated_data)
dim(simulated_data)




source("cHHMM_KD_methods.R")
clustered.structure_sim = cHHMM.cluster.features(simulated_data,method="Gap")
## Use plot for Gap method without Monte Carlo (“bootstrap”) samples (B)
#plot((clustered.structure_sim[["gap_stat"]]))
first_max_index <- which(diff(sign(diff(clustered.structure_sim[["gap_stat"]]$Gap))) == -2)[1]+1
# Plot the gap statistic using ggplot
ggplot(clustered.structure_sim[["gap_stat"]], aes(x = 1:length(Gap), y = Gap)) +
  geom_line(color = "blue") + 
  geom_point(color = "blue") +
  geom_vline(xintercept = first_max_index, color = "red", linetype = "dashed") +
  labs(x = "Number of Clusters", y = "Gap Statistic", title = "Gap Statistic vs. Number of Clusters")

####
fviz_cluster(clustered.structure_sim[["km_res"]], data = clustered.structure_sim[["data"]],
             ellipse.type = "convex",
             palette = "viridis",
             ggtheme = theme_minimal())



cross.sectional.obs.sim = cHHMM.cross.sectional(clustered.structure_sim)
phylogenetic.obs.sim = cHHMM.phylogenetic.estimation(clustered.structure_sim)



fit.cross.sectional.sim = HyperHMM(cross.sectional.obs.sim$cross_sectional_data)


## sIM
clustered.structure1 = cHHMM.cluster.features(simulated_data,method="NbClust")

sd_index <- data.frame((clustered.structure1[["gap_stat"]])$All.index)
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
fviz_cluster(clustered.structure1[["km_res"]], data = clustered.structure1[["data"]],
             ellipse.type = "convex",
             palette = "viridis",
             ggtheme = theme_minimal())

cross.sectional.obs = cHHMM.cross.sectional(clustered.structure1)
phylogenetic.obs = cHHMM.phylogenetic.estimation(clustered.structure1)
