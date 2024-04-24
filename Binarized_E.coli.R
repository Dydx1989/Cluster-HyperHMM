
# Read the data
E.data=read.csv("E.coli.data.csv",row.names = 1)
dim(E.data)
head(E.data[,1:8])


# Transpose the data frame
transposed_data <- t(E.data)
head(transposed_data[,1:8])
dim(transposed_data)

# Binarize the transposed data frame using median expression values
binarized_data <- apply(transposed_data, 2, function(x) ifelse(x > mean(x), 1, 0))
head(binarized_data [,1:8])
dim(binarized_data)
# Sample data frame (replace this with your actual data frame)

# Remove columns with all 1s or all 0s
#binarized_data_filtered <- binarized_data[, colSums(binarized_data) != 0 | colSums(binarized_data) != nrow(binarized_data)]

# Remove columns with all 0s or all 1s
binarized_data_filtered <- binarized_data[, apply(binarized_data, 2, function(col) !all(col == 0) & !all(col == 1))]

dim(binarized_data_filtered)

# Calculate the percentage of presence for each gene
presence_percent <- colMeans(binarized_data_filtered)

# Identify genes that are present in at least 70% of the samples
genes_to_keep <- colnames(binarized_data_filtered)[presence_percent >= 0.67]
# Subset the original data frame to retain only the selected genes
binarized_data_filtered <- binarized_data_filtered[, genes_to_keep]
dim(binarized_data_filtered)
head(binarized_data_filtered[,1:6])

# Add gene names as column names
write.csv(t(binarized_data_filtered),"Bi_Ecoli.data.csv")



