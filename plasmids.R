data <- read.csv("plasmids.tsv", sep='\t')
parse.percent <- \(x) as.numeric(substr(x, 1, nchar(x)-1)) / 100
data$Gene <- as.factor(data$Gene)
data$Cluster <- as.factor(data$Cluster)
data$NCBI.Chromosome <- parse.percent(data$NCBI.Chromosome)
data$NCBI.Plasmid <- parse.percent(data$NCBI.Plasmid)
data$NCBI.WGS <- parse.percent(data$NCBI.WGS)
data$NCBI.GI <- parse.percent(data$NCBI.GI)

library(ggplot2)
library(ggpubr)

ggarrange(ggplot(data) + 
  geom_text(aes(NCBI.Plasmid, NCBI.Chromosome, label=Gene, color=Cluster)),

ggplot(data) + 
  geom_text(aes(log(NCBI.Plasmid), log(NCBI.Chromosome), label=Gene, color=Cluster)), common.legend = TRUE)
