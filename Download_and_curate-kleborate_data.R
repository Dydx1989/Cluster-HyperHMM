## Kazeem, Olav, and Iain
library(R.utils)

if (.Platform$OS.type == "windows") {
  download.method = "wininet"
  
} else if (.Platform$OS.type == "unix") {
  download.method = "wget"
  
} else {
  stop("Unknown Operating System")
}

pathogen.watch.endpoint <- "https://pathogenwatch-public.ams3.cdn.digitaloceanspaces.com/"

files <- c("Klebsiella pneumoniae__kleborate.csv",
           "Klebsiella pneumoniae__cgmlst.csv",
           "Klebsiella pneumoniae__inctyper.csv",
           "Klebsiella pneumoniae__metadata.csv",
           "Klebsiella pneumoniae__metrics.csv",
           "Klebsiella pneumoniae__mlst.csv")

for (filename in files) {
  if (!file.exists(paste0("../raw/",filename))) {
    dir.create("../raw/", showWarnings=FALSE)
    download.file(url=paste0(pathogen.watch.endpoint,filename,".gz"),
                  destfile=paste0("../raw/",filename,".gz"),
                  method=download.method)
    gunzip(paste0("../raw/",filename,".gz"))
  }
}

## Preprocessing Kleborate Data to binary 

kleborate.df <- read.csv("../raw/Klebsiella pneumoniae__kleborate.csv")

selection <- colnames(kleborate.df)[endsWith(colnames(kleborate.df),
                                             "_acquired")]

resistance.df <- kleborate.df[, selection]
colnames(resistance.df) <- gsub("_acquired","",colnames(resistance.df))
selection <- colnames(resistance.df)

resistance.df <- as.data.frame(apply(resistance.df,
                                     c(1,2),
                                     \(x) ifelse(x=="-", 0, 1)))

dir.create("../clean/", showWarnings=FALSE)
write.csv(cbind(id=kleborate.df$Genome.Name, resistance.df),
          paste0("../clean/kleborate-ARGs-binary.csv"),
          row.names=FALSE)


## Clean the data
kleborate.df <- read.csv("../raw/Klebsiella pneumoniae__kleborate.csv")

# Pick only columns that have strong hits for antibiotic resistance genes
selection <- colnames(kleborate.df)[endsWith(colnames(kleborate.df),
                                             "_acquired")]

# Collapse columns to list of character vectors only containing present genes
obs <- apply(kleborate.df[, selection], 1, paste, collapse=';')
obs <- gsub("(-;|;-)","", obs)
obs <- strsplit(obs, ";")

# Remove " +13V" like syntax (unknown significance)
obs <- sapply(obs, \(x) gsub(" \\+[0-9][0-9][0-9]?[A-Z]$", "", x))

# Remove asterix, hats and questionmarks all denoting distance to known alleles
obs <- sapply(obs, \(x) gsub("[*^?]*$", "", x))

# Remove .v1, .v2 etc
obs <- sapply(obs, \(x) gsub("\\.v[1-9]$", "", x))

# Create a list of all unique columns
columns <- sort(unique(c(unlist(obs))))
columns <- columns[columns != '-']

# Create a matrix of 1s and 0s denoting presense/absense for all genes
m <- t(sapply(obs, \(x) as.numeric(columns %in% x)))

# Create the dataframe
df <- data.frame(m)
colnames(df) <- columns
df <- cbind(id = kleborate.df$Genome.Name, df)

# Write it to file
dir.create("../clean", showWarnings=FALSE)
write.csv(df, 
          file = "../clean/kleborate-wide-binary.csv",
          row.names=FALSE)
