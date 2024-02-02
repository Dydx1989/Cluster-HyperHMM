###   Data Preprocessing
### 
## Stage 1
library(readxl)
#### Input raw data such as K. pnuemonia
raw_dff=read.csv("BIGSdb_3708197_1272842014_68324.csv")
## Subset out everything with a nonzero entry in the resistance info column
data_dff_1= subset(raw_dff,resistance_info>0)

## Subset out everything that is exactly K. pneumoniae (our particular species of interest)
data_dff_2=subset(data_dff_1,taxonomic_designation=="K. pneumoniae")

## Assume (for now) that predictor genes are 0 (empty cell) or 1 (non-empty cell)
library(dplyr)
data_dff_new <- data_dff_2 %>% replace(is.na(.), 0)
head(data_dff_new[1:5])
#write.csv(data_dff_new,"data_dff_new.csv")
data_dff_new[data_dff_new==""] <- 0

data_dff_new_1=data_dff_new[1:14]
data_dff_new_2=data_dff_new[14:dim(data_dff_new)[2]]
dim(data_dff_new_2)

data_dff_new_3=data.frame(lapply(data_dff_new_2, function(x) { 
  if(is.factor(x)) {
    x <- as.character(x)
  }
  x[x!="0"]="1"
  x}))

Data_trans=data.frame(data_dff_new_1,data_dff_new_3)
## Extract Isolate and genes only 
df_isolate=Data_trans[2]
df_gene=Data_trans[,-c(1:13)]
Data_isolate_genes=data.frame(df_isolate,df_gene)

## Heatmap 1
Data_isolate_genes_1 = as.data.frame(sapply(Data_isolate_genes[,-1], as.numeric))
library(ggplot2)
library("gplots")
heatmap.2(t(Data_isolate_genes_1), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

###   Data Preprocessing: Stage 2
# In this stage, we remove all the isolates with completely empty cells and retain only isolates with nonempty cells
remove_zero_cols_and_rows <- function(df) {
  # Remove zero columns
  zero_cols <- sapply(df, function(col) all(col == 0))
  df <- df[, !zero_cols]
  
  # Remove zero rows (ignoring the first column)
  zero_rows <- rowSums(df[, -1] != 0, na.rm = TRUE) == 0
  df <- df[!zero_rows, , drop = FALSE]
  
  return(df)
}

Final_data=remove_zero_cols_and_rows(Data_isolate_genes)
dim(Final_data)

## Heatmap 2
Final_data_1 = as.data.frame(sapply(Final_data[,-1], as.numeric))
library(ggplot2)
library("gplots")
heatmap.2(t(Final_data_1), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
