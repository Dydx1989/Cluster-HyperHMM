#######################################################################
## ******************************************************##############
## **** Read in Continent_cluster_data for subsetting ******###########
## ******************************************************##############
#######################################################################
################ Continent_cluster_data.csv #############################
## the data consist of  continents and clusters (C1,....,C8)
Cont_Clsuter=read.csv("Continent_cluster_data.csv")
head(Cont_Clsuter)
dim(Cont_Clsuter)
 
## Now data subsetting 

##  Africa susbet

Africa=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="Africa")
#head(Africa)
#dim(Africa)
Africa=Africa[-1]
# Write data frame to a text file without column and row names
data_Africa <- as.matrix(Africa)
#Concatenate the values of each row without gaps
row_values <- apply(data_Africa, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "Africa.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


##  Asia susbet
Asia=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="Asia")
#head(Asia)
#dim(Asia)
Asia=Asia[-1]
# Write data frame to a text file without column and row names
data_Asia <- as.matrix(Asia)
#Concatenate the values of each row without gaps
row_values <- apply(data_Asia, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "Asia.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


##  Europe 
Europe=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="Europe")
#head(Europe)
#dim(Europe)
Europe=Europe[-1]
# Write data frame to a text file without column and row names
data_Europe <- as.matrix(Europe)
#Concatenate the values of each row without gaps
row_values <- apply(data_Europe, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "Europe.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)



##  North America
North_America=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="North America")
#head(North_America)
#dim(North_America)
North_America=North_America[-1]
# Write data frame to a text file without column and row names
data_North_America <- as.matrix(North_America)
#Concatenate the values of each row without gaps
row_values <- apply(data_North_America, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "North_America.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

##  Oceania
Oceania=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="Oceania")
#head(Oceania)
#dim(Oceania)
Oceania=Oceania[-1]
# Write data frame to a text file without column and row names
data_Oceania <- as.matrix(Oceania)
#Concatenate the values of each row without gaps
row_values <- apply(data_Oceania, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "Oceania.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)


##  South America
South_America=subset(Cont_Clsuter,  Cont_Clsuter$Continent=="South America")
#head(South_America)
#dim(South_America)
South_America=South_America[-1]
# Write data frame to a text file without column and row names
data_South_America <- as.matrix(South_America)
#Concatenate the values of each row without gaps
row_values <- apply(data_South_America, 1, function(x) paste0(x, collapse = ""))
write.table(row_values , "South_America.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)







