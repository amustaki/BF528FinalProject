#analyst section of project1 
library(tidyverse)
library(dplyr)
#load in the RMA normalized ComBat adjusted expression matrix 
expression <- read.csv('/projectnb/bf528/users/lush_2022/project_1/expression_data.csv', row.names = 1)
head(expression)
#filter 1 
#expressed in at least 20% of the samples must be >log2(15)
#so 20% of the columns must be >log2(15)
filter20 <- function(data){
  threshold <- log2(15)
  twenty <- .20*ncol(data)
  new <- mutate(data, thres = rowSums(data >threshold))
  above20 <- subset(new,thres >= twenty)
  #drop the thres column so it doesnt mess up the rest of the analysis
  above20 <- select(above20, -thres)
  return(above20)
}

filter1<- filter20(expression)
#head(filter1)
#have varinace slightly different then the median variance of all probe sets using a threshold of p > 0.01
#get the degrees of freedom for the chi squared test 

#filter2 is the chi sqaured filter 
chifilter <- function(data){
  degreefred <- ncol(data) -2 
  var <- apply(data,1,var)
  med <- median(var)
  chisquared <- data[(degreefred * (var/med)) > qchisq(0.99,degreefred,lower.tail = FALSE),]
  return(chisquared)
}
  
filter2 <- chifilter(filter1)

#filter 3 is the coefficent variance which is sd/mean 
coeff <- function(data){
  #calculate the coefficent of variance 
  meandat <- apply(data,1,mean)
  sddat <- apply(data, 1,sd)
  newdata <- mutate(data, coeff = (sddat/meandat))
  abovecoeff <- subset(newdata, coeff > 0.186)
  abovecoeff <- select(abovecoeff, -coeff)
}

filter3 <- coeff(filter2)

#write out csv for all the genes that pass all three of the filters 
write.csv(filter3, "pass3filtersgenes.csv")

#write out csv for the genes that pass the second filter, the chi squared filter 
write.csv(filter2, "chisquaredfiltergenes.csv")


#preform hierarchical clustering on the fully filtered dataset 
hier_cluster <- hclust(dist(t(filter3)))
#plot(hier_cluster, main = "Hierarchical Clustering Plot")

#cut the dendrogram into two clusters 
clust <- cutree(hier_cluster, 2)
#clust_vec <- as.vector(clust)
cluster1 <- sum(clust == 1)
cluster2 <- sum(clust ==2)
print(paste0("The number of samples in cluster 1 is ", cluster1))
print(paste0("The number of samples in cluster 2 is ", cluster2))

#creating the heatmap 
metadata <- read.csv("/project/bf528/project_1/doc/proj_metadata.csv")
metadata <- subset(metadata, select = c(geo_accession, cit.coloncancermolecularsubtype))
#do an ifelse statement for the colors for the heatmap for it to be red for C3 and blue for everything 
color <- ifelse(metadata$cit.coloncancermoleculars == "C3", "red", "blue")
#save the heatmap as a file jpeg 
jpeg("/usr4/bf527/amustaki/projectnb/bf528/students/amustaki/heatmap_allsamples.jpeg")
heatmap(as.matrix(filter3), ColSideColors = color, main = "Heatmap of Gene Expression across all Samples") 
legend(x = "topright", legend = c("C3","Other"), fill = c("red", "blue"), title = "Subtype")
dev.off()


#identify genes differentially expressed between the two clusters using welch T-test 
#used diff expressed matrix from 4.4, which is the list that past all 3 
# the clusters 1 and 2 were made earlier in the code called cluster1 and cluster2 
#make it a matrix 

clust1 <- filter3[, clust == 1]
clust2 <- filter3[, clust ==2]

welch_test <- apply(as.matrix(filter3), MARGIN = 1, function(x) t.test(x=x[clust==1], y=x[clust==2]))
welch_test


pvalue <- c(sapply(welch_test,function(x) x$p.value)) #getting the pvalue for the datafra,e
t_stat <- c(sapply(welch_test, function(x) x$statistic)) #getting the tstat for the dataframe
padj <- c(p.adjust(pvalue, method = "fdr")) #get the padj using fdr 
probeset_ID <- c(row.names(filter3)) #get the ids 
cluster_data <- data.frame(probeset_ID, t_stat, pvalue, padj)
head(cluster_data)

#filtering by pvalue and get the length to see how many genes that has 
diff_pvalue_expressed <- filter(cluster_data, pvalue < 0.05)
length(diff_pvalue_expressed$pvalue)
#how many genes are the most differentially expressed at adjusted p<0.05 
diff_expressed <- filter(cluster_data, padj < 0.05)
diff_expressed
#get the length to figure out how many genes
length(diff_expressed$padj)

#now repeat the same process again for filter 2 and save the results as a csv 
#add in log2fc bc it is needed for the biologist role 
welch_test_fil2 <- apply(as.matrix(filter2), MARGIN = 1, function(x) t.test(x=x[clust==1], y=x[clust==2]))
welch_test_fil2

pvaluefil2 <- c(sapply(welch_test_fil2, function(x) x$statistic))
tstatfil2 <- t_stat <- c(sapply(welch_test_fil2, function(x) x$statistic))
padjfil2 <- c(p.adjust(pvaluefil2, method = "fdr"))
probeset_IDfil2 <- c(row.names(filter2))
length(tstatfil2)
#adding in the log2FC in the dataset 

hier_cluster2 <- hclust(dist(t(filter2)))
#plot(hier_cluster, main = "Hierarchical Clustering Plot")

#cut the dendrogram into two clusters 
clustchi <- cutree(hier_cluster, 2)


clustchi1 <- filter2[, clustchi == 1]
clustchi2 <- filter2[, clustchi ==2]





log2FC <- apply(log2(clustchi2),1,mean) - apply(log2(clustchi1),1,mean)
length(log2FC)



cluster_datafil2 <- data.frame(probeset_IDfil2,tstatfil2,pvaluefil2,padjfil2,log2FC)
head(cluster_datafil2)


#differntial expressed genes with pvalue < 0.05
fil2_diff_pvalue <- filter(cluster_datafil2, pvaluefil2 < 0.05)
length(fil2_diff_pvalue$pvaluefil2)

#differential expressed genes with padj < 0.05 
fil2_diff <- filter(cluster_datafil2, padjfil2 < 0.05)
length(fil2_diff$padjfil2)
head(fil2_diff)

#wirte out the differnetially expressed genes for the padj < 0.05 into a csv 
write.csv(fil2_diff, row.names = F, "biologist.csv")

head(filter2)

