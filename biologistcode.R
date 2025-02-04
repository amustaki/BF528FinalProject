library("AnnotationDbi")
library("hgu133plus2.db")
library("GSEABase")
library("affy")
library("tidyverse")
library("dplyr")
#read in the expression data
data <- read.csv("biologist.csv", header = TRUE)


#using the select function from the package hgu133plus, map the porbeset ids to gene symbols specfiying by appriate jey and column 
match <- select(hgu133plus2.db, keys = as.character(data$probeset_IDfil2), columns = "SYMBOL", keyset = "PROBEID")
df <- merge(x=match, y=data, by.x = "PROBEID", by.y = "probeset_IDfil2")

#part 2 
#locate and download the KEGG, GO, and Hallmark gene sets and upload them to us later

#part 3
#get rid of the duplicated gene symbols 
#use x[!duplicated(x)]
df <- df[!duplicated(df$SYMBOL),]
#get rid of na and blank values 
df <- df[!(is.na(df$SYMBOL) | df$SYMBOL == ""),]
#select the top 1000 up regulated genes, log2Fc is postive 
posdf <- filter(df, log2FC > 0)
#put the - to order by descending 
#select the top 1000 of the positive log2FC
up1000 <- posdf[order(-posdf$log2FC),] %>% top_n(1000)

#now do the downregulared ones 
negdf <- filter(df, log2FC < 0)
down1000 <- negdf[order(-negdf$log2FC),] %>% top_n(1000)

#create a table of the top10 up and down regulated genes including the tstat, pvalue, padj, 
up10 <- up1000[order(-up1000$log2FC),] %>% top_n(10) 
down10 <- down1000[order(-down1000$log2FC),] %>% top_n(10)
total <- rbind(up10,down10)
total

#part 4
#use GSEABase 
#read in GMT files and find out how many gene sets they are in each of the collections 
hallmark <- getGmt('h.all.v7.5.1.symbols.gmt.txt')
kegg <- getGmt('c2.cp.kegg.v7.5.1.symbols.gmt.txt')
go <- getGmt('c5.go.v7.5.1.symbols.gmt.txt')

print(paste0("There are ", length(hallmark)," gene sets in the Hallmark collection"))
print(paste0("There are ", length(kegg)," gene sets in the Kegg collection"))
print(paste0("There are ", length(go)," gene sets in the GO collection"))

#part 5
#use the Fisher Exact Test to compute hypergeometric statistics and p-values comparing overlap for each gene set and each top 1000 increased 1000 decreased genes 

#get the genes that were not expressed
not_diff_up <- subset(df, !df$SYMBOL %in% up1000$SYMBOL)
length(not_diff_up$SYMBOL)
not_diff_down <- subset(df, !df$SYMBOL %in% down1000$SYMBOL)

#making the contingency table 
fishertest <- function(genelist, geneset, notdiff){
  diffexp_ingene <- length(intersect(genelist,geneset)) #diff exp gene that in geneset #6
  diffexp_notgene <- length(genelist) - diffexp_ingene #diff exp not in geneset
  notexp_ingene <- length(intersect(notdiff,geneset)) #not diff exp but in gene set 
  notexp_notgene <- length(notdiff) - notexp_ingene #not diff and not in gene set
  return(c(diffexp_ingene,diffexp_notgene,notexp_ingene,notexp_notgene))
}

#make a dataframe in order to store the results of the hallmark 
#gene set name, statistic est, pvalue, padj 

hallmark_results <- data.frame(geneset_name = character(), pvalue = numeric(), test_statistic = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(hallmark))
{
  geneid <- geneIds(hallmark[i])
  up_fisher <- fishertest(up1000$SYMBOL, geneid[[names(geneid)]], not_diff_up$SYMBOL)
  down_fisher <- fishertest(down1000$SYMBOL, geneid[[names(geneid)]], not_diff_down$SYMBOL)
  upper <- matrix(unlist(up_fisher),2)
  up <- fisher.test(upper)
  dowwn <- matrix(unlist(down_fisher),2)
  #down <- fisher.test(matrix(unlist(down_fisher), nrow=2))
  down <- fisher.test(dowwn)
  hallmark_results[nrow(hallmark_results) +1,] <- c(names(geneid), up$p.value, up$estimate, "UP")
  hallmark_results[nrow(hallmark_results) +1,] <- c(names(geneid), down$p.value, down$estimate, "DOWN")}
  

head(hallmark_results)
hallmark_results <- hallmark_results %>% mutate(pvalue = as.numeric(pvalue), test_statistic = as.numeric(test_statistic))

kegg_results <- data.frame(geneset_name = character(), pvalue = numeric(), test_statistic = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(kegg))
{
  geneid <- geneIds(kegg[i])
  up_fisher <- fishertest(up1000$SYMBOL, geneid[[names(geneid)]], not_diff_up$SYMBOL)
  down_fisher <- fishertest(down1000$SYMBOL, geneid[[names(geneid)]], not_diff_down$SYMBOL)
  upper <- matrix(unlist(up_fisher),2)
  up <- fisher.test(upper)
  dowwn <- matrix(unlist(down_fisher),2)
  #down <- fisher.test(matrix(unlist(down_fisher), nrow=2))
  down <- fisher.test(dowwn)
  kegg_results[nrow(kegg_results) +1,] <- c(names(geneid), up$p.value, up$estimate, "UP")
  kegg_results[nrow(kegg_results) +1,] <- c(names(geneid), down$p.value, down$estimate, "DOWN")}

head(kegg_results)
kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), test_statistic = as.numeric(test_statistic))

go_results <- data.frame(geneset_name = character(), pvalue = numeric(), test_statistic = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(go))
{
  geneid <- geneIds(go[i])
  up_fisher <- fishertest(up1000$SYMBOL, geneid[[names(geneid)]], not_diff_up$SYMBOL)
  down_fisher <- fishertest(down1000$SYMBOL, geneid[[names(geneid)]], not_diff_down$SYMBOL)
  upper <- matrix(unlist(up_fisher),2)
  up <- fisher.test(upper)
  dowwn <- matrix(unlist(down_fisher),2)
  #down <- fisher.test(matrix(unlist(down_fisher), nrow=2))
  down <- fisher.test(dowwn)
  go_results[nrow(go_results) +1,] <- c(names(geneid), up$p.value, up$estimate, "UP")
  go_results[nrow(go_results) +1,] <- c(names(geneid), down$p.value, down$estimate, "DOWN")}

head(go_results)
go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), test_statistic = as.numeric(test_statistic))


#add padj values for each of the dataframes: hallmark, kegg, go
hallmark_results$padj <- p.adjust(hallmark_results$pvalue, method = "BH")
kegg_results$padj <- p.adjust(kegg_results$pvalue, method = "BH")
go_results$padj <- p.adjust(go_results$pvalue, method = "BH")

#number of sig enriched genesets by padj < 0.05 

sig_enriched_hallmark <- filter(hallmark_results,padj < 0.05)
length(sig_enriched_hallmark$geneset_name)
sig_enriched_kegg <- filter(kegg_results,padj < 0.05)
length(sig_enriched_kegg$geneset_name)
sig_enriched_go <- filter(go_results,padj < 0.05)
length(sig_enriched_go$geneset_name)

print(paste0("The number of signifcantly enriched gene sets with padj < 0.05 in Hallmarks is ", length(sig_enriched_hallmark$geneset_name)))
print(paste0("The number of signifcantly enriched gene sets with padj < 0.05 in KEGG is ", length(sig_enriched_kegg$geneset_name)))
print(paste0("The number of signifcantly enriched gene sets with padj < 0.05 in GO is ", length(sig_enriched_go$geneset_name)))

#top 3 enriched gene sets for each gene set type 
most_enriched_hallmark <- sig_enriched_hallmark[order(-sig_enriched_hallmark$padj),] %>% top_n(3)
most_enriched_hallmark
most_enriched_kegg <- sig_enriched_kegg[order(-sig_enriched_kegg$padj),] %>% top_n(3)
most_enriched_kegg
most_enriched_go <- sig_enriched_go[order(-sig_enriched_go$padj),] %>% top_n(3)
most_enriched_go

write.csv(hallmark_results, "hallmark_results.csv")
write.csv(kegg_results, "kegg_results.csv")
write.csv(go_results, "go_result.csv")





