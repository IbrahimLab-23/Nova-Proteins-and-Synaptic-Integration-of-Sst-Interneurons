#All
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library("org.Mm.eg.db")
library(EnhancedVolcano)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

Counts <- read.delim("CountReads_csv.csv", header = TRUE, row.names = 1, sep = ",")
Counts <- Counts[,-1] #Reduce number of replocates 
library(tidyverse)
New_counts <-matrix(0, nrow= nrow(Counts), ncol=1)
for (i in 1:70) {
  for (u in 1:35){
    
    select_counts=rowMeans(Counts %>% select(i,i+1))
    
  }
  New_counts <- cbind(New_counts, select_counts)
  
}

New_counts= New_counts[,-1]
index= which(1:ncol(New_counts)%%2!=1)
Counts2=New_counts[,index]
Counts2= cbind(Counts2, Counts[,71:74])
condition <- factor(c("Nova1", "Nova1", "Nova1","Nova2", "Nova2", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO" ,"ECS_2KO" ,"ECS_2KO" , "N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1",  "Nova2", "Nova2","N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1", "Nova2", "Nova2", "Nova2", "Nova2", "ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT", "WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
colnames(Counts2)= condition
Counts= Counts2
coldata <- data.frame(row.names= colnames(Counts), condition )
adds <- DESeqDataSetFromMatrix(countData = round(Counts), colData = coldata, design = ~condition)
adds<- DESeq(adds)
vsdata <- vst(adds, blind = FALSE)
res <- as.matrix(results(adds))
sigs <- na.omit(res)
sigs<- as.data.frame(sigs)

#Plotting Volcano 
sigs.df <- as.data.frame(sigs)
sigs.df$symbol <- mapIds(org.Mm.eg.db, keys = rownames(sigs.df), keytype = "ENSEMBL", column = "SYMBOL")
EnhancedVolcano(sigs.df, x= "log2FoldChange", y = "padj", lab = sigs.df$symbol)

#Heatmap 
SynapticGenes= c( "Sema6a", "Sorbs1", "Hivep2", "Gnb1", "Cacna1c", "Wnk1", "Kif21a", "Birc6", "Rph3a", "Rogdi", "Dclk1", "Myo16", "Nrxn1", "Macf1", "Nrcam", "Ank2", "Stx3", "Hsf1", "Dnajc5", "Vars", "Hspa8" , "Sema6d", "Trpm7", "Arhgap26", "Fchsd2", "Cadm1", "Cask", "Tspyl2", "Nrxn3", "Tpd52l2", "Ank3", "Cacna1b", "Add1", "Rims1", "Smad2", "Pml", "Atf2", "Sptan1", "Mapk9", "Adam22", "Macf1", "Nrcam", "Gramd1a", "Nckap1", "Sh3glb1","Atxn2", "Bin1", "Hspa8", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8" , "Nrxn1", "Macf1", "Ercc5", "Setx", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8","Sema6d", "Trpm7", "Arhgap26", "Fchsd2", "Cadm1", "Cask", "Tspyl2", "Nrxn3", "Tpd52l2", "Ank3", "Cacna1b", "Add1", "Rims1", "Pml", "Atf2", "Sptan1", "Adam22", "Macf1", "Nrcam", "Gramd1a", "Sh3glb1","Atxn2", "Bin1", "Hspa8", "Vip", "SsT" )
SynapticGenes= unique(SynapticGenes)
sigs.df2=sigs.df[sigs.df$symbol %in% SynapticGenes ,]
mat<- counts(adds, normalized= T)[rownames(sigs.df2),]

#get z score for each row
mat.z<-t(apply(mat, 1, scale))
colnames(mat.z) <- rownames(coldata)
colnames(mat.z) <- condition
h<- Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score", row_labels = sigs.df[rownames(mat.z),]$symbol )
pheatmap(mat.z,labels_row =  sigs.df[rownames(mat.z),]$symbol)

#plotCounts
par(mfrow=c(2,3))
plotCounts(adds, gene="ENSMUSG00000021047", col= "lightblue", xlab = "Condition", main = "Nova1")#Nova1
plotCounts(adds, gene="ENSMUSG00000030411", col= "lightblue", xlab = "Condition", main = "Nova2")#Nova2
plotCounts(adds, gene="ENSMUSG00000019772", col= "lightblue", xlab = "Condition", main = "Vip")#VIP
plotCounts(adds, gene="ENSMUSG00000004366", col= "lightgreen", xlab = "Condition", main = "SST" )#SST
plotCounts(adds, gene="ENSMUSG00000010175", col= "lightgreen", xlab = "Condition", main = "Prox1" )#prox1
plotCounts(adds, gene="ENSMUSG00000033726", col= "lightgreen", xlab = "Condition", main = "Emx1" )#Emx1



#Combined plotCounts
Counts <- read.delim("CountReads_csv.csv", header = TRUE, row.names = 1, sep = ",")
Counts <- Counts[which(rowSums(Counts)>10),]
condition <- factor(c("Nova1", "Nova1", "Nova1", "Nova1", "Nova1", "Nova1", "Nova1" ,"Nova2", "Nova2", "Nova2", "Nova2", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO" ,"ECS_2KO" ,"ECS_2KO" ,"ECS_2KO" ,"ECS_2KO" ,"ECS_2KO", "ECS_2KO" ,"ECS_2KO", "ECS_2KO", "ECS_2KO", "ECS_2KO", "N1N2", "N1N2", "N1N2", "N1N2", "N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1", "Nova1", "Nova1",  "Nova2", "Nova2", "Nova2", "Nova2", "N1N2", "N1N2", "N1N2", "N1N2",  "N1N2", "N1N2", "N1N2", "N1N2", "Nova1", "Nova1", "Nova1", "Nova1", "Nova2", "Nova2", "Nova2", "Nova2", "Nova2", "Nova2", "Nova2", "Nova2","ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT", "ECS_WT","ECS_WT", "WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
coldata <- data.frame(row.names= colnames(Counts), condition )

adds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)
adds<- DESeq(adds)
par(mfrow=c(2,4))
plotCounts(adds, gene="ENSMUSG00000021047", col= "lightblue", xlab = "Condition", main = "Nova1")#Nova1
plotCounts(adds, gene="ENSMUSG00000030411", col= "lightblue", xlab = "Condition", main = "Nova2")#Nova2
plotCounts(adds, gene="ENSMUSG00000002718", col= "lightblue", xlab = "Condition", main = "Cse1l")#Nova2
plotCounts(adds, gene="ENSMUSG00000005973", col= "lightblue", xlab = "Condition", main = "Rcn1")#Nova2
plotCounts(adds, gene="ENSMUSG00000004366", col= "lightgreen", xlab = "Condition", main = "SST" )#SST
plotCounts(adds, gene="ENSMUSG00000010175", col= "lightgreen", xlab = "Condition", main = "Prox1" )#prox1
plotCounts(adds, gene="ENSMUSG00000033726", col= "lightgreen", xlab = "Condition", main = "Emx1" )#Emx1
plotCounts(adds, gene="ENSMUSG00000019772", col= "lightblue", xlab = "Condition", main = "Vip")#VIP

boxplot(assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021047"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000030411"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000004366"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000010175"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000033726"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000019772"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021250"], outline= FALSE, main = "Normlized Count", names= c("Nova1","Nova2","Sst", "prox1", "Emx1", "Vip", "Fos"), ylab= "Normlized Count" , col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00") )
#Without SST
boxplot(assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021047"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000030411"],  assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000010175"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000033726"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000019772"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021250"], outline= FALSE, main = "Normlized Count", names= c("Nova1","Nova2", "prox1", "Emx1", "Vip", "Fos"), ylab= "Normlized Count" , col=c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2") )
boxplot(assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000004366"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000010175"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000033726"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000019772"], outline= FALSE, main = "Normlized Count", names= c("Sst", "prox1", "Emx1", "Vip"), ylab= "Normlized Count" , col=c("#56B4E9","#009E73","#F0E442","#0072B2") )

#Scalling
dd=cbind(assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021047"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000030411"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000004366"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000010175"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000033726"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000019772"])
dd= as.data.frame(dd)
colnames(dd)= c("Sst", "Nova1","Nova2", "prox1", "Emx1", "Vip")
dd$Sst= rescale(dd$Sst)
dd$Nova1= rescale(dd$Nova1)
dd$Nova2= rescale(dd$Nova2)
dd$prox1= rescale(dd$prox1)
dd$Emx1= rescale(dd$Emx1)
dd$Vip= rescale(dd$Vip)

dd= melt(dd)
box_plot <- ggplot(dd, aes(x = variable, y = value, color = variable))
# Add the geometric object box plot
box_plot +
  geom_boxplot()


#ECS
dd=cbind(assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000004366"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000021250"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000038418"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000022602"], assays(adds)[[1]][rownames(assays(adds)[[1]]) %in% "ENSMUSG00000019772"])
dd= as.data.frame(dd)
colnames(dd)= c("Sst", "Fos","Egr1","Arc" , "Vip")
boxplot(dd)
dd$Sst= rescale(dd$Sst)
dd$Fos= rescale(dd$Fos)
dd$Egr1= rescale(dd$Egr1)

dd$Arc= rescale(dd$Arc)

dd$Vip= rescale(dd$Vip)


dd= melt(dd)
box_plot <- ggplot(dd, aes(x = variable, y = value, color = variable))
# Add the geometric object box plot
box_plot +
  geom_boxplot()

