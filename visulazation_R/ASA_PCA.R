library(splitstackshape )
library(matrixStats)

Nova1_AS <- read.delim("Nova1_AS.csv", header = TRUE, row.names = 1, sep = ",")
selected = c ("Acsl3", " Nrbp2", "Trem2", "Nrxn1", "Nrxn3", "Nrxn2", "Syngap1", "Nova1", "Nova2", "Emx1", "Vip", "Pvr", "Prox1", "Prox2", "Sst")
EnhancedVolcano(Nova1_AS, x= "IncLevelDifference", y = "PValue", lab = Nova1_AS$geneSymbol, selectLab = selected)
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]
dim(Nova1_AS)
rowCounts= data.frame(Nova1_AS$IC_SAMPLE_1, Nova1_AS$IC_SAMPLE_2)
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_1", sep=",")
rowCounts=cSplit(rowCounts, "Nova1_AS.IC_SAMPLE_2", sep=",")
rowCounts= cbind(rowCounts, Nova1_AS$geneSymbol)
rownames(rowCounts)= Nova1_AS$GeneID
genelist= rowCounts$V2
rowCounts=rowCounts[,-9]
condition <- factor(c("Nova1", "Nova1", "Nova1", "Nova1","WT_sham", "WT_sham", "WT_sham", "WT_sham"  ))
coldata <- data.frame(row.names= colnames(rowCounts), condition )
mat= data.frame(rowCounts)
mat[is.na(mat)] <- 0
#Calculate Z score 
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score", row_labels = genelist )

#significant genes 
sigs.df2 <- Nova1_AS[(Nova1_AS$PValue <0.05 ) & (abs(Nova1_AS$IncLevelDifference) > 0.1),]
mat= data.frame(rowCounts)
rownames(mat)=rownames(rowCounts)
mat[is.na(mat)] <- 0
mat<- mat[rownames(sigs.df2),]
mat= na.omit(mat)
mat.z=(mat-rowMeans(mat))/(rowSds(as.matrix(mat)))[row(mat)]
Heatmap(mat.z, cluster_rows = T, cluster_columns = T,  column_labels = colnames(mat.z), name = "Z-score" )
pheatmap(mat.z,labels_row =  sigs.df2$geneSymbol)

#ECS
Nova_AS <- read.delim("ECS_AS.csv", header = TRUE, row.names = 1, sep = ",")

#NovaDB
Nova_AS <- read.delim("NovaDB_AS.csv", header = TRUE, row.names = 1, sep = ",")

#ECS_DBL
Nova_AS <- read.delim("ECS_DB_AS.csv", header = TRUE, row.names = 1, sep = ",")
Nova1_AS=Nova1_AS[which(duplicated(Nova1_AS$GeneID)==FALSE),]
sigs.df <- sigs.df[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]
Nova1_AS <- Nova1_AS[which(sigs.df$symbol %in% Nova1_AS$geneSymbol),]
Combined_table<- data.frame(sigs.df$symbol,sigs.df$log2FoldChange, sigs.df$pvalue, Nova1_AS$IncLevelDifference, Nova1_AS$PValue)
rownames(Combined_table) <- rownames(sigs.df)

High<-Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) > 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) > 0.04),]
Low<- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) < 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) < 0.04),]
Diff_GE <- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) > 0.2 ) & (abs(Combined_table$Nova1_AS.IncLevelDifference) < 0.04),]
Diff_AS <- Combined_table[(abs(Combined_table$sigs.df.log2FoldChange) < 0.2) & (abs(Combined_table$Nova1_AS.IncLevelDifference) > 0.04),]

Scalled_table= Combined_table

aa= Scalled_table[, c(2,4)]
aa$GeneSymbol= Scalled_table$sigs.df.symbol
colnames(aa)= c("FC_GE", "FC_AS", "GeneSymbol")
library(reshape2)
reshape  <- melt(aa, id.var = "GeneSymbol")
library(ggplot2)
ggplot(reshape, aes(x = GeneSymbol, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 

Scalled_table= Combined_table
Scalled_table=Scalled_table[1:35,]

library(reshape2)
reshape  <- melt(aa, id.var = "GeneSymbol")
library(ggplot2)
ggplot(reshape, aes(x = GeneSymbol, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") 

#Sort in high to low 
ggplot(reshape, 
       aes(x=reorder(GeneSymbol, value),y= value), 
       fill = variable) + 
  geom_bar(stat = "identity", position = "dodge")  +
  ggtitle("Total ")

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable)) + 
  geom_bar(stat = "identity", position = "dodge")

reshape %>%
  mutate(class = fct_reorder(GeneSymbol, value)) %>%
  ggplot( aes(x=reorder(GeneSymbol, value), y=value, fill=variable)) + 
  geom_col()

ggplot(reshape, aes(x = reorder(GeneSymbol, -value), y = value, fill=variable)) + geom_bar(stat = "identity")

require(dplyr)
require(forcats)

reshape %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()


p= reshape %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()

require(scales) # for removing scientific notation
p <- p + scale_y_continuous(labels = comma)
# manually generate breaks/labels
labels <- seq(1901, 2000, length.out=10)
# and set breaks and labels
p <- p + scale_x_discrete(breaks=labels, labels=as.character(labels))

library(tidyverse)

ggplot(reshape %>% arrange(GeneSymbol, desc(value)) %>%
         mutate(GeneSymbol=factor(GeneSymbol, levels=unique(GeneSymbol))), 
       aes(x=GeneSymbol,y=value, fill = variable)) +
  geom_bar(stat="identity", show.legend=TRUE) +
  
  facet_grid(. ~ variable, scales="free_x", space="free_x") +
  scale_y_continuous(limits=c(-0.005, 1.05*max(reshape$value)), expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA))


p= ggplot(reshape %>% arrange(GeneSymbol, desc(value)) %>%
            mutate(GeneSymbol=factor(GeneSymbol, levels=unique(GeneSymbol))), 
          aes(x=GeneSymbol,y=value, fill = variable)) +
  geom_bar(stat="identity", show.legend=TRUE) +
  
  facet_grid(. ~ variable, scales="free_x", space="free_x") +
  scale_y_continuous(limits=c(-0.005, 1.05*max(reshape$value)), expand=c(0,0)) +
  theme_classic() +
  theme(panel.spacing=unit(0,"pt"), 
        panel.border=element_rect(colour="grey50", fill=NA))

require(scales) # for removing scientific notation
p <- p + scale_y_continuous(labels = comma)
# manually generate breaks/labels
labels <- seq(1901, 2000, length.out=10)
# and set breaks and labels
p <- p + scale_x_discrete(breaks=labels, labels=as.character(labels))
p

SynapticGenes= c( "Sema6a", "Sorbs1", "Hivep2", "Gnb1", "Cacna1c", "Wnk1", "Kif21a", "Birc6", "Rph3a", "Rogdi", "Dclk1", "Myo16", "Nrxn1", "Macf1", "Nrcam", "Ank2", "Stx3", "Hsf1", "Dnajc5", "Vars", "Hspa8" , "Sema6d", "Trpm7", "Arhgap26", "Fchsd2", "Cadm1", "Cask", "Tspyl2", "Nrxn3", "Tpd52l2", "Ank3", "Cacna1b", "Add1", "Rims1", "Smad2", "Pml", "Atf2", "Sptan1", "Mapk9", "Adam22", "Macf1", "Nrcam", "Gramd1a", "Nckap1", "Sh3glb1","Atxn2", "Bin1", "Hspa8", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8" , "Nrxn1", "Macf1", "Ercc5", "Setx", "Macf1", "Ercc5", "Setx", "Dnajc5", "Hspa8","Sema6d", "Trpm7", "Arhgap26", "Fchsd2", "Cadm1", "Cask", "Tspyl2", "Nrxn3", "Tpd52l2", "Ank3", "Cacna1b", "Add1", "Rims1", "Pml", "Atf2", "Sptan1", "Adam22", "Macf1", "Nrcam", "Gramd1a", "Sh3glb1","Atxn2", "Bin1", "Hspa8", "Vip", "SsT" )
SynapticGenes= unique(SynapticGenes)
AS_Genes=Combined_table$sigs.df.symbol[which(Combined_table$Nova1_AS.IncLevelDifference >0.1)]
AS_Genes= c(AS_Genes,SynapticGenes )
rr=reshape[which(reshape$GeneSymbol%in%AS_Genes),]
rr %>% 
  mutate(ordering = -as.numeric(variable) + value,
         GeneSymbol = fct_reorder(GeneSymbol, ordering, .desc = T)) %>% 
  ggplot(aes(GeneSymbol, value, fill = variable)) + geom_col()

#Barplot
GE= rr[which(rr$variable=="FC_GE"),]
AS= rr[which(rr$variable=="FC_AS"),]
rr2=rbind(GE[1,], AS[1,])
dd= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[2,], AS[2,])
dd2= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[3,], AS[3,])
dd3= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))


rr2=rbind(GE[101,], AS[101,])
dd4= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[30,], AS[30,])
dd5= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[67,], AS[67,])
dd6= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[18,], AS[18,])
dd7= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[19,], AS[19,])
dd8= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))

rr2=rbind(GE[54,], AS[54,])
dd9= ggplot(data=rr2, aes(x=GeneSymbol, y=value, fill=variable, color=variable, alpha=variable)) +
  geom_bar(stat="identity", position='dodge') +
  scale_colour_manual(values=c("lightblue4", "red")) +
  scale_fill_manual(values=c("lightblue", "pink")) +
  scale_alpha_manual(values=c(.3, .8))


library("ggpubr")

ggarrange(dd, dd2, dd3, dd4, dd5, dd6, dd7, dd8, dd9,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          ncol = 3, nrow = 3)

#Corrlation Pvalue
Nova1_coor= sigs.df

corrplot(cor(as.data.frame(Nova1_coor[1:12]), as.data.frame(Nova2_coor[1:12])),
         method = "number",
         type = "upper" # show only upper side
)
corrplot(cor(Cor_Events))
