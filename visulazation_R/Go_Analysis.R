if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
 
BiocManager::install("clusterProfiler")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)


genes_to_test <- Nova2_AS_All[Nova2_AS_All$PValue < 0.05,]

# for shared genes: 
#Nova2KO_WT

gene<- Nova2_GE[which(unique(Nova2_GE$gene) %in% unique(Nova2_AS_All$geneSymbol)),] # With only Nova2_GE
genes_to_test <- gene[gene$p_value < 0.05,]
colnames(genes_to_test)[2] = c("ENSEMBL")
GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=30)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "CC")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "MF")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

##NovaDB_WT
gene<- NovaDB_GE[which(unique(NovaDB_GE$gene) %in% unique(NovaDB_AS_All$geneSymbol)),]
genes_to_test <- gene[gene$p_value < 0.05,]
colnames(genes_to_test)[2] = c("ENSEMBL")
GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "BP")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=30)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "CC")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "MF")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)

#Most significant GE OF nova2

gene<- Nova2_GE[ -which(unique(Nova2_GE$gene) %in% unique(Nova2_AS_All$geneSymbol)),]
genes_to_test <- gene[gene$p_value < 0.05,]
colnames(genes_to_test)[2] = c("ENSEMBL")
GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 20 ))
dotplot(GO_results, showCategory=20)


#Most significant GE OF ECS

gene= ECS_GE_All[-which(unique(ECS_GE_All$gene) %in% unique(ECS_AS_All$geneSymbol)),]
#genes_to_test <- gene[gene$PValue < 0.05,]
#genes_to_test <- gene[gene$`log2(fold_change)`>0.1,]
genes_to_test<- gene
colnames(genes_to_test)[2] = c("ENSEMBL")
GO_results <- enrichGO(gene = genes_to_test$ENSEMBL, OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL", ont = "ALL")
head(summary(GO_results))
plot(barplot(GO_results, showCategory = 50 ))
dotplot(GO_results, showCategory=50)

###Other way 
genes_to_test<- gene
colnames(genes_to_test)[3] = c("SYMBOL")
ego <- enrichGO(gene = genes_to_test$SYMBOL, OrgDb = org.Mm.eg.db, ont= "ALL", keyType="SYMBOL")
plot(barplot(ego))
dotplot(ego)
dotplot(ego, showCategory = 15)
########################3

#Add new analysis in GO TERMS The data from GEA_PCA 
GO_terms=as.data.frame(GO_results)
#view(as.data.frame(grep("synap", GO_terms$Description, value = TRUE)))
length(grep("synap", GO_terms$Description, value = TRUE))
CT= c(grep("synap", GO_terms$Description, value = TRUE))

All_Synap_CT= GO_terms[ which(GO_terms$Description%in%CT),]
presynaptic_CT= All_Synap_CT[ which(All_Synap_CT$Description%in% c(grep("presynap", GO_terms$Description, value = TRUE))),]
postsynaptic_CT= All_Synap_CT[ which(All_Synap_CT$Description%in% c(grep("postsynap", GO_terms$Description, value = TRUE))),]
synaptic_CT= All_Synap_CT[  !(All_Synap_CT$Description %in% c(grep( "presynap", GO_terms$Description, value = TRUE))),]
synaptic_CT= synaptic_CT[  !(synaptic_CT$Description %in% c(grep( "postsynap", GO_terms$Description, value = TRUE))),]
