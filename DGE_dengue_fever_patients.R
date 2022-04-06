##Loading the required package
library(GEOquery)
library(tidyverse)
library(limma)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(EnhancedVolcano)
library(colorspace)
library(matrixStats)
library(biomaRt)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(pathview)
library(ReactomePA)
library(clusterProfiler)
library(RColorBrewer)


##Loading the GEO dataset (GEO ID: GDS5093) for differential gene expression analysis.
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

##Displaying the downloaded table containing quantification values of genes in samples.
X <- Table(gse)

##Conversion of the dataset into expression set object.  
eset <- GDS2eSet(gse, do.log2=TRUE)

##Conversion of the expression set object into table containing phenotypic information of the samples.
pData <- pData(eset)
pData$disease.state <- gsub(" ", "_", pData$disease.state)
head(pData)

##Retrieving the gene annotation
fData <- fData(eset)
fData
##Retrieving the expression data
exprs <- exprs(eset)
exprs

##Exprs get the expression levels as a data frame and get the distribution
summary(exprs(eset))

#Task 1: summarizing the dataset and making boxplot for the distribution of data in each sample.
##Getting summary of the dataset
summary(eset)
##Transformation of the exprs values and plotting a boxplot
exprs(eset) <- log2(exprs(eset))
boxplot(exprs(eset),outline=FALSE)

##Setting row names as gene symbols.
X <- Table(gse)
geneNames <- as.character(X$IDENTIFIER)
X <- exprs(eset)
rownames(X) <- geneNames

##Detection of the duplicate row names.
dim(X)
length(which(duplicated(rownames(X))))
##Removing duplicate row names.
X <- avereps(X)
dim(X)

#Task 2: Performing hierarchical clustering analysis 
##Condensing the data by transformation
X.t <- t(X)
##Calculating the distance and plotting the Hierarchical Clustering among the samples.
dist_samples <- dist(X.t, method = "euclidean")
plot(dist_samples)
samples_hclust1 <- hclust(dist_samples, method = "complete")
samples_hclust2 <- hclust(dist_samples, method = "average")
plot(samples_hclust1)
plot(samples_hclust2)

#Task 3: Plotting the dendogram with dendextend package and annotating the leaves.
dend <- as.dendrogram(samples_hclust1)
dend <- rotate(dend, 1:56)
dend <- color_branches(dend, k=4)
##Labelling the dendogram with the phenotype labels
labels_colors(dend) <-
  rainbow_hcl(4)[sort_levels_values(
    as.numeric(pData[,3])[order.dendrogram(dend)]
  )]
labels(dend) <- paste(as.character(pData[,3])[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")
dend <- hang.dendrogram(dend,hang_height=0.1)

##Plotting the dendogram
par(mar = c(3,3,3,7))
plot(dend, 
     horiz =  TRUE,  nodePar = list(cex = .001))

#Task 4: Plotting a dendogram of top 100 genes in the dataset based on variance.
##calculating variance and identification of top 100 genes
all_X <- apply(X, MARGIN=1, FUN=var, na.rm=TRUE)
all_X <- as.data.frame(all_X)
ID <- rownames(all_X)
rownames(all_X) <- NULL
all_X <- cbind(ID, all_X)
IDs <- all_X[order(all_X$all_X, decreasing = TRUE),]
IDs <- IDs[1:100,]

##Filtering ranked top 100 genes from the dataset.
X <- as.data.frame(X)
top100 <- subset(X, row.names(X) %in% IDs$ID)

##Plotting a dendogram of top 100 genes
dist100 <- dist(top100)
hclust <- hclust(dist100, method = "complete")
dend <- as.dendrogram(hclust)
plot(dend)

#Task 5: Heatmap of the top 100 genes from the dataset
## Get the rows corresponding to ids_of_interest and all columns
pheatmap(top100, labels_col = pData$disease.state)


#Task 6: Perfomring PCA of samples in the dataset and plotting a bar plot of the PCs percentage
##Transposing the matrix
pca <- prcomp(t(exprs(eset)), scale = T)
summary(pca)
pca_var <- pca$sdev^2
pca_var_perc <- round(pca_var/sum(pca_var) * 100, 1)
barplot(pca_var_perc, main = "Variation Plot", xlab = "PCs", ylab = "Percentage Variance", ylim = c(0, 100))

#Task 7: Plotting PCA plot of the samples in dataset
cbind(pData, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=disease.state, 
             label=paste(disease.state))) + 
  geom_point()

#Task 8: Plotting biplot
biplot(pca)

#Task 9: Producing the PCA scores plot without scaling the data and/or without transposing it first.
pcaw <- prcomp(exprs(eset), scale = T)
cbind(pData, pcaw$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=disease.state, 
             label=paste(disease.state))) + 
  geom_point()

#Task 10: Differential gene expression analysis
##Design for the study
design <- model.matrix(~0 + disease.state, pData)
design

##Grouping the the expression based on disease status
fit1 <- lmFit(X, design)
head(fit1$coefficients)
##Making contrasts between the infection status and controls 
contrasts1 <- makeContrasts(disease.stateConvalescent - disease.statehealthy_control, levels=design)
contrasts2 <- makeContrasts(disease.stateDengue_Fever - disease.statehealthy_control, levels=design)
contrasts3 <- makeContrasts(disease.stateDengue_Hemorrhagic_Fever - disease.statehealthy_control, levels=design)
contrasts4 <- makeContrasts(disease.stateConvalescent - disease.stateDengue_Fever, levels=design)
contrasts5 <- makeContrasts(disease.stateConvalescent - disease.stateDengue_Hemorrhagic_Fever, levels=design)
contrasts6 <- makeContrasts(disease.stateDengue_Fever - disease.stateDengue_Hemorrhagic_Fever, levels=design)
fit2_contrasts1 <- contrasts.fit(fit1, contrasts1)
fit2_contrasts2 <- contrasts.fit(fit1, contrasts2)
fit2_contrasts3 <- contrasts.fit(fit1, contrasts3)
fit2_contrasts4 <- contrasts.fit(fit1, contrasts4)
fit2_contrasts5 <- contrasts.fit(fit1, contrasts5)
fit2_contrasts6 <- contrasts.fit(fit1, contrasts6)
fit2_contrasts1 <- eBayes(fit2_contrasts1)
fit2_contrasts2 <- eBayes(fit2_contrasts2)
fit2_contrasts3 <- eBayes(fit2_contrasts3)
fit2_contrasts4 <- eBayes(fit2_contrasts4)
fit2_contrasts5 <- eBayes(fit2_contrasts5)
fit2_contrasts6 <- eBayes(fit2_contrasts6)

##Header of the table
topTable(fit2_contrasts1)
topTable(fit2_contrasts2)
topTable(fit2_contrasts3)
topTable(fit2_contrasts4)
topTable(fit2_contrasts5)
topTable(fit2_contrasts6)

##Creating as full results
full_results1 <- topTable(fit2_contrasts1, number=Inf)
full_results2 <- topTable(fit2_contrasts2, number=Inf)
full_results3 <- topTable(fit2_contrasts3, number=Inf)
full_results4 <- topTable(fit2_contrasts4, number=Inf)
full_results5 <- topTable(fit2_contrasts5, number=Inf)
full_results6 <- topTable(fit2_contrasts6, number=Inf)

##Exporting the statistically significant DEGs.
full_results1s <- subset(full_results1, adj.P.Val < 0.05 & abs(logFC) > 0.5)
full_results2s <- subset(full_results2, adj.P.Val < 0.05 & abs(logFC) > 0.5)
full_results3s <- subset(full_results3, adj.P.Val < 0.05 & abs(logFC) > 0.5)
full_results4s <- subset(full_results4, adj.P.Val < 0.05 & abs(logFC) > 0.5)
full_results5s <- subset(full_results5, adj.P.Val < 0.05 & abs(logFC) > 0.5)
full_results6s <- subset(full_results6, adj.P.Val < 0.05 & abs(logFC) > 0.5)

write.table(full_results1s, "C:/Users/name/Desktop/results1.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_results2s, "C:/Users/name/Desktop/results2.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_results3s, "C:/Users/name/Desktop/results3.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_results4s, "C:/Users/name/Desktop/results4.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_results5s, "C:/Users/name/Desktop/results5.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(full_results6s, "C:/Users/name/Desktop/results6.csv",
            sep = "\t", quote = FALSE, row.names = FALSE)

##Making volcanoplots
EnhancedVolcano(full_results1, lab = rownames(full_results1), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Convalescent vs healthy_control')

EnhancedVolcano(full_results2, lab = rownames(full_results2), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Dengue_Fever vs healthy_control')

EnhancedVolcano(full_results3, lab = rownames(full_results3), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Dengue_Hemorrhagic_Fever vs healthy_control')

EnhancedVolcano(full_results4, lab = rownames(full_results4), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Convalescent vs Dengue_Fever')

EnhancedVolcano(full_results5, lab = rownames(full_results5), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Convalescent vs Dengue_Hemorrhagic_Fever')

EnhancedVolcano(full_results6, lab = rownames(full_results6), x = 'logFC', y = 'adj.P.Val', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','logFC','adj.P',
                               'adj.P & logFC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, subtitle = 'Dengue_Fever vs Dengue_Hemorrhagic_Fever')

##Gene ID conversion of the result files
ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
DEGs1 <- getBM(attributes=c('hgnc_symbol', 'hgnc_symbol','entrezgene_id'), 
             filters = 'hgnc_symbol', values = row.names(full_results1), 
             mart = ensembl)

DEGs2 <- getBM(attributes=c('hgnc_symbol', 'hgnc_symbol','entrezgene_id'), 
               filters = 'hgnc_symbol', values = row.names(full_results2), 
               mart = ensembl)

DEGs3 <- getBM(attributes=c('hgnc_symbol', 'hgnc_symbol','entrezgene_id'), 
               filters = 'hgnc_symbol', values = row.names(full_results3), 
               mart = ensembl)

##Removing duplicate IDs
DEGs1 <- DEGs1[!duplicated(DEGs1$hgnc_symbol), ]
DEGs2 <- DEGs2[!duplicated(DEGs2$hgnc_symbol), ]
DEGs3 <- DEGs3[!duplicated(DEGs3$hgnc_symbol), ]

##Making uniques IDs list
full_results1 <- full_results1[!duplicated(row.names(full_results1)), ]
full_results2 <- full_results2[!duplicated(row.names(full_results2)), ]
full_results3 <- full_results3[!duplicated(row.names(full_results3)), ]

##Renaming the column names
DEGs1 <- DEGs1 %>% rename(ID = "hgnc_symbol")
DEGs2 <- DEGs2 %>% rename(ID = "hgnc_symbol")
DEGs3 <- DEGs3 %>% rename(ID = "hgnc_symbol")

##Joining the IDs with the data 
ID <- rownames(full_results1)
rownames(full_results1) <- NULL
full_results1 <- cbind(ID, full_results1)
full_results1 <- dplyr::inner_join(full_results1, DEGs1, by="ID")
ID <- rownames(full_results2)
rownames(full_results2) <- NULL
full_results2 <- cbind(ID, full_results2)
full_results2 <- dplyr::inner_join(full_results2, DEGs2, by="ID")
ID <- rownames(full_results3)
rownames(full_results3) <- NULL
full_results3 <- cbind(ID, full_results3)
full_results3 <- dplyr::inner_join(full_results3, DEGs3, by="ID")

####################################################
#Making the list ready for gene set enrichment analysis adn filtering for significantly expressed genes.
genelist1 <- cbind(full_results1$entrezgene_id, full_results1$logFC)
colnames(genelist1) <- c("ID", "logFC")
genelist1 <- as.data.frame(genelist1)
geneList1 <- as.numeric(genelist1[,2])
names(geneList1) <- as.character(genelist1[,1])
geneList1 <- sort(geneList1, decreasing = TRUE)
head(geneList1)
gene1 <- names(geneList1)[abs(geneList1) > 0.5]

genelist2 <- cbind(full_results2$entrezgene_id, full_results2$logFC)
colnames(genelist2) <- c("ID", "logFC")
genelist2 <- as.data.frame(genelist2)
geneList2 <- as.numeric(genelist2[,2])
names(geneList2) <- as.character(genelist2[,1])
geneList2 <- sort(geneList2, decreasing = TRUE)
head(geneList2)
gene2 <- names(geneList2)[abs(geneList2) > 0.5]

genelist3 <- cbind(full_results3$entrezgene_id, full_results3$logFC)
colnames(genelist3) <- c("ID", "logFC")
genelist3 <- as.data.frame(genelist3)
geneList3 <- as.numeric(genelist3[,2])
names(geneList3) <- as.character(genelist3[,1])
geneList3 <- sort(geneList3, decreasing = TRUE)
head(geneList3)
gene3 <- names(geneList3)[abs(geneList3) > 0.5]

#GSEA using msigdbr
msigdbr_species()
#Get all human gene sets
m_df <- msigdbr(species = "Homo sapiens")
head(m_df, 2) %>% as.data.frame

#Gene set names (gs_name) and gene identifiers (entrez_gene).
m_df_H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)
head(m_df_H_t2g)
collections <- msigdbr_collections()

#Select (subset) the Hallmark collection
m_df_H <- m_df[m_df$gs_cat=="H",]
head(m_df_H)

#For enrichemnt analysis we need a dataframe with two columns
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

#Enrichment by hypergeometric test implemented in the enricher function
res_enricher1 <- enricher(gene = gene1, TERM2GENE = msigdbr_t2g)
res_enricher2 <- enricher(gene = gene2, TERM2GENE = msigdbr_t2g)
res_enricher3 <- enricher(gene = gene3, TERM2GENE = msigdbr_t2g)

#Gene set enrichment analysis is implemented in the GSEA function
res_GSEA1 <- GSEA(geneList1, TERM2GENE = msigdbr_t2g)
res_GSEA2 <- GSEA(geneList2, TERM2GENE = msigdbr_t2g)
res_GSEA3 <- GSEA(geneList3, TERM2GENE = msigdbr_t2g)

#Convert it into a dataFrame:
res_GSEA_df1 <- as.data.frame(res_GSEA1)
res_GSEA_df2 <- as.data.frame(res_GSEA2)
res_GSEA_df3 <- as.data.frame(res_GSEA3)

#Visualization
dotplot(res_GSEA1, showCategory=20)
dotplot(res_GSEA2, showCategory=20)
dotplot(res_GSEA3, showCategory=20)

##KEGG analysis
kk1 <- enrichKEGG(gene = gene1, organism = 'hsa', pvalueCutoff = 0.05)
kk2 <- enrichKEGG(gene = gene2, organism = 'hsa', pvalueCutoff = 0.05)
kk3 <- enrichKEGG(gene = gene3, organism = 'hsa', pvalueCutoff = 0.05)
kk2_1 <- gseKEGG(geneList     = geneList1,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
kk2_2 <- gseKEGG(geneList     = geneList2,
                 organism     = 'hsa',
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
kk2_3 <- gseKEGG(geneList     = geneList3,
                 organism     = 'hsa',
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

#Reactome analysis
rt1 <- enrichPathway(gene=gene1, pvalueCutoff = 0.05, readable=TRUE)
rt1

rt2 <- enrichPathway(gene=gene2, pvalueCutoff = 0.05, readable=TRUE)
rt2

rt3 <- enrichPathway(gene=gene3, pvalueCutoff = 0.05, readable=TRUE)
rt3

#Plotting dot plots of the top KEGG pathways.
dotplot(kk1, showCategory=10)
dotplot(kk2, showCategory=10)
dotplot(kk3, showCategory=10)








