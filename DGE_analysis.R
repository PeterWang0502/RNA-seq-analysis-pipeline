#install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

library(limma)
library(edgeR)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(ggrepel)
library(tibble)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(tximport)
#install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)
#package.version('dbplyr')
#package.version('BiocFileCache')

mart = useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", host = "jul2023.archive.ensembl.org") # "may2017.archive.ensembl.org"
ens = getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = mart)

# tximport on RSEM output
dir <- '/Users/peter/Desktop/rsem_results'
files <- file.path(dir, list.files(dir))
names(files) <- lapply(files, function(x) substr(strsplit(x, '\\.')[[1]][1], 54, 100))
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
write.table(txi.rsem$counts, file="tximport-count-combined.txt", quote=F, sep='\t')

### Analysis starts

setwd("/Users/peter/Desktop/") # input your pipeline output directory 
count = read.delim("Quant.rsem.genes.age.matrix", stringsAsFactors = F, row.names = 1) # input your output file name here
rownames(count) = sub("\\.\\d+", "", rownames(count))
count = round(data.matrix(count))
meta = read.delim("MM_meta_age.txt", stringsAsFactors = F) # input your metadata here

meta = meta[meta$group %in% c('MA_OLD','MA_YOU'),] # input comparison group here
count = count[,colnames(count) %in% meta$ID]

genelist = DGEList(count, samples = meta)
plotMDS(genelist, col = brewer.pal(9,"Set1")[as.numeric(as.factor(meta$group))])
keep  = rowSums(cpm(genelist) >= 5) > ncol(genelist)/2
table(keep)
genelist.filtered = genelist[keep, keep.lib.sizes = FALSE]
genelist.norm = calcNormFactors(genelist.filtered, method = "TMM")
lcpm = cpm(genelist.norm, log = T)
mean(lcpm)
par(mar = c(10, 4, 2, 2) + 0.1)
boxplot(lcpm, ylab="Log2 counts per million",las=2)
abline(h=median(lcpm),col="red")
gene.var <- apply(lcpm, 1, var)
top.var <- names(sort(gene.var, decreasing=TRUE))[1:500]
highly_variable_lcpm <- lcpm[top.var,]
pheatmap(highly_variable_lcpm, show_rownames=F, cluster_cols=T, cluster_rows=T, scale="row")
pca <- prcomp(t(lcpm), scale = F)
summary(pca)
pca.proportionvariances <- ((pca$sdev^2) / (sum(pca$sdev^2)))*100
barplot(pca.proportionvariances, cex.names=1, xlab=paste("PC, 1-", length(pca$sdev)), ylab="Variation Explained (%)", main="Scree plot", ylim=c(0,100))

# EdgeR LRT
row.names(meta) = meta[,1]
pca.merge = merge(pca$x[,1:2], meta, by = "row.names")
qplot(PC1, PC2, data = pca.merge, colour = group, asp = 1)
design<-model.matrix(~0+meta$group)
colnames(design)<-levels(as.factor(meta$group))
genelist.disp = estimateDisp(genelist.norm, design, robust = TRUE)
fit <- glmFit(genelist.disp, design, robust=TRUE)
lrt <- glmLRT(fit, contrast = makeContrasts(MA_OLD-FE_OLD, levels = design))
resLRT = topTags(lrt, n = Inf)$table
sum(resLRT$FDR < 0.05, na.rm=TRUE)
resLRT = merge(resLRT, ens, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
resLRT = resLRT[order(resLRT$FDR),]
resLRT = rename(resLRT, Row.names = "geneID")
write.table(resLRT, "de_table_edger_lrt.txt", quote = F, sep = '\t', row.names = F)
rownames(resLRT) = resLRT$geneID
resLRT$geneID = NULL
resLRT$threshold = ifelse(abs(resLRT$logFC)>2 & resLRT$FDR<0.01, ifelse(abs(resLRT$logFC)>3 & resLRT$FDR<0.001, "A", "B"), "C")
table(resLRT$threshold)
volcanoData = cbind.data.frame(resLRT$logFC, -log10(resLRT$FDR), resLRT$external_gene_name, resLRT$threshold)
colnames(volcanoData) = c("logFC","negLogPval","annotation", "threshold")
ggplot(volcanoData, aes(x = logFC, y = negLogPval)) + ggtitle("Volcano Plot") +
  geom_point(aes(color = threshold), size = 2) + 
  scale_color_manual(values = c("A" = "red", "B" = "orange", "C" = "black")) +
  geom_text_repel(aes(label = ifelse(threshold == "A", as.character(annotation), '')), hjust = 0, vjust = 0)

# EdgeR QLF
qlfit = glmQLFit(genelist.disp, design, robust = TRUE)
qlf = glmQLFTest(qlfit, contrast = makeContrasts(MA_OLD-FE_OLD, levels = design))
resQLF = topTags(qlf, n = Inf)$table
sum(resQLF$FDR < 0.05, na.rm=TRUE)
resQLF = merge(resQLF, ens, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
resQLF = resQLF[order(resQLF$FDR),]
resQLF = rename(resQLF, Row.names = "geneID")
write.table(resQLF, "de_table_edger_qlf.txt", quote = F, sep = '\t', row.names = F)
rownames(resQLF) = resQLF$geneID
resQLF$geneID = NULL
resQLF$threshold = ifelse(abs(resQLF$logFC)>2 & resQLF$FDR<0.01, ifelse(abs(resQLF$logFC)>3 & resQLF$FDR<0.001, "A", "B"), "C")
table(resQLF$threshold)
volcanoData = cbind.data.frame(resQLF$logFC, -log10(resQLF$FDR), resQLF$external_gene_name, resQLF$threshold)
colnames(volcanoData) = c("logFC","negLogPval","annotation", "threshold")
ggplot(volcanoData, aes(x = logFC, y = negLogPval)) + ggtitle("Volcano Plot") +
  geom_point(aes(color = threshold), size = 2) + 
  scale_color_manual(values = c("A" = "red", "B" = "orange", "C" = "black")) +
  geom_text_repel(aes(label = ifelse(threshold == "A", as.character(annotation), '')), hjust = 0, vjust = 0)

# DeSeq2
meta.level = meta
meta.level$group = factor(meta.level$group, levels=c("MA_OLD","MA_YOU")) # input groups of comparison here
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta.level, design = ~group)
dds <- DESeq(dds)
rld <- rlog(dds, blind=TRUE)

#pdf(file="testPCA.pdf", height=15, width=15)
plotPCA(rld, intgroup=c("group")) + # before normalization
  geom_text(label=meta.level$ID, check_overlap=F, hjust="inward", vjust="inward")
#dev.off()

rld_cor <- cor(assay(rld))
pheatmap(rld_cor)
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:100]
pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=F, cluster_cols=T, scale = "row")
resultsNames(dds)
plotDispEsts(dds)
resLFC <- lfcShrink(dds, coef = 2, type="apeglm")
summary(resLFC)
resLFC <- resLFC[order(resLFC$padj),]
sum(resLFC$padj < 0.05, na.rm=TRUE)
plotCounts(dds, gene=which.min(resLFC$padj), intgroup="group", pch = 16)
plotMA(resLFC, ylim=c(-2,2))
hist(resLFC$pvalue[resLFC$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white")
res = as.data.frame(resLFC) #[1:3000,]
res = res[!is.na(res$padj),]
res = merge(res, ens, by.x = "row.names", by.y = "ensembl_gene_id", all.x = T, sort = F)
colnames(res)[1] = 'geneID'
write.table(res, "de_table_deseq2.txt", quote = F, sep = '\t', row.names = F)
res = column_to_rownames(res, var = "geneID")
res$threshold = ifelse(abs(res$log2FoldChange)>2 & res$padj<0.01 , ifelse(abs(res$log2FoldChange)>3 & res$padj<0.001, "A", "B"), "C")
table(res$threshold)
volcanoData = cbind.data.frame(res$log2FoldChange, -log10(res$padj), res$external_gene_name, res$threshold)
colnames(volcanoData) = c("logFC","negLogPval","annotation", "threshold")
rownames(volcanoData) = rownames(res)
volcanoData$label = ifelse(volcanoData$threshold=='A', volcanoData$annotation, "")

ggplot(volcanoData, aes(x = logFC, y = negLogPval)) + ggtitle("Volcano Plot") +
  geom_point(aes(color = threshold), size = 2.5) + 
  scale_color_manual(values = c("A" = "red", "B" = "orange", "C" = "green", "D" = "black")) +
  geom_text(label = volcanoData$label, nudge_x = 0.25, nudge_y = 0.25, check_overlap = F, vjust="inward",hjust='inward')

rld <- rlog(dds, blind=FALSE)
de <- rownames(res[res$padj<0.05 & !is.na(res$padj),])
pheatmap(assay(rld)[de,], cluster_rows=T, show_rownames=F, cluster_cols=T, scale = "row")

gene_lrt = rownames(resLRT)[resLRT$FDR < 0.05]
gene_qlf = rownames(resQLF)[resQLF$FDR < 0.05]
gene_deseq2 = rownames(res)[res$padj < 0.05]
deseq_red = rownames(volcanoData)[volcanoData$threshold=="A"]

paste(sum(deseq_red %in% gene_lrt)/length(deseq_red),sum(deseq_red %in% gene_qlf)/length(deseq_red),sum(deseq_red %in% gene_deseq2)/length(deseq_red))

a1 = length(gene_lrt)
a2 = length(gene_qlf)
a3 = length(gene_deseq2)
n12 = table(gene_lrt %in% gene_qlf); n12=n12[rownames(n12)=='TRUE']
n13 = table(gene_lrt %in% gene_deseq2); n13=n13[rownames(n13)=='TRUE']
n23 = table(gene_qlf %in% gene_deseq2); n23=n23[rownames(n23)=='TRUE']
n123 = table((gene_lrt %in% gene_qlf) & (gene_lrt %in% gene_deseq2)); n123=n123[rownames(n123)=='TRUE']
grid.newpage()
draw.triple.venn(area1 = a1, area2 = a2, area3 = a3, n12 = n12, n23 = n23, n13 = n13, n123 = n123, 
                 category = c("EdgeR LRT", "EdgeR QLF", "Deseq2"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))