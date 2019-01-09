library(dplyr)
library(tximport)
library(rjson)
library(DESeq2)
library(readr)
library(magrittr)
library('AnnotationDbi')
library('calibrate')
library("RColorBrewer")
library(ggplot2)
library(ggrepel)
library(gplots)

my_concat <- function(x){paste(x, sep="|", collapse="|")}

# Organism
library('org.Hs.eg.db')
orgDB <- org.Hs.eg.db

# Input and Output

directory<-"~/depot/projects/Ma/July_2018/"
outDir <- "~/depot/projects/Ma/July_2018/deseq"
outPrefix <- 'Group_Early_Late'

setwd(directory)

# The Stats

PCA_Group <- 'Group'
design =~ Group 
contrast <- c('Group','Early','Late')




tx2gene <- read.table('hg38/IDs', header=T, sep="\t", stringsAsFactors = F)

metadata <- read.table('meta', header = T, sep="\t", stringsAsFactors = T)

meta <- metadata %>% dplyr::filter(Group %in% c('Late', 'Early'))




timePoint <- c('0','120d','14d','30d','3d','3h','60d','6d','8h','9d','wo')
for(idx1 in 1:10){
	for (idx2 in (idx1+1):11){
		G1 <- timePoint[idx1]
		G2 <- timePoint[idx2]
		print(paste("G1 is",G1,"and G2 is ",G2))
		meta <- metadata %>% filter(Time %in% c(G1,G2))
		
		contrast <- c('Time',G1,G2)
		outPrefix <- paste0("Time_",G1,"_",G2)
		
		
		doItAll()		
		
	}
}



#G1 <- '0'
#G2 <- '120d'


doItAll <- function(){


#meta <- metadata %>% filter(Time == '0' | Time =='120d')
#meta <- metadata
meta$ID <- meta$Sample
samples <- meta$ID

files <- paste0(directory, 'salmon/', samples, '/quant.sf')

txi <- tximport(files, type='salmon', tx2gene = tx2gene)

dds <- DESeqDataSetFromTximport(txi, meta, design)
#dds$Type %<>% relevel('C')


vsd <- vst(dds, blind=F)

# How many genes, out of those with at least a single count, have three samples with a count of 10 or more
dds <- dds[rowSums(counts(dds)) > 0,]
keep <- rowSums(counts(dds) >= 10) >= 3
#table(keep)
dds <- dds[keep,] # filter them out


dds <- DESeq(dds)

res<-results(dds, contrast=contrast)
res<-res[order(res$padj),]
res <- as.data.frame(res)
head(res)

# Save dds
saveRDS(dds, paste0(outDir, '/', outPrefix, '_dds.rds'))



my_concat <- function(x){paste(x, sep="|", collapse="|")}

# Organism
library('org.Hs.eg.db')
orgDB <- org.Hs.eg.db


# Get gene names
res$Gene <- mapIds(orgDB, keys=row.names(res), column='SYMBOL', keytype='ENSEMBL', multiVals=my_concat)
res$ID <- row.names(res)

# Use Gene ID in place where there is no gene name
idx <- is.na(res$Gene)
res$Gene[idx] <- res$ID[idx]
idx <- which(res$Gene == 'NA')  # mapIDs with my_concat can return "NA", not NA
res$Gene[idx] <- res$ID[idx]


# Write Results
outResults <- data.frame(GeneID=res$ID, Gene=res$Gene, baseMean=res$baseMean, stat=res$stat, log2FoldChange=res$log2FoldChange, pvalue=res$pvalue, padj=res$padj)
name <- paste(outDir, '/', outPrefix, '_results.txt', sep="") 
write.table(outResults, file=name, sep="\t", quote=F, row.names=F)

# Significant genes
r2 <- res[!(is.na(res$padj)),]
resSig <- r2[ r2$padj < 0.05, ]
resTable <- data.frame(GeneID=row.names(resSig), Gene=resSig$Gene, baseMean=resSig$baseMean, stat=resSig$stat, log2FoldChange=resSig$log2FoldChange, pvalue=resSig$pvalue, padj=resSig$padj)
write.table(resTable,file=paste(outDir, "/", outPrefix, "_significant.txt", sep=""), sep="\t", quote=F, row.names=F)



##########  Sanity Check
# Plot counts of most significant, to check if fold change is right
png(paste0(outDir, '/',outPrefix,'_sanity.check.png'))
plotCounts(dds, gene=res[1,]$ID, intgroup = PCA_Group, main=res[1,]$Gene, pch=19)
dev.off()

#########  MA Plot   #########

name <- paste(outDir, '/', outPrefix, '_MAplot.png', sep="") 
png(name)
plotMA(dds)
dev.off()

name <- paste(outDir, '/', outPrefix, '_MA_20.png', sep="") 
png(name)
res[1:20,] %>% ggplot(aes(x=baseMean,y=log2FoldChange)) + geom_point() + 
			geom_text_repel(aes(label=Genus)) + 
			scale_x_log10() + geom_hline(yintercept=0, color='black') +  theme_bw()
dev.off()


#########  Heatmap   #########

select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(vsd)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(Time,ID , sep=" : "))

name <- paste(outDir, '/', outPrefix, '_heatmap.png', sep="") 
png(name)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16), density.info = 'none')
dev.off()


###########  Cluster   ###########

name <- paste(outDir, '/', outPrefix, '_cluster.png', sep="") 
png(name)
plot(hclust(dist(t(assay(vsd)))), label=with(colData(dds), paste(Time,ID, sep=" : ")), main='Dendrogram', xlab='', sub='')
dev.off()



############  PCA    ###########s

name <- paste(outDir, '/', outPrefix, '_PCA.png', sep="") 
png(name)
print(plotPCA(vsd, intgroup=c(PCA_Group)))
dev.off()

#tiff(file=name, width=1800, height=1200, units='px', res=300)

name <- paste(outDir, '/', outPrefix, '_PCA_names.png', sep="") 
png(name)
p <- plotPCA(vsd, intgroup=c(PCA_Group))
p <- p + geom_text_repel(aes_string(x="PC1", y="PC2", label=colData(dds)$ID), point.padding = unit(2,"points"))
print(p)
dev.off()



pcaData <- plotPCA(vsd, intgroup=PCA_Group, returnData=TRUE)
pcaData$f <- colData(dds)$File
pcaData$Location <- colData(dds)$Location
pcaData$BMI <- colData(dds)$BMI
pcaData$Name <- colData(dds)$PS
pcaData$Patient <- colData(dds)$Patient

percentVar <- round(100 * attr(pcaData, "percentVar"))

name <- paste(outDir, '/', outPrefix, '_PCA_shape.png', sep="") 
png(name)
ggplot(pcaData, aes(PC1, PC2, color=Patient, shape=Location)) +
	geom_point(aes(size=BMI)) +
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	geom_text_repel(aes(x=PC1, y=PC2, label=Patient), point.padding = unit(2,"points")) +
		coord_fixed()
dev.off()


########### Heatmap

# This picks the top 50 genes, ranked by absolute fold change, and does a heatmap of the normalized expression

d <- as.data.frame(assay(vsd))
names(d) <- paste(colData(vsd)$Time, colData(vsd)$ID, sep=":")

d$Gene <- mapIds(orgDB, keys=row.names(d), column='SYMBOL', keytype='ENSEMBL', multiVals=my_concat)
d$ID <- row.names(d)

# Use Gene ID in place where there is no gene name
idx <- is.na(d$Gene)
d$Gene[idx] <- d$ID[idx]
idx <- which(d$Gene == 'NA')  # mapIDs with my_concat can return "NA", not NA
d$Gene[idx] <- d$ID[idx]

best <- res[order(abs(res$log2FoldChange),decreasing = T)[1:50],]

m <- d[row.names(d) %in% row.names(best),]
m <- m[!duplicated(m$Gene),]   # Make sure there are only unique gene names


row.names(m) <- m$Gene
m$Gene <- NULL
m$ID <- NULL
m <- m[order(m[,1], decreasing=T),]
m <- as.matrix(m)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

name <- paste(outDir, '/', outPrefix, '_heat.png', sep="") 
png(name, height = 850, width=1200)
#tiff(file=name, width=1500, height=2100, units='px', res=300)
heatmap.2(m, col=hmcol, dendrogram='column', trace='none', margin=c(10,6), density.info='none', Colv=T, Rowv=F)
dev.off()


####  Volcano

name <- paste(outDir, '/', outPrefix, '_volcano.png', sep="") 
png(name)

par(pch = 16)
with(res, plot(log2FoldChange, -log10(pvalue), main = "Volcano plot"))
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), col = "red"))
with(subset(res, abs(log2FoldChange) > 2), points(log2FoldChange, -log10(pvalue),  col = "orange"))

with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), points(log2FoldChange,  -log10(pvalue), col = "green"))

# Add legend
legend("topleft", legend = c("FDR<0.05", "|LFC|>2", "both"), pch = 16, col = c("red", "orange", "green"))

# Label Extra significant points
#with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 1))

# Label all significant
#with(subset(res, padj < 0.05), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 1))

dev.off()

}
