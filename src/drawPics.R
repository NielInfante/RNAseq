# Pics for Linda

library(tidyverse)
library(RColorBrewer)

setwd('~/depot/projects/Davis/Nov_2018/')


#d <- read.table('FromJim/normalized_counts.csv', header=T)

# Normalized counts from DESeq
d <- read_csv('FromJim/normalized_counts.csv')

head(d)



bmiSigList <- c('ENSG00000183148','ENSG00000172014','ENSG00000011132','ENSG00000130612','ENSG00000215374','ENSG00000171049','ENSG00000260026','ENSG00000234111',
'ENSG00000214237','ENSG00000224698','ENSG00000240445','ENSG00000251914','ENSG00000252945','ENSG00000254497','ENSG00000259048','ENSG00000275928',
'ENSG00000265043','ENSG00000259057','ENSG00000263154','ENSG00000253153','ENSG00000236508','ENSG00000254707','ENSG00000276107','ENSG00000235058','ENSG00000180974','ENSG00000069011','ENSG00000205176')

dat <- filter(d, GeneID %in% bmiSigList )
row.names(dat) <- dat$GeneID
dat$GeneID <- NULL


select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(dat))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(Sample))

name <- paste(outDir, '/', outPrefix, '_heatmap.png', sep="") 
png(name)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(6, 6), density.info = 'none')
dev.off()


###########  Cluster   ###########

name <- paste(outDir, '/', outPrefix, '_cluster.png', sep="") 
png(name)
plot(hclust(dist(t(assay(vsd)))), label=with(colData(dds), paste0(Sample,':',Age_Day)), main=outPrefix, xlab='', sub='')
dev.off()


###### Heatmap

d <- as.data.frame(assay(vsd))
names(d) <- paste(colData(vsd)$Status, colData(vsd)$Location, sep=":")

d$Gene <- mapIds(orgDB, keys=row.names(d), column='SYMBOL', keytype='ENSEMBL', multiVals='first')
d$ID <- row.names(d)

# Use Gene ID in place where there is no gene name
idx <- is.na(d$Gene)
d$Gene[idx] <- d$ID[idx]
idx <- which(d$Gene == 'NA')  # mapIDs with my_concat can return "NA", not NA
d$Gene[idx] <- d$ID[idx]

#best <- res[order(abs(res$log2FoldChange),decreasing = T)[1:50],]

m <- d[row.names(d) %in% bmiSigList,]
m <- m[!duplicated(m$Gene),]   # Make sure there are only unique gene names


row.names(m) <- m$Gene
m$Gene <- NULL
m$ID <- NULL
m <- m[order(m[,1], decreasing=T),]
m <- as.matrix(m)

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

name <- paste(outDir, '/', outPrefix, '_heat.png', sep="") 
png(name, height = 650)
#tiff(file=name, width=1500, height=2100, units='px', res=300)
heatmap.2(m, col=hmcol, dendrogram='column', trace='none', margin=c(10,6), density.info='none', Colv=T, Rowv=F)
dev.off()



#### Pheatmap

library(genefilter)
library(DESeq2)
library(tidyverse)
library(pheatmap)

dds <- readRDS('~/depot/projects/Kovinich/Kovinich_2019_01/deseq/H2O_MYB29A2_vs_pGWB2_dds.rds')
rld <- rlogTransformation(dds)

topVarGenes <- head(order(-rowVars(assay(rld))),20)


mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Group","Treatment")])

# Need the same names in both the df and the matrix
rownames(df) <- colData(rld)[,c('SampleID')]
colnames(mat) <- colData(rld)[,c('SampleID')]

pheatmap(mat, annotation_col=df, cutree_rows = 3, cutree_cols = 2)




