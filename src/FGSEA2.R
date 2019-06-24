
library(fgsea)


# Trying new pathway stuff

# Trying it with expression levels directly.
# Probably not right.


d <- as.data.frame(assay(vsd))
#names(d) <- paste(colData(vsd)$Patient, colData(vsd)$Status, colData(vsd)$Location, sep=":")

meta <- read.table('meta', header=T)
meta$Patient <- as.factor(meta$Patient)
meta$Num <- as.factor(meta$Num)
meta <- filter(meta, Status == 'ovob')
meta$FN <- translate( paste0('X', meta$File), '-', '.')
d <- dplyr::select(d, meta$FN)
#head(meta)

outPrefix <- "Ovob"

d$RM <- rowMeans(d)

d$Gene <- mapIds(orgDB, keys=row.names(d), column='SYMBOL', keytype='ENSEMBL', multiVals='first')



#head(res)
#res$row <- row.names(res)

#library(org.Hs.eg.db)

#ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
#																		key=res$row, 
#																		columns="ENTREZ",
#																		keytype="ENSEMBL")
#ens2symbol <- as_tibble(ens2symbol)
#ens2symbol

#res <- inner_join(res, ens2symbol, by=c("row"="ENSEMBL"))
#head(res)

#res$Gene <- mapIds(orgDB, keys=row.names(res), column='ENTREZID', keytype='ENSEMBL', multiVals=my_concat)




res2 <-  d %>% 
	dplyr::select(Gene, RM) %>% 
	na.omit() %>% 
	distinct() %>% 
	group_by(Gene) %>% 
	dplyr::summarize(stat=mean(RM))
#res2

ranks <- deframe(res2)
head(ranks, 20)

#doPathways('~/depot/projects/Davis/Nov_2018/path',prefix = 'Status', foldchanges = ranks)


pathways.hallmark <- gmtPathways("data/h.all.v6.2.symbols.gmt")

# Show the first few pathways, and within those, show only the first few genes. 
#pathways.hallmark %>% 
#	head() %>% 
#	lapply(head)


fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=1000)



fgseaResTidy <- fgseaRes %>%
	as_tibble() %>%
	arrange(desc(NES))

fgseaOut <- dplyr::select(fgseaResTidy, pathway, pval, padj, NES, nMoreExtreme, size)

#write.table(fgseaRes, file=paste0(outDir, '/../gsea/', outPrefix, '_hallmark.txt'), sep="\t", quote=F, row.names=F)

write_tsv(fgseaOut, paste0(outDir, '/../gsea/', outPrefix, '_hallmark.txt'), col_names = T)


# Show in a nice table:
#fgseaResTidy %>% 
#	dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
#	arrange(padj) %>% 
#	DT::datatable()


png(name <- paste0(outDir, '/../gsea/', outPrefix, '_hallmark.png'))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill=padj<0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score",
			 title="Hallmark pathways NES from GSEA") + 
	theme_minimal()
dev.off()

#pathways.hallmark %>% 
#	enframe("pathway", "SYMBOL") %>% 
#	unnest() %>% 
#	inner_join(res, by="SYMBOL")



## Do KEGG
kegg <- fgsea(pathways=gmtPathways("data/c2.cp.kegg.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
	as_tibble() %>% 
	arrange(padj)

# Sort it and write
kegg <- kegg %>%
	as_tibble() %>%
	arrange(desc(NES)) %>% 
	dplyr::select(pathway, pval, padj, NES, nMoreExtreme, size)

write_tsv(kegg, paste0(outDir, '/../gsea/', outPrefix, '_kegg.txt'), col_names = T)


png(name <- paste0(outDir, '/../gsea/', outPrefix, '_kegg.png'))
ggplot(kegg, aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill=padj<0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score",
			 title="KEGG pathways NES from GSEA") + 
	theme_minimal()
dev.off()



## Do Biocarta

biocarta <- fgsea(pathways=gmtPathways("data/c2.cp.biocarta.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
	as_tibble() %>% 
	arrange(padj)

biocarta <- biocarta %>%
	as_tibble() %>%
	arrange(desc(NES)) %>% 
	dplyr::select(pathway, pval, padj, NES, nMoreExtreme, size)

write_tsv(biocarta, paste0(outDir, '/../gsea/', outPrefix, '_biocarta.txt'), col_names = T)


png(name <- paste0(outDir, '/../gsea/', outPrefix, '_biocarta.png'))
ggplot(kegg, aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill=padj<0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score",
			 title="Biocarta pathways NES from GSEA") + 
	theme_minimal()
dev.off()


##  Do onco genes

onco <- fgsea(pathways=gmtPathways("data/c6.all.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
	as_tibble() %>% 
	arrange(padj)

onco <- onco %>%
	as_tibble() %>%
	arrange(desc(NES)) %>% 
	dplyr::select(pathway, pval, padj, NES, nMoreExtreme, size)

write_tsv(onco, paste0(outDir, '/../gsea/', outPrefix, '_onco.txt'), col_names = T)


png(name <- paste0(outDir, '/../gsea/', outPrefix, '_onco.png'))
ggplot(kegg, aes(reorder(pathway, NES), NES)) +
	geom_col(aes(fill=padj<0.05)) +
	coord_flip() +
	labs(x="Pathway", y="Normalized Enrichment Score",
			 title="Onco pathways NES from GSEA") + 
	theme_minimal()
dev.off()


## Do GO

GO <- fgsea(pathways=gmtPathways("data/c5.all.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
	as_tibble() %>% 
	arrange(padj)

GO <- GO %>%
	as_tibble() %>%
	arrange(desc(NES)) %>% 
	dplyr::select(pathway, pval, padj, NES, nMoreExtreme, size)

write_tsv(GO, paste0(outDir, '/../gsea/', outPrefix, '_GO.txt'), col_names = T)

head(GO)




#all <- fgsea(pathways=gmtPathways("data/msigdb.v6.2.symbols.gmt"), ranks, nperm=1000) %>% 
#	as_tibble() %>% 
#	arrange(padj)

#all



