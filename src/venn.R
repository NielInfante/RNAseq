library(tidyverse)

na <- read_tsv('gsea/Norm_adj_hallmark.txt')
nd <- read_tsv("gsea/Norm_dist_hallmark.txt")
oa <- read_tsv('gsea/Ovov_adj_hallmark.txt')
od <- read_tsv("gsea/Ovov_dist_hallmark.txt")

norm <- read_tsv("gsea/Norm_hallmark.txt")
ovob <- read_tsv("gsea/Ovob_hallmark.txt")



na <- filter(na, padj <= 0.1 ) %>% dplyr::select(pathway)
nd <- filter(nd, padj <= 0.1 ) %>% dplyr::select(pathway)
oa <- filter(oa, padj <= 0.1 ) %>% dplyr::select(pathway)
od <- filter(od, padj <= 0.1 ) %>% dplyr::select(pathway)

norm <- filter(norm, padj <= 0.1 )
ovob <- filter(ovob, padj <= 0.1 )

norm <- dplyr::select(norm, pathway, norm_padj=padj)
ovob <- dplyr::select(ovob, pathway, ovob_padj=padj)

comb <- inner_join(ovob, norm)
comb <- filter(comb, ovob_padj <= 0.1)

filter(comb, ovob_padj <= 0.05)




#library(gplots)
itemList <- venn(list(na=na$pathway, nd=nd$pathway, oa=oa$pathway, od=od$pathway), show.plot=F)

library(VennDiagram)
png('pics/hallmark_venn.png')
vd <- venn.diagram(list("Norm Adjacent"=na$pathway, "Norm Distant"=nd$pathway, 
												"Obese Adjacent"=oa$pathway, "Obese Distant"=od$pathway), 
									 fill=2:5, alpha=0.4, filename=NULL, cex=1.3, cat.cex=1.4,
									 main="Hallmark Pathways")
grid.newpage()  ;  grid.draw(vd)
dev.off()


png('pics/hallmark_Norm_Ob.png')
vd <- venn.diagram(list("Normal"=norm$pathway, "Obese"=ovob$pathway), 
									 fill=4:5, alpha=0.4, filename=NULL, cex=1.3, cat.cex=1.4,
									 main="Hallmark Pathways")
grid.newpage()  ;  grid.draw(vd)
dev.off()






