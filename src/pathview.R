# Pathways for Raj

library(pathview)
library(gage)
library(gageData)
library(biomaRt)
library(dplyr)
library("AnnotationDbi")
library("org.Rn.eg.db")
#columns(org.Hs.eg.db)

data(kegg.sets.rn)
data(sigmet.idx.rn)
kegg.sets.rn = kegg.sets.rn[sigmet.idx.rn]
head(kegg.sets.rn, 3)

#doPathways('PD-HOBvLTMC')


doPathways <- function(dir, prefix, foldchanges){

setwd(dir)

newdir <- paste0(prefix)  
  
dir.create(newdir, showWarnings=F)
setwd(newdir)

# Do the stats
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)


# Write out the up pathways
tmp <- as.data.frame(keggres$greater)
sigPaths <- na.omit(tmp[tmp$q.val < 0.05,])
pn <- row.names(sigPaths)
ids <- substr(pn, start=1, stop=8)
names <- substr(pn, start=10, stop=nchar(pn))

out <- data.frame(ID=ids, Pathway=names, FDR=sigPaths$q.val, set.size=sigPaths$set.size)
write.table(out, file=paste0(prefix,'up.txt'), col.names=T, row.names=F, sep="\t", quote=F) 

# Write the down pathways
tmp <- as.data.frame(keggres$less)
sigPaths <- na.omit(tmp[tmp$q.val < 0.05,])
pn <- row.names(sigPaths)
ids <- substr(pn, start=1, stop=8)
names <- substr(pn, start=10, stop=nchar(pn))

out <- data.frame(ID=ids, Pathway=names, FDR=sigPaths$q.val, set.size=sigPaths$set.size)
write.table(out, file=paste0(prefix,'down.txt'), col.names=T, row.names=F, sep="\t", quote=F) 

# Make pretty pictures

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="rno", new.signature=FALSE)


# do Up
dir.create('Up')
setwd('Up')

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(.$q.val < 0.05 ) %>%       # Get significant pathways
  .$id %>% 
  as.character()

keggresids = substr(keggrespathways, start=1, stop=8)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
#detach("package:dplyr", unload=TRUE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="rno"))


# do Down
setwd('..')
dir.create('Down')
setwd('Down')

# Get the pathways
keggrespathways = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(.$q.val < 0.05 ) %>%       # Get significant pathways
  .$id %>% 
  as.character()

keggresids = substr(keggrespathways, start=1, stop=8)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
#detach("package:dplyr", unload=TRUE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="rno"))

}







