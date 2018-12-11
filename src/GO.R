library(gage)
library(gageData)
library(biomaRt)
library(dplyr)
library("AnnotationDbi")
library("org.Hs.eg.db")
#columns(org.Hs.eg.db)

#data(go.sets.hs)
#data(go.subs.hs)
#gobpsets = go.sets.hs[go.subs.hs$BP]       # BP - Biological Processes  
                                           # Also available - CC : Cellular Component 
    
setwd('~/depot/projects/Ma/July_2018')


                                       # and MF : Molecular Function

doAllGo <- function(foldchanges, prefix){
	for (st in c('BP','CC','MF')){
	  doGo(foldchanges, prefix, st)
	}
}

 
doGo <- function(foldchanges, prefix, set){



data(go.sets.hs)
data(go.subs.hs)
#gobpsets = go.sets.hs[go.subs.hs$BP]
gobpsets = go.sets.hs[go.subs.hs[[set]]]
    
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=T)

#lapply(gobpres, head)

gt <- as.data.frame(gobpres$greater)
#dim(gt)
gt <- na.omit(gt[gt$q.val < 0.05,])

terms <- row.names(gt)
    term <- substr(terms, start=1, stop=10)
    name <- substr(terms, start=12, stop=nchar(terms))
    
out <- data.frame(Term=term, Name=name, FDR=gt$q.val, set.size=gt$set.size)    
    

write.table(out, file=paste0( prefix, "_sig.", set, ".gt.txt"), sep="\t", quote=F, col.names=T, row.names = F)

lt <- as.data.frame(gobpres$less)
lt <- na.omit(lt[lt$q.val < 0.05,])
#dim(lt)

terms <- row.names(lt)
term <- substr(terms, start=1, stop=10)
name <- substr(terms, start=12, stop=nchar(terms))
    
out <- data.frame(Term=term, Name=name, FDR=lt$q.val, set.size=lt$set.size)    

write.table(out, file=paste0( prefix, "_sig.", set, ".lt.txt"), sep="\t", quote=F, col.names=T, row.names = F)


gobpres = gage(foldchanges, gsets=gobpsets, same.dir=F)
gt <- as.data.frame(gobpres$greater)
gt <- na.omit(gt[gt$q.val < 0.05,])

terms <- row.names(gt)
term <- substr(terms, start=1, stop=10)
name <- substr(terms, start=12, stop=nchar(terms))
    
out <- data.frame(Term=term, Name=name, FDR=gt$q.val, set.size=gt$set.size)    

write.table(out, file=paste0( prefix, "_sig.", set, ".both.txt"), sep="\t", quote=F, col.names=T, row.names = F)

}
