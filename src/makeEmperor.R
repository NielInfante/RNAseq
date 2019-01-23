<<<<<<< HEAD
library(DESeq2)
dds <- readRDS(file='~/depot/projects/Pistilli/Pistilli_2019_01/deseq/All_dds.rds')

vsd <- vst(dds, blind=F)
res<-results(dds)
res <- as.data.frame(res)



=======
>>>>>>> 2fa6b21439b8173f68a9a7bd1dcfa752967ab92b

# Get most expressed genes
topExp <- res[order(res$baseMean, decreasing=T)[1:500],]$ID

vsd$SampleID
a <- as.data.frame(assay(vsd))

a <- a[topExp,]
dim(a)
names(a) <- vsd$SampleID
a <- t(a)


t1 <- prcomp(a, center = T, scale. = T)


summary(t1)


str(t1)
dim(t1$x)

coords <- as.data.frame(t1$x[,1:3])



write.table(coords, file='~/depot/projects/Pistilli/Pistilli_2019_01/emperor/pca.dat',
						quote=F, sep="\t",row.names = T)

# add "pc vector number" to header manually



varexp <- t1$sdev^2 / sum(t1$sdev^2)
varexp <- c('% variation explained', varexp[1:3])


eigenvals <- (t1$sdev ^ 2)[1:3] 

eig <- c('eigvals',eigenvals)


meta <- data.frame(SampleID=vsd$SampleID, Surgeon=vsd$Surgeon, Group=vsd$Group, 
									 Run=vsd$Run, Cancer=vsd$Cancer, EdName=vsd$EdName)

write.table(meta, file='~/depot/projects/Pistilli/Pistilli_2019_01/emperor/meta.dat',
						quote=F, sep="\t",row.names = F)


# Add # to header manually



# On the command line,
# conda activate emperor
# make_emperor.py -i pca.dat -m meta.dat
# conda deactivate




