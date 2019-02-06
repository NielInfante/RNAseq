library(DESeq2)

# Files
inputFile <- ''  # dds as rds
outputFileMain <- ''   # Will append _pca.dat 

# Number of Principal Components you want in the final output
numberPCs <- 5

# Read dds
dds <- readRDS(file=inputFile)

vsd <- vst(dds, blind=F)
res<-results(dds)
res <- as.data.frame(res)
res$ID <- row.names(res)

# Get most expressed genes
topExp <- res[order(res$baseMean, decreasing=T)[1:500],]$ID

# Get normalized counts, arrange data
a <- as.data.frame(assay(vsd))
a <- a[topExp,]
dim(a)
names(a) <- vsd$SampleID
a <- t(a)

# Do the PCA
t1 <- prcomp(a, center = T, scale. = T)


# Get coordinates that you want
coords <- as.data.frame(t1$x[,1:numberPCs])


coords$Name <- row.names(coords)
coords <- coords[c(numberPCs+1,1:numberPCs)]
names(coords)[1] <- 'pc vector number'

write.table(coords, file=paste0(outputFileMain, '_pca.dat'), quote=F, sep="\t",row.names = F)


# Get variation explained and eigenvalues
varexp <- t1$sdev^2 / sum(t1$sdev^2)
varexp <- c('% variation explained', varexp[1:numberPCs])

eigenvals <- (t1$sdev ^ 2)[1:numberPCs] 

eig <- c('eigvals',eigenvals)

# Add data to table
cat('\n\n', 
		paste0(eig, collapse='\t'), '\n',
		paste0(varexp, collapse='\t'),
		file=paste0(outputFileMain, '_pca.dat'),'\n', append=T)


####  Step to be done manually  #####

# Copy metadata file, add # to header
# The first column should be #SampleID



# On the command line,
# conda activate emperor
# make_emperor.py -i pca.dat -m meta.dat -o outdir
# conda deactivate




