# RNAseq
Template for basic RNA seq analysis


## Steps:

Steps that have been done are marked with a  :white_check_mark: while those that have yet to be started are marking with a  :x:. Steps in progress are marked by a :large_orange_diamond:.


### Quality Control
In the reads folder, create a fastqc folder. Then run:

```
for f in *.gz; do echo $f; fastq -o fastqc -t 40 $f; done
```

This makes QC reports for each file. To get a single report, cd to the fastqc folder and use:
```
multiqc -o multiqc .
```

Usually quality is good enough to proceed without any trimming, but do make note of the total read numbers, to see if any samples are too low, or there are outliers.

### Quantification

Download the appropriate transcripts and build a salmon index

In the reads folder, use salmon to quantify each read

### Differential expression

Use the R script doDESeq.R to find differentially expressed genes, and produce some figures

### GO Analysis

### Pathway Analysis

