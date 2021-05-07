#!/usr/bin/env Rscript

library(cn.mops)
load(file = "wgs_cnv_controls.RData")

genes <- read.table("hg19_genes.bed")
colnames(genes) <- c('chrom', 'start', 'stop', 'gene')

args = commandArgs(trailingOnly = TRUE)
sample_id = args[1]

bam <- paste(sample_id, "-sort.bam", sep = '')

bamDataRanges <- getReadCountsFromBAM(bam)

XandCB <- bamDataRanges
elementMetadata(XandCB) <- cbind(
    elementMetadata(XandCB),
    elementMetadata(controls)
)

result <- calcIntegerCopyNumbers(cn.mops(XandCB[1:length(XandCB)-1]))

segm <- as.data.frame(segmentation(result))
sampleSEGM <- subset(segm, sampleName == bam)
sampleSEGM$genes <- "NA"

for(i in 1:nrow(sampleSEGM)) {
    cnv <- sampleSEGM[i, ]
    genes_for_cnv <- c()
    rows_to_remove <- c()
    for(j in 1:nrow(genes)) {
        geneRow <- genes[j, ]
        if(geneRow$chrom == cnv$seqnames && geneRow$start >= cnv$start && geneRow$stop <= cnv$end) {
            genes_for_cnv <- c(genes_for_cnv, geneRow$gene)
            rows_to_remove <- c(rows_to_remove, j)
        }
    }
    if(length(rows_to_remove) != 0) {
        genes <- genes[-c(rows_to_remove), ]
    }
    sampleSEGM[i, ]$genes <- paste(genes_for_cnv, collapse = ",")
}

outFile <- paste(sample_id, "_cnvs.csv", sep = '')
write.csv(sampleSEGM, file = outFile)
