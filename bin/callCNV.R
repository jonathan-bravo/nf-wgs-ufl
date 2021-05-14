#!/usr/bin/env Rscript

library(cn.mops)
load(file = "wgs_cnv_controls.RData")

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

outFile <- paste(sample_id, "_cnvs.csv", sep = '')
write.csv(sampleSEGM, file = outFile)
