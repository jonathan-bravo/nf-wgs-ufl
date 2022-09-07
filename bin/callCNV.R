#!/usr/bin/env Rscript

library(cn.mops)

args = commandArgs(trailingOnly = TRUE)
cnv_control = args[1]
bam_file = args[2]
out_file = args[3]

load(file = cnv_control)

bam <- paste(bam_file, sep = '')

bamDataRanges <- getReadCountsFromBAM(bam, WL = 3500)

XandCB <- bamDataRanges
elementMetadata(XandCB) <- cbind(
    elementMetadata(XandCB),
    elementMetadata(controls)
)

seqlevels(XandCB, pruning.mode="coarse") <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

result <- calcIntegerCopyNumbers(cn.mops(XandCB))

segm <- as.data.frame(segmentation(result))
sampleSEGM <- subset(segm, sampleName == bam)

outFile <- paste(out_file, sep = '')
write.csv(sampleSEGM, file = outFile)
