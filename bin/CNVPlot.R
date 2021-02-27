#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
sample_id = args[1]

library(CopyNumberPlots)

calls <- scan(paste(sample_id,"_cnv_calls.csv",sep=''), what="", sep="\n")
cnclass <- scan(paste(sample_id,"_cnv_class.csv",sep=''), what="", sep="\n")

cn.data <- toGRanges(calls)
cn.data$CopyNumberInteger <- cnclass
cn.data$LossHetero <- cn.data$CopyNumberInteger<2
cn.calls <- loadCopyNumberCalls(cn.data)

pdf(paste(sample_id,"_cnv.pdf",sep=''))
plotCopyNumberCalls(plotKaryotype(
    genome = "hg19",
    chromosomes= c(
        "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
        "chr9","chr10","chr11","chr12","chr13","chr14","chr15",
        "chr16","chr17","chr18","chr19","chr20","chr21","chr22")),
    cn.calls = cn.calls,
    labels = sample_id,
    label.cex = 0.5)
dev.off()
