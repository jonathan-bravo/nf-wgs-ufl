#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
sample_id = args[1]

library(CopyNumberPlots)

calls <- scan(
    paste(sample_id, "_cnv_calls.csv", sep = ''),
    what = "",
    sep = "\n"
)
cnclass <- scan(
    paste(sample_id, "_cnv_class.csv", sep = ''),
    what = "",
    sep = "\n"
)

for(i in 1:length(cnclass)){
    if(cnclass[i] <= 2){
        cnclass[i] <- 0
    }
    else if(cnclass[i] > 2){
        cnclass[i] <- 4
    }
}

cn.data <- toGRanges(calls)
cn.data$CopyNumberInteger <- cnclass
cn.data$LossHetero <- cn.data$CopyNumberInteger<2
cn.calls <- loadCopyNumberCalls(cn.data)

pdf(paste(sample_id, "_cnv.pdf", sep = ''))
plotCopyNumberCalls(
    plotKaryotype(
        genome = "hg19",
        chromosomes= c(
            seqlevels(cn.data)
        )
    ),
    cn.calls = cn.calls,
    labels = sample_id,
    label.cex = 0.5
)
cn.cols <- getCopyNumberColors()[1]
cn.cols <- c(cn.cols, getCopyNumberColors()[5])
legend("bottomright", legend=c('DEL', 'DUP'), fill = cn.cols, ncol=length(cn.cols), cex=0.75)
dev.off()
