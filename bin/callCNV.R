#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
sample_id = args[1]
control_data = args[2]

library(panelcn.mops)

load(file=control_data)

testBam <- paste(sample_id,"-sort.bam",sep='')

test <- countBamListInGRanges(
countWindows = countWindows,
bam.files    = testBam,
read.width   = 150)

XandCB <- test
elementMetadata(XandCB) <- cbind(
    elementMetadata(XandCB),
    elementMetadata(control))

resultsList <- runPanelcnMops(
    XandCB,
    testiv        = 1:ncol(elementMetadata(test)),
    countWindows  = countWindows,
    selectedGenes = NULL,
    maxControls   = 25)

sampleNames <- colnames(elementMetadata(test))
resulttable <- createResultTable(
    resultlist    = resultsList,
    XandCB        = XandCB,
    countWindows  = countWindows,
    selectedGenes = NULL,
    sampleNames   = sampleNames)

write.table(resulttable[[1]], paste(sample_id,'_cnv_table.csv',sep=''), append = TRUE)
