#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)
sample_id = args[1]
male_control_data = args[2]
female_control_data = args[3]

library(panelcn.mops)

m_or_f <- paste(sample_id, "_m_or_f.txt", sep = '')
sample_sex <- read.table(m_or_f, header = FALSE, sep = "\n")

if(sample_sex[1, 1] == "Female"){
    sex <- "Female"
    load(file = female_control_data)
} else {
    sex <- "Male"
    load(file = male_control_data)
}

testBam <- paste(sample_id, "-sort.bam", sep = '')

test <- countBamListInGRanges(
    countWindows = countWindows,
    bam.files    = testBam,
    read.width   = 150
)

XandCB <- test
elementMetadata(XandCB) <- cbind(
    elementMetadata(XandCB),
    elementMetadata(control)
)

resultsList <- runPanelcnMops(
    XandCB,
    testiv        = 1:ncol(elementMetadata(test)),
    countWindows  = countWindows,
    selectedGenes = NULL,
    maxControls   = 25,
    sex           = sex
)

sampleNames <- colnames(elementMetadata(test))
resulttable <- createResultTable(
    resultlist    = resultsList,
    XandCB        = XandCB,
    countWindows  = countWindows,
    selectedGenes = NULL,
    sampleNames   = sampleNames
)

write.table(
    resulttable[[1]],
    paste(sample_id, '_cnv_table.csv', sep = ''),
    append = TRUE
)
