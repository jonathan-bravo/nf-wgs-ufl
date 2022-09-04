#!/usr/bin/env Rscript

library(cn.mops)

controlBam <- list.files(pattern="*.bam$")

sample_id <- unlist(
    strsplit(
        tools::file_path_sans_ext(controlBam),
        "_"
    )
)[1]

assign(
    paste("control_", sample_id, sep=""),
    getReadCountsFromBAM(controlBam, WL = 3500)
)

save(
    list = paste("control_", sample_id, sep=""),
    file = paste(sample_id, "_cnv_control.RData", sep="")
)
