#!/usr/bin/env Rscript

library(cn.mops)

controlRData <- list.files(pattern="*.RData$")

getSampleId <- function(x) {
    sample_id <- unlist(
        strsplit(
            tools::file_path_sans_ext(x),
            "_"
        )
    )[1]
    return(sample_id)
}

for (i in 1:length(controlRData)) {
    load(file = controlRData[i])
    sample_id <- getSampleId(controlRData[i])
    if(i == 1) {
        controls <- get(paste("control_", sample_id, sep=""))
    } else {
        elementMetadata(controls) <- cbind(
            elementMetadata(controls),
            elementMetadata(get(paste("control_", sample_id, sep="")))
        )
    }
}

save(controls, file = "wgs_cnv_controls.RData")
