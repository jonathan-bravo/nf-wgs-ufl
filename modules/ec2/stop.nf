#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process STOP_EC2 {

    label 'local'

    input:
    file "${params.run_id}.html"

    script:
    """
    sudo shutdown -P +1
    """
}