#!/usr/bin/env nextflow

params.data = "${launchDir}"
params.studies = "${params.data}/studies.csv"
params.visualize = true

workflow {
    // Read the studies table
    studies = Channel.fromPath(params.studies)
        .splitCsv(header: true, sep: ",")

    // Do some work on the studies

    studies | import_data | step1

    // Visualize if desired
    if (params.visualize) {
        visualize(import_data.out)
    }

}

process import_data {
    publishDir "${params.data}/imports", mode: 'copy', overwrite: true

    cpus 4
    memory "8GB"
    time "2h"

    input:
    val(study)

    output:
    tuple val(study), path("*.txt")

    script:
    if (study.layout == "paired") {
        """
        echo "importing from ${params.data}/${study.location}/manifest.tsv in paired-end layout." > '${study.id}.txt'
        """
    } else if (study.layout == "single") {
        """
        echo "importing from ${params.data}/${study.location}/manifest.tsv in single-end layout." > '${study.id}.txt'
        """
    } else {
        error "Invalid libray layout specified. Must be 'paired' or 'single' :("
    }
}

process step1 {
    publishDir "${params.data}/step1", mode: 'copy', overwrite: true

    cpus 1
    memory "2GB"
    time "1h"

    input:
    tuple val(study), path(imported)

    output:
    tuple val(study), path("*.result")

    script:
    """
    echo "processed ${study.id} data with truncations of ${study.trunc_f},${study.trunc_r}" > '${study.id}.result'
    """
}

process visualize {
    publishDir "${params.data}/viz", mode: 'copy', overwrite: true

    cpus 1
    memory "2GB"
    time "1h"

    input:
    tuple val(study), path(imported)

    output:
    tuple val(study), path("${study.id}.viz")

    script:
    """
    echo "visualized ${study.id} import" > '${study.id}.viz'
    """
}