#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.data_dir = "${launchDir}/data"
params.runtable = "${params.data_dir}/links.txt"

workflow {
    Channel.fromPath(params.runtable)
        .splitCsv()
        .map{ row -> row[0] }
        .filter{ s -> s=~/RawData.*\.fq\.gz$/ }
        .set{runs}

    runs | download
}


process download {
    cpus 6
    publishDir "${params.data_dir}/raw", mode: 'copy'

    input:
    val(link)

    output:
    path("*.fq.gz")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${link}
    """
}
