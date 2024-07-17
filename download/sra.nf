nextflow.enable.dsl=2

params.runtable = "${launchDir}/data/runtable.csv"

process download {
    cpus 4
    memory "2GB"
    time "2h"
    maxRetries 3
    errorStrategy { (task.attempt <= maxRetries)  ? 'retry' : 'ignore' }
    publishDir "${launchDir}/data/raw", mode: 'copy'

    input:
    val(run)

    output:
    path("*.fastq.gz")

    """
    aws s3 sync s3://sra-pub-run-odp/sra/${run} sra_${run} --no-sign-request
    fasterq-dump -e ${task.cpus} -f -3 ./sra_${run}/${run} -t -x && rm -rf ./sra_${run}
    pigz -p ${task.cpus} *.fastq
    """
}

// Use this to merge files by groups
process mergeFastq {
    publishDir "${launchDir}/data/merged"

    cpus 1
    memory "1GB"
    time "1h"

    input:
    tuple val(id), path(runs)

    output:
    path("${id}.fastq.gz")

    """
    cat ${runs} > ${id}.fastq.gz
    """
}

workflow {
    Channel.fromPath(params.runtable)
        .splitCsv(header: true)
        .map{row -> row.Run}
        .set{runs}

    runs | download
}
