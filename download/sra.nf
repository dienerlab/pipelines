nextflow.enable.dsl=2

params.runtable = "${launchDir}/data/runtable.csv"
params.out = "${launchDir}/data"

process download {
    cpus 4
    memory "2GB"
    time "2h"
    maxRetries 3
    errorStrategy { (task.attempt <= task.maxRetries)  ? 'retry' : 'ignore' }
    publishDir "${params.out}/raw", mode: 'copy', overwrite: true

    input:
    val(run)

    output:
    path("*.fastq.gz")

    script:
    """
    aws s3 sync s3://sra-pub-run-odp/sra/${run} sra_${run} --no-sign-request
    fasterq-dump -e ${task.cpus} -f -3 ./sra_${run}/${run} -t -x && rm -rf ./sra_${run}
    pigz -p ${task.cpus} *.fastq
    """
}

// Use this to merge files by groups
process mergeFastq {
    publishDir "${params.out}/merged"

    cpus 1
    memory "1GB"
    time "1h"

    input:
    tuple val(id), path(runs)

    output:
    path("${id}.fastq.gz")

    script:
    """
    cat ${runs} > ${id}.fastq.gz
    """
}
