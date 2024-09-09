#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.data_dir = "${launchDir}/data"
params.single_end = false
params.min_contig_length = 5000
params.min_bin_size = 100000
params.gtdb = "${launchDir}/refs/gtdb"
params.checkm = "${launchDir}/refs/checkm2"
params.maxcpus = 12


process contig_align {
    cpus 8
    memory "16 GB"
    time "2h"

    input:
    tuple val(id), path(contigs), path(reads)

    output:
    tuple val(id), path("${id}.bam"), path("${id}.bai")

    """
    minimap2 -ax sr -N 100 -t ${task.cpus} ${contigs} ${reads} | \
    samtools sort -@${task.cpus} -o ${id}.bam && \
    samtools index ${id}.bam ${id}.bai
    """
}

process coverage {
    cpus 1
    memory "32 GB"
    time "12h"
    publishDir "${params.data_dir}"

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path("${id}_coverage.txt")

    """
    jgi_summarize_bam_contig_depths --outputDepth ${id}_coverage.txt ${bam}
    """
}

process metabat {
    cpus 8
    memory "16 GB"
    time "12h"
    publishDir "${params.data_dir}"

    input:
    tuple val(id), path(contigs), path(coverage)

    output:
    tuple val(id), path("bins/${id}_*.fa.gz")

    """
    metabat2 -i ${contigs} -a ${coverage} -o bins/${id}_ \
        -t ${task.cpus} -m ${params.min_contig_length} -s 100000
    pigz -p ${task.cpus} bins/*.fa
    """
}

process checkm {
    cpus params.maxcpus
    memory "32 GB"
    time "2h"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(bins)

    output:
    path("checkm")

    """
    checkm2 predict --threads ${task.cpus} --input bins --output-directory checkm
    """
}

process gtdb_classify {
    cpus params.maxcpus
    memory "64 GB"
    time "8h"

    publishDir "${params.data_dir}"

    input:
    path(bins)

    output:
    path("gtdb")

    """
    GTDBTK_DATA_PATH=${params.gtdb} gtdbtk classify_wf \
        --genome_dir bins/ --prefix bins \
        --cpus ${task.cpus} --out_dir gtdb
    """

}


workflow {
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/preprocessed/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{reads}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/preprocessed/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/preprocessed/*_{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
            .set{reads}
    }

    Channel.
        fromPath("${params.data_dir}/assembled/contigs/*.fa")
        .map{row -> tuple(row.baseName.split("\\.contigs\\.fa")[0], row)}
        .set{assemblies}

    contig_align(assemblies.join(reads)) | coverage
    binned = metabat(assemblies.join(coverage.out))
    all_bins = binned.out.map{it -> it[1]}.collect()
    checkm(all_bins)
    gtdb_classify(all_bins)
}
