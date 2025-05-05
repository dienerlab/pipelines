#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.assembly_levels = "Complete Genome"
params.dbs = "archaea,bacteria,vertebrate_mammalian/Homo_sapiens"
params.url = "rsync://ftp.ncbi.nlm.nih.gov/genomes"
params.out = "${launchDir}/data"


workflow {
    good_levels = params.assembly_levels.split(",")
    cats = params.dbs.split(",")

    Channel.from(cats).set{summaries}
    summaries | assembly_summaries

    assembly_summaries.out
        .map{it.text}
        .splitCsv(header: true, sep: "\t", skip: 1)
        .filter{ good_levels.contains(it.assembly_level) }
        .map{
            row -> tuple(
                row.ftp_path.split("/").last(),
                row.species_taxid,
                row.ftp_path.replaceFirst("https:", "rsync:")
            )
        }
        .set{files}

    files | download
    download.out.groupTuple() | merge_species
}


process assembly_summaries {
    maxRetries 4
    errorStrategy { (task.attempt <= task.maxRetries)  ? 'retry' : 'ignore' }
    cpus 2
    memory "2GB"
    time "2h"

    publishDir params.out, mode: "copy", overwrite: true

    input:
    val(category)

    output:
    path("assembly_summary_*.tsv")

    script:
    """
    rsync --no-motd ${params.url}/refseq/${category}/assembly_summary.txt \
        ./assembly_summary_${category.replaceAll("/", "_")}.tsv
    """
}

process download {
    maxRetries 4
    errorStrategy { (task.attempt <= task.maxRetries)  ? 'retry' : 'ignore' }
    cpus 2
    memory "4 GB"
    time "2h"

    input:
    tuple val(acc), val(sid), val(url)

    output:
    tuple val(sid), path("*.fna.gz")

    script:
    """
    rsync --no-motd ${url}/${acc}_genomic.fna.gz .
    """
}

process merge_species {
    publishDir "${params.out}/background_genomes", mode: "copy", overwrite: true
    cpus 1
    memory "4 GB"
    time "1h"

    input:
    tuple val(sid), path(fasta)

    output:
    path("${sid}.fna.gz")

    script:
    """
    cat ${fasta} > ${sid}.fna.gz
    """

}
