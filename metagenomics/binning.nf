#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.data_dir = "${launchDir}/data"
params.single_end = false
params.min_contig_length = 5000
params.min_bin_size = 100000
params.gtdb = "${launchDir}/refs/gtdb"
params.checkm = "${launchDir}/refs/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
params.ani = 0.99
params.maxcpus = 12
params.preset = "illumina"
params.conda_path = "$HOME/miniforge3/envs"
params.manifest = null


process contig_align {
    cpus 8
    memory "8 GB"
    time "2h"
    conda "${params.conda_path}/metagenomics"

    input:
    tuple val(id), path(contigs), path(reads)

    output:
    tuple val(id), path("${id}.bam"), path("${id}.bai")

    script:
    def mode = params.preset == "nanopore" ? "map-ont" : "sr"
    """
    minimap2 -ax ${mode} -N 100 -t ${task.cpus} ${contigs} ${reads} | \
    samtools sort -@${task.cpus} -o ${id}.bam && \
    samtools index ${id}.bam ${id}.bai
    """
}

process coverage {
    cpus 1
    memory "1 GB"
    time "30m"
    publishDir "${params.data_dir}"
    conda "${params.conda_path}/metagenomics"

    input:
    tuple val(id), path(bam), path(bai)

    output:
    tuple val(id), path("${id}_coverage.txt")

    """
    jgi_summarize_bam_contig_depths --outputDepth ${id}_coverage.txt ${bam}
    """
}

process metabat {
    cpus 4
    memory "2 GB"
    time "12h"
    conda "${params.conda_path}/metagenomics"

    input:
    tuple val(id), path(contigs), path(reads), path(coverage)

    output:
    tuple val(id), path("bins/${id}.*.fa.gz")

    """
    metabat2 -i ${contigs} -a ${coverage} -o bins/${id} \
        -t ${task.cpus} -m ${params.min_contig_length} -s 100000
    pigz -p ${task.cpus} bins/*.fa
    """
}

process dereplicate {
    cpus 4
    memory "8 GB"
    time "2h"
    conda "${params.conda_path}/binchecks"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(checkm_report)
    tuple path(gtdb_report), path(raw)

    output:
    tuple path("dereplicated/*.fna"), path("dereplicated/figures")

    """
    dRep dereplicate ./dereplicated -g ${raw} \
        --genomeInfo ${checkm_report} \
        --S_algorithm skani --S_ani 0.99 \
        --processors ${task.cpus}
    mv dereplicated/dereplicated_genomes/*.fa.gz dereplicated/
    """
 }

 process checkm {
    cpus params.maxcpus
    memory "32 GB"
    time "2h"
    conda "${params.conda_path}/binchecks"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(bins)

    output:
    path("checkm/quality_report.tsv")

    """
    mkdir bins && mv ${bins} bins
    checkm2 predict \
        --threads ${task.cpus} \
        --database_path ${params.checkm} \
        --extension .fa.gz \
        --input bins --output-directory checkm
    """
 }

 process format_report {
    cpus 1
    memory "250 MB"
    time "10m"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    conda "${params.conda_path}/binchecks"

    input:
    path(report)

    output:
    path("checkm2_report.csv")

    """
    #!/usr/bin/env python

    import pandas as pd

    report = pd.read_csv("${report}", sep="\\t").rename(columns={"Name": "genome"})
    report.columns = report.columns.str.lower()
    report.genome = report.genome + ".gz"
    report.to_csv("checkm2_report.csv", index=False)
    """
}


process gtdb_classify {
    cpus params.maxcpus
    memory "80 GB"
    time "8h"
    conda "${params.conda_path}/binchecks"

    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    tuple path(bins)

    output:
    path("bins.*.summary.tsv")

    """
    mkdir bins && mv ${bins} bins
    GTDBTK_DATA_PATH=${params.gtdb} gtdbtk classify_wf \
        --genome_dir bins --prefix bins --extension fa.gz \
        --mash_db mash --extension fa.gz\
        --cpus ${task.cpus} --out_dir gtdb
    cp -L gtdb/bins.*.summary.tsv .
    """
}

process rename {
    cpus 1
    memory "100 MB"
    time "30 m"
    conda "${params.conda_path}/binchecks"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(bins)
    path(reports)

    output:
    tuple path("gtdb_report.csv"), path("bins/*.fna")

    """
    #!/usr/bin/env python

    import gzip
    import shutil
    import os
    import re
    import pandas as pd

    def first_name(s):
        names = [s for s in s.split(";") if re.search("\\\\w__\$", s) is None]
        return names[-1].replace(" ", "_")

    def extract(row):
        with gzip.open(row["user_genome"] + ".fa.gz", "rb") as gz:
            with open(f"bins/{row['id']}.fna", "wb") as out:
                shutil.copyfileobj(gz, out)

    report = []
    for r in "${reports}".split():
        report.append(pd.read_csv(r, sep="\\t"))
    report = pd.concat(report)
    report["id"] = report.classification.apply(first_name)
    report["id"] = report["id"] + "_" + report["user_genome"]
    
    os.mkdir("bins")
    report.apply(extract, axis=1)
    report.to_csv("gtdb_report.csv", index=False)
    """
}


workflow {
    if (params.manifest == null) {
        if (params.single_end) {
            Channel
                .fromPath("${params.data_dir}/preprocessed/*.fastq.gz")
                .map{row -> tuple(row.baseName, tuple(row))}
                .set{reads}
        } else {
            Channel
                .fromFilePairs([
                    "${params.data_dir}/preprocessed/*_R{1,2}_001.fastq.gz",
                    "${params.data_dir}/preprocessed/*_{1,2}.fastq.gz",
                    "${params.data_dir}/preprocessed/*_R{1,2}.fastq.gz"
                ])
                .ifEmpty { error "Cannot find any read files in ${params.data_dir}!" }
                .set{reads}
        }

        clean = reads.map{tuple it[0].replace("_filtered", ""), it[1]}

        Channel.
            fromPath("${params.data_dir}/assembled/contigs/*.contigs.fa")
            .map{row -> tuple(row.baseName.split("\\.contigs")[0], row)}
            .set{assemblies}

        merged = assemblies.join(clean)
    } else {
        log.info("Using manifest file ${params.manifest}.")
        Channel
            .fromPath("${launchDir}/${params.manifest}")
            .splitCsv(header: true)
            .set{rows}
        if (params.single_end) {
            rows.map{row -> tuple(
                    row.id,
                    file("${launchDir}/${row.contigs}"),
                    file("${launchDir}/${row.forward}")
                    )
                }
                .set{merged}
        } else {
            rows.map{row -> tuple(
                    row.id,
                    file("${launchDir}/${row.contigs}"),
                    [file("${launchDir}/${row.forward}"), file("${launchDir}/${row.reverse}")]
                    )
                }
                .set{merged}
        }
    }

    merged | contig_align | coverage
    binned = metabat(merged.join(coverage.out))
    all_bins = binned.map{it -> it[1]}.collect()
    all_bins | checkm | format_report
    all_bins | gtdb_classify
    rename(all_bins, gtdb_classify.out)
    dereplicate(format_report.out, rename.out)

}
