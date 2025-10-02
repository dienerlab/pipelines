#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.read_length = 150
params.out = "${launchDir}/data"
params.genomes = "${params.out}/genomes"
params.mutation_rate = 0.0
params.manifest = "${params.out}/manifest.csv"
params.preset = "illumina"
params.fullAmpliconSize = 4000


def relative_depth(abundance, depth) {
    Math.round(Float.parseFloat(abundance) * Float.parseFloat(depth))
}

workflow {
    Channel
        .fromPath("${params.manifest}")
        .splitCsv(header: true)
        .map{row -> tuple(
                row.sample_id, row.id, "${params.genomes}/${row.file}",
                relative_depth(row.relative_abundance, row.depth)
            )
        }
        .set{manifest}

    manifest | sample

    if (params.preset == "illumina") {
        sample.out
            .map{ it -> tuple(it[0], it[2][0]) }
            .groupTuple()
            .set{forward_groups}
        sample.out
            .map{ it -> tuple(it[0], it[2][1]) }
            .groupTuple()
            .set{reverse_groups}
        forward_groups
            .combine(reverse_groups, by: 0)
        .set{fastq_groups}
        fastq_groups | merge_paired | random_order
    } else if (params.preset == "nanopore") {
        sample.out
            .map{ it -> tuple(it[0], it[2]) }
            .groupTuple()
            .set{fastq_groups}
        fastq_groups | merge_single | random_order
    } else {
        error "Unsupported preset: "${params.preset}""
    }

}


process sample {
    cpus 1
    memory "4 GB"
    time "1h"
    beforeScript "ulimit -Sf unlimited"

    input:
    tuple val(sample_id), val(id), path(genome), val(n)

    output:
    tuple val(sample_id), val(id), path("*.fastq.gz")

    script:
    if (params.preset == "illumina")
        """
        trap "rm ref.fna" EXIT

        zcat ${genome} > ref.fna

        dwgsim -N ${n} -H \
            -c 0 -S 0 -A 0 -e "0.001-0.005" -E "0.005-0.01" \
            -1 ${params.read_length} -2 ${params.read_length} -n ${params.read_length} \
            -d ${params.read_length * 2} \
            -r ${params.mutation_rate} -y 0 \
            -P ${id} -o 1 ref.fna sim

        mv sim.bwa.read1.fastq.gz ${sample_id}_${id}_R1.fastq.gz
        mv sim.bwa.read2.fastq.gz ${sample_id}_${id}_R2.fastq.gz
        """
    else if (params.preset == "nanopore")
        """
        badread simulate --reference ${genome} \
            --quantity ${n * params.fullAmpliconSize} \
            | head -n ${4 * n} | gzip > ${sample_id}_${id}_R1.fastq.gz
        """
}

process merge_paired {
    cpus 1
    memory "4 GB"
    time "1h"

    input:
    tuple val(sample_id), path(forward), path(reverse)

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.gz")

    script:
    """
    cat ${forward} > ${sample_id}_R1.fastq.gz
    cat ${reverse} > ${sample_id}_R2.fastq.gz
    """

}

process merge_single {
    cpus 1
    memory "4 GB"
    time "1h"

    input:
    tuple val(sample_id), path(forward)

    output:
    tuple val(sample_id), path("${sample_id}_*.fastq.gz")

    script:
    """
    cat ${forward} > ${sample_id}_R1.fastq.gz
    """

}

process random_order {
    cpus 4
    memory "16 GB"
    time "4h"
    publishDir params.out, mode: "copy", overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("sampled/*.fastq.gz")

    script:
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import gzip
    import os
    from random import shuffle
    from pathlib import Path
    import shutil

    base = Path("sampled")
    forward = "${reads[0]}" if "${params.preset}" == "illumina" else "${reads}"

    with gzip.open(forward, "rt") as infile, open("fwd.fastq", "w") as outfile:
        shutil.copyfileobj(infile, outfile)
    fidx = SeqIO.index("fwd.fastq", "fastq")
    fnames = list(fidx.keys())

    if "${params.preset}" == "illumina":
        with gzip.open("${reads[1]}", "rt") as infile, open("rev.fastq", "w") as outfile:
            shutil.copyfileobj(infile, outfile)
        ridx = SeqIO.index("rev.fastq", "fastq")
        rnames = list(ridx.keys())

    print(f"Finished indexing {len(fidx)} records...")
    ord = list(range(len(fidx)))
    shuffle(ord)

    os.mkdir(base)
    print("Sampling records...", flush=True)

    try:
        fwd = gzip.open(base / "${sample_id}_R1.fastq.gz", "wb")
        if "${params.preset}" == "illumina":
            rev = gzip.open(base / "${sample_id}_R2.fastq.gz", "wb")

        for i, rec in enumerate(ord):
            fwd.write(fidx.get_raw(fnames[rec]))
            if "${params.preset}" == "illumina":
                rev.write(ridx.get_raw(rnames[rec]))
            if i % 10000 == 0:
                print(f"Processed {i} records...", flush=True)

    finally:
        fwd.close()
        os.remove("fwd.fastq")
        if "${params.preset}" == "illumina":
            rev.close()
            os.remove("rev.fastq")
    print("Done.")
    """

}
