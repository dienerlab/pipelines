#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${launchDir}/data"
params.raw_data = "raw"
params.refs = env("DLP") ? "/home/isilon/dienerlab/refs" : "${launchDir}/refs"
params.eggnog_refs = "${params.refs}/eggnog"
params.metapackage = "${params.refs}/GlobDB_r232.metapackage_v4.smpkg"
params.viralpackage = "${params.refs}/lyrebird_v0.3.1_phrog_v4.1_metapackage_20250720.smpkg.zb"

params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threshold = 10
params.contig_length = 500
params.overlap = 0.8
params.identity = 0.99
params.threads = 12
params.method = "illumina"


def helpMessage() {
    log.info"""
    ~~~ Diener Lab Metagenomics Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf --resume

    An example run could look like:
    > nextflow run main.nf -with-conda /my/envs/metagenomics -resume \
                            --data_dir=./data --single_end=false --refs=/my/references \
                            --read_length=150

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --read_length [str]           The length of the reads.
      --single_end [bool]           Specifies that the input is single-end reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
      --method [str]                What sequencing technology was used. Can be "illumina", "nanopore", or "pacbio".
      --raw_data [str]             The folder inside data_dir containing the raw read files.
      --refs [str]                  Folder in which to find references DBs.
    Reference DBs:
      --refs [str]                  Folder in which to find references DBs.
      --eggnogg_refs [str]          Where to find EGGNOG references. Defaults to <refs>/eggnog.
      --metapackage [str]           Where to find the singleM metapackage.
      --viralpackage [str]          Where to find the lyrebird metapackage for viral detection.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
      --threshold [str]             Smallest abundance threshold used by Kraken.

    Assembly:
      --contig_length [int]         Minimum length of a contig.
      --identity [double]           Minimum average nucleotide identity.
      --overlap [double]            Minimum required overlap between contigs.
    """.stripIndent()
}



params.help = false


workflow {
    main:

    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    // find files
    if (params.single_end) {
        channel
            .fromPath("${params.data_dir}/${params.raw_data}/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}/${params.raw_data}!" }
            .set{raw}
    }

    // quality filtering
    preprocess(raw)

    // quantify taxa abundances
    singleM(preprocess.out)
    lyrebird(preprocess.out)
    summarizeProfiles(
        singleM.out.map{it -> it[1]}.collect(),
        lyrebird.out.map{it -> it[1]}.collect()
    )
    summarizeOTUs(
        singleM.out.map{it -> it[2]}.collect(),
        lyrebird.out.map{it -> it[2]}.collect()
    )

    mergePF(
        singleM.out.map{it -> it[3]}.collect()
    )

    // quality overview
    multiqc(
        preprocess.out
        .map{sa -> sa[2]}
        .collect()
    )

    // assemble de novo
    assemble(preprocess.out)

    // find ORFs and count them
    find_genes(assemble.out)
    preprocess.out.combine(find_genes.out, by: 0) | map_and_count
    merge_counts(map_and_count.out.collect())

    // cluster proteins, collapse mapping counts, and annotate clusters
    find_genes.out.map{sample -> sample[2]}.collect() | cluster_proteins
    filter_transcripts(
        find_genes.out.map{sample -> sample[1]}.collect(),
        cluster_proteins.out.map{sample -> sample[0]}
    )
    cluster_counts(merge_counts.out, cluster_proteins.out)
    annotate(cluster_proteins.out)

    publish:

    preprocessed = preprocess.out
    taxonomic_profiles = summarizeProfiles.out
    otus = summarizeOTUs.out
    multiqc_report = multiqc.out
    assemblies = assemble.out
    txns = filter_transcripts.out
    clusters = cluster_proteins.out
    counts = cluster_counts.out
    annotations = annotate.out
}

output {
    preprocessed {
        path "preprocessed"
    }

    taxonomic_profiles {
        mode "copy"
        overwrite true
    }

    otus {
        mode "copy"
        overwrite true
    }

    multiqc_report {
        mode "copy"
        overwrite true
    }

    assemblies {
        path "assemblies"
    }

    txns {
        mode "copy"
        overwrite true
    }

    clusters {
        mode "copy"
        overwrite true
    }

    counts {
        mode "copy"
        overwrite true
    }

    annotations {
        mode "copy"
        overwrite true
    }
}

process preprocess {
    cpus 3
    memory "4GB"
    time "30m"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), path("${id}_fastp.json"), path("${id}.html")

    script:
    if (params.single_end && params.method == "illumina")
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -w ${task.cpus} \
            --max_len1 ${params.read_length}
        """

    else if (!params.single_end && params.method == "illumina")
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -w ${task.cpus} \
            --max_len1 ${params.read_length} --max_len2 ${params.read_length}
        """
    else if (params.method == "nanopore" || params.method == "pacbio")
        """
        fastplong -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
             -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -w ${task.cpus}
        """
    else
        error "Unsupported method: ${params.method}"
}

process singleM {
    cpus 3
    memory 8.GB
    time 2.h

    input:
    tuple val(id), path(fastqs), path(json), path(html)

    output:
    tuple val(id), path("${id}_profile.tsv"), path("${id}_otus.tsv"), path("${id}_spf.tsv")

    script:
    if (params.single_end)
        """
        singlem pipe -1 ${fastqs} --threads ${task.cpus} \
            --metapackage ${params.metapackage} \
            -p ${id}_profile.tsv \
            --otu-table ${id}_otus.tsv

        singlem prokaryotic_fraction -1 ${fastqs} \
            -p ${id}_profile.tsv > ${id}_spf.tsv
        """
    else
        """
        singlem pipe -1 ${fastqs[0]} -2 ${fastqs[1]} --threads ${task.cpus} \
            --metapackage ${params.metapackage} \
            -p ${id}_profile.tsv \
            --otu-table ${id}_otus.tsv

        singlem prokaryotic_fraction -1 ${fastqs[0]} -2 ${fastqs[1]} \
            -p ${id}_profile.tsv > ${id}_spf.tsv
        """
}

process lyrebird {
    cpus 3
    memory 8.GB
    time 2.h

    input:
    tuple val(id), path(fastqs), path(json), path(html)

    output:
    tuple val(id), path("${id}_profile.tsv"), path("${id}_otus.tsv")

    script:
    if (params.single_end)
        """
        lyrebird pipe -1 ${fastqs} --threads ${task.cpus} \
            --metapackage ${params.viralpackage} \
            -p ${id}_profile.tsv \
            --otu-table ${id}_otus.tsv
        """
    else
        """
        lyrebird pipe -1 ${fastqs[0]} -2 ${fastqs[1]} --threads ${task.cpus} \
            --metapackage ${params.viralpackage} \
            -p ${id}_profile.tsv \
            --otu-table ${id}_otus.tsv
        """
}

process summarizeProfiles {
    cpus 1
    memory 16.GB
    time 2.h

    input:
    path(microbes)
    path(viral)

    output:
    path("*_abundances.tsv")

    script:
    """
    singlem summarise --input-taxonomic-profile ${microbes} \
        --output-taxonomic-profile-with-extras microbial_abundances.tsv

    singlem summarise --input-taxonomic-profile ${viral} \
        --output-taxonomic-profile-with-extras viral_abundances.tsv
    """
}

process mergePF {
    cpus 1
    memory 16.GB
    time 2.h

    input:
    path(spfs)

    output:
    path("prokaryotic_fractions.csv")

    script:
    """
    #!/usr/bin/env python

    import polars as pl

    df = pl.scan_csv("*_spf.tsv", separator="\\t").collect()
    df.write_csv("prokaryotic_fractions.csv")
    """
}

process summarizeOTUs {
    cpus 1
    memory 16.GB
    time 2.h

    input:
    path(microbial)
    path(viral)

    output:
    path("*_otu_abundances.tsv")

    script:
    """
    singlem summarise --input-otu-tables ${microbial} --cluster \
        --output-otu-table microbial_otu_abundances.tsv

    singlem summarise --input-otu-tables ${viral} --cluster \
        --output-otu-table viral_otu_abundances.tsv \
        --cluster-id 0.99
    """
}


process multiqc {
    cpus 1
    memory "8GB"
    time "1h"

    input:
    path(jsons)

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc ./
    """
}


process assemble {
    cpus 4
    memory "16GB"
    time "12h"

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("spades/contigs.fasta")

    script:
    if (params.single_end && params.method == "illumina")
        """
        spades.py -s ${reads} \
            -o spades -t ${task.cpus} -m ${task.memory.toGiga()}
        sed -i -e "s/^>/>${id}_/" spades/contigs.fasta
        """
    else if (!params.single_end && params.method == "illumina")
        """
        spades.py -1 ${reads[0]} -2 ${reads[1]} \
            -o spades -t ${task.cpus} -m ${task.memory.toGiga()} \
            --meta
        sed -i -e "s/^>/>${id}_/" spades/contigs.fasta
        """
    else if (params.method == "nanopore" || params.method == "pacbio")
        """
        metaMDBG asm --out-dir ./contigs --in-ont ${reads} --threads ${task.cpus}
        zcat contigs/contigs.fasta.gz | sed -e "s/^>/>${id}_/" > contigs/${id}.contigs.fa
        """
}

process find_genes {
    cpus 1
    memory "2GB"
    time "1h"

    input:
    tuple val(id), path(assembly)

    output:
    tuple val(id), path("${id}.ffn"), path("${id}.faa")

    script:
    """
    if grep -q ">" ${assembly}; then
        pyrodigal -p anon -i ${assembly} -o ${id}.gff -d ${id}.ffn -a ${id}.faa
    else
        touch ${id}.faa
        touch ${id}.ffn
    fi
    """
}

process cluster_proteins {
    cpus params.threads/2
    memory "40GB"
    time "2h"

    input:
    path(proteins)

    output:
    tuple path("proteins.faa"), path("proteins_cluster.tsv")

    script:
    """
    mmseqs easy-linclust ${proteins} proteins tmp \
        --cov-mode 0 -c ${params.overlap} \
        --min-seq-id ${params.identity} \
        --split-memory-limit 32G --threads ${task.cpus}
    rm -rf proteins_all_seqs.fna tmp
    mv proteins_rep_seq.fasta proteins.faa
    """
}

process filter_transcripts {
    cpus 1
    memory "8GB"
    time "8h"

    input:
    path(transcripts)
    path(proteins)

    output:
    path("transcripts.fna.gz")

    script:
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import os
    import gzip

    os.system("cat ${transcripts} > merged.fna")
    print("Reading transcript indices...")
    transcripts = SeqIO.index("merged.fna", "fasta")
    print("Reading protein indices...")
    proteins_idx = set(SeqIO.index("${proteins}", "fasta"))
    print("Writing filtered transcripts...")

    with gzip.open("transcripts.fna.gz", "wb") as out:
        for i, id in enumerate(proteins_idx, start=1):
            out.write(transcripts.get_raw(id))
            if (i % 100000) == 0:
                print(f"Processed {i} proteins.")
    os.system("rm merged.fna")
    """
}

process map_and_count {
    cpus 2
    memory "32 GB"
    time "4h"

    input:
    tuple val(id), path(reads), path(json), path(html), path(genes), path(proteins)

    output:
    path("${id}.sf")

    script:
    if (params.single_end && params.method == "illumina")
        """
        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -r ${reads} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        rm -rf ${id}_index
        """
    else if (!params.single_end && params.method == "illumina")
        """
        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -1 ${reads[0]} -2 ${reads[1]} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        rm -rf ${id}_index
        """
    else if (params.method == "nanopore" || params.method == "pacbio")
        """
        minimap2 -ax map-ont -p 1.0 -N 100 -t ${task.cpus} ${genes} ${reads} | samtools view -bS > ${id}.bam
        salmon quant -t ${genes} -q --ont --meta-l U -a ${id}.bam -o ${id} -p ${task.cpus} &&
            mv ${id}_salmon/quant.sf ${id}.sf || touch ${id}.sf
        rm ${id}.bam
        """
}

process merge_counts {
    cpus 1
    memory "8GB"
    time "4h"

    input:
    path(salmon_quants)

    output:
    path("gene_counts.csv.gz")

    script:
    """
    #!/usr/bin/env python

    from sys import stdin
    from os import path
    import pandas as pd
    import gzip

    paths = "${salmon_quants}"
    paths = paths.split(" ")
    nread = 0
    with gzip.open("gene_counts.csv.gz", "ab") as gzf:
        for p in paths:
            sample = path.splitext(path.basename(p))[0]
            print(f"Processing sample {sample}...")
            try:
                counts = pd.read_csv(p, sep="\t").query("NumReads > 0.1")
            except Exception:
                continue
            nread += 1

            counts.columns = [
                "locus_tag", "length", "effective_length", "tpm", "reads"]
            counts["sample_id"] = sample
            print(f"writing compressed output for sample {sample}...")
            counts.to_csv(gzf, header=(nread==1),
                          index=False)
    """
}

process cluster_counts {
    cpus 1
    memory "8GB"
    time "4h"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(gene_counts)
    tuple path(proteins), path(clusters)

    output:
    path("cluster_counts.csv.gz")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd

    counts = pd.read_csv("${gene_counts}")
    clusters = pd.read_csv("${clusters}", sep="\t")
    clusters.columns = ["representative", "contig"]
    clusters.set_index("contig", inplace=True)
    found = counts.locus_tag.isin(clusters.index)
    if (~found).any():
        not_clustered = counts.locus_tag[~found].unique()
        print(
            f"The following {len(not_clustered)} genes were omitted"
            f" because they were only observed once: {', '.join(not_clustered)}"
        )
        counts = counts[found]
    counts["cluster"] = clusters.representative[counts.locus_tag].values
    collapsed = counts.groupby(["sample_id", "cluster"]).agg({"tpm": "sum", "reads": "sum"}).reset_index()
    collapsed.to_csv("cluster_counts.csv.gz", index=False)
    """
}

process annotate {
    cpus params.threads
    memory "64GB"
    time "2d"

    input:
    tuple path(proteins), path(clusters)

    output:
    path("proteins.emapper.annotations")

    script:
    """
    EMTMP=\$(mktemp -d -t eggnog_results_XXXXXXXXXX)
    emapper.py -i ${proteins} --output proteins -m diamond \
        --data_dir ${params.eggnog_refs} --scratch_dir \$EMTMP --temp_dir \$TMPDIR \
        --cpu ${task.cpus}
    rm -rf \$EMTMP
    """
}
