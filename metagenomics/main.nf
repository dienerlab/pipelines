#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${launchDir}/data"
params.raw_data = "raw"
params.refs = "${launchDir}/refs"
params.eggnog_refs = "${params.refs}/eggnog"
params.kraken2_db = "${params.refs}/kraken2"
params.kraken2_mem = null

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
params.confidence = 0.3
params.ranks = "D,P,G,S"
params.batchsize = 50


// Helper to calculate the required RAM for the Kraken2 database
def estimate_db_size(hash, extra) {
    def db_size = null

    // Calculate db memory requirement
    if (params.dbmem) {
        db_size = MemoryUnit.of("${params.dbmem} GB")
    } else {
        db_size = MemoryUnit.of(file(hash).size()) + extra
        log.info("Based on the hash size I am reserving ${db_size.toGiga()}GB of memory for Kraken2.")
    }

    return db_size
}


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
    Reference DBs:
      --refs [str]                  Folder in which to find references DBs.
      --eggnogg_refs [str]          Where to find EGGNOG references. Defaults to <refs>/eggnog.
      --kraken2_db [str]            Where to find the Kraken2 reference. Defaults to <refs>/kraken2_default.
      --kraken2_mem [str]           The maximum amount of memory for Kraken2. If not set will choose this automatically
                                    based on the database size. Thus, only use to limit Kraken2 to less memory.
    Quality filter:
      --trim_front [str]            How many bases to trim from the 5' end of each read.
      --min_length [str]            Minimum accepted length for a read.
      --quality_threshold [str]     Smallest acceptable average quality.
      --threshold [str]             Smallest abundance threshold used by Kraken.

    Assembly:
      --contig_length [int]         Minimum length of a contig.
      --identity [double]           Minimum average nucleotide identity.
      --overlap [double]            Minimum required overlap between contigs.

    Taxonomic classification:
      --batchsize [int]             The batch size for Kraken2 jobs. See documentation
                                    for more info. Should be 1 on single machine setups
                                    and much larger than one on HPC setups.
      --kraken2_mem [int]           Maximum memory in GB to use for Kraken2. If not set
                                    this will be determined automatically from the database.
                                    So, only set this if you want to overwrite the automatic
                                    detection.
    """.stripIndent()
}



params.help = false


workflow {
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    Channel
    .of(params.ranks.split(","))
    .set{levels}

    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/${params.raw_data}/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
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

    // buffer the samples into batches
    batched = preprocess.out
        .collate(params.batchsize)
        .map{it -> tuple it.collect{a -> a[0]}, it.collect{a -> a[1]}.flatten()}
    // run Kraken2
    kraken(batched)
    reports = kraken.out
        .flatMap{it[1]}
        .map{tuple it.baseName.split(".tsv")[0], it}

    count_taxa(reports.combine(levels))
    count_taxa.out.map{s -> tuple(s[1], s[3])}
        .groupTuple()
        .set{merge_groups}
    merge_taxonomy(merge_groups)

    // quality overview
    multiqc(merge_taxonomy.out.collect())

    // assemble de novo
    megahit(preprocess.out)

    // find ORFs and count them
    find_genes(megahit.out)
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
}

process preprocess {
    cpus 3
    memory "4GB"
    time "30m"

    publishDir "${params.data_dir}/preprocessed"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id), path("${id}_filtered_R*.fastq.gz"), path("${id}_fastp.json"), path("${id}.html")

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus} \
            --max_len1 ${params.read_length}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus} \
            --max_len1 ${params.read_length} --max_len2 ${params.read_length}
        """
}


process kraken {
    cpus params.threads
    memory { estimate_db_size("${params.kraken2_db}/hash.k2d", 16.GB) }
    time { 2.h + ids.size() * 0.5.h }
    scratch false
    publishDir "${params.data_dir}/kraken2"

    input:
    tuple val(ids), path(reads)

    output:
    tuple path("*.k2"), path("*.tsv")

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    from subprocess import run

    base_args = [
        "kraken2", "--db", "${params.kraken2_db}",
        "--confidence", "${params.confidence}",
        "--threads", "${task.cpus}", "--gzip-compressed"
    ]

    ids = "${ids.join(' ')}".split()
    reads = "${reads}".split()

    se = ${params.single_end ? "True" : "False"}
    if se:
        fwd = reads
    else:
        fwd = reads[0::2]
        rev = reads[1::2]
        base_args += ["--paired"]
        assert len(fwd) == len(rev)

    assert len(ids) == len(fwd)

    for i, idx in enumerate(ids):
        args = base_args + [
            "--output", f"{idx}.k2",
            "--report", f"{idx}.tsv",
            "--memory-mapping", fwd[i]
        ]
        if not se:
            args.append(rev[i])
        res = run(args)
        if res.returncode != 0:
            if os.path.exists(f"{idx}.k2"):
                os.remove(f"{idx}.k2")
            sys.exit(res.returncode)
    """
}

process count_taxa {
    cpus 4
    memory "8GB"
    time "30m"

    input:
    tuple val(id), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${id}.b2"), path("${lev}/${id}_bracken_mpa.tsv")

    script:
    """
    mkdir ${lev} && \
        sed 's/\\tR1\\t/\\tD\\t/g' ${report} > ${lev}/${report} && \
        bracken -d ${params.kraken2_db} -i ${lev}/${report} \
        -l ${lev} -o ${lev}/${id}.b2 -r ${params.read_length} \
        -t ${params.threshold} -w ${lev}/${id}_bracken.tsv && \
        kreport2mpa.py -r ${lev}/${id}_bracken.tsv -o ${lev}/${id}_bracken_mpa.tsv \
        --no-intermediate-ranks
    """
}

process merge_taxonomy {
    cpus 1
    memory "4GB"
    time "4h"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    tuple val(lev), path(reports)

    output:
    path("${lev}_counts.csv")

    script:
    """
    #!/usr/bin/env python

    from sys import stdin
    from os import path
    import pandas as pd
    import re

    ranks = pd.Series({
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species"
    })

    def str_to_taxa(taxon):
        taxon = taxon.split("|")
        taxa = pd.Series({ranks[t.split("_", 1)[0]]: t.split("_", 1)[1] for t in taxon})
        return taxa

    read = []
    lev = "${lev}"
    input = "${reports.join(",")}"
    paths = input.split(",")

    for p in paths:
        print(f"processing {p}...")
        id = re.findall("(.+)_bracken_mpa.tsv", p)[0]
        try:
            counts = pd.read_csv(p, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            continue
        counts = counts[counts.iloc[:, 0].str.contains(
            str("d" if lev == "D" else lev).lower() + "_")]
        taxa = counts.iloc[:, 0].apply(str_to_taxa)
        taxa["reads"] = counts.iloc[:, 1]
        taxa["sample"] = id
        read.append(taxa.dropna())
    pd.concat(read, sort=False).to_csv("%s_counts.csv" % lev, index=False)
    """
}

process multiqc {
    cpus 1
    memory "8GB"
    time "1h"

    publishDir "${params.data_dir}", mode: "copy", overwrite: true

    input:
    path(taxonomy)

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc ${params.data_dir}/preprocessed ${params.data_dir}/kraken2
    """
}


process megahit {
    cpus 4
    memory "16GB"
    time "12h"

    publishDir "${params.data_dir}/assembled"

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("contigs/${id}.contigs.fa")

    script:
    if (params.single_end)
        """
        megahit -r ${reads} -o contigs -t ${task.cpus} -m ${task.memory.toBytes()} \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        """
    else
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -o contigs -t ${task.cpus} -m ${task.memory.toBytes()} \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        """
}

process find_genes {
    cpus 1
    memory "2GB"
    time "1h"
    publishDir "${params.data_dir}/genes"

    input:
    tuple val(id), path(assembly)

    output:
    tuple val(id), path("${id}.ffn"), path("${id}.faa")

    script:
    """
    if grep -q ">" ${assembly}; then
        prodigal -p meta -i ${assembly} -o ${id}.gff -d ${id}.ffn -a ${id}.faa
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

    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
    memory "16GB"
    time "2h"

    input:
    tuple val(id), path(reads), path(json), path(html), path(genes), path(proteins)

    output:
    path("${id}.sf")

    script:
    if (params.single_end)
        """
        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -r ${reads} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        rm -rf ${id}_index
        """
    else
        """
        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -1 ${reads[0]} -2 ${reads[1]} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        rm -rf ${id}_index
        """
}

process merge_counts {
    cpus 1
    memory "8GB"
    time "4h"
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
    publishDir "${params.data_dir}", mode: "copy", overwrite: true

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
