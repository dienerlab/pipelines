#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.data_dir = "${launchDir}/data"
params.raw_data = "raw"
params.refs = env("DLP") ? "/home/isilon/dienerlab/refs" : "${launchDir}/refs"
params.eggnog_refs = "${params.refs}/eggnog"
params.metapackage = "${params.refs}/GlobDB_r232.metapackage_v4.smpkg"
params.viralpackage = "${params.refs}/lyrebird_v0.3.1_phrog_v4.1_metapackage_20250720.smpkg.zb"

params.single_end = false
params.trim_front = 3
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threshold = 10
params.contig_length = 1000
params.overlap = 0.8
params.identity = 0.97
params.method = "illumina"

params.assemblyPreset = "default"

params.neighbors = 8
params.metric = "braycurtis"
params.neighborDomain = "microbial"


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

    // Prepare sample sheets for mag
    neighborhoods(summarizeProfiles.out)
    sample_sheet(preprocess.out.collect{it -> it[1]}, neighborhoods.out)
    assembly_sheet(assemble.out.collect{it -> it[1]}, neighborhoods.out)

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
    sheets = sample_sheet.out.mix(assembly_sheet.out)
    pf = mergePF.out
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

    sheets {
        mode "copy"
        overwrite true
    }

    pf {
        mode "copy"
        overwrite true
    }
}

process preprocess {
    cpus 3
    memory "4GB"
    time "30m"
    tag { id }

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
    tag { id }

    input:
    tuple val(id), path(fastqs), path(json), path(html)

    output:
    tuple val(id), path("${id}_microbial_profile.tsv"), path("${id}_microbial_otus.tsv"), path("${id}_spf.tsv")

    script:
    if (params.single_end)
        """
        singlem pipe -1 ${fastqs} --threads ${task.cpus} \
            --metapackage ${params.metapackage} \
            -p ${id}_microbial_profile.tsv \
            --otu-table ${id}_microbial_otus.tsv

        singlem prokaryotic_fraction -1 ${fastqs} \
            -p ${id}_microbial_profile.tsv > ${id}_spf.tsv
        """
    else
        """
        singlem pipe -1 ${fastqs[0]} -2 ${fastqs[1]} --threads ${task.cpus} \
            --metapackage ${params.metapackage} \
            -p ${id}_microbial_profile.tsv \
            --otu-table ${id}_microbial_otus.tsv

        singlem prokaryotic_fraction -1 ${fastqs[0]} -2 ${fastqs[1]} \
            --metapackage ${params.metapackage} \
            -p ${id}_microbial_profile.tsv > ${id}_spf.tsv
        """
}

process lyrebird {
    cpus 3
    memory 8.GB
    time 4.h
    tag { id }

    input:
    tuple val(id), path(fastqs), path(json), path(html)

    output:
    tuple val(id), path("${id}_viral_profile.tsv"), path("${id}_viral_otus.tsv")

    script:
    if (params.single_end)
        """
        lyrebird pipe -1 ${fastqs} --threads ${task.cpus} \
            --metapackage ${params.viralpackage} \
            -p ${id}_viral_profile.tsv \
            --otu-table ${id}_viral_otus.tsv
        """
    else
        """
        lyrebird pipe -1 ${fastqs[0]} -2 ${fastqs[1]} --threads ${task.cpus} \
            --metapackage ${params.viralpackage} \
            -p ${id}_viral_profile.tsv \
            --otu-table ${id}_viral_otus.tsv
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
    memory 16.GB
    time 12.h
    tag { id }

    input:
    tuple val(id), path(reads), path(json), path(report)

    output:
    tuple val(id), path("contigs/${id}.contigs.fa")

    script:
    def args = params.assemblyPreset == "default" ? "" : " --preset ${params.assemblyPreset}"
    if (params.single_end && params.method == "illumina")
        """
        megahit -r ${reads} -o contigs -t ${task.cpus} -m ${task.memory.toBytes()} \
                ${args} \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
        """
    else if (!params.single_end && params.method == "illumina")
        """
        megahit -1 ${reads[0]} -2 ${reads[1]} -o contigs -t ${task.cpus} -m ${task.memory.toBytes()} \
                ${args} \
                --min-contig-len ${params.contig_length} --out-prefix ${id}
        sed -i -e "s/^>/>${id}_/" contigs/${id}.contigs.fa
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
    tag { id }

    input:
    tuple val(id), path(assembly)

    output:
    tuple val(id), path("${id}.ffn"), path("${id}.faa")

    script:
    """
    if grep -q ">" ${assembly}; then
        pyrodigal -p meta -i ${assembly} -o ${id}.gff -d ${id}.ffn -a ${id}.faa
    else
        touch ${id}.faa
        touch ${id}.ffn
    fi
    """
}

process cluster_proteins {
    cpus 12
    memory "40GB"
    time "2h"

    input:
    path(proteins)

    output:
    tuple path("proteins.faa"), path("proteins_cluster.tsv")

    script:
    """
    trap "rm -rf all.faa" EXIT

    cat ${proteins} > all.faa

    diamond cluster -d all.faa -o proteins_cluster.tsv \
        --id ${params.identity} --member-cover ${params.coverage} \
        -M ${task.memory.toGiga()}G -p ${task.cpus}

    seqkit grep -f <(cut -f 1 proteins_cluster.tsv) all.faa > proteins.faa
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
    seqkit grep -f <(seqkit seq -n -i ${proteins}) \
        ${transcripts} -o transcripts.fna.gz
    """
}

process map_and_count {
    cpus 2
    memory "32 GB"
    time "4h"
    tag { id }

    input:
    tuple val(id), path(reads), path(json), path(html), path(genes), path(proteins)

    output:
    path("${id}.sf")

    script:
    if (params.single_end && params.method == "illumina")
        """
        trap "rm -rf ${id}_index" EXIT

        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -r ${reads} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        """
    else if (!params.single_end && params.method == "illumina")
        """
        trap "rm -rf ${id}_index" EXIT

        salmon index -p ${task.cpus} -t ${genes} -i ${id}_index || touch ${id}_index
        salmon quant --meta -p ${task.cpus} -l A -i ${id}_index -1 ${reads[0]} -2 ${reads[1]} -o ${id} &&
            mv ${id}/quant.sf ${id}.sf || touch ${id}.sf
        """
    else if (params.method == "nanopore" || params.method == "pacbio")
        """
        trap "rm -rf ${id}.bam" EXIT

        minimap2 -ax map-ont -p 1.0 -N 100 -t ${task.cpus} ${genes} ${reads} | samtools view -bS > ${id}.bam
        salmon quant -t ${genes} -q --ont --meta-l U -a ${id}.bam -o ${id} -p ${task.cpus} &&
            mv ${id}_salmon/quant.sf ${id}.sf || touch ${id}.sf
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
    cpus 12
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

process sample_sheet {
    cpus 1
    memory "1GB"
    time "1h"

    input:
    path(reads)
    path(neighbors)

    output:
    path("samplesheet.csv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from pathlib import Path

    reads = "${reads}".split()
    forward = sorted(["preprocessed/" + r for r in reads if "_R1" in reads])
    reverse = sorted(["preprocessed/" + r for r in reads if "_R2" in reads])
    ids = [r.split("_filtered_R")[0] for r in forward]
    if len(reverse) == len(forward):
        df = pd.DataFrame({
            "sample": ids,
            "short_reads_1": forward,
            "short_reads_2": reverse,
            "short_reads_platform": "${params.method}".upper()
        })
    else:
        df = pd.DataFrame({
            "sample": ids,
        })

        if "${params.method}" == "illumina":
            df["short_reads_1"] = forward
            df["short_reads_platform"] = "${params.method}".upper()
        else:
            df["long_reads"] = forward
            df["long_reads_platform"] = "OXFORD_NANOPORE_HQ" if "${params.method}" == "nanopore" else "PACBIO_HFI"

    neighbors = pd.read_csv("${neighbors}")
    df = df.merge(neighbors, on="sample", how="inner")

    df.to_csv("samplesheet.csv", index=False)
    """
}

process assembly_sheet {
    cpus 1
    memory "1GB"
    time "1h"

    input:
    path(assemblies)
    path(neighbors)

    output:
    path("assembly_sheet.csv")

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from pathlib import Path

    assemblies = sorted("${assemblies}".split())
    ids = ["assemblies/" + a.split(".contigs")[0] for a in assemblies]
    df = pd.DataFrame({
        "id": ids,
        "group": range(len(ids)),
        "assembler": "megahit" if "${params.method}" == "illumina" else "metaMDBG",
        "fasta": assemblies
    })

    neighbors = pd.read_csv("${neighbors}")
    df = df.merge(neighbors, left_on="id", right_on="sample", how="inner")

    df.to_csv("assembly_sheet.csv", index=False)
    """
}

process neighborhoods {
    cpus 1
    memory "8GB"
    time "2h"

    input:
    path(abundances)

    output:
    path("${params.neighborDomain}_neighborhood.csv")

    script:
    """
    #!/usr/bin/env python

    import numpy as np
    import pandas as pd
    from scipy.spatial.distance import pdist, squareform

    files = "${abundances}".split()
    domains = {s.split("_")[0]: s for s in files}
    abundances = domains.get("${params.neighborDomain}", None)
    if abundances is None:
        raise ValueError(f"Domain ${params.neighborDomain} not found in abundances: {list(domains.keys())}")

    df = pd.read_csv(abundances, sep="\\t")
    df = df[df["level"] == "species"]
    mat = df.pivot_table(
        index="sample", columns="taxonomy",
        values="relative_abundance", fill_value=0
    )

    def group_samples_fixed_size(df, k, metric="braycurtis"):
        '''Groups a pandas DataFrame of samples into clusters of fixed size k

        minimizing the pairwise distance within each group.

        Parameters:
        - df: pd.DataFrame with samples as index and relative abundances as columns.
        - k: int, desired size of each group.
        - metric: str or callable, distance metric to pass to scipy.spatial.distance.pdist
                (e.g., 'braycurtis', 'euclidean', 'cosine', 'jensenshannon').

        Returns:
        - pd.DataFrame with columns ['sample', 'group', 'dist']
        '''
        # 1. Compute pairwise distance matrix using the specified metric
        dist_array = pdist(df.values, metric=metric)
        dist_matrix = squareform(dist_array)

        samples = df.index.tolist()
        n_samples = len(samples)

        unassigned = set(range(n_samples))
        group_records = []
        group_id = 0

        # 2. Greedy allocation loop
        while len(unassigned) >= k:
            unassigned_list = list(unassigned)

            # Find the pair of unassigned samples with the minimum distance to seed the group
            sub_dist = dist_matrix[np.ix_(unassigned_list, unassigned_list)]
            np.fill_diagonal(sub_dist, np.inf)
            min_idx = np.unravel_index(np.argmin(sub_dist), sub_dist.shape)

            seed1 = unassigned_list[min_idx[0]]
            seed2 = unassigned_list[min_idx[1]]

            current_group = [seed1, seed2]
            unassigned.remove(seed1)
            unassigned.remove(seed2)

            # Greedily add the closest remaining samples until group reaches size k
            while len(current_group) < k and unassigned:
                unassigned_list = list(unassigned)
                mean_distances = np.mean(
                    dist_matrix[unassigned_list][:, current_group], axis=1
                )
                closest_idx = unassigned_list[np.argmin(mean_distances)]
                current_group.append(closest_idx)
                unassigned.remove(closest_idx)

            # Calculate mean distance for each sample to all *other* members of its group
            for idx in current_group:
                others = [o for o in current_group if o != idx]
                if others:
                    mean_dist = np.mean(dist_matrix[idx, others])
                else:
                    mean_dist = 0.0  # Case where k=1

                group_records.append(
                    {
                        "sample": samples[idx],
                        "group": group_id,
                        "dist": mean_dist,
                    }
                )

            group_id += 1

        # 3. Handle leftover samples if total samples is not a multiple of k
        if unassigned:
            unassigned_list = list(unassigned)
            for idx in unassigned_list:
                # If the final group has other members, compute mean distance to them
                others = [o for o in unassigned_list if o != idx]
                if others:
                    mean_dist = np.mean(dist_matrix[idx, others])
                if group_members:
                    mean_dist = np.mean(dist_matrix[idx, group_members])
                else:
                    mean_dist = 0.0

                group_records.append(
                    {
                        "sample": samples[idx],
                        "group": group_id,
                        "dist": mean_dist,
                    }
                )

        # Convert to DataFrame and sort to match original index order
        result_df = pd.DataFrame(group_records)

        return result_df

    res = group_samples_fixed_size(mat, k=${params.neighbors}, metric='${params.metric}')
    res.to_csv("${params.neighborDomain}_neighborhood.csv", index=False)
    """
}
