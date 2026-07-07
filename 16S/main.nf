params.trim_left = 3
params.read_length = 250
params.cut_forward = 5
params.cut_reverse = 20
params.maxEE = 2
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.data_dir = "${launchDir}/data"
params.refs = env("DLP") ? "/home/isilon/dienerlab/refs" : "${launchDir}/refs"
params.taxa_db = "${params.refs}/silva_nr99_v138.2_toGenus_trainset.fa.gz"
params.species_db = "${params.refs}/silva_v138.2_assignSpecies.fa.gz"
params.threads = 16
params.manifest = null
params.pattern = "illumina"

include { find_files; quality_control; trim; denoise; tables; tree } from "./modules/16S.nf"

def helpMessage() {
    log.info"""
    ~~~ Diener Lab 16S Workflow ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf -resume

    General options:
      --data_dir [str]              The main data directory for the analysis (must contain `raw`).
      --read_length [int]           The length of the reads.
      --forward-only [bool]         Run analysis only on forward reads.
      --threads [int]               The maximum number of threads a single process can use.
                                    This is not the same as the maximum number of total threads used.
      --manifest [str]              A manifest file listing the files to be processed. Should be a CSV file with
                                    columns "id", "forward", "reverse" (optional), and "run" (optional). Listing the
                                    sample IDs and read files. If samples were sequenced in different runs indicate
                                    this with the run column.
      --pattern [str]               The file pattern for the FASTQ files. Options are illumina, sra, and simple.
                                    Only used if no manuscript was provided.

    Reference DBs:
      --taxa_db [str]               Path to the default taxonomy database.
      --species_db [str]            Path to species database to perform exact matching to ASVs.

    Quality filter:
      --trim_left [int]             How many bases to trim from the 5' end of each read.
      --cut_forward [int]           How many bases to cut from forward reads at 3' end. Default 5.
      --cut_reverse [int]           How many bases to cut from reverse reads at 3' end. Default 20.
      --maxEE [int]                 Maximum number of expected errors per read.

    Denoising:
      --min_overlap [int]           Minimum overlap between reverse and forward ASVs to merge them.
      --merge [bool]                Whether to merge several runs into a single output.
    """.stripIndent()
}

workflow {
    main:

    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }
    println(params)

    if (params.manifest) {
        manifest = tuple(
            channel.fromPath("${params.manifest}"),
            channel.fromPath("raw")
        )
    } else {
        manifest = find_files(channel.fromPath("${params.data_dir}/raw"))
    }

    manifest | quality_control | trim | denoise | tables
    denoise.out | tree

    publish:
    results = quality_control.out
        .mix(trim.out.map{it -> tuple(it[1], it[2])})
        .mix(denoise.out)
        .mix(tables.out)
        .mix(tree.out)
        .flatten()
}

output {
    results {
        path { file ->
            if (file.extension == "png") {
                return "${params.data_dir}/figures/"
            }
            else if (file.extension == "rds") {
                return "${params.data_dir}/r_data/"
            }
            else if (file.extension == "log") {
                return "${params.data_dir}/logs/"
            }
            else {
                return "${params.data_dir}"
            }
        }
        mode "copy"
        overwrite true
    }
}