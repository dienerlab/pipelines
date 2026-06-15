params.trim_left = 3
params.read_length = 280
params.trunc_forward = params.read_length - 5
params.trunc_reverse = params.read_length - 20
params.maxEE = 8
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.refs = env("DLP") ? "/home/isilon/dienerlab/refs" : "${launchDir}/refs"
params.taxa_db = "${params.refs}/silva_nr99_v138.2_toGenus_trainset.fa.gz"
params.species_db = "${params.refs}/silva_v138.2_assignSpecies.fa.gz"
params.threads = 16
params.manifest = null
params.pattern = "patho"
params.run = "latest"
params.data_dir = "${launchDir}/${params.run.replaceAll('-', '').split('__')[0]}"


include { find_files; quality_control; trim; denoise; tables; tree } from "./modules/16S.nf"

def helpMessage() {
    log.info"""
    ~~~ Diener Lab 16S Workflow for Pathology Sequencing ~~~

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
      --trunc_forward [int]         Where to truncate forward reads. Default length - 5
      --trunc_reverse [int]         Where to truncate reverse reads. Default length - 20.
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

    if (!params.run) {
        log.error "No run specified. Please provide a run name with --run."
        exit 1
    }


    log.info "Will save results to '${params.data_dir}'."

    manifest = download_raw_files | find_files
    manifest | quality_control | trim | denoise | tables
    denoise.out | tree

    report(
        channel.fromPath("${projectDir}/report.qmd")
        .mix(denoise.out.map{it -> it[1]})
        .mix(tree.out.map{it -> it[1]})
        .mix(quality_control.out.map{it -> it[2]})
        .mix(download_raw_files.out)
        .collect()
    )

    merged = quality_control.out.map{it -> tuple(it[2], it[3])}
        .mix(trim.out.map{it -> tuple(it[1], it[2])})
        .mix(denoise.out)
        .mix(tables.out)
        .mix(tree.out)
        .mix(report.out)
        .flatten().collect()
   
    upload(merged)


    publish:
    results = merged
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
            else if (file.extension == "tree") {
                return "${params.data_dir}/trees/"
            }
            else {
                return "${params.data_dir}/"
            }
        }
        mode "copy"
        overwrite true
    }
}

process download_raw_files {
    cpus 1
    memory "4 GB"
    time "1h"

    output:
    path("raw")

    script:
    """
    mkdir raw
    pre=\$(Rscript -e "paste0(format(Sys.Date(), '%Y'), quarters(Sys.Date())) |> cat()")
    rclone copy -P \
        "nextcloud:/Analysisresult_Sequenzierung_Hygiene_16s_Diener/input_\$pre/${params.run}" \
        raw
    """
}


process report {
    cpus 1
    memory "8 GB"
    time "1h"

    input:
    tuple path(template), path(denoised), path(ps_with_tree), path(qc), path(raw)

    output:
    tuple path("report.html"), path("figures/*.*"), path("tables/*.*")

    script:
    """
    mkdir r_data && mv *.rds r_data && mkdir figures && mkdir tables
    quarto render ${template} --execute --to html --output report.html
    """
}

process upload {
    cpus 1
    memory 4.GB
    time 1.h

    input:
    path(files)

    script:
    """
    mkdir out && mv ${files} out && cd out
    mkdir r_data && mv *.rds r_data
    mkdir tables && mv *.csv tables
    mkdir figures && mv *.png figures
    mkdir trees && mv *.tree trees
    mkdir logs && mv *.log logs
    rclone copy -P -L . \
        "nextcloud:/Patho 16S sequencing/${params.run.replaceAll('-', '').split('__')[0]}"
    """
}


