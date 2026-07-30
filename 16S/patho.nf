params.trim_left = 3
params.read_length = 280
params.cut_forward = 5
params.cut_reverse = 20
params.min_score = 20
params.maxEE = 8
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.refs = env("DLP") ? "/home/isilon/dienerlab/refs" : "${launchDir}/refs"
params.taxa_db = "${params.refs}/silva_nr99_v138.2_toGenus_trainset.fa.gz"
params.species_db = "${params.refs}/silva_v138.2_assignSpecies.fa.gz"
params.manifest = null
params.pattern = "patho"
params.run = "latest"
params.upload = false


include { find_files; quality_control; trim; denoise; tables; tree } from "./modules/16S.nf"

def helpMessage() {
    log.info"""
    ~~~ Diener Lab 16S Workflow for Pathology Sequencing ~~~

    Usage:
    A run using all,default parameters can be started with:
    > nextflow run main.nf -resume

    General options:
      --read_length [int]           The length of the reads.
      --forward-only [bool]         Run analysis only on forward reads.
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

    if (!params.run) {
        log.error "No run specified. Please provide a run name with --run."
        exit 1
    }


    runID = channel.of(params.run)
    manifest = download_raw_files(runID) | find_files | annotate_samples
    length_check(manifest)
    manifest | quality_control | trim | denoise | tables
    denoise.out | tree

    report(
        channel.fromPath("${projectDir}/report.qmd")
        .mix(denoise.out.map{it -> it[1]})
        .mix(tree.out.map{it -> it[1]})
        .mix(quality_control.out.map{it -> it[2]})
        .mix(length_check.out)
        .mix(download_raw_files.out)
        .collect()
    )

    merged = quality_control.out.map{it -> tuple(it[2], it[3])}
        .mix(trim.out.map{it -> tuple(it[1], it[2])})
        .mix(length_check.out)
        .mix(denoise.out)
        .mix(tables.out)
        .mix(tree.out)
        .mix(report.out)
        .flatten()

//    upload(merged.collect())


    publish:
    results = merged
    reports = report.out.flatten()
}

output {
    results {
        path { file ->
            if (file.extension == "png") {
                return "figures/"
            }
            else if (file.extension == "rds") {
                return "r_data/"
            }
            else if (file.extension == "log") {
                return "logs/"
            }
            else if (file.extension == "tree") {
                return "trees/"
            }
            else if (file.extension == "csv") {
                return "tables"
            }
            else {
                return ""
            }
        }
        mode "copy"
        overwrite true
    }

    reports {
        mode "copy"
        overwrite true
    }
}

process download_raw_files {
    cpus 1
    memory 512.MB
    time "1h"

    input:
    val runID

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

process annotate_samples {
    cpus 1
    memory 512.MB
    time "1h"

    input:
    tuple path(manifest), path(raw_dir)

    output:
    tuple path("manifest_annotated.csv"), path(raw_dir)

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    files <- read_csv("${manifest}") |> mutate(Barcode = as.character(id)) |> select(!id)
    man <- readxl::read_excel(Sys.glob("raw/*.xlsx")[1], skip=9) |>
        mutate(Barcode = str_split_i(Barcode, " ", 2), id = `Externe ID`) |>
        drop_na(id) |>
        mutate(type = c("sample", "control")[str_detect(id, "pos|neg") + 1])
    merged <- man |> inner_join(files, by="Barcode")
    write_csv(merged, "manifest_annotated.csv")
    """
}

process length_check {
    cpus 1
    memory 512.MB
    time "1h"

    input:
    tuple path(manifest), path(raw_dir)

    output:
    tuple path("amplicon_types.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(ShortRead)
    library(data.table)

    AMPLEN <- ${params.read_length} + 1

    man <- fread("${manifest}")
    res <- list()
    for (fq in man[["forward"]]) {
        reads <- readFastq(fq)
        lengths <- width(sread(reads))
        target_amplicons = sum(between(lengths, AMPLEN - 10, AMPLEN + 10))
        mitochondrial_amplicons = sum(between(lengths, 197, 217))
        daisy_chains = sum((lengths[lengths > AMPLEN + 100] %% AMPLEN) <= 10)
        res[[fq]] <- data.table(
            id = man[forward == fq, id],
            total_reads = length(reads),
            avg_length = mean(lengths),
            target = target_amplicons,
            mitochondrial = mitochondrial_amplicons,
            daisy_chains = daisy_chains,
            other = length(reads) - target_amplicons - mitochondrial_amplicons - daisy_chains,
            target_fraction = target / length(reads),
            mitochondrial_fraction = mitochondrial_amplicons / length(reads),
            daisy_chain_fraction = daisy_chains / length(reads)
        )
    }
    rbindlist(res) |> fwrite("amplicon_types.csv")
    """
}

process report {
    cpus 1
    memory 4.GB
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
    memory 512.MB
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


