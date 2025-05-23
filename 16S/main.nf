params.trim_left = 3
params.read_length = 250
params.trunc_forward = params.read_length - 5
params.trunc_reverse = params.read_length - 20
params.maxEE = 2
params.merge = true
params.min_overlap = 8
params.forward_only = false
params.data_dir = "${launchDir}/data"
params.taxa_db = "${launchDir}/refs/GTDB_bac120_arc53_ssu_r220_genus.fa.gz"
params.species_db = "${launchDir}/refs/GTDB_bac120_arc53_ssu_r220_species.fa.gz"
params.threads = 16
params.manifest = null
params.pattern = "illumina"

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
      --trunc_forward [int]         Where to truncate forward reads. Default length - 5
      --trunc_reverse [int]         Where to truncate reverse reads. Default length - 20.
      --maxEE [int]                 Maximum number of expected errors per read.

    Denoising:
      --min_overlap [int]           Minimum overlap between reverse and forward ASVs to merge them.
      --merge [bool]                Whether to merge several runs into a single output.
    """.stripIndent()
}

params.help = false

workflow {
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (params.manifest) {
        manifest = Channel.fromPath("${params.manifest}")
    } else {
        manifest = find_files()
    }
    manifest | quality_control | trim | denoise | tables
    denoise.output | tree
}


process find_files {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    cpus 1
    memory "4 GB"
    time "10 m"

    output:
    path("manifest_autogen.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(miso)

    files <- find_read_files(
        "${params.data_dir}/raw",
        format=file_formats[["${params.pattern}"]],
        dirs_are_runs = T
    )

    if (nrow(files) == 0) {
        flog.error("Could not find any read files matching the ${params.pattern} pattern.")
        stop("No files to process.")
    }

    fwrite(files, "manifest_autogen.csv")
    """
}

process quality_control {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    cpus params.threads
    memory "32 GB"
    time "8h"

    input:
    path(manifest)

    output:
    tuple path(manifest), path("qc.rds"), path("*.png"), path("qc.log")

    script:
    """
    #!/usr/bin/env Rscript

    options(mc.cores = ${task.cpus})

    library(Biostrings)
    library(miso)
    library(futile.logger)

    PREFIX = "${params.data_dir}/raw"

    flog.appender(appender.file("qc.log"))

    files <- fread("${manifest}")
    if (!file.exists(files[["forward"]][1])) {
        flog.warn("Files do not seem to exist, looking in %s.", PREFIX)
        files[, "forward" := file.path(PREFIX, forward)]
        if ("reverse" %in% names(files)) {
            files[, "reverse" := file.path(PREFIX, `reverse`)]
        }
    }

    if (${params.forward_only ? "T" : "F"}) {
        files[, "reverse" := NULL]
    }
    if ("reverse" %in% names(files)) {
        good <- !(is.na(files[["forward"]]) | is.na(files[["reverse"]]))
        bad <- files[!file.exists(forward) | !file.exists(`reverse`)]
        files <- files[good]
    } else {
        files <- files[!is.na(forward)]
        bad <- files[!file.exists(forward)]
    }

    # Check if all files are there
    if (nrow(bad) > 0) {
        flog.error("The following IDs have missing files: %s", paste(bad[, id], sep=","))
        print("The following files do not exist at the specified location:")
        print(bad)
        stop("Can't continue with missing files :(")
    }

    qc <- quality_control(files, min_score = 20)
    saveRDS(qc, "qc.rds")
    ggsave("qualities.png", pl = qc[["quality_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    ggsave("entropy.png", pl = qc[["entropy_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    ggsave("lengths.png", pl = qc[["length_plot"]] + theme_minimal(),
           width = 8, height = 4, dpi = 300)
    """
}

process trim {
    publishDir "${params.data_dir}"
    cpus params.threads
    memory "16 GB"
    time "24h"

    input:
    tuple path(manifest), path(qc), path(pl), path(log)

    output:
    tuple path("preprocessed"), path("preprocessed.rds"), path("trim.log")

    script:
    """
    #!/usr/bin/env Rscript
    library(miso)
    library(futile.logger)

    flog.appender(appender.file("trim.log"))

    qc <- readRDS("${qc}")
    manifest <- fread("${manifest}")

    if ("reverse" %in% names(manifest)) {
        trunc <- c(${params.trunc_forward}, ${params.trunc_reverse})
    } else {
        trunc <- ${params.trunc_forward}
    }

    procced <- preprocess(
        qc,
        trimLeft = ${params.trim_left},
        truncLen = trunc,
        maxEE = ${params.maxEE},
        out_dir = "preprocessed",
        threads = ${task.cpus}
    )
    saveRDS(procced, "preprocessed.rds")
    """
}

process denoise {
    publishDir "${params.data_dir}/denoise"
    cpus params.threads
    memory "64 GB"
    time "2d"

    input:
    tuple path(procced), path(artifact), path(log)

    output:
    tuple path("read_stats.csv"), path("denoised.rds"), path("phyloseq.rds"), path("denoise.log")

    script:
    """
    #!/usr/bin/env Rscript
    library(miso)
    library(futile.logger)

    flog.appender(appender.file("denoise.log"))

    procced <- readRDS("${artifact}")
    passed_reads <- procced[["passed"]][["preprocessed"]]
    procced[["files"]] <- procced[["files"]][passed_reads > 0]
    denoised <- denoise(
        procced,
        hash = T,
        threads = ${task.cpus},
        merge = ${params.merge ? "T" : "F"},
        min_overlap = ${params.min_overlap},
        taxa_db = "${params.taxa_db}",
        species_db = "${params.species_db}",
        nbases=1e8
    )
    saveRDS(denoised, "denoised.rds")
    fwrite(denoised[["passed_reads"]], "read_stats.csv")
    sdata <- as.data.frame(procced[["files"]])
    rownames(sdata) <- sdata[, "id"]
    ps <- as_phyloseq(denoised, sdata)
    saveRDS(ps, "phyloseq.rds")
    """
}

process tree {
    publishDir "${params.data_dir}", mode: 'copy', overwrite: true

    cpus params.threads
    memory "64 GB"
    time "24h"

    input:
    tuple path(stats), path(denoised), path(ps), path(log)

    output:
    tuple path("asvs.tree"), path("phyloseq_with_tree.rds"), path("tree.log")

    script:
    """
    #!/usr/bin/env Rscript

    library(futile.logger)
    library(miso)

    flog.appender(appender.file("tree.log"))

    denoised <- readRDS("${denoised}")

    seqs <- denoised[["taxonomy"]][, "sequence"]
    alignments <- DECIPHER::AlignSeqs(
        Biostrings::DNAStringSet(seqs),
        processors = ${task.cpus}
    )

    am <- as(alignments, "matrix")
    freqs <- t(apply(am, 2, function(x)
        table(factor(x, levels = c("-", "A", "C", "G", "T"))) / length(x)))
    absent <- freqs[, "-"] > 0.5
    flog.info(paste("%d alignment columns are absent in >50%% of alignments,",
                    "removing them."), sum(absent))
    alignments <- am[, !absent] %>% apply(1, paste0, collapse = "") %>%
                  Biostrings::DNAStringSet()
    flog.info("Final alignment has %d positions. Starting FastTree...",
              width(alignments)[1])
    Biostrings::writeXStringSet(alignments, "asvs_aligned.fna")
    args <- c("-fastest", "-gtr", "-gamma",
              "-nt", "asvs_aligned.fna", ">", "asvs.tree")
    out <- system2("FastTreeMP",
                   args = args,
                   env = "OMP_NUM_THREADS=${task.cpus}"
            )
    tree <- read_tree("asvs.tree")

    ps <- readRDS("${ps}")
    phy_tree(ps) <- tree

    saveRDS(ps, "phyloseq_with_tree.rds")
    """
}


process tables {
    publishDir "${params.data_dir}", mode: "copy", overwrite: true
    cpus 1
    memory "16 GB"
    time "1h"

    input:
    tuple path(stats), path(arti), path(ps), path(log)

    output:
    tuple path("asvs.csv"), path("taxonomy.csv")

    script:
    """
    #!/usr/bin/env Rscript
    library(miso)

    denoised <- readRDS("${arti}")
    ids <- rownames(denoised[["feature_table"]])
    asvs <- as.data.table(denoised[["feature_table"]])[, "id" := ids]
    asvs <- melt(asvs, id.vars="id", variable.name="hash",
                 value.name="count")[count > 0]
    fwrite(asvs, "asvs.csv")
    ids <- rownames(denoised[["taxonomy"]])
    tax <- as.data.table(denoised[["taxonomy"]])[, "id" := ids]
    fwrite(tax, "taxonomy.csv")
    """
}
