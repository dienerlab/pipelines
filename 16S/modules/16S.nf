
process find_files {
    cpus 1
    memory "4 GB"
    time "10 m"

    input:
    path("raw")

    output:
    tuple path("manifest_autogen.csv"), path("raw")

    script:
    """
    #!/usr/bin/env Rscript

    library(miso)

    formats <- file_formats
    formats[["patho"]] <- list(
        pattern="([A-Za-z0-9\\\\-]+)__IonXpress_(\\\\d+)_.+.basecaller.fastq.gz",
        annotations=c("patho_id", "id")
    )

    files <- find_read_files(
        "raw/",
        format=formats[["${params.pattern}"]],
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
    cpus params.threads
    memory "32 GB"
    time "8h"

    input:
    tuple path(manifest), path(raw_dir)

    output:
    tuple path(manifest), path(raw_dir), path("qc.rds"), path("*.png")

    script:
    """
    #!/usr/bin/env Rscript

    options(mc.cores = ${task.cpus})

    library(Biostrings)
    library(miso)
    library(futile.logger)

    PREFIX = "raw_dir"

    flog.appender(appender.file("qc.log"))

    files <- fread("${manifest}")[, "id" := as.character(id)]
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
    cpus params.threads
    memory "16 GB"
    time "24h"

    input:
    tuple path(manifest), path(raw_dir), path(qc), path(pl)

    output:
    tuple path("preprocessed"), path("preprocessed.rds"), path("trim.log")

    script:
    """
    #!/usr/bin/env Rscript
    library(miso)
    library(futile.logger)

    flog.appender(appender.file("trim.log"))

    qc <- readRDS("${qc}")
    manifest <- fread("${manifest}")[, "id" := as.character(id)]

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
    cpus params.threads
    memory "64 GB"
    time "2d"

    input:
    tuple path(procced), path(artifact), path(log)

    output:
    tuple path("read_stats.csv"), path("denoised.rds"), path("phyloseq.rds"), path("denoise.log"), path("denoise_error_model_*.png")

    script:
    """
    #!/usr/bin/env Rscript
    library(miso)
    library(futile.logger)
    library(ggplot2)

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

    for (run in names(denoised[["errors"]])) {
        p <- denoised[["errors"]][[run]]
        ggsave(
            paste0("denoise_error_model_", run, ".png"),
            pl = p + theme_minimal(),
            width = 8, height = 4, dpi = 300
        )
    }
    """
}

process tree {
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
