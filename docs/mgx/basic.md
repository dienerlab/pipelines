# Metagenomics

## Overview

The metagenomics workflow(s) consist of a basic gene cluster-based workflow followed by
add-ons for binning an replication rate inference.

The basic workflow is a protein/gene-centric workflow. It will quantify taxon abundances
as well as abundances for protein-coding orthologous gene clusters (groups of genes with
putatively similar function across organisms). It is feasible to answer the following
questions:

- Are the taxa whose abundances are related to a phenotype?
- Are there proteins whose genes are more prevalent/abundant in a specific phenotype?
- What functional potential is present in particular metagenomes?

### Feasible data:

- paired or single end metagenomic shotgun sequencing
- any depth (there is no separate shallow pipeline anymore)

!!! WARNING

    The basic workflow currently does not work with Nanopore data due to the salmon step.
    We will add compatibility soon.
    You can still run all the previous steps which are fully compatible.


### Steps:

<div class="annotate" markdown>

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp) (1)
2. Read annotation with [Kraken2](https://github.com/DerrickWood/kraken2) with some custom HPC optimization
3. Taxon counting using [Bracken](https://github.com/jenniferlu717/Bracken)
4. Assembly with [MegaHit](https://github.com/voutcn/megahit)
5. *De novo* gene prediction with [prodigal](https://github.com/hyattpd/Prodigal)
6. Clustering of all genes on protein identity using [mmseqs2 linclust](https://github.com/soedinglab/MMseqs2)
7. Pufferfish mapping index creation (needed for next step)
8. Gene quantification (mapping + counting) with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
9. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)

</div>

1.  *quality reports in HTML and JSON are provided for each file*

## Setup

For the basic workflow set up the `metagenomics` environment.

``` bash
conda env create -f pipelines/metagenomics/conda.yml
```


## Workflow options

``` text
~~~ Diener Lab Metagenomics Workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run main.nf --resume

An exampl erun could look like:
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
```

Additional workflows can be run *after* the basic workflow has finished.

!!! TIP

    On HPC systems we recommend to set the `--batchsize` option to something close to
    the number of samples unless you either have lots (>1000) or problematic samples.

See [concepts](../concepts/batching.md) for more info on batches.