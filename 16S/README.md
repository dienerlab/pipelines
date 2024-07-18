# 16S amplicon sequencing workflow

**Feasible data:**

- single-end or paired-end 16S amplicon sequencing data
- decent depth (>10K reads per sample)

## Installation

### Option 1: Everything in a conda env

```bash
mamba env create -f conda.yml
conda activate 16S
Rscript -e "remotes::install_github('dienerlab/miso')"
```

### Options 2: Separate R everything else in conda


> [!TIP]
> This option will give you the fastest pipeline because natively compiled DADA2 is
> at least twice as fast as the conda one. The installation will take a while though.

```bash
mamba env create -f conda-no-r.yml
Rscript -e "remotes::install_github('dienerlab/miso')"
```

### Download taxonomy databases

`cd` into your woirk directory, then

```bash
mkdir -p refs
wget https://zenodo.org/records/10403693/files/GTDB_bac120_arc53_ssu_r214_genus.fa.gz?download=1 -O GTDB_bac120_arc53_ssu_r214_genus.fa.gz
wget https://zenodo.org/records/10403693/files/GTDB_bac120_arc53_ssu_r214_species.fa.gz?download=1 -O GTDB_bac120_arc53_ssu_r214_species.fa.gz
```

## Basic workflow

**Definition**: `main.nf`

```text
~~~ Diener Lab 16S Workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run main.nf -resume

General options:
    --data_dir [str]              The main data directory for the analysis (must contain `raw`).
    --read_length [str]           The length of the reads.
    --forward-only [bool]         Run analysis only on forward reads.
    --threads [int]               The maximum number of threads a single process can use.
                                This is not the same as the maximum number of total threads used.
    --pattern [str]               The file pattern for the FASTQ files. Options are illumina, sra, and simple.
Reference DBs:
    --taxa_db [str]               Path to the default taxonomy database.
    --species_db [str]            Path to species database to perform exact matching to ASVs.
    --eggnogg_refs [str]          Where to find EGGNOG references. Defaults to <refs>/eggnog.
    --kraken2_db [str]            Where to find the Kraken2 reference. Defaults to <refs>/kraken2_default.
    --kraken2_mem [str]           The maximum amount of memory for Kraken2. If not set will choose this automatically
                                based on the database size. Thus, only use to limit Kraken2 to less memory.
Quality filter:
    --trim_left [str]             How many bases to trim from the 5' end of each read.
    --trunc_forward [int]         Where to truncate forward reads. Default length - 5
    --trunc_reverse [int]         Where to truncate reverse reads. Default length - 20.
    --maxEE                       Maximum number of expected errors per read.
    --threshold [str]             Smallest abundance threshold used by Kraken.
Denoising:
    --min_overlap [int]           Minimum overlap between reverse and forward ASVs to merge them.
    --merge [bool]                Whether to merge several runs into a single output.
    --overlap [double]            Minimum required overlap between contigs.
```

Most likely run using:

```bash
nextflow run main.nf -resume -with-conda [PATH/TO/CONDA]/envs/16S[-minimal]
```

### Folder structure

```
project/
    > [WORKFLOW].nf
    > data
        > raw
            > run1
                | sample1_S1_L001_R1_001.fastq.gz
                | ...
            > run2
                | sample75_S1_L001_R1_001.fastq.gz
                | ...
        > step1
        > step2
        > asvs.csv
        > ...
        > qc.png
        > ...
```

If you have only one run put all raw FASTQ files under `raw`.

### Steps:

1. Automatic detection of Illumina read files and management of multiple runs
2. Trimming, filtering and quality metrics (base qualities, entropy, lengths) with DADA2 and mbtools
3. Denoising using DADA2
4. Taxonomy assignment using DADA2 and GTDB 214
5. Alignment of ASVs with DECIPHER and building of a phylogenetic tree from core alignment (FastTree).
6. Creation of tabular and phyloseq output

