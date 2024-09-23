# Metagenomics

**Feasible data:**

- paired or single end metagenomic shotgun sequencing
- any depth (there is no separate shallow pipeline anymore)

The metagenomics workflow(s) consist of a basic gene cluster-based workflow followed by
add-ons for binning an replication rate inference.

## General file structure

```
project root
├─ data
│  ├─ raw
│  │  ├─ sample1_R1.fastq.gz
│  │  ├─ sample1_R2.fastq.gz
│  │  └─ ...
│  └─ ...
└─ refs
   ├─ eggnog
   └─ kraken2
```

All generated data will be placed into the `data` folder. Raw sequencing data is expected
in `data/raw`. Refernce databases can be placed anywhere and set with command line arguments.
However, for readability we recommend to create a `refs` folder and place symlinks to the
databases there. In that case the pipeline will pick those up automatically.

## Setup

For the basic workflow set up the `metagenomics` environment.

```bash
conda env create -f conda.yml
```

If you want to run metagenomic binning as well also set up the `binchecks` environment.

```bash
conda env create binchecks.yml
```

## Basic (functional) workflow

**Definition**: `main.nf`

> [!IMPORTANT]
> The basic workflow currently does not work with Nanopore data due to the salmon step.
> We will add compatibility soon.
> You can still run all the previous steps which are fully compatible.

### Steps:

1. Adapter and quality trimming with [fastp](https://github.com/OpenGene/fastp)
   *quality reports in HTML and JSON are provided for each file*
2. Read annotation with [Kraken2](https://github.com/DerrickWood/kraken2) with some custom HPC optimization
3. Taxon counting using [Bracken](https://github.com/jenniferlu717/Bracken)
4. Assembly with [MegaHit](https://github.com/voutcn/megahit)
5. *De novo* gene prediction with [prodigal](https://github.com/hyattpd/Prodigal)
6. Clustering of all genes on protein identity using [mmseqs2 linclust](https://github.com/soedinglab/MMseqs2)
7. Pufferfish mapping index creation (needed for next step)
8. Gene quantification (mapping + counting) with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)
9. Protein annotation using the [EGGNoG mapper](https://github.com/eggnogdb/eggnog-mapper)

### Options:

```
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

> [!TIP]
> On HPC systems we recommend to set the `--batchsize` option to something close to
> the number of samples unless you either have lots (>1000) or problematic samples.

See [concepts](../docs/concepts.md) for more info on batches.

## Replication rates

**Definition**: replication.nf

1. Alignment to [~3K high quality assemblies](https://www.nature.com/articles/s41586-019-1058-x) from the gut microbiome with bowtie2
2. extraction of coverage maps
3. Quanitifaction of peak-to-trough ratios (PTRs) with [coptr](https://github.com/tyjo/coptr)

This will use the preprocessed data from the basic workflow. The output will be a single
CSV file called `data/annotated_rates.csv` that contains the replication rates and
species annotations.

### Options:

```
~~~ Diener Lab Replication Rate workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run replication.nf --resume

A run with all parametrs set would look like:
> nextflow run main.nf --data_dir=./data --single_end true --threads 12 --min_reads 2500 --IGG /refs/IGG

General options:
  --data_dir [str]              The main data directory for the analysis (must contain `raw`).
  --single_end [bool]           Whether the data is single-end sequencing data.

COPTR options:
  --min_reads [int]             Minimum number of reads for a genome to calculate PTRs.
  --IGG [path]                  Location of the IGG bowtie2 reference.
```

## Binning workflow

**Definition**: `binning.nf`


### Steps:

1. metagenomic binning with [Metabat2](https://bitbucket.org/berkeleylab/metabat/)
2. Quality checks with [checkM2](https://ecogenomics.github.io/CheckM/)
3. taxonomy assignment using [GTDB-TK](https://github.com/Ecogenomics/GTDBTk)
4. Dereplication with [dRep](https://github.com/MrOlm/drep) using the checkM2 metrics and skANI

If no manifest is provided this will use the contigs and preprocessed reads from the basic
workflow. Otherwise, you can provide a CSV manifest with the columns "id", "contigs", "forward" and
(optionally) "reverse" to specify the sample id, contigs, and read files.

> [!NOTE]
> When specifying a manifest please note that the file paths are assumed to be relative
> to the current project directory.

### Options:

```
~~~ Diener Lab Binning workflow ~~~

Usage:
A run using all,default parameters can be started with:
> nextflow run binning.nf --resume

An example run would look like
> nextflow run binning.nf -resume \
                          --conda_path $HOME/miniconda3/envs \
                          --min_contig_length 5000 --ani 0.95

General options:
  --data_dir [str]              The main data directory for the analysis (must contain `raw`).
  --single_end [bool]           Whether the data is single-end sequencing data.
  --preset [str]                What sequencing technology was used. Can be "illumina" or
                                "nanopore".
  ---manifest [path]            Location of a manifest containing paths to contigs and reads.
                                This location is relative to the project directory.

Binning options:
  --min_contig_length [int]     Minimum length of the contigs to include them.
  --min_bin_size [int]          Miinimum length of the binned genome.

Dereplication options:
  --ani [float]                 On what ANI to dereplicate the MAGs. The default is
                                adequate for strain level analyses.

Reference DBs:
  --checkm [path]                Location of the checkM2 uniref DB.
  --gtdb [path]                  Location of the GTDB-TK database.
```
