If you want to run metagenomic binning as well also set up the `binchecks` environment.

```bash
conda env create binchecks.yml
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
