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