# :hammer: :triangular_ruler: Pipelines

This repo contains various analysis pipelines for the lab. Here are the basic
rules:

- each folder includes pipelines for a particular analysis - data type combination
- pipelines are [nextflow](https://www.nextflow.io/) workflows
- each pipeline comes with a list of conda environment files that manage the required software

## Data layout

Pipelines will usually operate from a top level project
directory structured in the following way:

```
[project root]
├─ [pipeline].nf
├─ data
│  ├─ raw
│  │  ├─ sample1_R1.fastq.gz
│  │  ├─ sample1_R2.fastq.gz
│  │  └─ ...
│  ├─ figures
│  ├─ fig1.png
│  │  └─ ...
│  └─ ...
└─ refs
   ├─ eggnog
   └─ kraken2
```

The initial raw data lives in `data/raw` and all analysis artifacts should
be written into `data/` as well. Figures go into `figures/`.

## Setup

The first step is to copy or symlink the pipeline files into the top project
directory. After that you can set up a conda environment that includes all software
for the pipeline (please see individual pipelines for variations on that).

```bash
conda env create -f conda.yml
```

Either activate the environment (usualy named after the pipeline):

```bash
conda activate metagenomics
```

or run the pipeline with the `-with-conda /my/envs/metagenomics` option (required for HPC).

## Nextflow Configuration

You may also create a [nextflow config](https://www.nextflow.io/docs/latest/config.html) either in the project
directory as `nextflow.config` or in your user HOME as `~/.nextflow/config`. A template config is
[included in this repo](nextflow.config). If you are a lab member please use the [optimized
version from the wiki](https://github.com/dienerlab/internal/wiki/configs).

To install it as a global configuration:

On the server run the following to create the config directory

```bash
mkdir ~/.nextflow
```

After that edit and copy the config:

```bash
cp /path/to/pipelines/nextlow.config ~/.nextflow/config
```

Add in your token if you want to use [Nextflow Tower](https://tower.nf) to track your pipeline.

For slurm substitute the partition name `default` with the SLURM partition.

## Run the pipeline

After setup you can test the pipeline with

```bash
nextflow run [WORKFLOW].nf -profile local -resume
```

By default this will use all available 12 CPUs and 128 GB RAM unless specified otherwise in the personal [netxflow config](https://www.nextflow.io/docs/latest/config.html#scope-executor).
