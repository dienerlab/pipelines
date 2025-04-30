# :hammer: :triangular_ruler: Pipelines

This repo contains various analysis pipelines for the lab. Here are the basic
rules:

- each folder includes pipelines for a particular analysis - data type combination
- pipelines are [nextflow](https://www.nextflow.io/) workflows
- each pipeline comes with a list of conda environment files that manage the required software
  or a docker image that packages the required software

## Data layout

Pipelines will usually operate from a top level project
directory structured in the following way:

```
[project root]
├─ [pipeline].nf              # optional, see setup
├─ data                       # anything that is not code
│  ├─ raw
│  │  ├─ sample1_R1.fastq.gz
│  │  ├─ sample1_R2.fastq.gz
│  │  └─ ...
│  └─ ...
├─ figures
│  ├─ fig1.png
│  └─ ...
└─ refs                       # often a symbolic link
   ├─ eggnog
   └─ kraken2
```

The initial raw data lives in `data/raw` and all analysis artifacts should
be written into `data/` as well. Figures go into `figures/`.