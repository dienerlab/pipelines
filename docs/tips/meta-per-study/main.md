# Simple meta analysis

This is a strategy for a simple meta-analysis nextflow workflow. Here "simple" means
that all you individual processing steps are already expected to run on all the data
in a single study. This would be the case for Qiime2 analyses for instance because each
Qiime2 always operates on all data from a study.

Here we first have a main table that lists the studies with the required parameters/settings
for each study.

``` csv title="main table"
--8<-- "docs/tips/meta-per-study/studies.csv"
```

The nextflow pipeline then reads the CSV and passed the study information as value
in all the steps. The initial channel has one entry for each study.

``` groovy title="pipeline example"
--8<-- "docs/tips/meta-per-study/main.nf"
```

!!! note

    The individuals processng scripts don't make much sense here. Those just serve
    as an example how to inject the study parameters.

Running this will then distribute the row for each study.

``` bash
$ nextflow run main.nf

 N E X T F L O W   ~  version 25.03.1-edge

Launching `main.nf` [stoic_murdock] DSL2 - revision: 0e2f97c09b

executor >  local (6)
[ec/8110b0] process > import_data (1) [100%] 2 of 2 ✔
[12/497c3c] process > step1 (1)       [100%] 2 of 2 ✔
[1e/c9f3fe] process > visualize (2)   [100%] 2 of 2 ✔
```

This also shows an example how to globally disable a part of the pipeline (visualization)
while still retaining the cache.

``` bash
$ nextflow run main.nf -resume --visualize false

 N E X T F L O W   ~  version 25.03.1-edge

Launching `main.nf` [intergalactic_goldwasser] DSL2 - revision: 0e2f97c09b

[f9/53dc3a] process > import_data (2) [100%] 2 of 2, cached: 2 ✔
[04/179b27] process > step1 (1)       [100%] 2 of 2, cached: 2 ✔
```

You can see the injection of the library layout and the truncation parameters from the
main table.

``` bash
cat imports/*.txt
```

``` text
importing from /home/cdiener/code/pipelines/docs/tips/meta-per-study/study_1/manifest.tsv in paired-end layout.
importing from /home/cdiener/code/pipelines/docs/tips/meta-per-study/study_2/manifest.tsv in single-end layout.
```

```bash
cat steps/*.result
```

``` text
processed study 1 data with truncations of 220,200
processed study 2 data with truncations of 220,0
```