# Review/Test of binning pipelines

!!! note

    This is a very opinionated assessment. Others might disagree with some details of course.

## Criteria

One of the plans for the lab in 2025 was to add a binning pipeline to our normally
gene-centric workflows. Going through the literature I mostly focused on the insights from
some recent benchmarks:

<div class="annotate" markdown>

- [Kim et al.](https://www.nature.com/articles/s41467-026-71521-w)
- [Han et al.](https://www.nature.com/articles/s41467-025-57957-6)
- [Coleman et al.](https://www.biorxiv.org/content/10.64898/2026.04.06.712906v1)
- [Yazhini et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12636519/)
- [CAMI II benchmark](https://www.nature.com/articles/s41592-022-01431-4) (1)

</div>

1.  This is somehwhat older (2022).


Based on those my criteria for a pipeline were:

1. must include the newer deep learning binners, especially SemiBin2 and ComeBin (seem to work well in general)
2. must provide multi-sample binning and allow choosing small groups
3. should include MIMAG evaluations (so completeness, contamination, rRNAs, tRNAs)
4. should reuse mapping
5. should include chimera detection like GUNC
6. Newer bin refiners
7. should include bin dereplication
8. scalable to 1000+ samples
9. support for long reads

!!! info

    Some additional insights: The Coleman paper suggests that CheckM2 estimates are only correct for non-chimeric genomes. DasTool was shown to sometimes makes bins worse, especially the deep learning ones. The Kim paper suggests that multi-sample has diminishing returns that start with quite small sample groups (~10).

Back in 2025 there was not a single pipeline that would implement more than 2 of those, so the
idea was to make a new pipeline from scratch. But, I never got to that.

## New pipelines from 2026

In 2026 there were two pipelines that got closer to the criteria. [Tofu-MAAPO](https://github.com/ikmb/TOFU-MAaPO)
was published and [nf-core/mag](https://nf-co.re/mag/5.5.0) got a major update.

| Pipeline    | [bus factor](https://en.wikipedia.org/wiki/Bus_factor) | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|-------------|--------------------------------------------------------|---|---|---|---|---|---|---|---|---|
| Tofu-MAAPO  | 1  | ✅ | ❌ | ❌ | 😐 | ❌ | ✅ | ❌ | ✅ | ❌ |
| nf-core/mag | 3+ | ✅ | ✅ | ❌ | ✅ | ✅ | ❌ | ❌ | 😐 | ✅ |

[Bus Factor](https://en.wikipedia.org/wiki/Bus_factor) was calculated as active contributors in the last 3 months with more than 10 lines of code.

So both close. One thing that can be noticed that all the missing features are mostly post-binning (dereplication, annotations, bin refinement, GUNC). So that makes things easier.

## Pipeline testing

I tested both of the pipelines with 10 samples (paired-end) from our "larger" bechmark data set.

Data is from infant gut microbiome time series, 400+ samples from [this paper](https://doi.org/10.1038/s41591-018-0216-2). I use this because:

- fairly low depth (~10M reads per sample)
- time series so can check samples from same individual
- high amount of host reads (typical for infant microbiome)
- some ground truth about common taxa (Bifidos, low Bacteroides)
- increasing bacterial diversity across age

Config had to be adapted for Tofu-MAAPO. For nf-core/mag I added some of our standard config but I suspect it does not overwrite
the standard profile as aggressively as Tofu-MAAPO so I suspect the normal nextflow setup might have worked without adjustements.
The configs can be found in the `configs/` directory in the repo.

### Cool things and issues

#### Tofu-MAAPO

**Cool**

- has a single container for large parts of the pipeline (less disk space)
- can download from SRA directly
- supports GPUs for deep learning binners including VAMB
- quite fast
- uses containers → no setup or installation

**Issues**

<div class="annotate" markdown>

- needed several adjustments to config
- vamb want to run in exclusive mode which causes massive delays (1)
- parts of it (MagScot) still rely on GTDB 207
- databases need to be downloaded manually
- docs are lacking
- multiqc report not that informative (only read QC, no mapping)

</div>

1. They claim this is because of numpy taking up all threads but you can just limit this with
   `OMP_NUM_THREADS`. I never observed VAMB using more than the allocated CPUs. I disabled this
   in the config without any impact.

#### nf-core/mag

**Cool**

- lot's of features and parameters to tune
- supports long reads and hybrid setups
- uses containers → no setup or installation
- great community and docs (Slack, issues)
- some basic summary tables

**Issues**

<div class="annotate" markdown>

- some bugs with parameters being ignored or not working (1)
- somehwat slower because of the many binners and multi-sample binning
- no GPU support
- BUSCO could not be turned off and ran forever (24h+ for a *single* genome)
- generates a lot of redundant bins
- multiqc report not that informative (only read QC, no mapping)

</div>

1.  This is why BUSCO could not be turned off. Also had to download the GUNC DB manually as their
    parameter to do so was broken.

The retry strategy for both pipelines can be sluggish if you have actual failures.

### Summary

Both recovered around 50 dereplicated MAGs that were at least MQ with several NC MAGs. There
were only very few HQ MAGs by MIMAG criteria because the 23S gene could only be found in 5 MAGs (not
surprising because the ITS is hard to assemble). Overall pretty good results for this low depth data
set with a lot of host contamination.

Overall I liked nf-core/mag a bit more. It covered more of our features though it was also
a bit slower. However, since you only do this once this should be fine.

Neither of the tool gave that great reports over the acquired bins. Tofu-MAAPO ends with
MagScot and does no annotation of the bins. None of the tools pools the final bins in a single
directory. nf-core/mag does provide a merged bin summary from QUAST and CheckM but does not
include GUNC.

So I think the best strategy would be the following:

1. Adjust our MGX pipeline to output sample and assembly sheets for nf-core/mag
2. Use nf-core/mag for binning
3. Add a small summary pipeline on top of nf-core/mag output to dereplicate and annotate
   according to MIMAG

It might also be fun to assign groups on similarity of SingleM OTU profiles to improve
multi-sample binning.
